#!/usr/bin/perl

use strict;
use Switch;
use warnings;
use Getopt::Long;
use File::stat;
use File::Copy;
use Time::localtime;
use Pod::Usage;
use Cwd qw(abs_path);
use File::Basename;

# This BEGIN block opens the config file
# to find where uSeq's custom perllib is,
# then adds it to @INC so that any custom
# modules can be loaded
BEGIN{
	my $configFile = sprintf("%s/config.txt",dirname(abs_path($0)));
	if(defined $ENV{'JOB_ID'}) {
		# script was launched using SGE
		# config.txt must be obtained from $JOB_SCRIPT
		$configFile = $ENV{'USEQ_CONFIG'};
	}
	my $config = do($configFile);

	die "Could not parse config file: $@" if $@;
	die "Could not read config file: $!" 	unless defined $config;
	die "Could not run config file"				unless $config;
	unshift @INC, $config->{uSeqPerllib};
}

use MicroSeqFunc qw(:all);
use JobManager::SGE qw(:all);
use DirOps qw(:all);
use Certificate qw(:all);
use Timer;
use Logger;

sub loadParams;
sub checkParams;
sub tetrascanParamComments;
sub scanAndAlign;
sub joinMsiFiles;
sub markDups;
sub setLevel;

sub checkGenOpts;
sub setSgeOpts;
sub loadSeqDir;
sub setTetrascanOpts;
sub prepareOutDirs;
sub getScanAndAlignCmd;
sub getProfilerCmd;
sub getMarkDupCmd;

my $prog_paths = getProgramPaths();
my $virtual_free_limit = 16; #in GB, should be set from config file in future

#Load and check parameters
my ($gen_params,$sge_params,$tet_params) = loadParams();
my $files = checkParams($gen_params,$sge_params,$tet_params);

#start log file and print params to log
my $log = Logger->new($gen_params->{target_dir},"MicroSeq.pl");
my $timer = Timer->new();
print "A MicroSeq log file for this run can be found at ",$log->fname(),"\n";

my ($tOptExplanations,$tChoiceExplanations) = tetrascanParamComments();
$log->recordParams($gen_params,"General");
$log->recordParams($tet_params,"Tetrascan",$tOptExplanations,$tChoiceExplanations);
if(defined $gen_params->{distributed}){
	$log->recordParams($sge_params,"SGE") if $gen_params->{distributed} eq "sge";
}

#print MicroSeq mode
if(defined $gen_params->{local})
	{
		print "Running MicroSeq in local mode\n";
		$log->record("Running MicroSeq in local mode");
	}
else
	{
		printf("Running MicroSeq in distributed mode on %s\n",$gen_params->{distributed});
		$log->record(sprintf("Running MicroSeq in distributed mode on %s",$gen_params->{distributed}));
	}

#Prepare directories for STDOUT and STDERR from MicroSeq subprocesses	
#this may fail if run using pe1/pe2 flags, should be tested and fixed
prepareOutDirs($gen_params,$log) if defined $gen_params->{dir};

$gen_params->{counts_dir} = sprintf("%s/counts/",$gen_params->{target_dir});
$log->mkDir($gen_params->{counts_dir});
#mkdir $gen_params->{counts_dir};
#my $action = sprintf("Created directory %s",$gen_params->{counts_dir});
# $log->record($action);


#send files to scanAndAlign, combine them afterwards, create microsatellite profiles and cleanup directory
if($gen_params->{level} <= 3)
{
	scanAndAlign($gen_params,$tet_params,$sge_params,$files,$prog_paths,$log);
}
markDups($gen_params,$sge_params,$prog_paths,$log);

printf("Completed MicroSeq in %s\n", $timer->elapsed());

##################################SUBROUTINES###############################################

sub scanAndAlign
{
	my ($gen_params,$tet_params,$sge_params,$files,$prog_paths,$log) = @_;

	my $cmd;
	my $out;
	my $jobs;
	my $job_id;
	my $lib_name;

	#scans sequence files for microsatellites, align them to reference genome,
	#and partition by chromosome
	foreach my $library (keys %{$files})
	{
		$lib_name = sprintf($files->{$library}->{local_pair_format},$files->{$library}->{file_id});
		print "Processing $lib_name";
		if($sge_params->{lb_scheme} == 1)
			{
				$cmd = getScanAndAlignCmd($sge_params,$tet_params,$gen_params,$files->{$library},$prog_paths);
				$out = $log->runCommand($cmd);
				if(defined $gen_params->{distributed})
					{
						if($gen_params->{distributed} eq "sge")
							{
								($jobs, $job_id) = addSgeJob($jobs,$out,$lib_name);
								print " as $job_id";
							}
					}
			}
		elsif($sge_params->{lb_scheme} == 2)
			{
				my $job_dependency_list;
				
				#scan step
				$cmd = getScanCmd($sge_params,$tet_params,$gen_params,$files->{$library},$prog_paths);
				$out = $log->runCommand($cmd);
				#add jobs from scan step to job dependency list for downstream 
				#jobs to wait until the scan jobs finish
				if($gen_params->{distributed} eq "sge")
					{
						($jobs, $job_id) = addSgeJob($jobs,$out,$lib_name);
						push @{$job_dependency_list}, $job_id;
					}

				#aln step
				$cmd = getAlnCmd($sge_params,$gen_params,$files->{$library},$prog_paths,$job_dependency_list);
				$out = $log->runCommand($cmd);
				if($gen_params->{distributed} eq "sge")
					{
						($jobs, $job_id) = addSgeJob($jobs,$out,$lib_name);
						push @{$job_dependency_list}, $job_id;
					}

				#post-processing step
				$cmd = getPProCmd($sge_params,$gen_params,$files->{$library},$prog_paths,$job_dependency_list);
				$out = $log->runCommand($cmd);
				if($gen_params->{distributed} eq "sge")
					{
						($jobs, $job_id) = addSgeJob($jobs,$out,$lib_name);
						push @{$job_dependency_list}, $job_id;
					}
				printf(" as jobs %s", join ", ", @{$job_dependency_list});
				push @{$gen_params->{dependencies}}, @{$job_dependency_list};
			}
		print "\n";

	}
	sgeWait($jobs) if $gen_params->{distributed} eq "sge";
	printSgeCompletionInfo($jobs,$log);
	$jobs = undef;
	print "All fastq lanes have been scanned, aligned and reindexed\n";
}

sub markDups
{
	my $gen_params = shift;
	my $sge_params = shift;
	my $prog_paths = shift;
	my $log = shift;
	
	my $cmd;
	my $out;
	my $jobs;
	my $job_id;
	
	print "Profiling microsatellites in ",$gen_params->{base_project_id};
	$cmd = getMarkDupCmd($sge_params,$gen_params,$prog_paths);
	$out = $log->runCommand($cmd);
	if(defined $gen_params->{distributed})
		{
			if($gen_params->{distributed} eq "sge")
				{
					($jobs, $job_id) = addSgeJob($jobs,$out,$gen_params->{base_project_id});
					print " as $job_id";
				}
		}
	print "\n";			

	if(defined $gen_params->{distributed})
		{
			sgeWait($jobs) if $gen_params->{distributed} eq "sge";
			printSgeCompletionInfo($jobs,$log);
			$jobs = undef;
		}
	
	print "Completed marking duplicates\n";
}

sub checkParams
{
	my $gen_params = shift;
	my $sge_params = shift;
	my $tet_params = shift;
	
	my $files;
	
	#check options and load directories
	checkGenOpts($gen_params);
	setLevel($gen_params);

	$tet_params = setTetrascanOpts($tet_params);

	if(defined $gen_params->{dir})
	{
		$files = loadSeqDir($gen_params, $tet_params,$log);
	}
	$sge_params = setSgeOpts($sge_params, $gen_params);

	return $files;
}

#determines what steps of the pipeline need to be run based on the stage
#requested by the user
sub setLevel
{
	my $params = shift;
	
	switch($params->{stage})
		{
			case "scan"			{$params->{level} = 0}
			case "align"		{$params->{level} = 1}
			case "reindex"	{$params->{level} = 2}
			case "split"		{$params->{level} = 3}
			case "sort"			{$params->{level} = 4}
			case "profile"	{$params->{level} = 5}
			case "stats"		{$params->{level} = 6}
			else						{$params->{level} = 0}
		}
}

sub loadParams
{
	my ($gen_params,$sge_params,$tet_params);

	GetOptions
	(
		"dir:s"													=> \$gen_params->{dir},
		"I|casava_18_input|new"					=> sub {$gen_params->{casava_version} = "new"},
		"i|pre_casava_18_input|old"			=> sub {$gen_params->{casava_version} = "old"},
		"pe1|1|file:s"									=> \$gen_params->{pe1},
		"pe2|2:s"												=> \$gen_params->{pe2},
		"genome:s"											=> \$gen_params->{genome},
		"target:s"											=> \$gen_params->{target_dir},
		"offset:s"											=> \$gen_params->{offset_index}, #optional--will be inferred from --genome
		"msdb|db|ms_db:s"						 		=> \$gen_params->{ms_db}, #optional--will be inferred from --genome
		"split" 												=> \$gen_params->{split},
		"sge"					 									=> sub {$gen_params->{distributed} = "sge"},
		"local" 												=> \$gen_params->{local},
		"stage=s"												=> \$gen_params->{stage},
		"clean"													=> \$gen_params->{clean},
		"cleanBAM|clean_bam"						=> \$gen_params->{clean_bam},
		"zip"														=> \$gen_params->{zip_output},
		"keep_all"											=> \$gen_params->{keep_all},
		"project_id:s"									=> \$gen_params->{project_id},
		"decentralized"									=> \$gen_params->{decentralized},
		"keep_duplicates"								=> \$gen_params->{keep_duplicates}, #profiles will include duplicate reads

		"virtual_free:s"								=> \$sge_params->{virtual_free},
		"threads:i"											=> \$sge_params->{threads},
		"node|q:s"											=> \$sge_params->{node},
		"queue_threshold_type:s"				=> \$sge_params->{threshold_type},
		"queue_threshold:f"							=> \$sge_params->{threshold},
		"lb_scheme:i"										=> \$sge_params->{lb_scheme},
		
		"allms|allMS|a"									=> \$tet_params->{print_all_MS},
		"onlyms|onlyMS|M"								=> sub {$tet_params->{print_all_ms} = 0},
		"maxN|N:i"											=> \$tet_params->{max_n},
		"msUnit|unit|u:i"								=> \$tet_params->{min_unit},
		"msLength|length|l:i"						=> \$tet_params->{min_length},
		"motifMin|min|m:i"							=> \$tet_params->{min_usize},
		"motifMax|max|x:i"							=> \$tet_params->{max_usize},
		"minQual|min_qual|Q:i"					=> \$tet_params->{min_qual},
		"barcodeLength|t:i"							=> \$tet_params->{barcodeLength},
		"sanger|a"											=> sub {$tet_params->{qual_format} = "a"},
		"solexa|s"											=> sub {$tet_params->{qual_format} = "s"},
		#set the trim type and the minqual with the same command
		"BWA|B:i"												=> sub {$tet_params->{trim_type} = "B";$tet_params->{min_qual} = $_[1]},
		"lastGood|g:i"									=> sub {$tet_params->{trim_type} = "g";$tet_params->{min_qual} = $_[1]},
		"firstBad|b:i"									=> sub {$tet_params->{trim_type} = "b";$tet_params->{min_qual} = $_[1]},
		"noTrim|n"											=> sub {$tet_params->{trim_type} = "n"},
	);
	
	return ($gen_params,$sge_params,$tet_params);
}

sub  tetrascanParamComments
{
	#Explanations for different option names
	my $optExplanations = {
		"M" => "print only reads with microsatellites",
		"N" => "max number of Ns in read",
		"u" => "min repeats of microsatellite motif",
		"l" => "min microsatellite tract length",
		"m" => "min microsatellite motif length",
		"x" => "max microsatellite motif length",
		"Q" => "min basecall qual (independent of Phred score)",
		"t" => "length of barcode sequence at the beginning of every read"
	};

	#Explanations for different option choices (e.g. which qual format was specified)
	my $choiceExplanations = {
		"a" => "Sanger format (Phred + 33)",
		"s" => "Illumina format (Phred + 64)",
		"B" => "BWA-style trimming",
		"g" => "trim to last base pair greater than Q",
		"b" => "trim to first base pair less than Q",
		"n" => "do not trim",
	};

	return ($optExplanations,$choiceExplanations);
}

sub getMarkDupCmd
{
	my ($sge_params,$gen_params,$paths) = @_;
	my $err_file = sprintf("%s%s.md.err",$gen_params->{err_dir},$gen_params->{base_project_id});
	my $out_file = sprintf("%s%s.md.out",$gen_params->{out_dir},$gen_params->{base_project_id});
	my $command = "";

	if(defined $gen_params->{distributed})
		{
			if($gen_params->{distributed} eq "sge")
				{
					$command = sprintf("qsub -S /usr/bin/perl -l virtual_free=%s -o %s -e %s -v USEQ_CONFIG=%s ", 
															$sge_params->{virtual_free},$out_file,$err_file, $paths->{config});
					if(defined $sge_params->{node})
						{
							$command .= sprintf("-q %s ", $sge_params->{node}) if defined $sge_params->{node};
						}
					elsif(`hostname` =~ m/wigclust/)
						{
							if($gen_params->{level} <= 4)
							{
								my $msi_dir = sprintf("%s/MSI/",$gen_params->{target_dir});
								opendir(MSI,$msi_dir) or $log->confess("Could not open ".$msi_dir.": $!");
								my $node_resident_count;
								FILE:while(my $file = readdir MSI)
								{
									next FILE if $file =~ m/^\.{1,2}$/;
									$file = sprintf("%s%s",$msi_dir,$file);
									my $linked_file = readlink($file);
									my $key = `hostname`;
									chomp $key;
									if(defined $linked_file)
										{
											$key = $1 if $linked_file =~ m/\/mnt\/(wigclust\d+)\//;
											$key .= ".cshl.edu";
										}
									$node_resident_count->{$key}++;
								}

								my $max_count = 0;
								my $max_node;
								foreach(keys %{$node_resident_count})
								{
									if($node_resident_count->{$_} > $max_count)
										{
											$max_count = $node_resident_count->{$_};
											$max_node = $_;
										}
								}
							
								my $local_queue = sprintf("all.q@%s",$max_node);
								$command .= sprintf("-q %s ",$local_queue);
							}
							else
							{
								my $individual_name = [split /\//, $gen_params->{target_dir}];
								$individual_name = $individual_name->[$#{$individual_name}];

								my $msi_file = sprintf("%s%s.mdup.msi.bam",$gen_params->{target_dir},$individual_name);
								my $linked_file = readlink($msi_file);
								
								my $msi_node = $1 if $linked_file =~ m/\/mnt\/(wigclust\d+)\//;
								my $local_queue = sprintf("all.q@%s",$msi_node);
								$command .= sprintf("-q %s ",$local_queue);
							}
						}
					elsif(defined $sge_params->{queue_threshold})
						{
							my $queues = getGoodQueues($sge_params->{threshold_type},$sge_params->{threshold});
							$command .= sprintf("-q %s ", $queues) if defined $queues;
						}
				}
		}

	$command .= sprintf("%smsMarkDups2.pl --target_dir %s --db %s --level %d",
											$paths->{uSeq},$gen_params->{target_dir},$gen_params->{ms_db},$gen_params->{level});
	
	if(defined $sge_params->{virtual_free} and $gen_params->{distributed} eq "sge")
	{
		$command .= sprintf(" --virtual_free %s",$sge_params->{virtual_free});
	}
		
	foreach(sort qw(zip_output keep_all decentralized keep_duplicates))
	{
		$command .= sprintf(" --%s",$_) if defined $gen_params->{$_};
	}
	$command .= sprintf(" --project_id %s",$gen_params->{project_id}) if $gen_params->{project_id};

	if(defined $gen_params->{local} and !defined $gen_params->{distributed})
		{
			$command .= sprintf(" 1>%s 2>%s",$out_file,$err_file);
		}
	return $command;
}

sub getScanCmd
{
	my ($sge_params,$tet_params,$gen_params,$file,$paths) = @_;

	my $err_file = sprintf("%s%s_%s.scan.err",$gen_params->{err_dir}, $gen_params->{base_project_id},$file->{file_id});
	my $out_file = sprintf("%s%s_%s.scan.out",$gen_params->{out_dir}, $gen_params->{base_project_id},$file->{file_id});

	my $command = "";
	if(defined $gen_params->{distributed})
		{
			#microsatellite scanning is always single-threaded
			if($gen_params->{distributed} eq "sge")
				{
					$command = sprintf("qsub -S /usr/bin/perl -l virtual_free=%s -o %s -e %s -v USEQ_CONFIG=%s ", 
															$sge_params->{virtual_free},$out_file,$err_file, $paths->{config});
					if(defined $sge_params->{node})
						{
							$command .= sprintf("-q %s ", $sge_params->{node}) if defined $sge_params->{node};
						}
					elsif(`hostname` =~ m/wigclust/)
						{
							#Assumes that both paired end files are on the same node, which should always be true
							my $file_prefix;
							if($file->{suffix} =~ m/\.bam$/)
								{
									$file_prefix = sprintf($file->{local_pair_format}, $file->{file_id});
								}
							else
								{
									$file_prefix = sprintf($file->{local_solo_format}, $file->{file_id},1);
								}
							my $linked_file = sprintf("%s%s%s",$gen_params->{dir},$file_prefix,$file->{suffix});
							$linked_file .= ".gz" if $file->{zipped} and $file->{suffix} !~ m/\.bam$/;
							my $local_file = readlink($linked_file);

							my $local_node;
							my $local_queue;
							if(defined $local_file)
								{
									$local_node = $1 if $local_file =~ m/\/mnt\/(wigclust\d+)\//;
									$local_node .= ".cshl.edu";
								}
							else
								{
									$local_node = `hostname`;
									chomp $local_node;
								}
							$local_queue = sprintf("all.q@%s",$local_node);
							$command .= sprintf("-q %s ",$local_queue);
						}
					elsif(defined $sge_params->{queue_threshold})
						{
							my $queues = getGoodQueues($sge_params->{threshold_type},$sge_params->{threshold});
							$command .= sprintf("-q %s ", $queues) if defined $queues;
						}
				}
		}
	$command .= sprintf("%smsScan.pl --fastq %s", $paths->{uSeq}, $file->{file_id});
								
	#general params
	foreach(sort qw(dir casava_version target_dir project_id))
		{
			$command .= sprintf(" --%s %s", $_, $gen_params->{$_}) if defined $gen_params->{$_};
		}
	foreach(sort qw(clean zip_output keep_all decentralized))
		{
			$command .= sprintf(" --%s",$_) if defined $gen_params->{$_};
		}
		
	#file params
	foreach(sort qw(count solo_format pair_format qual_format suffix))
		{
					$command .= sprintf(" --%s %s", $_, $file->{$_}) if defined $file->{$_};
		}
	foreach(sort qw(paired zipped isDir))
		{
			$command .= sprintf(" --%s",$_) if defined $file->{$_};
		}

	#tetrascan params
	foreach(keys %{$tet_params})
		{
			$command .= sprintf(" --%s %s", $_, $tet_params->{$_});
		}

	if(defined $gen_params->{local} and !defined $gen_params->{distributed})
		{
			$command .= sprintf(" 1>%s 2>%s",$out_file,$err_file);
		}
	return $command;
}

sub getAlnCmd
{
	my ($sge_params,$gen_params,$file,$paths,$job_dependencies) = @_;

	my $err_file = sprintf("%s%s_%s.aln.err",$gen_params->{err_dir}, $gen_params->{base_project_id},$file->{file_id});
	my $out_file = sprintf("%s%s_%s.aln.out",$gen_params->{out_dir}, $gen_params->{base_project_id},$file->{file_id});
	my $command = "";
	my $threads = $sge_params->{threads};
	$threads /= 2 if $file->{paired};

	if(defined $gen_params->{distributed})
		{
			if($gen_params->{distributed} eq "sge")
				{
					$command = sprintf("qsub -S /usr/bin/perl -l virtual_free=%s -pe threads %d -o %s -e %s -v USEQ_CONFIG=%s ", 
															$sge_params->{virtual_free},$threads,$out_file,$err_file, $paths->{config});
					$command .= sprintf("-hold_jid %s ", join ",", @{$job_dependencies});
					if(defined $sge_params->{node})
						{
							$command .= sprintf("-q %s ", $sge_params->{node}) if defined $sge_params->{node};
						}
					elsif(`hostname` =~ m/wigclust/)
						{
						#Assumes that both paired end files are on the same node, which should always be true
							my $file_prefix;
							if($file->{suffix} =~ m/\.bam$/)
								{
									$file_prefix = sprintf($file->{local_pair_format}, $file->{file_id});
								}
							else
								{
									$file_prefix = sprintf($file->{local_solo_format}, $file->{file_id},1);
								}
							my $linked_file = sprintf("%s%s%s",$gen_params->{dir},$file_prefix,$file->{suffix});
							$linked_file .= ".gz" if $file->{zipped} and $file->{suffix} !~ m/\.bam$/;
							my $local_file = readlink($linked_file);

							my $local_node;
							my $local_queue;
							if(defined $local_file)
								{
									$local_node = $1 if $local_file =~ m/\/mnt\/(wigclust\d+)\//;
									$local_node .= ".cshl.edu";
								}
							else
								{
									$local_node = `hostname`;
									chomp $local_node;
								}
							$local_queue = sprintf("all.q@%s",$local_node);
							$command .= sprintf("-q %s ",$local_queue);
						}
					elsif(defined $sge_params->{queue_threshold})
						{
							my $queues = getGoodQueues($sge_params->{threshold_type},$sge_params->{threshold});
							$command .= sprintf("-q %s ", $queues) if defined $queues;
						}
				}
		}
		
	$command .= sprintf("%smsAlign.pl --fastq %s --threads %d", $paths->{uSeq}, $file->{file_id}, $threads);
		
	#general params
	foreach(sort qw(genome dir target_dir project_id casava_version))
		{
			$command .= sprintf(" --%s %s", $_, $gen_params->{$_}) if defined $gen_params->{$_};
		}
	foreach(sort qw(clean zip_output keep_all decentralized))
		{
			$command .= sprintf(" --%s",$_) if defined $gen_params->{$_};
		}
		
	#file params
	foreach(sort qw(count solo_format pair_format suffix))
		{
					$command .= sprintf(" --%s %s", $_, $file->{$_}) if defined $file->{$_};
		}
	foreach(sort qw(paired zipped isDir))
		{
			$command .= sprintf(" --%s",$_) if defined $file->{$_};
		}

	if(defined $gen_params->{local} and !defined $gen_params->{distributed})
		{
			$command .= sprintf(" 1>%s 2>%s",$out_file,$err_file);
		}
	return $command;
}

sub getPProCmd
{
	my ($sge_params,$gen_params,$file,$paths,$job_dependencies) = @_;
	my $err_file = sprintf("%s%s_%s.ppro.err",$gen_params->{err_dir}, $gen_params->{base_project_id},$file->{file_id});
	my $out_file = sprintf("%s%s_%s.ppro.out",$gen_params->{out_dir}, $gen_params->{base_project_id},$file->{file_id});
	my $command = "";

	if(defined $gen_params->{distributed})
		{
			if($gen_params->{distributed} eq "sge")
				{
					$command = sprintf("qsub -S /usr/bin/perl -l virtual_free=%s -o %s -e %s -v USEQ_CONFIG=%s ", 
															$sge_params->{virtual_free},$out_file,$err_file, $paths->{config});
				#this should only wait for scan jobs running on the same lane
					$command .= sprintf("-hold_jid %s ", join ",", @{$job_dependencies});
					if(defined $sge_params->{node})
						{
							$command .= sprintf("-q %s ", $sge_params->{node}) if defined $sge_params->{node};
						}
					elsif(`hostname` =~ m/wigclust/)
						{
						#Assumes that both paired end files are on the same node, which should always be true
							my $file_prefix;
							if($file->{suffix} =~ m/\.bam$/)
								{
									$file_prefix = sprintf($file->{local_pair_format}, $file->{file_id});
								}
							else
								{
									$file_prefix = sprintf($file->{local_solo_format}, $file->{file_id},1);
								}
							my $linked_file = sprintf("%s%s%s",$gen_params->{dir},$file_prefix,$file->{suffix});
							$linked_file .= ".gz" if $file->{zipped} and $file->{suffix} !~ m/\.bam$/;
							my $local_file = readlink($linked_file);

							my $local_node;
							my $local_queue;
							if(defined $local_file)
								{
									$local_node = $1 if $local_file =~ m/\/mnt\/(wigclust\d+)\//;
									$local_node .= ".cshl.edu";
								}
							else
								{
									$local_node = `hostname`;
									chomp $local_node;
								}
							$local_queue = sprintf("all.q@%s",$local_node);
							$command .= sprintf("-q %s ",$local_queue);
						}
					elsif(defined $sge_params->{queue_threshold})
						{
							my $queues = getGoodQueues($sge_params->{threshold_type},$sge_params->{threshold});
							$command .= sprintf("-q %s ", $queues) if defined $queues;
						}
				}
		}
			
	$command .= sprintf("%smsPostProcess2.pl --fastq %s", $paths->{uSeq}, $file->{file_id});
	if(defined $sge_params->{virtual_free} and $gen_params->{distributed} eq "sge") {
		$command .= sprintf(" --virtual_free %s",$sge_params->{virtual_free});
	}
		
	#general params
	foreach(sort qw(genome dir target_dir project_id offset_index)) {
			$command .= sprintf(" --%s %s", $_, $gen_params->{$_}) if defined $gen_params->{$_};
	}
	foreach(sort qw(split clean zip_output keep_all decentralized)) {
			$command .= sprintf(" --%s",$_) if defined $gen_params->{$_};
	}
		
	#file params
	foreach(sort qw(count solo_format pair_format suffix)) {
					$command .= sprintf(" --%s %s", $_, $file->{$_}) if defined $file->{$_};
	}
	foreach(sort qw(paired zipped isDir)) {
			$command .= sprintf(" --%s",$_) if defined $file->{$_};
	}
		
	if(defined $gen_params->{local} and !defined $gen_params->{distributed}) {
			$command .= sprintf(" 1>%s 2>%s",$out_file,$err_file);
	}
	return $command;
}

sub getScanAndAlignCmd
{
	my ($sge_params,$tet_params,$gen_params,$file,$paths) = @_;
	my $file_prefix = sprintf($file->{local_pair_format},$file->{file_id});
	my $err_file = sprintf("%s%s.sanda.err",$gen_params->{err_dir},$file_prefix);
	my $out_file = sprintf("%s%s.sanda.out",$gen_params->{out_dir},$file_prefix);
	my $command = "";
	my $threads = $sge_params->{threads};
	$threads /= 2 if $file->{paired} and !$file->{isDir};

	if(defined $gen_params->{distributed})
		{
			if($gen_params->{distributed} eq "sge")
				{
					$command = sprintf("qsub -S /usr/bin/perl -l virtual_free=%s -pe threads %d -o %s -e %s -v USEQ_CONFIG=%s", 
															$sge_params->{virtual_free},$threads,$out_file,$err_file, $paths->{config});
					if(defined $sge_params->{node})
						{
							$command .= sprintf("-q %s ", $sge_params->{node}) if defined $sge_params->{node};
						}
					elsif(`hostname` =~ m/wigclust/)
						{
						#Assumes that both paired end files are on the same node, which should always be true
							my $linked_file;
							my $local_file;
							if(!$file->{isDir})
								{
									$file_prefix = sprintf($file->{local_solo_format}, $file->{file_id},1) if !$file->{isDir};
									$linked_file = sprintf("%s%s%s",$gen_params->{dir},$file_prefix,$file->{suffix});
									$linked_file .= ".gz" if $file->{zipped} and !$file->{isDir};
									$local_file = readlink($linked_file);
								}
							else
								{
									$linked_file = sprintf("%s%s",$gen_params->{dir},$file->{file_id});
									$local_file = readlink($linked_file);
								}
								
							my $local_node;
							my $local_queue;
							if(defined $local_file)
								{
									$local_node = $1 if $local_file =~ m/\/mnt\/(wigclust\d+)\//;
									$local_node .= ".cshl.edu";
								}
							else
								{
									$local_node = `hostname`;
									chomp $local_node;
								}
							$local_queue = sprintf("all.q@%s",$local_node);
							$command .= sprintf("-q %s ",$local_queue);
						}
					elsif(defined $sge_params->{queue_threshold})
						{
							my $queues = getGoodQueues($sge_params->{threshold_type},$sge_params->{threshold});
							$command .= sprintf("-q %s ", $queues) if defined $queues;
						}
				}
		}
	
	my $sanda_script = "msScanAndAlign.pl";
	$sanda_script = "msScanAndAlignDir.pl" if $file->{isDir};
	
	print $sanda_script,"\n";
	print $command,"\n";
	
	$command .= sprintf("%s%s --fastq %s --threads %d",$paths->{uSeq},$sanda_script,$file->{file_id},$threads);
						
	foreach(keys %{$tet_params})
		{
			$command .= sprintf(" --%s %s", $_, $tet_params->{$_});
		}
		
	#general params
	foreach(sort qw(genome dir casava_version project_id offset_index level target_dir))
		{
					$command .= sprintf(" --%s %s", $_, $gen_params->{$_}) if defined $gen_params->{$_};
		}
	foreach(sort qw(split clean zip_output clean_bam keep_all decentralized))
		{
			$command .= sprintf(" --%s",$_) if defined $gen_params->{$_};
		}
		
	#file params
	foreach(sort qw(count solo_format pair_format qual_format suffix))
		{
					$command .= sprintf(" --%s %s", $_, $file->{$_}) if defined $file->{$_};
		}
	foreach(sort qw(paired zipped))
		{
			$command .= sprintf(" --%s",$_) if defined $file->{$_};
		}

	if(defined $gen_params->{local} and !defined $gen_params->{distributed})
		{
			$command .= sprintf(" 1>%s 2>%s",$out_file,$err_file);
		}
	
	return $command;
}

sub prepareOutDirs
{
	my $params = shift;
	my $log = shift;
	
	if(defined $params->{distributed})
		{
			if($params->{distributed} eq "sge")
				{
					$params->{out_dir} = sprintf("%sSGEOUT/",$params->{target_dir});
					$params->{err_dir} = sprintf("%sSGEERR/",$params->{target_dir});
				}
		}
	else
		{
			$params->{out_dir} = sprintf("%sSTDOUT/",$params->{target_dir});
			$params->{err_dir} = sprintf("%sSTDERR/",$params->{target_dir});
		}

	foreach(qw(out_dir err_dir))
		{
			mkdir $params->{$_};
			my $action = sprintf("Created directory %s",$params->{$_});
			$log->record($action);
		}
}

sub setTetrascanOpts
{
	my $params = shift;
	
	#set defaults
	$params->{print_all_MS} ||= 0; #default is to only print reads with microsatellites and their mate pairs
	$params->{max_n}			  ||= 3; #default max number of Ns in a read is 3
	$params->{min_unit}		  ||= 3; #default minimum motif number is 3
	$params->{min_length}	  ||= 8; #default minimum microsatellite length is 8
	$params->{min_usize}	  ||= 1; #default minimum motif length is 1
	$params->{max_usize}	  ||= 6; #default maximum motif length is 6
	$params->{barcodeLength} ||= 0; #default barcode length is 0
	$params->{min_qual}			||= 20; #default minimum trimming quality is 20
	$params->{trim_type}    ||= "B"; #default is BWA-style trimming
	
	#warnings
	if($params->{min_usize} > $params->{max_usize})
	{
		warn "min unit size is greater than max unit size, switching them around and trying again\n";
		my $temp = $params->{min_usize};
		$params->{min_usize} = $params->{max_usize};
		$params->{max_usize} = $temp;
	}

	#fatal error messages
	die "Error! maxN cannot be a negative number\n" 		if $params->{max_n} < 0;
	die "Error! msUnit must be greater than 1\n"				if $params->{min_unit} <= 1;
	die "Error! msLength must be greater than 1\n"	 		if $params->{min_length} <= 1;
	die "Error! motifMin must be at least 1\n"					if $params->{min_usize} < 1;
	die "Error! motifMax cannot be larger than 100\n"		if $params->{max_usize} > 100; #need to change this so that max_usize is limited to largest integer possible on computer	

	return $params;
}

sub setSgeOpts
{
	my $params = shift;
	my $genParams = shift;
	my $log = shift;
	
	#set defaults
	$params->{lb_scheme} ||= 1;
	$params->{virtual_free} ||= "6G";
	$params->{threads} = 8;
		
	#check load balance scheme
	if($params->{lb_scheme} > 2)
		{
			warn "Unknown load balance scheme ID, using default scheme\n";
			$params->{lb_scheme} = 1;
		}

	#set special sge options for wigclust
	if(defined $gen_params->{distributed})
		{
			if(`hostname` =~ m/wigclust/ and $gen_params->{distributed} eq 'sge')  #if running sge on wigclust, by default, make sure nodes have at least 10% df
				{
					if(!defined $params->{threshold_type} and !defined $params->{threshold})
						{
							$params->{threshold_type} = "percent";
							$params->{threshold} = 90;
# log file not initialized before setSgeOpts
#							my $action = sprintf("Running on %s, maximum of %d %s of space may be occupied on the node",$host,$params->{threshold},$params->{threshold_type});
#							logAction($action,$log);
						}

					if((defined $params->{threshold_type} and !defined $params->{threshold}) or (!defined $params->{threshold_type} and defined $params->{threshold}))
						{
							die "You must specify a threshold and value to threshold on using --queue_threshold and --queue_threshold_type for MicroSeq to limit nodes based on free disk space\n";
						}
					$gen_params->{decentralized} = 1;
				}
		}

	#check to see if virtual_free ends with M or G, if it doesn't, guess one
	my $temp;
	if($params->{virtual_free} =~ /^\d+$/)
		{
			($params->{virtual_free} <= $virtual_free_limit) ? ($params->{virtual_free} .="G") : ($params->{virtual_free} .= "M");
		}
	elsif($params->{virtual_free} =~ m/G$/)
		{
			$temp = $1 if $params->{virtual_free} =~ m/^(.*)G$/;
			die "virtual_free is too large (virtual_free limit = $virtual_free_limit GB)\n" if $temp > $virtual_free_limit;
		}
	elsif($params->{virtual_free} =~ m/M$/)
		{
			$temp = $1 if $params->{virtual_free} =~ m/^(.*)M$/;
			die "virtual_free is too large (virtual_free limit = $virtual_free_limit GB)\n" if ($temp / 1000) > $virtual_free_limit;
		}
	else #if character is unrecognized, die
		{
			die "unrecgonized virtual free specification ",$params->{virtual_free},"\n";
		}
	return $params
}

sub loadSeqDir
{
	my $params = shift;
	my $tetrascan_params = shift;
	my $log = shift;
	my $files;
	
	#get sequence files if sequence is from earlier versions of casava pipeline or are grouped in files by lane
	if(!defined $params->{casava_version} or $params->{casava_version} eq "old")
		{
			opendir SEQDIR, $params->{dir} or die "Could not open ",$params->{dir},": $!\n";
			$files = {map {$_ => undef} grep {!/^\.{1,2}$/ && /(fastq|txt|fq|bam)/ && !/(mod|ms|relationship|hdr|dis|seq\.noise|dups|md|sam)/ } readdir SEQDIR};
			closedir SEQDIR;
		}
	elsif($params->{casava_version} eq "new")
		{
			#if output is from illumina's casava 1.8, sequences are stored in directory per lane instead of file per lane
			opendir SEQDIR, $params->{dir} or die "Could not open ",$params->{dir},": $!\n";
			foreach(@{[grep { !/^\.{1,2}/ && !/certificates/ && !/logs/ & !/counts/ && !/SGE*/} readdir SEQDIR]})
				{
					my $obj = sprintf("%s%s",$params->{dir},$_);
					if(-l $obj)
						{
							if(-d readlink($obj))
								{
									$files->{$_}->{isDir} = 1;
								}
						}
					elsif(-d $obj)
						{
							$files->{$_}->{isDir} = 1 if($obj !~ m/(STD|SGE)(OUT|ERR)/ and $obj !~ m/logs|stats|ModMs|MSI|counts|certificates/);
						}
				}
			closedir SEQDIR;
		}
	
	#guess file naming format and get basenames
	my $file_ids;
	my $fname_format;
	
	#file name has format s_%d_%d_sequence where first %d is lane number and second %d is pair
	#number, e.g. s_1_1_sequence.txt and s_1_2_sequence.txt or s_1_1_sequence.txt.gz
	foreach my $fname (keys %{$files})
		{
			if($fname =~ m/^s_(\d+)_(1|2)_sequence(\.\w+)+$/)
				{
					$files->{$fname}->{local_solo_format} = "s_%d_%d_sequence";
					$files->{$fname}->{local_pair_format} = "s_%d_sequence";
					$files->{$fname}->{solo_format} = "s_\@d_\@d_sequence";
					$files->{$fname}->{pair_format} = "s_\@d_sequence";
					$files->{$fname}->{fname_format} = qr/^s_(\d+)_(1|2)_sequence(\.\w+)+$/;
				}
			if($fname =~ m/^.+(1|2)(\.\w+)+$/ and $fname !~ m/\.bam$/) #BAM files will contain both reads in a read pair in the same file, so won't have the flowcell_(1/2) file naming convention
				{
					$files->{$fname}->{local_solo_format} = "%s%d";
					$files->{$fname}->{local_pair_format} = "%s";
					$files->{$fname}->{solo_format} = "\@s\@d";
					$files->{$fname}->{pair_format} = "\@s";
					$files->{$fname}->{fname_format} = qr/^(.+)(1|2)(\.\w+)+/;
				}
			#file name solo and pair formats are the same for file name formats sequence_(1|2).fq
			#and sequence2.bam
			else
				{
					$files->{$fname}->{local_solo_format} = "%s_%d";
					$files->{$fname}->{local_pair_format} = "%s";
					$files->{$fname}->{solo_format} = "\@s_\@d";
					$files->{$fname}->{pair_format} = "\@s";

					#file name has format %s_%d where %s is sequence name and %d is pair number, e.g. 
					#sequence_1.fq and sequence_2.fq or sequence_1.fq.gz
					if($fname =~ m/^(.+)_(1|2)(\.\w+)+$/)
						{
							$files->{$fname}->{fname_format} = qr/^(.+)_(1|2)(\.\w+)+$/;
						}
					#files are in BAM format, with paired reads in the same file, e.g. flowcell2.bam
					#N.B. THIS WILL NOT MATCH A FILE WITH TWO SUFFIXES
					#the previous two file formats will....can be changed...
					elsif($fname =~ m/^(.+)\.\w+$/)
						{
							$files->{$fname}->{fname_format} = qr/^(.+)\.\w+$/;
						}
					#files are in casava 1.8+ directory format
					elsif($fname =~ m/^(.+)$/)
						{
							$files->{$fname}->{fname_format} = qr/^(.+)$/;
						}
				}
		}
	foreach my $fname (keys %{$files})
		{
			my $format = $files->{$fname}->{fname_format};
			$files->{$fname}->{file_id} = $1 if $fname =~ m/$format/;
		}

	#get suffix of fastq files (tested, picks up .bam and .fq.gz suffixes, should be able to pick up any other suffix as well
	foreach my $fname (keys %{$files})
		{
			my $file = $fname;
			if($files->{$fname}->{isDir})
				{
					#CLEAN UP
					my $targetDir = sprintf("%s%s",$params->{dir},$fname);
					$targetDir = readlink($targetDir) if -l $targetDir;
					
					opendir DIR,$targetDir or die "Failed to open $targetDir: $!";
					my $filesInDir = [grep {!/^\.{1,2}$/} readdir DIR];
					closedir DIR;
					
					$files->{$fname}->{sample_file} = sprintf("%s%s/%s",$params->{dir},$fname,$filesInDir->[0]);
					$file = $files->{$fname}->{sample_file};
					$files->{$fname}->{count} = scalar grep {/(fastq|fq|txt|bam)/} readdir DIR;
					$files->{$fname}->{count}++;				
					closedir DIR;
				}
			$files->{$fname}->{suffix} = $1 if $file =~ m/.*(\.\w+)$/;
			
			if($files->{$fname}->{suffix} eq ".gz")
				{
					$files->{$fname}->{zipped} = 1;
					$files->{$fname}->{suffix} = $1 if $file =~ m/.*(\.\w+)\.gz$/;
				}
		}


	#check fastq files for correct format
	foreach my $fname (keys %{$files})
		{
			if($files->{$fname}->{suffix} =~ m/(fastq|txt|fq)$/)
				{
					my $name = sprintf("%s%s",$params->{dir},$fname);
					$name = $files->{$fname}->{sample_file} if $files->{$fname}->{isDir};
					
					if(!checkFastq($name, "quick", $files->{$fname}->{zipped}))
						{
							warn "$name is not in FASTQ format, skipping\n" if !$files->{$fname}->{isDir};
							warn "$name from directory $fname does not appear to be in FASTQ format, skipping directory\n" if $files->{$fname}->{isDir};
							delete $files->{$fname};
						}
				}
		}


	#if quality score isn't defined, guess it here from one of the files
	if(!defined $tetrascan_params->{qual_format})
		{
			print "No quality score format specified, guessing quality format for each sequence library in ".$params->{dir}."\n";
			foreach my $fname (keys %{$files})
				{
					my $file = sprintf("%s%s",$params->{dir},$fname);
					$file = $files->{$fname}->{sample_file} if $files->{$fname}->{isDir};
					$files->{$fname}->{qual_format} = guessQualScore($file,$files->{$fname}->{zipped});
				}
		}
		
	my $paired_files;
	PAIR:foreach my $fname (keys %{$files})
		{
			if($files->{$fname}->{isDir})
				{
					DIR:foreach(keys %{$files->{$fname}})
						{
							next DIR if $_ =~ m/(fname_format|sample_file)/;
							$paired_files->{$files->{$fname}->{file_id}}->{$_} = $files->{$fname}->{$_};
						}
					$paired_files->{$files->{$fname}->{file_id}}->{paired} = 1;
					if($files->{$fname}->{suffix} ne "bam")
						{
							$paired_files->{$files->{$fname}->{file_id}}->{paired} = 0 if ($paired_files->{$fname}->{count} % 2) == 1;
						}
				}
			else
				{
					if(!defined $paired_files->{$files->{$fname}->{file_id}})
						{
							FILE:foreach(keys %{$files->{$fname}})
								{
									next FILE if $_ =~ m/(fname_format)/;
									$paired_files->{$files->{$fname}->{file_id}}->{$_} = $files->{$fname}->{$_};
								}
							$paired_files->{$files->{$fname}->{file_id}}->{paired} = 1 if $files->{$fname}->{suffix} =~ m/bam/;

						}
					else
						{
							$paired_files->{$files->{$fname}->{file_id}}->{paired} = 1;
						}
				}
		}
				
	return ($paired_files);
}

sub checkGenOpts
{
	my $params = shift;
	
	#set default parameters
	$params->{local} = 1 if !defined $params->{distributed};
	$params->{split} = 1 if !defined $params->{split};
	$params->{zip_output} = 1 if !defined $params->{zip_output};
	$params->{clean_bam} = 1 if !defined $params->{clean_bam};

	#check that either a directory of sequencing file or 1 or 2 sequencing files are specified
	if(!defined $params->{dir} and !defined $params->{pe1} and !defined $params->{pe2})
		{
			pod2usage(-message => "You must specify a directory of sequence files (BAM or fastq format) using the --dir flag or a FASTQ file with --pe1",
									-verbose => 2, -output => \*STDERR)
		}

	#check if both paired ends are present if paired flag is defined and using pe1/pe2 tags
	if(defined $params->{pe1} and !defined $params->{pe2})
		{
			$params->{isPaired} == 1 ? die "Error: isPaired true; pe2 not defined\n" : print "Processing ", $params->{pe1}, " as a single-end read library\n";
		}
	
	#UNRESOLVED: what to do if no project_id is specified and using pe1/pe2?
	#infer project_id if not specified
	if(!defined $params->{project_id} and defined $params->{dir})
		{
			$params->{project_id} = $2 if $params->{dir} =~ m/(.*\/)?([^\/]+)(\/)*?$/;
			print "Inferring project id as ",$params->{project_id},"\n";
		}
	$params->{base_project_id} = [split /\//, $params->{project_id}];
	$params->{base_project_id} = $params->{base_project_id}->[$#{$params->{base_project_id}}];
	
	die "Error!  You must specify a genome for the alignment step\n" if !defined $params->{genome};
	die "Error: genome ",$params->{genome}," does not exist\n" if ! -s $params->{genome};

	#get genome root for figuring out for inferring offset index and usat db files, genome input format
	#must look like /path/to/genome/genome.fa
	if($params->{genome} !~ m/(.*\/)?.+\..*/)
		{
			die "Incorrect genome specification.  Please specify genome as /path/to/genome/genome.mod.fa\n";
		}
	my $genome_root = [split /\./, $params->{genome}]; #root is everything before first dot in name
	$genome_root = $genome_root->[0];
	my $genome_dir = $1 if $params->{genome} =~ m/(.*)\/.*/;
	$genome_dir = "." if !defined $genome_dir;

	my %file_suffix_regex = (offset_index => qr/\.osi$/, ms_db => qr/\.ms.*\.fa/);
	my %file_suffixes = (offset_index => ".osi", ms_db => ".ms.fa");
	my %file_tags = (offset_index => "--offset", ms_db => "--ms_db");
	my %file_names = (offset_index => 'an offset index', ms_db => 'a microsatellite database');
	#check if offset_index and ms_db exist from inferred file paths or specified --ms_db/--offset file paths
	foreach(qw(offset_index ms_db))
		{
			if(defined $params->{$_})
				{
					die "Error: $_ ",$params->{$_}," is invalid. Create ",$file_names{$_}," using  MicroSeq --complete or tetrascan\n" if $params->{$_} !~ m/$file_suffix_regex{$_}/;
					die "Error: $_ ",$params->{$_}," does not exist\n" if ! -s $params->{$_};
				}
			else
				{
					$params->{$_} = sprintf("%s%s", $genome_root,$file_suffixes{$_});
					if(! -s $params->{$_})
						{
							my $errstring = sprintf("Error: inferred %s %s does not exist.\nPlease specify %s with %s, move %s that exists to %s, or run MicroSeq --complete/tetrascan to generate %s.\n",
																	 $_,$params->{$_},$file_names{$_},$file_tags{$_},$file_names{$_},$genome_dir,$file_names{$_});
							die $errstring;
						}
				}
		}
		
	my %dir_desc = (dir => "sequencing library", target_dir => "target");
	foreach(qw(dir target_dir))
		{
			if(defined $params->{$_})
				{
					checkDirExist($params->{$_},$dir_desc{$_});
					$params->{$_} = checkTrailingSlash($params->{$_});
				}
		}
	$params->{target_dir} = $params->{dir} if !defined $params->{target_dir};
}
