#!/usr/bin/perl

use lib "/mnt/wigclust4/home/bekritsk/tools/microsatellite/perllib/"; #need to change this to config file variable

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path qw(mkpath);
use POSIX qw(:sys_wait_h);

use PipelineLogger qw(:all);
use MicroSeqFunc qw(forkWait createLocalMicroSeqDir getProgramPaths); #change to config file var
use Timer qw(:all);

sub findMicrosatellites;
sub alignReads;
sub getSamFile;
sub reindexReads;
sub splitReadsByChr;
sub changeFormatStrings;
sub convertToBam;
sub zip;

sub printParamsToLogFile;
sub filesToString;
sub makeTetrascanCmd;
sub makeAlnCmd;
sub makeSamCmd;
sub makeReindexCmd;
sub makeSplitCmd;
sub makeBamConvertCmd;
sub generateFileNames;

my $tet_params;
my $gen_params;
my $files;

my $created;
my $action;

GetOptions
(
	"pair_format=s"				=> \$files->{pair_f},
	"solo_format=s"				=> \$files->{solo_f},
	"fastq=s"							=> \$files->{fastq_base},
	"dir=s"								=> \$files->{dir},
	"target_dir:s"				=> \$files->{target_dir},
	"genome=s"						=> \$files->{genome},
	"offset_index=s"			=> \$files->{offset_index},
	"paired"							=> \$files->{paired},
	"suffix=s"						=> \$files->{suffix},
	"zipped"							=> \$files->{zipped},
	"count:i"							=> \$files->{count},
	"isDir"								=> \$files->{isDir},
	
	"level=i"							=> \$gen_params->{level},
	"threads=i"						=> \$gen_params->{threads},
	"split"								=> \$gen_params->{split},
	"clean"								=> \$gen_params->{clean},
	"zip_output"					=> \$gen_params->{zip_output},
	"clean_bam"						=> \$gen_params->{clean_bam},
	"keep_all"						=> \$gen_params->{keep_all},
	"decentralized"				=> \$gen_params->{decentralized},
	"project_id:s"				=> \$gen_params->{project_id},
	
	"print_all_MS=i"			=> \$tet_params->{Mi},
	"max_n=i"							=> \$tet_params->{N},
	"min_unit=i"					=> \$tet_params->{u},
	"min_length=i"				=> \$tet_params->{l},
	"min_usize=i"					=> \$tet_params->{m},
	"max_usize=i"					=> \$tet_params->{x},
	"min_qual=i"					=> \$tet_params->{Q},
	"qual_format=s"				=> \$tet_params->{qual_format},
	"trim_type=s"					=> \$tet_params->{trim_type},
);

my $hostname = `hostname`;
my $nodename = [split /\./,$hostname];
$nodename = $nodename->[0];
my $prog_paths = getProgramPaths();

print "Running msScanAndAlignDir; level ", $gen_params->{level}," on ",$hostname,"\n";

changeFormatStrings($files);
my ($log_file, $start_time) = xInit($files->{target_dir},"msScanAndAlignDir.pl",$files->{fastq_base});
print "Log file for this run can be found at ", $log_file->{filename},"\n";

$files->{local_dir} = $files->{target_dir};
if($gen_params->{decentralized})
	{
		($files->{local_dir},$created) = createLocalMicroSeqDir($gen_params->{project_id},$nodename);
		if($created)
			{
				$action = sprintf("Created local MicroSeq directory %s",$files->{local_dir});
			}
		else
			{
				$action = sprintf("%s already exists",$files->{local_dir});
			}
		xLog($action,$log_file);
	}
	
prepareOutDirs($files,$log_file);
generateFileNames($files,$gen_params,$log_file);
printParamsToLogFile($log_file,$files, $tet_params, $gen_params);

printf("msScanAndAlignDir processing %s\n", $files->{fastq_base});
my $mini_sched = {threads => $gen_params->{threads}, jobs => scalar keys %{$files->{fastq}}};

foreach my $index (keys %{$files->{fastq}})
	{
		submitToScheduler(getDFileCmd($gen_params,$tet_params,$files,$index,$prog_paths),$mini_sched);
	}
finishScheduler($mini_sched);

my $fhs;
for(my $i = 0; $i < 2; $i++)
	{
		$files->{master_mod_fastq}->[$i] = sprintf("%s%s_%d.mod.ms.txt",$files->{local_dir},$files->{fastq_base},$i + 1);
		open my $fh, '>', $files->{master_mod_fastq}->[$i] or die "Could not open ",$files->{master_mod_fastq}->[$i],": $!\n";
		push @{$fhs}, $fh;
	}

foreach my $index (keys %{$files->{mod_fastq}})
	{
		for(my $i = 0; $i < scalar @{$files->{mod_fastq}->{$index}}; $i++)
			{
				xCat($files->{mod_fastq}->{$index}->[$i],$files->{master_mod_fastq}->[$i],$fhs->[$i]);
				xUnlink($files->{mod_fastq}->{$index}->[$i]);
			}
	}
close $_ foreach(@{$fhs});
zipModFastq($files,$log_file);

$files->{master_split_dir} = sprintf("%s/.%ssplit/",$files->{local_dir},$files->{fastq_base});
$files->{target_split_dir} = sprintf("%s/.%ssplit",$files->{target_dir},$files->{fastq_base});

xUnlink($files->{target_split_dir},"old link") if -l $files->{target_split_dir};
xRmdir($files->{target_split_dir},1,"old directory") if -d $files->{target_split_dir};
xUnlink($files->{target_split_dir},"old file") if -f $files->{target_split_dir};

xUnlink($files->{master_split_dir},"old link") if -l $files->{master_split_dir};
xRmdir($files->{master_split_dir},1,"old directory") if -d $files->{master_split_dir};
xUnlink($files->{master_split_dir},"old file") if -f $files->{master_split_dir};
xMkdir($files->{master_split_dir});

my $ofh;
foreach my $index (keys %{$files->{msi}})
	{
		my $coord_time = time;
		my $split_dir = sprintf("%s/%s/.%s_%0*dsplit/",$files->{local_dir},$files->{fastq_base},$files->{prefix},$files->{index_length},$index);
		print "Adding $split_dir to ",$files->{master_split_dir},"\n";
		opendir(IDIR,$split_dir) or die "Could not open $split_dir: $!\n";
		foreach my $file (grep {!/^\.{1,2}$/} readdir IDIR)
			{
				my $fname = sprintf("%s%s",$split_dir,$file);
				open my $ifh, '<', $fname or die "Could not open $fname: $!\n";
				my $coord = @{[split /\./, $file]}[0];
				if(!defined $ofh->{$coord})
					{
						$ofh->{$coord}->{name} = sprintf("%s%s.%s_%0*d.msi.sam",$files->{master_split_dir},$coord,$files->{prefix},$files->{index_length},$index);
						$ofh->{$coord}->{handle} = xCreateFile($ofh->{$coord}->{name});
					}
				
				$action = sprintf("Adding %s to %s",$fname,$ofh->{$coord}->{name});
				xLog($action,$log_file);
				HDR:while(my $line = <$ifh>)
					{
						if($line =~ m/^@/ and !defined $ofh->{$coord}->{hdr})
							{
								print {$ofh->{$coord}->{handle}} $line;
							}
						elsif($line !~ m/^@/)
							{
								print {$ofh->{$coord}->{handle}} $line;
								$ofh->{$coord}->{hdr} = 1;
								last HDR;
							}
					}
				print {$ofh->{$coord}->{handle}} $_ while <$ifh>;
				close $ifh;
				print {$log_file->{ptr}} "\t",getElapsedTime($coord_time),"\n";
				xUnlink($fname);
			}
		xRmdir($split_dir);
	}

foreach(keys %{$ofh})
	{
		close $ofh->{$_}->{handle};
	}

my $basename;
my $link_name;

#link files in local work directory back to project target directory
#if not writing files back to target directory
if($gen_params->{decentralized})
	{
		foreach(@{$files->{master_mod_fastq}})
			{
				$basename = [split /\//, $_];
				$basename = $basename->[$#{$basename}];
				$link_name = sprintf("%s%s",$files->{target_dir},$basename);
				my $source_name = $_;
				if($gen_params->{zip_output})
					{
						$source_name .= ".gz";
						$link_name .= ".gz";
					}
				xUnlink($link_name,"old file") if -f $link_name;
				xUnlink($link_name,"old link") if -l $link_name;
				xSymlink($source_name,$link_name) if $gen_params->{decentralized};
			}
		$link_name = sprintf("%s.%ssplit",$files->{target_dir},$files->{fastq_base});

		xSymlink($files->{master_split_dir},$link_name);
	}
#xRmdir(); delete parent directory

#################################SUBROUTINES###############################################

#scheduler functions derived from celera assembler source code
sub submitToScheduler
{
	my $command = shift;
	my $scheduler_info = shift;
	
	push @{$scheduler_info->{queue}}, $command;
}

sub schedulerForkProcess
{
	my $command = shift;
	my $pid;
	
	FORK:
		{
			if($pid = fork)
				{
					return $pid;
				}
			elsif(defined $pid)
				{
					xRun($command);
					exit 0;
				}
			elsif($! =~ m/No more processes/)
				{
					sleep 1;
					redo FORK;
				}
			else
				{
					die "Can't fork: $!\n";
				}
		}
}

sub reapProcess
{
	my $pid = shift;
	
	if(waitpid($pid,WNOHANG) > 0)
		{
			return 1;
		}
	else
		{
			return 0;
		}
}

sub runScheduler
{
	my $scheduler_info = shift;
	$scheduler_info->{new} = [];
	
	foreach my $curr_pid (@{$scheduler_info->{running}})
		{
			#if current pid is still running, add it to list of new processes
			push @{$scheduler_info->{new}}, $curr_pid if !reapProcess($curr_pid);
		}
	$scheduler_info->{running} = [];
	$scheduler_info->{running} = $scheduler_info->{new};
	
	my $num_running;
	defined $scheduler_info->{running} ? $num_running = scalar @{$scheduler_info->{running}} : $num_running = 0;
		
	my $num_added = 0;
	while($num_running < $scheduler_info->{threads} 
		and scalar @{$scheduler_info->{queue}} > 0)
		{
			my $command = shift @{$scheduler_info->{queue}};
			push @{$scheduler_info->{running}},schedulerForkProcess($command);
			$num_running = scalar @{$scheduler_info->{running}};
		}
}

sub finishScheduler
{
	my $scheduler_info = shift;
	my $child;
	my $remaining_jobs;
	
	$remaining_jobs = scalar @{$scheduler_info->{queue}};
	while($remaining_jobs > 0)
		{
			runScheduler($scheduler_info);
			$remaining_jobs = scalar @{$scheduler_info->{queue}};
			
			if($remaining_jobs > 0)
				{		
					$child = waitpid(-1,0); #waits for a child process to die
					$scheduler_info->{new} = [];
					foreach my $curr_pid (@{$scheduler_info->{running}})
						{
							push @{$scheduler_info->{new}}, $curr_pid if $curr_pid != $child;
						}
					$scheduler_info->{running} = [];
					$scheduler_info->{running} = $scheduler_info->{new};
				}
		}

	my $wait_job;
	while(scalar @{$scheduler_info->{running}} > 0)
		{
			waitpid(shift @{$scheduler_info->{running}}, 0);
		}		
}

sub changeFormatStrings
{
	my $params = shift;
	
	$params->{pair_f} =~ s/\@/%/g;
	$params->{solo_f} =~ s/\@/%/g;
}

sub prepareOutDirs
{
	my $params = shift;
	my $log = shift;
	$params->{out_dir} = sprintf("%sSTDOUT/",$params->{target_dir});
	$params->{err_dir} = sprintf("%sSTDERR/",$params->{target_dir});

	foreach(qw(out_dir err_dir))
		{
			xMkdir($params->{$_}) if !-d $params->{$_};
		}
}

sub getDFileCmd
{
	my ($gen_params,$tet_params,$files,$index,$paths) = @_;

	my $file_prefix = sprintf("%s_%0*d",$files->{prefix},$files->{index_length},$index);
	my $err_file = sprintf("%s%s.sanda.err",$files->{err_dir},$file_prefix);
	my $out_file = sprintf("%s%s.sanda.out",$files->{out_dir},$file_prefix);

	my $command = undef;
			
	$command = sprintf("%smsScanAndAlignDFile.pl --index %d", $paths->{MicroSeq},$index);
						
	foreach(keys %{$tet_params})
		{
			$command .= sprintf(" --%s %s", $_, $tet_params->{$_});
		}
		
	#general params
	foreach(sort qw(project_id level))
		{
					$command .= sprintf(" --%s %s", $_, $gen_params->{$_}) if defined $gen_params->{$_};
		}
	foreach(sort qw(split clean zip_output clean_bam keep_all decentralized))
		{
			$command .= sprintf(" --%s",$_) if defined $gen_params->{$_};
		}
		
	#file params
	foreach(sort qw(count solo_format pair_format suffix fastq_base dir target_dir genome offset_index prefix index_length))
		{
					$command .= sprintf(" --%s %s", $_, $files->{$_}) if defined $files->{$_};
		}
	foreach(sort qw(paired zipped isDir))
		{
			$command .= sprintf(" --%s",$_) if defined $files->{$_};
		}

	$command .= sprintf(" 1>%s 2>%s",$out_file,$err_file);
	
	return $command;
}

sub zipModFastq
{
	my $files = shift;
	my $log = shift;
	
	my $action;
	my $children;
	
	for(my $i = 0; $i < scalar @{$files->{master_mod_fastq}}; $i++)
		{
			my $pid = fork();
			if($pid)
				{
					push @{$children}, $pid;
				}
			else
				{
					$action = sprintf("Zipping %s\n",$files->{master_mod_fastq}->[$i]);
					xLog($action,$log);
					my $command = sprintf("gzip -f %s",$files->{master_mod_fastq}->[$i]);
					my $zipOut = xRun($command);
					printf("Zipped %s\n",$files->{master_mod_fastq}->[$i]); #this will overwrite an existing *mod.ms.* file
					exit 0;
				}
		}
	forkWait($children);
}

sub printParamsToLogFile
{
	my $log_file = shift;
	my $files = shift;
	my $tetrascan_params = shift;
	my $gen_params = shift;
	
	my $opt_explanations = {
		"M" => "#print only reads with microsatellites",
		"N" => "#max number of Ns in read",
		"u" => "#min repeats of microsatellite motif",
		"l" => "#min microsatellite tract length",
		"m" => "#min microsatellite motif length",
		"x" => "#max microsatellite motif length",
		"Q" => "#min basecall qual (independent of Phred score)",
	};

	my $qual_format_explanations = {
		"a" => "#sanger format (Phred + 33)",
		"s" => "#illumina format (Phred + 64)",
	};

	my $trim_explanations = {
		"B" => "#BWA-style trimming",
		"g" => "#trim to last base pair greater than Q",
		"b" => "#trim to first base pair less than Q",
		"n" => "#do not trim",
	};

	print {$log_file->{ptr}} "File parameters:\n";
	foreach(sort keys %{$files})
		{
			printf {$log_file->{ptr}} ("\t%s: %s",$_, $files->{$_}) if defined $files->{$_} and ref($files->{$_}) ne 'ARRAY';
			printf {$log_file->{ptr}} ("\t%s: %s",$_, join "; ", @{$files->{$_}}) if ref($files->{$_}) eq 'ARRAY';
			print {$log_file->{ptr}} "\n";
		}

	print {$log_file->{ptr}} "Tetrascan parameters:\n";
	foreach(sort keys %{$tetrascan_params})
		{
			next if !defined $tetrascan_params;
			print STDERR $_,"\n";
			printf {$log_file->{ptr}} ("\t%s: %s",$_, $tetrascan_params->{$_});
			printf {$log_file->{ptr}} (" %s", $qual_format_explanations->{$tetrascan_params->{$_}}) if $_ eq "qual_format";
			printf {$log_file->{ptr}} (" %s", $trim_explanations->{$tetrascan_params->{$_}}) if $_ eq "trim_type";
			printf {$log_file->{ptr}} (" %s", $opt_explanations->{$_}) if defined $opt_explanations->{$_};
			print {$log_file->{ptr}} "\n";
		}

	print {$log_file->{ptr}} "General parameters:\n";
	foreach(sort keys %{$gen_params})
		{
			printf {$log_file->{ptr}} ("\t%s: %s\n",$_, $gen_params->{$_}) if defined $gen_params->{$_};
		}
	
	xBreak($log_file);
}

sub generateFileNames
{
	my $files = shift;
	my $gen_params = shift;
	my $log = shift;
	my $file_prefix;
	my $file_stub;

	my $dir = sprintf("%s%s/",$files->{dir},$files->{fastq_base});
	my $local_dir = sprintf("%s%s/",$files->{local_dir},$files->{fastq_base});
	$files->{local_subdir} = sprintf("%s%s/",$files->{local_dir},$files->{fastq_base});
	opendir(SEQ,$dir) or die "Could not open $dir: $!\n";
	my $filenames = [readdir SEQ];
	closedir SEQ;
	
	my ($prefixes,$read_indices,$indices);
	my ($prefix,$read_index,$index);
	my $file_info;

	if($files->{suffix} ne '.bam' and $files->{paired})
		{
			FILE:foreach(@{$filenames})
				{
					next FILE if $_ =~ m/^\.{1,2}$/ or $_ !~ m/$files->{suffix}/;
					($prefix,$read_index,$index) = ($1,$2,$3) if $_ =~ m/(.+)_(R.+)_(\d+)\..+(\..+)*$/; #assumes that all output has an R1/R2 flag
					$read_indices->{$read_index}++;
					$prefixes->{$prefix}++;
					$indices->{$index}++;
				}
		if(scalar keys %{$read_indices} != 2)
			{
				warn "Files have incorrect mate pair identifiers in ",$files->{dir},": ";
				print STDERR map {"$_,"} keys %{$read_indices} if defined $read_indices;
				print STDERR "\n";
				exit 1;
			}
			foreach my $read1 (keys %{$read_indices})
				{
					foreach my $read2 (keys %{$read_indices})
						{
							if($read_indices->{$read1} != $read_indices->{$read2})
								{
									printf STDERR ("Unbalanced number of mate pairs: %s => %d; %s => %d",$read1,$read_indices->{$read1},$read2,$read_indices->{$read2});
									exit 1;
								}
						}
				}
		}
	else
		{
			FILE:foreach(@{$filenames})
				{
					next FILE if $_ =~ m/^\.{1,2}$/ or $_ !~ m/$files->{suffix}/;
					($prefix,$index) = ($1,$2) if $_ =~ m/(.+)_(\d+)\..+$/; #assumes that all files have BAM suffix
					$prefixes->{$prefix}++;
					$indices->{$index}++;
				}
		}	
	if(scalar keys %{$prefixes} > 1)
		{
			warn "Too many prefixes in ",$files->{dir},": ";
			print STDERR map {"$_,"} keys %{$prefixes};
			print STDERR "\n";
			exit 1;
		}
	my $index_length;
	$index_length->{length($_)}++ foreach(keys %{$indices});
	if(scalar keys %{$index_length} > 1)
		{
			warn "Inconsistent index padding length in ",$files->{dir},": ";
			print STDERR map {"length $_ =>".$index_length->{$_}.", "} keys %{$index_length};
			print STDERR "\n";
			exit 1;
		}
	$files->{index_length} = @{[keys %{$index_length}]}[0];
	$files->{prefix} = @{[keys %{$prefixes}]}[0];
	$index_length = @{[keys %{$index_length}]}[0];
	$prefix = @{[keys %{$prefixes}]}[0];

	my $fastq_file;
	#Create names for output files
	if(!$files->{paired})
		{
			foreach $index (keys %{$indices})
				{
					$fastq_file = sprintf("%s%s_%0*d%s",					$dir,$prefix,$index_length,$index,$files->{suffix});
					$fastq_file .= ".gz" if $files->{zipped};
					push @{$files->{fastq}->{$index}},$fastq_file;

					push @{$files->{mod_fastq}->{$index}}, sprintf("%s%s_%0*d.mod.ms.txt", $local_dir,$prefix,$index_length,$index);
					push @{$files->{ms_fastq}->{$index}},  sprintf("%s%s_%0*d.ms.txt", $local_dir,$prefix,$index_length,$index);
					push @{$files->{hdr}->{$index}},			 sprintf("%s%s_%0*d.hdr.txt", $local_dir,$prefix,$index_length,$index);
					push @{$files->{sai}->{$index}}, 			 sprintf("%s%s_%0*d.mod.sai", $local_dir,$prefix,$index_length,$index);
					$files->{sam}->{$index} = sprintf("%s%s_%0*d.mod.sam",$local_dir,$prefix,$index_length,$index);
					$files->{msi}->{$index} = sprintf("%s%s_%0*d.msi.sam",$local_dir,$prefix,$index_length,$index);
				}
		}
	else
		{
			foreach $index (keys %{$indices})
				{
					foreach $read_index (qw(1 2))
						{
							if($files->{suffix} !~ m/\.bam$/)
								{
									$fastq_file = sprintf("%s%s_R%d_%0*d%s",					$dir,$prefix,$read_index,$index_length,$index,$files->{suffix});
									$fastq_file .= ".gz" if $files->{zipped};
									push @{$files->{fastq}->{$index}},$fastq_file;
								}
							else
								{
									$files->{isBam} = 1;
									$file_stub = sprintf($files->{pair_f}, $files->{fastq_base});
									$files->{fastq}->{$index}->[0] = sprintf("%s%s_%0*d%s",$dir,$prefix,$index_length,$index,$files->{suffix});
								}
							push @{$files->{mod_fastq}->{$index}}, sprintf("%s%s_R%d_%0*d.mod.ms.txt", $local_dir,$prefix,$read_index,$index_length,$index);
							push @{$files->{ms_fastq}->{$index}},  sprintf("%s%s_R%d_%0*d.ms.txt", $local_dir,$prefix,$read_index,$index_length,$index);
							push @{$files->{hdr}->{$index}},			 sprintf("%s%s_R%d_%0*d.hdr.txt", $local_dir,$prefix,$read_index,$index_length,$index);
							push @{$files->{sai}->{$index}}, 			 sprintf("%s%s_R%d_%0*d.mod.sai", $local_dir,$prefix,$read_index,$index_length,$index);
						}
					$files->{sam}->{$index} = sprintf("%s%s_%0*d.mod.sam",$local_dir,$prefix,$index_length,$index);
					$files->{msi}->{$index} = sprintf("%s%s_%0*d.msi.sam",$local_dir,$prefix,$index_length,$index);
				}
		}
}

sub filesToString
{
	my $files = shift;

	return $files->[0] if scalar @{$files} == 1;
	return ((join ",", @{$files}[0..(scalar @{$files} - 2)])." and ".($files->[(scalar @{$files} - 1)]));
}
