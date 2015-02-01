#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Path qw(mkpath);
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

	print $config->{uSeqPerllib}, "\n";
	die "Could not parse config file: $@" if $@;
	die "Could not read config file: $!" 	unless defined $config;
	die "Could not run config file"				unless $config;
	unshift @INC, $config->{uSeqPerllib};
}

use MicroSeqFunc qw(createLocalMicroSeqDir);
use Certificate;
use Logger;

sub findMicrosatellites;
sub alignReads;
sub getSamFile;
sub reindexReads;
sub splitReadsByChr;
sub changeFormatStrings;
sub convertToBam;
sub zip;

sub filesToString;
sub makeTetrascanCmd;
sub makeAlnCmd;
sub makeSamCmd;
sub makeReindexCmd;
sub makeSplitCmd;
sub makeBamConvertCmd;
sub generateFileNames;
sub deleteFiles;

my $tetrascan_params;
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
	"paired"							=> \$files->{isPaired},
	"suffix=s"						=> \$files->{suffix},
	"zipped"							=> \$files->{zipped},
	
	"level=i"							=> \$gen_params->{level},
	"threads=i"						=> \$gen_params->{threads},
	"split"								=> \$gen_params->{split},
	"clean"								=> \$gen_params->{clean},
	"zip_output"					=> \$gen_params->{zip_output},
	"clean_bam"						=> \$gen_params->{clean_bam},
	"keep_all"						=> \$gen_params->{keep_all},
	"decentralized"				=> \$gen_params->{decentralized},
	"project_id:s"				=> \$gen_params->{project_id},
	
	"print_all_MS=i"			=> \$tetrascan_params->{M},
	"max_n=i"							=> \$tetrascan_params->{N},
	"min_unit=i"					=> \$tetrascan_params->{u},
	"min_length=i"				=> \$tetrascan_params->{l},
	"min_usize=i"					=> \$tetrascan_params->{m},
	"max_usize=i"					=> \$tetrascan_params->{x},
	"min_qual=i"					=> \$tetrascan_params->{Q},
	"barcodeLength=i"			=> \$tetrascan_params->{t},
	"qual_format=s"				=> \$tetrascan_params->{qual_format},
	"trim_type=s"					=> \$tetrascan_params->{trim_type},
);

my $hostname = `hostname`;
my $nodename = [split /\./,$hostname];
$nodename = $nodename->[0];

print "Running msScanAndAlign; level ", $gen_params->{level}," on ",$hostname,"\n";

changeFormatStrings($files);
my $log = Logger->new($files->{target_dir},"msScanAndAlign.pl",sprintf($files->{pair_f}, $files->{fastq_base}));
my $certifier = Certificate->new($files->{target_dir})

$files->{local_dir} = $files->{target_dir};
if($gen_params->{decentralized})
	{
		($files->{local_dir},$created) = createLocalMicroSeqDir($gen_params->{project_id},$nodename,$log);
		print "Local directory: ",$files->{local_dir},"\n";
	}
	
#get filenames
generateFileNames($files,$gen_params);
print "Log file for this run can be found at ", $log->fname(),"\n";
$log->recordParams($files,"File");
$log->recordParams($gen_params,"General");
$log->recordParams($tetrascan_params,"Tetrascan");

my $lane_name = sprintf($files->{pair_f}, $files->{fastq_base});
printf("msScanAndAlign Processing %s\n", $lane_name);

if($gen_params->{level}==0)
	{
		findMicrosatellites($MicroSeq_path,$tetrascan_params, $files,$log);
		print STDERR "Found microsatellites\n";
		$log->rm($files->{ms_fastq});
	}
if($gen_params->{level} <= 1)
	{
		alignReads($bwa_path,$files, $gen_params);
		print STDERR "Aligned reads\n";
		getSamFile($bwa_path,$files);
		print STDERR "Got SAM file\n";
		deleteFiles($files->{sai},$log_file);
		if($gen_params->{zip_output})
			{
				zipModFastq($files,$log_file);
			}
		elsif($gen_params->{clean})
			{
				deleteFiles($files->{mod_fastq},$log_file);
			}
	}
if($gen_params->{level} <= 2)
	{
		reindexReads($MicroSeq_path, $files);
		print STDERR "Reindexed reads\n";
		if(!$gen_params->{keep_all})
			{
				deleteFiles($files->{hdr},$log_file);
				deleteFiles($files->{sam},$log_file);
			}
	}

if($gen_params->{level} <= 3)
	{
	#split files by chromosome if instructed to (generally, this is necessary
	#to build profiles if one individual has multiple lanes of sequence data to
	#make sure that all the reads that map to the same locus are in the same file
		if($gen_params->{split})
			{
				splitReadsByChr($MicroSeq_path, $files, $log_file);
				print STDERR "Split reads\n";
				if($gen_params->{zip_output} and ! $gen_params->{clean_bam})
					{
						convertToBam($samtools_path,$files->{msi},$log_file);
					}
				if(!$gen_params->{keep_all})
					{
						deleteFiles($files->{msi_link},$log_file);
						deleteFiles($files->{msi},$log_file);
					}
			}
	}
	
my $basename;
my $link_name;

#link files in local work directory back to project target directory
#if not writing files back to target directory
if($gen_params->{decentralized})
	{
		foreach(@{$files->{mod_fastq}})
			{
				$basename = [split /\//, $_];
				$basename = $basename->[$#{$basename}];
				$link_name = sprintf("%s/%s",$files->{target_dir},$basename);
				my $origin_name = $_;
				if($gen_params->{zip_output})
					{
						$origin_name .= ".gz";
						$link_name .= ".gz";
					}
				if(-f $link_name or -l $link_name)
					{
						print "$link_name already exists!  Removing...\n";
						$action = sprintf("%s already exists, removing...\n",$link_name);
						logAction($action,$log_file);
						unlink $link_name;			
					}
				if($gen_params->{decentralized})
					{
						print "Linking $origin_name to $link_name\n";
						$action = sprintf("Linking %s to %s",$origin_name,$link_name);
						logAction($action,$log_file);
						symlink($origin_name,$link_name);
					}
			}

		$basename = [split /\//, $files->{split_dir}];
		$basename = $basename->[$#{$basename}];
		$link_name = sprintf("%s/%s",$files->{target_dir},$basename);
		print "Loooking to see if $link_name already exists\n";
		if(-d $link_name and ! -l $link_name)
			{
				print "$link_name is already a directory, removing\n";
				opendir LINK, $link_name;
				FILE:foreach(readdir LINK)
					{
						next FILE if $_ =~ m/^\.{1,2}$/;
						unlink sprintf("%s/%s",$link_name,$_) or die "Could not delete $link_name/$_: $!\n";
					}
				closedir LINK;
				$action = sprintf("Removing old copy of directory %s",$link_name);
				logAction($action,$log_file);
				rmdir $link_name;
			}
		elsif(-l $link_name)
			{
				my $linked_dir = readlink($link_name);
				if(-d $linked_dir)
					{
						opendir LINK, $linked_dir or warn "Could not open linked directory: $!\n";
						DEL:foreach(readdir LINK)
							{
								next DEL if $_ =~ m/^\.{1,2}$/;
								unlink sprintf("%s/%s",$linked_dir,$_) or die "Could not delete $_ in $linked_dir: $!\n";
							}
						closedir LINK;
						$action = sprintf("Deleted old linked directory %s",$linked_dir);
						logAction($action,$log_file);
					}
				elsif(-f $linked_dir)
					{
						$action = sprintf("Deleting linked files %s",$linked_dir);
						logAction($action,$log_file);
						unlink $linked_dir;
					}
				print "$link_name already exists!  Removing...\n";
				$action = sprintf("Deleting old link %s",$link_name);
				logAction($action,$log_file);
				unlink $link_name;
			}
		elsif(-f $link_name)
			{
				print "$link_name is already a file! Removing...\n";
				$action = sprintf("%s is already a file, removing...\n",$link_name);
				deleteFiles($link_name,$log_file);
			}
	 print "Linking ",$files->{split_dir}," to $link_name\n";
	 $action = sprintf("Linking %s to %s",$files->{split_dir},$link_name);
	 logAction($action,$log_file);
	 symlink($files->{split_dir},$link_name);
	}

certifyComplete($files->{target_dir},"sanda",$files->{fastq_base});

#################################SUBROUTINES###############################################

sub deleteFiles
{
	my $to_delete = shift;
	my $log = shift;
	my $action;
	if(ref($to_delete) eq "ARRAY")
		{
			foreach(@{$to_delete})
				{
					$action = sprintf("Deleting %s",$_);
					logAction($action,$log);
					unlink $_;
				}
		}
	elsif(!defined ref($to_delete))
		{
			$action = sprintf("Deleting %s",$to_delete);
			logAction($action,$log);
			unlink $to_delete;
		}
}

sub changeFormatStrings
{
	my $params = shift;
	
	$params->{pair_f} =~ s/\@/%/g;
	$params->{solo_f} =~ s/\@/%/g;
}

sub splitReadsByChr
{
	my $MicroSeq_path = shift;
	my $files = shift;
	my $log = shift;
	my $action;
	
	$files->{split_dir} = sprintf("%s.%ssplit",$files->{local_dir}, $files->{fastq_base});
	#remove split directory if it already exists
	if(-l $files->{split_dir})
		{
			my $linked_object = readlink($files->{split_dir});
			if(-d $linked_object)
				{
					$action = sprintf("Deleting old linked directory %s",$linked_object);
					logAction($action,$log);
					opendir(DIR,$linked_object);
					foreach(readdir DIR)
						{
							next FILE if $_ =~ m/^\.{1,2}$/;
							unlink sprintf("%s/%s",$linked_object,$_) or die "Could not delete $_ in $linked_object: $!\n";
						}
					closedir DIR;
					rmdir $linked_object or die "Could not remove directory $linked_object: $!\n";
				}
			else
				{
					$action = sprintf("Deleting old symbolic link %s",$linked_object);
					logAction($action,$log);
					unlink $linked_object;
				}
			$action = sprintf("Deleting old version of %s",$files->{split_dir});
			logAction($action,$log);
			unlink $files->{split_dir};
		}
	elsif(-d $files->{split_dir})
		{
			$action = sprintf("Deleting old directory %s",$files->{split_dir});
			logAction($action,$log);
			opendir(DIR,$files->{split_dir});
			FILE:foreach(readdir DIR)
				{
					next FILE if $_ =~ m/^\.{1,2}$/;
					unlink sprintf("%s/%s",$files->{split_dir},$_) or die "Could not delete $_ in ",$files->{split_dir},": $!\n";
				}
			closedir DIR;
			rmdir $files->{split_dir} or die "Could not remove directory ",$files->{split_dir},": $!\n";
		}
	elsif(-f $files->{split_dir})
		{
			$action = sprintf("Removing old file %s",$files->{split_dir});
			logAction($action,$log);
			unlink $files->{split_dir};
		}
	$files->{split_dir} .= '/';
	mkdir $files->{split_dir} or die "Could not create ",$files->{split_dir},": $!\n";
	$action = sprintf("Created directory %s",$files->{split_dir});
	logAction($action,$log);
	
	$files->{msi_link} = sprintf("%s%s.msi.sam",$files->{split_dir},$files->{fastq_base});
	print "Linked ",$files->{msi}, " to ", $files->{msi_link},"\n";
	$action = sprintf("Linked %s to %s",$files->{msi},$files->{msi_link});
	logAction($action,$log);
	symlink($files->{msi}, $files->{msi_link});
	my $split_cmd = makeSplitCmd($MicroSeq_path,$files->{msi_link});
	my $split_out = runCommand($split_cmd);
	printf("%s split\n",$files->{fastq_base});
	printf("\tstdout:\n%s\n",join "\t", @{$split_out->{stdout}}) if defined $split_out->{stdout};
	printf("\tstderr:\n%s\n",join "\t", @{$split_out->{stderr}}) if defined $split_out->{stderr};

}

sub reindexReads
{
	my $MicroSeq_path = shift;
	my $files = shift;

	#reindex mapped reads to unmodified genome
	#N.B.  ReindexSam needs to be retooled to work with se reads
	my $reindex_cmd = makeReindexCmd($MicroSeq_path,$files);
	my $reindex_out = runCommand($reindex_cmd);
	printf("%s reindex\n",$files->{fastq_base});
	printf("\tstdout:\n%s\n",join "\t", @{$reindex_out->{stdout}}) if defined $reindex_out->{stdout};
	printf("\tstderr:\n%s\n",join "\t", @{$reindex_out->{stderr}}) if defined $reindex_out->{stderr};
}

sub getSamFile
{
	my $bwa_path = shift;
	my $files = shift;
	#Generate sam file
	my $bwa_sam_cmd = makeSamCmd($bwa_path,$files);
	my $bwa_sam_out = runCommand($bwa_sam_cmd);
	print "Completed converting to SAM format\n";
}

sub alignReads
{
	my $bwa_path = shift;
	my $files = shift;
	my $gen_params = shift;
	
	#Align modified reads to modified genome
	my $children;

	for(my $i = 0; $i < ($files->{isPaired} + 1); $i++)
		{
			my $pid = fork();
			if($pid)
				{
					push @{$children}, $pid;
				}
			else
				{
					my $bwa_aln_cmd = makeAlnCmd($bwa_path, $gen_params, $files, $i);
					my $bwa_aln_out = runCommand($bwa_aln_cmd);
					printf("  completed %s -> %s alignment successfully\n",$files->{mod_fastq}->[$i], $files->{genome});
					exit 0;
				}
		}

	forkWait($children);

	print "Completed alignments\n";
}

sub zipModFastq
{
	my $files = shift;
	my $log = shift;
	
	my $action;
	my $children;
	
	for(my $i = 0; $i < scalar @{$files->{mod_fastq}}; $i++)
		{
			my $pid = fork();
			if($pid)
				{
					push @{$children}, $pid;
				}
			else
				{
					$action = sprintf("Zipping %s\n",$files->{mod_fastq}->[$i]);
					logAction($action,$log);
					my $command = sprintf("gzip -f %s",$files->{mod_fastq}->[$i]);
					my $zipOut = runCommand($command);
					printf("Zipped %s\n",$files->{mod_fastq}->[$i]); #this will overwrite an existing *mod.ms.* file
					exit 0;
				}
		}
	forkWait($children);
}

sub convertToBam
{
	my $samtools_path = shift;
	my $file = shift;
	my $log = shift;
	
	my $action = sprintf("Converting %s to BAM format",$file);
	logAction($action,$log);
	
	my $convert_cmd = makeBamConvertCmd($samtools_path,$file);
	my $convert_out = runCommand($convert_cmd);
}

sub makeBamConvertCmd
{
	my $samtools_path = shift;
	my $sam_file = shift;
		
	my $chr = $1 if $sam_file =~ m/(.*)\.msi\.sam/;
	my $bam_file = sprintf("%s.msi", $chr);
	
	my $command = sprintf("bash -c \'%s/samtools sort <(%s/samtools view -bS %s) %s\'",$samtools_path,$samtools_path,$sam_file,$bam_file);
	return $command;
}

sub findMicrosatellites
{
	my $MicroSeq_path = shift;
	my $tetrascan_params = shift;
	my $files = shift;
	
	my $tetrascan_cmd = makeTetrascanCmd($MicroSeq_path,$tetrascan_params,$files);
	my $tetrascan_out = runCommand($tetrascan_cmd);
	printf("%s tetrascan\n",$files->{fastq_base});
	printf("\tstdout:\n%s\n",join "\t", @{$tetrascan_out->{stdout}}) if defined $tetrascan_out->{stdout};
	printf("\tstderr:\n%s\n",join "\t", @{$tetrascan_out->{stderr}}) if defined $tetrascan_out->{stderr};
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
		"t" => "length of barcode sequence at the beginning of every read"
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
	
	printLogSectionBreak($log_file);
}

sub generateFileNames
{
	my $files = shift;
	my $gen_params = shift;
	my $file_prefix;
	my $file_stub;

	print "local dir is now ",$files->{local_dir},"\n";

	#Create names for output files
	if(! $files->{isPaired})
	{
		if(!$files->{zipped})
			{
				push @{$files->{fastq}}, sprintf("%s%s%s",					$files->{dir},$files->{fastq_base},$files->{suffix});
			}
		else
			{
				push @{$files->{fastq}}, sprintf("%s%s%s.gz",				$files->{dir},$files->{fastq_base},$files->{suffix});
			}
		push @{$files->{mod_fastq}}, sprintf("%s%s.mod.ms.txt", $files->{local_dir},$files->{fastq_base});
		push @{$files->{ms_fastq}},	 sprintf("%s%s.ms.txt",			$files->{local_dir},$files->{fastq_base});
		push @{$files->{hdr}}, 			 sprintf("%s%s.hdr.txt",		$files->{local_dir},$files->{fastq_base});
		push @{$files->{sai}}, 			 sprintf("%s%s.mod.sai",		$files->{local_dir},$files->{fastq_base});
	}
	else
	{
		for (my $i = 1; $i < 3; $i++)
		{
			$file_prefix = sprintf($files->{solo_f}, $files->{fastq_base},$i);
			if($files->{suffix} !~ m/\.bam$/)
				{
					if(!$files->{zipped})
						{
							push @{$files->{fastq}}, sprintf("%s%s%s",					$files->{dir},$file_prefix,$files->{suffix});
						}
					else
						{
							push @{$files->{fastq}}, sprintf("%s%s%s.gz",				$files->{dir},$file_prefix,$files->{suffix});
						}
				}
			else
				{
					$files->{isBam} = 1;
					$file_stub = sprintf($files->{pair_f}, $files->{fastq_base});
					$files->{fastq}->[0] = sprintf("%s%s%s",$files->{dir},$file_stub,$files->{suffix});
				}
			push @{$files->{mod_fastq}}, sprintf("%s%s.mod.ms.txt", $files->{local_dir},$file_prefix);
			push @{$files->{ms_fastq}},	 sprintf("%s%s.ms.txt",			$files->{local_dir},$file_prefix);
			push @{$files->{hdr}}, 			 sprintf("%s%s.hdr.txt",	  $files->{local_dir},$file_prefix);
			push @{$files->{sai}}, 			 sprintf("%s%s.mod.sai",	  $files->{local_dir},$file_prefix);
		}
	}
	
	for(my $i = 0; $i < scalar @{$files->{fastq}}; $i++)
		{
			if(-l $files->{fastq}->[$i])
				{
					my $file_path = [split /\//, $files->{fastq}->[$i]];
					my $new_file_name = sprintf("%s/%s",$files->{local_dir},$file_path->[$#{$file_path}]);
					print "linking ", readlink($files->{fastq}->[$i])," to $new_file_name\n";
					symlink(readlink($files->{fastq}->[$i]),$new_file_name);
					$files->{fastq}->[$i] = $new_file_name;
				}
			#if the file is a symbolic link in original dir, link the original file to the local directory
		}

	$file_prefix = sprintf($files->{pair_f}, $files->{fastq_base});
	$files->{sam} = sprintf("%s%s.mod.sam",$files->{local_dir},$file_prefix);
	$files->{msi} = sprintf("%s%s.msi.sam",$files->{local_dir},$file_prefix);
}

sub makeTetrascanCmd
{
	my $path = shift;
	my $params = shift;
	my $files = shift;

	my $command = undef;


	$command = "bash -c \'"	if $files->{zipped};

	$command .= sprintf("%s/tetrascan ", $path);

	foreach(keys %{$params})
	{
		if($_ !~ m/(M|qual_format|trim_type)/)
		{
			$command .= sprintf(" -%s %d", $_, $params->{$_});
		}
	}

	$command .= " -M" if $params->{M} == 0;
	$command .= sprintf(" -%s", $params->{$_}) foreach(qw(qual_format trim_type));
	if(!$files->{isBam})
		{
			if(!$files->{zipped})
				{
					$command .= sprintf(" -q %s", $files->{fastq}->[0]);
					$command .= sprintf(" -p %s", $files->{fastq}->[1]) if $files->{isPaired};
				}
			else
				{
					$command .= sprintf(" -q <(gunzip -c %s)", $files->{fastq}->[0]);
					if($files->{isPaired})
						{
							#this will not work for both file naming conventions
							my $output_prefix = sprintf($files->{pair_f},$files->{fastq_base});
							$output_prefix = $files->{local_dir}.$output_prefix;
							$command .= sprintf(" -p <(gunzip -c %s) -O %s", $files->{fastq}->[1],$output_prefix);
						}
					#Need to come up with alternative if files are unpaired
					$command .= "\'";
				}
		}
	else
		{	
			my $output_prefix = sprintf($files->{pair_f},$files->{fastq_base});
			$output_prefix = $files->{local_dir}.$output_prefix;
			$command .= sprintf(" -b %s -O %s",$files->{fastq}->[0], $output_prefix);
		}
	return $command;
}

sub makeAlnCmd
{
	my $bwa_path = shift;
	my $params = shift;
	my $files = shift;
	my $index = shift;
	
	my $command = sprintf("%s/bwa aln -t %d %s %s 1> %s", $bwa_path, $params->{threads}/($files->{isPaired} + 1), 
												$files->{genome}, $files->{mod_fastq}->[$index], $files->{sai}->[$index]);
	return $command;
}

sub makeSamCmd
{
	my $bwa_path = shift;
	my $files = shift;
	my $sam_prog = "samse";
	$sam_prog = "sampe" if $files->{isPaired};

	my $command = sprintf("%s/bwa %s %s %s %s 1> %s", $bwa_path, $sam_prog, $files->{genome},(join " ", @{$files->{sai}}),
												(join " ", @{$files->{mod_fastq}}), $files->{sam});
	return $command;
}

sub makeReindexCmd
{
	my $MicroSeq_path = shift;
	my $files = shift;
	
	my $command = sprintf("%s/reindexSam -S %s -I %s -1 %s", $MicroSeq_path, $files->{sam}, $files->{offset_index},
													$files->{hdr}->[0]);
	$command .= sprintf(" -2 %s", $files->{hdr}->[1]) if $files->{isPaired};
	
	return $command;
}

sub makeSplitCmd
{
	my $MicroSeq_path = shift;
	my $link = shift;
	
	my $command = sprintf("%s/splitByChr -S %s", $MicroSeq_path, $link);
	return $command;
}

sub filesToString
{
	my $files = shift;

	return $files->[0] if scalar @{$files} == 1;
	return ((join ",", @{$files}[0..(scalar @{$files} - 2)])." and ".($files->[(scalar @{$files} - 1)]));
}
