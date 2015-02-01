#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Switch;
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


use MicroSeqFunc qw(createLocalMicroSeqDir virtualFreeString2Int getProgramPaths);
use Logger;
use Certificate;

sub findMicrosatellites;
sub alignReads;
sub getSamFile;
sub reindexReads;
sub splitReadsByChr;
sub changeFormatStrings;
sub convertToBam;

sub filesToString;
sub makeTetrascanCmd;
sub makeAlnCmd;
sub makeSamCmd;
sub makeReindexCmd;
sub makeSplitCmd;
sub makeBamConvertCmd;
sub generateFileNames;
sub virtualFreeString2Int;

my $prog_paths = getProgramPaths();
if(defined $ENV{'JOB_ID'}) {
	$prog_paths = getProgramPaths($ENV{'USEQ_CONFIG'});
}

my $program_name = "msPostProcess.pl";
my $step_name = "msPostProcess";

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
	"paired"							=> \$files->{isPaired},
	"suffix=s"						=> \$files->{suffix},
	"zipped"							=> \$files->{zipped},
	"genome=s"						=> \$files->{genome},
	"offset_index=s"			=> \$files->{offset_index},
	
	"clean"								=> \$gen_params->{clean},
	"level:i"							=> \$gen_params->{level},
	"split"								=> \$gen_params->{split},
	"zip_output"					=> \$gen_params->{zip_output},
	"keep_all"						=> \$gen_params->{keep_all},
	"decentralized"				=> \$gen_params->{decentralized},
	"project_id:s"				=> \$gen_params->{project_id},
	"virtual_free=s"			=> \$gen_params->{virtual_free},
);

my $hostname = `hostname`;
my $nodename = [split /\./,$hostname];
$nodename = $nodename->[0];

$gen_params->{level} = 0 if !defined $gen_params->{level};
if(defined $gen_params->{virtual_free} and $gen_params->{virtual_free} =~ m/[A-Za-z]$/)
{
	$gen_params->{virtual_free_int} = virtualFreeString2Int($gen_params->{virtual_free});
}

print "Running $program_name on $hostname\n";

changeFormatStrings($files);
my $log = Logger->new($files->{target_dir},$program_name,sprintf($files->{pair_f},$files->{fastq_base}));
my $certifier = Certificate->new($files->{target_dir});

$files->{local_dir} = $files->{target_dir};
$files->{individual_name} = [split /\//, $files->{target_dir}];
$files->{individual_name} = $files->{individual_name}->[$#{$files->{individual_name}}];

if($gen_params->{decentralized})
	{
		$files->{local_dir} = createLocalMicroSeqDir($gen_params->{project_id},$nodename,$log);
		print "Local directory: ",$files->{local_dir},"\n";
	}
	
#get filenames
generateFileNames($files,$gen_params,$log);
print "Log file for this run can be found at ", $log->fname(),"\n";

#print parameters to log file
$log->recordParams($files,"File");
$log->recordParams($gen_params,"General");
my $lane_name = sprintf($files->{pair_f}, $files->{fastq_base});
printf("%s processing %s\n", $program_name,$lane_name);

#convert sai files to sorted, reindexed BAM
getSortedReindexedBam($prog_paths,$files,$gen_params->{virtual_free_int},$gen_params->{virtual_free},$log);
print STDERR "Got sorted and reindexed BAM file\n";

#finalize mod fastq files--either delete them if MicroSeq is being run in clean mode,
#or zip them if it's being told to zip output
if(!$gen_params->{keep_all})
{
	$log->rm($files->{sai});
	$log->rm($files->{hdr});
}

#move count files from local directory to home directory
my $fastq_stub=sprintf($files->{pair_f},$files->{fastq_base});
my $fullSourcePath = sprintf("%s/%s.reindex.count",$files->{local_counts_dir},$fastq_stub);
my $fullTargetPath = sprintf("%s/%s.reindex.count",$files->{target_counts_dir},$fastq_stub);
$log->mv($fullSourcePath,$fullTargetPath) if $gen_params->{decentralized};

my $basename;
my $link_name;

#link mod fastq files to home directory
if(!$gen_params->{keep_all})
{
	$log->mkDir($files->{linked_modMs});
	foreach(@{$files->{mod_fastq}})
	{
		$basename = [split /\//, $_];
		$basename = $basename->[$#{$basename}];
		$link_name = sprintf("%s/%s",$files->{linked_modMs},$basename);
		my $origin_name = $_;
		if($gen_params->{zip_output})
			{
				$log->zip($_);
				$origin_name .= ".gz";
				$link_name .= ".gz";
			}
		if(-e $link_name)
			{
				$log->record("$link_name already exists!  Removing...");
				$log->rm($link_name);
			}
		$log->symlink({source => $origin_name,target => $link_name,overwrite => 1});
	}
}
else
{
	$log->rm($files->{mod_fastq});
}

#link sorted, reindexed BAM files back to home directory
 $log->mkDir($files->{sbam_dir});
 $log->symlink({source => $files->{sbam}, target => $files->{linked_sbam}, overwrite => 1});

 $certifier->completed($step_name,$lane_name);

###############################################

sub getSortedReindexedBam
{
	my $paths = shift;
	my $files = shift;
	my $virtual_free_int = shift;
	my $virtual_free = shift;
	my $log = shift;
	
	#Generate sam file
	my $bwa_sam_cmd = makeSortedReindexedBamCmd($paths->{bwa},$paths->{samtools},$paths->{picard},
																							$paths->{uSeq},$files,$virtual_free_int,$virtual_free);
	my $bwa_sam_out = $log->runCommand($bwa_sam_cmd);
	print "Converted to SAM format\n";
}

sub makeSortedReindexedBamCmd
{
	my $bwa_path = shift;
	my $samtools_path = shift;
	my $picard_path = shift;
	my $uSeq_path = shift;
	my $files = shift;
	my $virtual_free_int = shift;
	my $virtual_free = shift;
	my $sam_prog = "samse";
	$sam_prog = "sampe" if $files->{isPaired};
	
	$virtual_free = "6G" if !defined $virtual_free;
	# virtual free is used here to cap the amount of memory available
	# to Picard when running AddOrReplaceGroups.jar.  We set a default
	# so that the JVM doesn't go out of control

	#there are occasionally INVALID_INDEXING_BIN errors when running Picard tools, but they seem to have no downstream issues.  Instead of troubleshooting, we're choosing
	#to set VALIDATION_STRINGENCY to LENIENT
	my $command = sprintf("%s/bwa %s %s %s %s | %s/samtools view -Shu - | %s/reindexBam -B stdin -I %s -1 %s -2 %s -T %s -P %s -o | %s/samtools sort", 
												$bwa_path, $sam_prog, $files->{genome},(join " ", @{$files->{sai}}),(join " ", @{$files->{mod_fastq}}), $samtools_path, $uSeq_path, 
												$files->{offset_index}, $files->{hdr}->[0], $files->{hdr}->[1], $files->{local_dir}, $files->{file_prefix}, $samtools_path);

	if(defined $virtual_free_int) {
	# set the maximum memory used by samtools sort, given in bytes
		$command .= sprintf(" -m %d",$virtual_free_int)
	}
	$command .= sprintf(" -o - %s | java -jar -Xmx%s -XX:+UseSerialGC %s/AddOrReplaceReadGroups.jar INPUT=/dev/stdin OUTPUT=%s ID=%s LB=%s PL=ILLUMINA SM=%s PU=%s VALIDATION_STRINGENCY=LENIENT",
												$files->{sbam_prefix},$virtual_free,$picard_path,$files->{sbam},sprintf("%s_%s",$files->{individual_name},$files->{file_prefix}),substr($files->{file_prefix},0,-1),$files->{individual_name},
												substr($files->{file_prefix},-1,1));
	return $command;
}

sub changeFormatStrings
{
	my $params = shift;
	
	$params->{pair_f} =~ s/\@/%/g;
	$params->{solo_f} =~ s/\@/%/g;
}

sub generateFileNames
{
	my $files = shift;
	my $gen_params = shift;
	my $log = shift;
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
							push @{$files->{fastq}}, sprintf("%s%s%s",		$files->{dir},$file_prefix,$files->{suffix});
						}
					else
						{
							push @{$files->{fastq}}, sprintf("%s%s%s.gz",	$files->{dir},$file_prefix,$files->{suffix});
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
					if( (!-e $new_file_name) and readlink($new_file_name) eq readlink($files->{fastq}->[$i]))
					{
						print "linking ", readlink($files->{fastq}->[$i])," to $new_file_name\n";
						$log->symlink({source => readlink($files->{fastq}->[$i]), target => $new_file_name});
					}
					$files->{fastq}->[$i] = $new_file_name;
				}
			#if the file is a symbolic link in original dir, link the original file to the local directory
		}

	$files->{file_prefix} = sprintf($files->{pair_f}, $files->{fastq_base});
	$files->{sbam_prefix} = sprintf("%s%s.msi",$files->{local_dir},$files->{file_prefix});
	$files->{sbam} = sprintf("%s%s.msi.bam",$files->{local_dir},$files->{file_prefix});
	$files->{sbam_dir} = sprintf("%s/MSI",$files->{target_dir});	
	$files->{linked_sbam} = sprintf("%s/%s.msi.bam",$files->{sbam_dir},$files->{file_prefix});

	$files->{linked_modMs} = sprintf("%s/ModMs",$files->{target_dir});

	$files->{target_counts_dir} = sprintf("%s/counts",$files->{target_dir});
	$files->{local_counts_dir}	= sprintf("%s/counts",$files->{local_dir});
}
