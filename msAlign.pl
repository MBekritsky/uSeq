#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
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
	die "Could not read config file $configFile: $!" 	unless defined $config;
	die "Could not run config file"				unless $config;
	unshift @INC, $config->{uSeqPerllib};
}

use Logger;
use Certificate;
use MicroSeqFunc qw(forkWait createLocalMicroSeqDir getProgramPaths); #change to config file var

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

my $prog_paths = getProgramPaths();
if(defined $ENV{'JOB_ID'}) {
	$prog_paths = getProgramPaths($ENV{'USEQ_CONFIG'});
}

my $gen_params;
my $files;

my $prog_name = "msAlign.pl";
my $step_name = "msAlign";

my $created;
my $action;

$files->{isPaired} = 0; #default is that files are not paired

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
	
	"clean"								=> \$gen_params->{clean},
	"threads:i"						=> \$gen_params->{threads},
	"zip_output"					=> \$gen_params->{zip_output},
	"keep_all"						=> \$gen_params->{keep_all},
	"decentralized"				=> \$gen_params->{decentralized},
	"project_id:s"				=> \$gen_params->{project_id},
);

my $hostname = `hostname`;
my $nodename = [split /\./,$hostname];
$nodename = $nodename->[0];

##WHAT IS THREADS DEFAULT?

print "Running $prog_name on $hostname\n";

changeFormatStrings($files);
my $log = Logger->new($files->{target_dir},$prog_name,sprintf($files->{pair_f}, $files->{fastq_base})); #need explicit name here in case job is submitted to SGE
my $certifier = Certificate->new($files->{target_dir});

$files->{local_dir} = $files->{target_dir};
if($gen_params->{decentralized})
	{
		$files->{local_dir} = createLocalMicroSeqDir($gen_params->{project_id},$nodename,$log);
		print "Local directory: ",$files->{local_dir},"\n";
	}
	
#get filenames
generateFileNames($files,$gen_params,$log);
print "Log file for this run can be found at ", $log->fname(),"\n";
$log->recordParams($files,"File");
$log->recordParams($gen_params,"General");

my $lane_name = sprintf($files->{pair_f}, $files->{fastq_base});
printf("msAlign Processing %s\n", $lane_name);

alignReads($prog_paths->{bwa},$files, $gen_params,$log);
print STDERR "Aligned reads\n";

$certifier->completed($step_name,$lane_name);

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
					if(! (-l $new_file_name and readlink($new_file_name) eq readlink($files->{fastq}->[$i])))
					{
						print "linking ", readlink($files->{fastq}->[$i])," to $new_file_name\n";
						$log->symlink({source => readlink($files->{fastq}->[$i]),target => $new_file_name, overwrite => 1});
					}
					$files->{fastq}->[$i] = $new_file_name;
				}
			#if the file is a symbolic link in original dir, link the original file to the local directory
		}

	$file_prefix = sprintf($files->{pair_f}, $files->{fastq_base});
	$files->{sam} = sprintf("%s%s.mod.sam",$files->{local_dir},$file_prefix);
	$files->{msi} = sprintf("%s%s.msi.sam",$files->{local_dir},$file_prefix);
}

sub alignReads
{
	my $bwa_path = shift;
	my $files = shift;
	my $gen_params = shift;
	my $log = shift;
	
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
					my $bwa_aln_out = $log->runCommand($bwa_aln_cmd);
					printf("  completed %s -> %s alignment successfully\n",$files->{mod_fastq}->[$i], $files->{genome});
					exit 0;
				}
		}

	forkWait($children);

	print "Completed alignments\n";
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
