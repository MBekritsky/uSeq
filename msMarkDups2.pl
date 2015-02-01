#!/usr/bin/perl
use strict;
use Getopt::Long;
use warnings;
use Switch;
use Carp qw(confess croak);
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
use Certificate;
use Logger;

sub getMergeBamCmd;
sub getMarkDupCmd;
sub getProfilerCmd;
sub cleanOldCountsDirs;

my $prog_name = "msMarkDups.pl";
my $step_name = "msMarkDups";
my $prog_paths = getProgramPaths();
if(defined $ENV{'JOB_ID'}) {
	$prog_paths = getProgramPaths($ENV{'USEQ_CONFIG'});
}

my $params;
my $action;
my $created;

GetOptions
(
	"target_dir=s"		=> \$params->{target_dir},
	"sam_basename=s"	=> \@{$params->{sam_basename}}, #this relies on the assumption that all the basenames are within target_dir...this can be fixed, but not yet
	"decentralized"		=> \$params->{decentralized},
	"keep_duplicates" => \$params->{keep_duplicates},
	"keep_all"				=> \$params->{keep_all},
	"zip_output"			=> \$params->{zip_output},
	"project_id:s"		=> \$params->{project_id},
	"db=s"						=> \$params->{msdb},
	"level=i"					=> \$params->{level},
	"virtual_free:s"		=> \$params->{virtual_free},
);

$params->{target_dir} = sprintf("%s/", $params->{target_dir}) if $params->{target_dir} !~ m/\/$/;

my $hostname = `hostname`;
my $nodename = [split /\./,$hostname];
$nodename = $nodename->[0];

print "Running $prog_name on $hostname\n";
my $log = Logger->new($params->{target_dir},$prog_name);
my $certifier = Certificate->new($params->{target_dir});
printf("Log file for this run of %s can be found at %s\n",$prog_name,$log->fname());

$params->{local_dir} = $params->{target_dir};
if($params->{decentralized})
{
	$params->{local_dir} = createLocalMicroSeqDir($params->{project_id},$nodename,$log);
	print "Local directory: ",$params->{local_dir},"\n";
}

if(defined $params->{virtual_free} and $params->{virtual_free} =~ m/[A-Za-z]$/)
{
	$params->{virtual_free_int} = virtualFreeString2Int($params->{virtual_free});
}

my $merge_cmd;
my $merge_out;
my $markDups_cmd;
my $markDups_out;
my $profiler_cmd;
my $profiler_out;

$params->{individual_name} = [split /\//, $params->{target_dir}];
$params->{individual_name} = $params->{individual_name}->[$#{$params->{individual_name}}];

$params->{msi_dir} = sprintf("%sMSI",$params->{target_dir});
$params->{merged_unmarked_bam} = sprintf("%s/%s.msi.bam",$params->{local_dir},$params->{individual_name});


if($params->{level} <= 4)
{
	opendir(MSI,$params->{msi_dir}) or confess "Could not open ",$params->{msi_dir},": $!";
	FILE:while(my $lane_file = readdir MSI)
	{
		next FILE if $lane_file =~ m/^\.{1,2}$/;
		$lane_file = sprintf("%s/%s",$params->{msi_dir},$lane_file);
		push @{$params->{target_dir_msi}}, $lane_file;
		if(-l $lane_file)
		{
			push @{$params->{local_msi}}, readlink($lane_file);
		}
		else
		{
			push @{$params->{local_msi}}, $lane_file;
		}
	}
	cleanOldCountsDirs($params,$log);
	
	$merge_cmd = getMergeBamCmd($prog_paths,$params);
	$merge_out = $log->runCommand($merge_cmd);

	$params->{local_mdup_file} = sprintf("%s%s.mdup.msi.bam",$params->{local_dir},$params->{individual_name});
	$params->{target_mdup_file} = sprintf("%s%s.mdup.msi.bam",$params->{target_dir},$params->{individual_name});

	$params->{mdup_metrics} = sprintf("%s%s.mdup.metrics.txt",$params->{target_dir},$params->{individual_name});
	$markDups_cmd = getMarkDupsCmd($prog_paths->{picard},$params);
	$markDups_out = $log->runCommand($markDups_cmd);
	$log->symlink({source => $params->{local_mdup_file}, target => $params->{target_mdup_file}, overwrite => 1});

	#delete the merged BAM file that was used by Picard MarkDuplicates
	$log->rm($params->{merged_unmarked_bam});
	#delete the BAM files from each individual lane
	$log->rm($params->{msi_dir},undef,1);

	$params->{local_msp} 				= sprintf("%s%s.msp",									$params->{local_dir},$params->{individual_name});
	$params->{local_dup} 				= sprintf("%s%s.dups.txt",						$params->{local_dir},$params->{individual_name});
	$params->{local_noise} 			= sprintf("%s%s.seq.noise.bam",				$params->{local_dir},$params->{individual_name});
	$params->{local_prof_count} = sprintf("%scounts/%s.profile.count",$params->{local_dir},$params->{individual_name});
}
else
{
	$params->{local_msp} 				= sprintf("%s%s.msp",									$params->{local_dir},$params->{individual_name});
	$params->{local_dup}				= sprintf("%s%s.dups.txt",						$params->{local_dir},$params->{individual_name});
	$params->{local_noise} 			= sprintf("%s%s.seq.noise.bam",				$params->{local_dir},$params->{individual_name});
	$params->{local_prof_count} = sprintf("%scounts/%s.profile.count",$params->{local_dir},$params->{individual_name});
	$params->{local_mdup_file} = sprintf("%s%s.mdup.msi.bam",$params->{local_dir},$params->{individual_name});
	$params->{target_mdup_file} = sprintf("%s%s.mdup.msi.bam",$params->{target_dir},$params->{individual_name});
}

$profiler_cmd = getProfilerCmd($prog_paths->{uSeq},$params);
$profiler_out = $log->runCommand($profiler_cmd);

$params->{master_msp} = sprintf("%s%s.msp",$params->{target_dir},$params->{individual_name});
$params->{master_dup} = sprintf("%s%s.dups.txt",$params->{target_dir},$params->{individual_name});
$params->{overlap_noise} = sprintf("%s%s.seq.noise.bam",$params->{target_dir},$params->{individual_name});
$params->{prof_count}	= sprintf("%scounts/%s.profile.count",$params->{target_dir},$params->{individual_name});

if($params->{zip_output})
{
	$log->zip($params->{local_msp});
	$log->zip($params->{local_dup});
	
	$params->{local_dup} = sprintf("%s.gz",$params->{local_dup});
	$params->{local_msp} = sprintf("%s.gz",$params->{local_msp});
	$params->{master_dup} = sprintf("%s.gz",$params->{master_dup});
	$params->{master_msp} = sprintf("%s.gz",$params->{master_msp});
}

if($nodename eq "wigclust10") #THIS SHOULD BE FIXED SO THAT THE NODE NAME CAN BE FLEXIBLE
{
	$log->mv({source => $params->{local_msp}, target => $params->{master_msp}});
	$log->mv({source => $params->{local_prof_count}, target => $params->{prof_count}});
	$log->mv({source => $params->{local_dup}, target => $params->{master_dup}}) if -f $params->{local_dup};
	$log->mv({source => $params->{local_noise}, target => $params->{overlap_noise}}) if -f $params->{local_noise};
}
else #if not on the home node, use scp to move files from one node to the other
{
	$log->smv({source => $params->{local_msp}, target => $params->{master_msp}, targetHost => "wigclust10"});
	$log->smv({source => $params->{local_prof_count}, target => $params->{prof_count}, targetHost => "wigclust10"});
	$log->smv({source => $params->{local_dup}, target => $params->{master_dup}, targetHost => "wigclust10"}) if -f $params->{local_dup};
	$log->smv({source => $params->{local_noise}, target => $params->{overlap_noise}, targetHost => "wigclust10"}) if -f $params->{local_noise};
}

$certifier->completed($step_name);	
#################################SUBROUTINES###############################################
sub cleanOldCountsDirs
{
	my $params = shift;
	my $log = shift;
	my $localCountDir;
	
	foreach(@{$params->{local_msi}})
	{
		$localCountDir = substr($_,0,rindex($_,"/")+1)."counts/";
		$log->rmDir($localCountDir,undef,1) if -d $localCountDir;
	}
}

sub getMergeBamCmd
{
	my $paths = shift;
	my $params = shift;
	my $virtual_free = $params->{virtual_free};
	$virtual_free = "6G" if !defined $virtual_free;
	# virtual free is used here to cap the amount of memory available
	# to Picard when running MergeSamFiles.jar.  We set a default
	# so that the JVM doesn't go out of control

	my $command = sprintf("java -jar -Xmx%s -XX:+UseSerialGC %s/MergeSamFiles.jar OUTPUT=%s AS=true MSD=true VALIDATION_STRINGENCY=LENIENT"
												,$virtual_free,$paths->{picard},$params->{merged_unmarked_bam});
	for(my $i = 0; $i < scalar @{$params->{local_msi}}; $i++)
	{
		$command .= sprintf(" INPUT=%s",$params->{local_msi}->[$i]);
	}
	return $command;
}

sub getMarkDupsCmd
{
	my $picard_path = shift;
	my $params = shift;
	my $virtual_free = $params->{virtual_free};
	$virtual_free = "6G" if !defined $virtual_free;
	# virtual free is used here to cap the amount of memory available
	# to Picard when running MarkDuplicates.jar.  We set a default
	# so that the JVM doesn't go out of control
	
	my $command = sprintf("java -jar -Xmx%s -XX:+UseSerialGC %s/MarkDuplicates.jar INPUT=%s OUTPUT=%s METRICS_FILE=%s AS=true VALIDATION_STRINGENCY=LENIENT"
												,$virtual_free,$picard_path,$params->{merged_unmarked_bam},$params->{local_mdup_file},$params->{mdup_metrics});
	return $command;
}

sub getProfilerCmd
{
	my $uSeq_path = shift;
	my $params = shift;
	
	my $command = sprintf("%s/profiler -bam %s -msdb %s",$uSeq_path,$params->{local_mdup_file},$params->{msdb});
	$command  .= " -k" if $params->{keep_duplicates};
	return $command;
}
