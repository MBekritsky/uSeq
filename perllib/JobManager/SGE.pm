package JobManager::SGE;

use lib "/mnt/wigclust4/home/bekritsk/tools/microsatellite/perllib";

use strict;
use warnings;
use Exporter;
use MicroSeqFunc qw(runCommand);
use Timer;
use Logger;
use Carp;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION = 1.20;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(sgeWait addSgeJob printSgeCompletionInfo);
%EXPORT_TAGS = ( All 		 => [@EXPORT_OK],
								 all 		 => [@EXPORT_OK],
								 DEFAULT => [qw(sgeWait addSgeJob printSgeCompletionInfo)],
							 );

sub printSgeCompletionInfo
{
	my $jobs = shift;
	my $log = shift;
	
	my $success = 1;
	
#	$fh = *STDOUT if !defined $fh;
	
		$log->record("Time until job completion:");
		foreach(sort keys %{$jobs})
			{
				$log->record(sprintf("%s: %s\n", $jobs->{$_}->{name}, $jobs->{$_}->{elapsed})) if $jobs->{$_}->{status} eq "c";
				if($jobs->{$_}->{status} eq "Eqw")
					{
						$log->record(sprintf("%s: killed due to %s\n", $jobs->{$_}->{name}, $jobs->{$_}->{error_reason}));
						$success = 0;
					}
			}
		if(!$success)
		{
			die "Some jobs did not complete successfully.  Please check the errors and attempt to fix them, then restart MicroSeq\n";
		}
}

sub addSgeJob
{
	my $jobs = shift;
	my $output = shift;
	my $description = shift;
	
	my $job_id;
	
	$job_id = $1 if $output->{stdout}->[0] =~ m/Your job ([0-9]+).*/;
	$jobs->{$job_id}->{start} = Timer->new();
	$jobs->{$job_id}->{name} = $description;
	$jobs->{$job_id}->{status} = "u";
		
	return ($jobs,$job_id);
}

sub setJobStates
{
	my $jobs = shift;

	my $command;
	my $job_status;
	my $error_reason;
	
	my $qstat = MicroSeqFunc::runCommand("qstat",[256],0);
	my $current_job_states = [ grep { /^\d+/ } @{$qstat->{stdout}}];

	my $current_job_hash;
	my $job_array;
	foreach(@{$current_job_states})
		{
			$job_array = [split /\s+/, $_];
			$current_job_hash->{$job_array->[0]}->{state} = $job_array->[4];
			$current_job_hash->{$job_array->[0]}->{node} = $job_array->[7];
		}
		
	foreach	my $job_id (keys %{$jobs})
		{
			if(defined $current_job_hash->{$job_id})
				{
					if($current_job_hash->{$job_id}->{state} eq "hqw")
						{
							$jobs->{$job_id}->{status} = "hqw";
						}
					elsif($current_job_hash->{$job_id}->{state} eq "qw")
						{
							$jobs->{$job_id}->{status} = "qw";
						}
					elsif($current_job_hash->{$job_id}->{state} eq "r")
						{
							$jobs->{$job_id}->{status} = "r";
							$jobs->{$job_id}->{node} = $current_job_hash->{$job_id}->{node};
						}
					elsif($current_job_hash->{$job_id}->{state} eq "Eqw")
						{
							$jobs->{$job_id}->{status} = "Eqw";
							$error_reason = jobInErrorState($job_id, $jobs);
							warn $error_reason,"\n";
						}
				}
			else
				{
					$jobs->{$job_id}->{status}  = 'c';
					$jobs->{$job_id}->{elapsed} = $jobs->{$job_id}->{start}->elapsed() or Carp::confess("Couldn't get elapsed time using Timer");
				}
		}
}

sub jobOrJobs
{
	my $count = shift;
	$count > 1 ? return "jobs" : return "job";
}

sub jobsInState
{
	my $jobs = shift;
	my $state = shift;
	
	my $count = 0;
	
	foreach(keys %{$jobs})
		{
			$count++ if $jobs->{$_}->{status} eq $state;
		}
	return $count;
}

sub jobsNotDone
{
	my $jobs = shift;
	
	return (jobsInState($jobs,"r") + jobsInState($jobs,"qw") + jobsInState($jobs,"hqw"));
}

sub jobInErrorState
{
	my $job_id = shift;
	my $jobs = shift;
	
	my $kill_status;
	my $reason;
	
	my $command = sprintf("qstat -j %d",$job_id);
	my $status = MicroSeqFunc::runCommand($command,[256],0);
	
	foreach(@{$status->{stdout}})
		{
			if($_ =~ m/^error reason/)
				{
					warn "Job $job_id is in an error state, killing job\n";
					
					$reason = $1 if $_ =~ m/(error:.*)$/;
					$jobs->{$job_id}->{error_reason} = $reason;

					$command = sprintf("qdel %s", $job_id);
					$kill_status = MicroSeqFunc::runCommand($command);
					if(defined $kill_status->{stderr})
						{
							warn "Could not delete job $job_id:\n";
							warn join @{$kill_status->{stderr}};
						}
					return $reason;
				}
		}
	return undef;
}

sub jobIsWaiting
{
	my $status = shift;
	
	foreach(@{$status->{stdout}})
		{
			return 0 if $_ =~ m/^usage/;
		}
	return 1;
}

sub sgeWait
{
	my $jobs = shift;
	my $report_interval = shift; #time between reportings of how many jobs remain
	my $command;
	my $job_status;
	my $jobStatusCounts;

	$report_interval ||= 120;
	my $sleep_time = 240;
	
	my $checks_per_interval = ($report_interval / $sleep_time);
	my $checks_since_last_report = 0;

	printf("Monitoring %d job", scalar keys %{$jobs});
	print "s" if (scalar keys %{$jobs}) > 1;
	print "\n";

	setJobStates($jobs);

	while(jobsNotDone($jobs) > 0)
	{
		sleep $sleep_time;

		setJobStates($jobs);

		$jobStatusCounts->{Eqw} = jobsInState($jobs,"Eqw");
		$jobStatusCounts->{hqw} = jobsInState($jobs,"hqw");
		$jobStatusCounts->{qw} = jobsInState($jobs,"qw");
		$jobStatusCounts->{r} = jobsInState($jobs,"r");
		$jobStatusCounts->{c} = jobsInState($jobs,"c");

		$checks_since_last_report++;

		if($checks_since_last_report == $checks_per_interval)
			{
				if(jobsNotDone($jobs) > 0)
					{
						printf("%d %s running", $jobStatusCounts->{r},jobOrJobs($jobStatusCounts->{r}));
						printf(", %d %s waiting", $jobStatusCounts->{qw},jobOrJobs($jobStatusCounts->{qw})) if $jobStatusCounts->{qw} > 0;
						printf(", %d %s holding", $jobStatusCounts->{hqw},jobOrJobs($jobStatusCounts->{hqw})) if $jobStatusCounts->{hqw} > 0;
						printf(", %d %s completed", $jobStatusCounts->{c},jobOrJobs($jobStatusCounts->{c})) if $jobStatusCounts->{c} > 0;
						printf(", %d %s in error state", $jobStatusCounts->{Eqw}, jobOrJobs($jobStatusCounts->{Eqw})) if $jobStatusCounts->{Eqw} > 0;
						print "\n";
						$checks_since_last_report = 0;
					}
			}
	}
}
