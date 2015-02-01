package MicroSeqFunc;

use strict;
use warnings;
use Sys::Hostname;
use Switch;
use Exporter;
use Cwd qw(abs_path);
use File::Basename;

BEGIN{
	# This BEGIN block opens the config file
	# to find where uSeq's custom perllib is,
	# then adds it to @INC so that any custom
	# modules can be loaded
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

use Timer;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION = 1.20;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(runCommand guessQualScore checkFastq getScriptName forkWait getGoodQueues createLocalMicroSeqDir getProgramPaths virtualFreeString2Int);
%EXPORT_TAGS = ( All 		 => [@EXPORT_OK],
								 all 		 => [@EXPORT_OK],
								 DEFAULT => [qw(runCommand getScriptName getProgramPaths virtualFreeString2Int)],
							 );

sub getProgramPaths {
	# loads program paths for BWA, SamTools, Picard
	# and uSeq from config.txt
	my $configFile = shift;
	
	$configFile = sprintf("%s/config.txt",dirname(abs_path($0))) if !defined $configFile;
	if(defined $ENV{'JOB_ID'}) {
		# script was launched using SGE
		# config.txt must be obtained from $JOB_SCRIPT
		# overrides any previous config file specification
		$configFile = $ENV{'USEQ_CONFIG'};
	}
	my $programPaths = do($configFile);

	die "Could not parse config file: $@" if $@;
	die "Could not read config file $configFile: $!" 	unless defined $programPaths;
	die "Could not run config file"				unless $programPaths;

	return $programPaths;
}

sub virtualFreeString2Int {
	my $vfString = shift;
	
	my $power  = substr($vfString,-1,1);
	my $number = substr($vfString,0,-1);
	
	switch($power) {
		case /G/i		{return $number * (10 ** 9)}
		case /M/i		{return $number * (10 ** 6)}
		case /K/i		{return $number * (10 ** 3)}
		else				{warn "$power is not recognized, allowing program to use its default\n"; 
									return undef;}
	}
}

sub createLocalMicroSeqDir {
	my $projectID = shift;
	my $nodename  = shift;
	my $log       = shift;
	my $created   = 0;
	my $uname     = $ENV{LOGNAME} || $ENV{USER};
	
	if(!defined $projectID) {
		die "You must specify a project ID with the --project_id tag if you are using MicroSeq in decentralized mode\n";
	}
	
	my $localDir = sprintf("/mnt/%s/data/unsafe/%s/MicroSeq/%s/",$nodename,$uname,$projectID);
	$log->mkDir($localDir);
	
	return $localDir;
}

sub getGoodQueues {
	# this function will currently only work on wigclust
	my $threshold_type = shift;
	my $threshold      = shift;
		
	$threshold_type = "percent" if !defined $threshold_type;
	$threshold      = 90        if !defined $threshold;
	
	my $mounted_fs_dir = "/mnt/";
	my $mnt_fh;

	opendir($mnt_fh, $mounted_fs_dir) or die "Could not open $mounted_fs_dir: $!\n";
	my $mounted_fs = [readdir $mnt_fh];
	closedir $mnt_fh;

	my $command;
	my $output;
	my $wigclust_nodes = undef;
	my $info;

	foreach(@{[grep { /wigclust/ } @{$mounted_fs}]}) {
		$command = sprintf("df -h /mnt/%s/data/ | grep -v Filesystem", $_);
		$output  = runCommand($command, undef, 0);

		chomp $output->{stdout}->[0];
		$info = [split / +/, $output->{stdout}->[0]];
		$wigclust_nodes->{$_}->{size}      = $info->[1];
		$wigclust_nodes->{$_}->{used}      = $info->[2];
		$wigclust_nodes->{$_}->{available} = $info->[3];
		$wigclust_nodes->{$_}->{percent}   = $1 if $info->[4] =~ m/(\d+)%/;
	}

	return undef if !defined $wigclust_nodes;

	my $allowed_queues;
	foreach(keys %{$wigclust_nodes}) {
		if($wigclust_nodes->{$_}->{$threshold_type} < $threshold){
			push @{$allowed_queues}, sprintf("all.q@%s.cshl.edu",$_);
		}
	}

	return undef if keys %{$wigclust_nodes} == scalar @{$allowed_queues};
	# if all queues are valid, returns undef
	return join ",", @{$allowed_queues};
}

sub forkWait {
# wait for forked processes to complete
	my $children = shift;
	
	foreach(@{$children}) {
		waitpid($_, 0);
	}
}

sub getScriptName{
	my $full_name = shift;
	my $prog_name = basename($full_name);
	
	#trim any suffix
	$prog_name =~ s/\.[^.]+$//;
	
	return $prog_name;
}

sub guessQualScore {
	my $seq_file = shift;
	my $zipped   = shift;
	my $isBam    = 0;
	
	if($seq_file =~ m/\.bam$/) {
		$isBam = 1;
		open SEQ, "/data/software/samtools/default/samtools view $seq_file |" or die "Could not open samtools view pipe to $seq_file: $!\n";
	}
	else {
	# assumes it's a BAM file if it's not a FASTQ file
		if($zipped){
			open SEQ, "gunzip -c $seq_file |" or die "Could not open gunzip pipe to file $seq_file: $!\n";
		}
		else {
			open SEQ, '<', $seq_file or die "guessQualScore: could not open $seq_file: $!\n";
		}
	}
		
	my $counter = 0;
	my $fields;
	if($isBam) {
		GUESS:while(my $line = <SEQ>) {
			chomp $line;
			$fields = [split /\t/, $line];

			if($fields->[10] =~ m/[0123456789:!"#\$%&'()*+,-.\/]/) {
				# qualities are in Sanger format (Phred + 33)
				return "a";
			}
			if($fields->[10] =~ m/^B+$/) {
			# qualities are in old Illumina format (Phred + 64)
				return "s";
			}
			if($fields->[10] =~ m/[KLMNOPQRSTUVWXYZ[\\\]\^_`abcdefgh]/) {
			# qualities are in old Illumina format (Phred + 64)
				return "s";
			}

			#At this point, I'm not trying to distinguish solexa + 64 (log-odds score) from phred + 64
			$counter++;
			last GUESS if $counter > 1000;
		}
	}
	else{
		GUESS:while(my $line = <SEQ>) {
			#skip to quality line
			$line = <SEQ>;
			$line = <SEQ>;
			die "$seq_file does not have the proper format\n" if $line !~ /^\+/;
			$line = <SEQ>;
			chomp $line;

			if($line =~ m/[0123456789:!"#\$%&'()*+,-.\/]/) {
			# qualities are in Sanger format (Phred + 33)
				return "a";
			}
			if($line =~ m/[KLMNOPQRSTUVWXYZ[\\\]\^_`abcdefgh]/) {
				return "s";
			}
			#At this point, I'm not trying to distinguish solexa + 64 (log-odds score) from phred + 64
			$counter++;
			last GUESS if $counter > 1000;
		}
	}
	warn "After looking through 1000 reads in $seq_file, still cannot guess the quality score encoding\n";
	return undef;
}

sub checkFastq {
	my $fastq_file = shift;
	my $mode       = shift;
	my $zipped     = shift;
		
	$mode ||= "quick";
	if($mode !~ /(brief|quick|full|complete)/) {
		$mode = "quick";
		print "Unrecognized option for checkFastq, setting mode to quick\n";
	}
	
	if($zipped) {
		open FASTQ, "gunzip -c $fastq_file |" or die "Could not open pipe to file $fastq_file: $!\n";
	}
	else {
		open FASTQ, '<', $fastq_file or die "Could not open $fastq_file: $!\n";
	}
	my $num_lines = 0;

	# check if first line of file starts with '@'
	my $line = <FASTQ>;
	if($line !~ m/^@/) {
		warn "Seq header is invalid in $fastq_file\n";
		return 0;
	}
	
	# check if second line only has ACGTNacgtn
	$line = <FASTQ>;
	chomp $line;
	if($line !~ m/^[ACGTNacgtn]+$/) {
		warn "Sequence contains unrecognized characters in $fastq_file\n";
		return 0;
	}
	
	# check that third line begins with '+'
	$line = <FASTQ>;
	if($line !~ m/^\+/) {
		warn "Quality header is invalid in $fastq_file\n";
		return 0;
	}
	
	if($mode eq "complete" or $mode eq "full") {
		$line = <FASTQ>;
		$num_lines += 4;
		while($line = <FASTQ>) {
			
			if($line !~ m/^@/) {
				warn "Seq header is invalid\n";
				return 0;
			}
			
			$line = <FASTQ>;
			chomp $line;
			if($line !~ m/^[ACGTNacgtn]$/) {
				warn "Sequence contains unrecognized characters\n";
				return 0;
			}

			$line = <FASTQ>;
			if($line !~ m/^\+/) {
				warn "Quality header is invalid\n";
				return 0;
			}
			
			$line = <FASTQ>; # no checks done on base quality line
			$num_lines += 4;
		}
		if(($num_lines % 4) != 0) {
		# check to make sure that number of lines in FASTQ file is
		# divisible by 4
			warn "$fastq_file does not have the correct amount of lines\n" ;
			return 0;
		}
	}
	close FASTQ;
	return 1;
}

sub runCommand {
	# deprecated, should use runCommand in Logger.pm instead
	my $command        = shift;
	my $allowed_errors = shift;
	my $print_command  = shift;

	$print_command = 1 if !defined $print_command;
	
	my $prog_name = @{[split / /, $command]}[0];
	$prog_name    = $1 if $prog_name =~ m/\/(.+)$/; 
	# splits off any leading slashes (e.g. for bwa)
	
	my $read_out = 0;
	my $read_err = 0;

	my $output;
	my $prog_out = sprintf(".%d.out", $$);
	my $prog_err = sprintf(".%d.err", $$);

	# check to see if command already redirects STDOUT or STDERR	
	if(($command !~ m/(1|12|21)>/ or $command =~ m/ *[^12]>/) and $command !~ m/>{1,2}/) {
		$command .= " 1>$prog_out";
		$read_out++;
	}
	if($command !~ m/(2|12|21)>/) {
		$command .= " 2>$prog_err";
		$read_err++;
	}
	
	my $action = sprintf("Running command %s\n", $command);
	
	#run command
	system($command);

	# load $prog_out from file if redirect wasn't specified
	if($read_out) {
		if(-s $prog_out) {
			open OUT, '<', $prog_out;
			$output->{stdout} = [<OUT>];
			close OUT;
		}
		else {
			$output->{stdout} = undef;
		}
		unlink $prog_out or warn "Could not delete $prog_out: $!\n";
	}
	
	#load $prog_err from file if redirect wasn't specified
	if($read_err) {
		if(-s $prog_err) {
			open ERR, '<', $prog_err;
			$output->{stderr} = [<ERR>];
			close ERR;
		}
		else {
			$output->{stderr} = undef;
		}
		unlink $prog_err or warn "Could not delete $prog_err: $!\n";
	}

	# if exit value is non-zero, print exit value from command, and 
	# any other output to STDERR if exit value isn't allowed
	if($? != 0 and !isAllowedError($?, $allowed_errors)) {
		my $err_out = join "", @{$output->{stderr}} if defined $output->{stderr};

		warn ">> Error running $prog_name\n\texit val: $?\n";
		warn "command attempted: $command\n";
		warn $err_out,"\n" if defined $err_out;
		exit 1;
	}
		
	return $output;
}

sub isAllowedError {
	my $error          = shift;
	my $allowed_errors = shift;
	
	return 0 if ! defined $allowed_errors->[0];
	foreach(@{$allowed_errors}) {
		return 1 if $_ == $error;
	}
	return 0;
}

1;
