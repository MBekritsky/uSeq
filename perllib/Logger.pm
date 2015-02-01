package Logger;

use strict;
use warnings;
use Carp ();
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Handle;
use Sys::Hostname;
use File::Copy;
use File::Path;
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
use MicroSeqFunc qw(getScriptName);

my $VERSION = 1.05;

sub new {
	my $class     = shift;
	my $targetDir = shift;
	my $progname  = shift;
	my $filename  = shift;
	
	my $logDir = sprintf("%s/logs",$targetDir);
	if(! -d $logDir) {
		#avoids mkdir race condition, thanks to http://www.perlmonks.org/?node_id=897280
		my $error;
		mkdir $logDir or $error = $!;
		unless( -d $logDir) {
			Carp::confess("Error: Failed to create log directory $logDir: $!");
		}
	}
	
	$progname = getScriptName($progname);

	my $timer = Timer->new();
	my $fname = sprintf("%s/%s_%s",$logDir, $progname, $timer->fstring());
	$fname    = sprintf("%s_%s", $fname, $filename) if defined $filename;
	$fname    = sprintf("%s.log", $fname);
	
	open my $fptr, '>', $fname or Carp::confess("Error: Failed to initialize logfile $fname: $!");
	$fptr->autoflush(1);
	my $self = {
		_fname => $fname,
		_ptr   => $fptr,
		_timer => $timer,
		_pname => $progname,
		_init  => 0,
		_hostname => undef,
	};
	bless $self, $class;
	
	$self->init();
	return $self;
}

sub DESTROY {
	my $self = shift;
	close $self->{_ptr};
}

sub init {
# prints a line at beginning of log file with program name, hostname
# and start time
	my $self = shift;
	$self->{_hostname} = hostname;
	chomp $self->{_hostname};
	
	if(!$self->{_init}) {
		printf {$self->{_ptr}} ("%s launched on %s at %s\n", $self->{_pname},
		$self->{_hostname}, $self->{_timer}->string());
		$self->{_init} = 1;
	}
	else {
		Carp::cluck("Warning: ", $self->{_fname}, " has already been initialized");
	}
}

sub sectionBreak {
# prints a section break in the log file
	my $self = shift;
	print {$self->{_ptr}} "-" x 100;
	print {$self->{_ptr}} "\n";
}

sub longmessToLog {
# When calling cluck or confess, prints stack trace to log file as well
	my $self 		= shift;
	my $message = shift;
	
	print {$self->{_ptr}} Carp::longmess($message);
}

sub confess {
# a tweaked version of Carp's confess
	my $self    = shift;
	my $message = shift;
	
	$self->longmessToLog("Error: ".$message);
	Carp::confess($message);
}

sub cluck {
# a tweaked version of Carp's cluck
	my $self    = shift;
	my $message = shift;
	
	$self->longmessToLog("Warning: ".$message);
	Carp::cluck($message);
}

sub record {
# prints a message to the log file with a time stamp
	my $self        = shift;
	my $action      = shift;
	my $actionTimer = shift;
	
	if(defined $action) {
		printf {$self->{_ptr}} ("%s %s", $self->{_timer}->stamp(1), $action);
		printf {$self->{_ptr}} (" in %s", $actionTimer->elapsed()) if defined $actionTimer;
		print  {$self->{_ptr}} "\n";
	}
	else {
		my $errMessage = "Logger::record called without an action to record";
		$self->cluck($errMessage);
	}
}

sub recordParams {
# prints a list of parameters to the log file, with comments explaining 
# the parameters if they are provided
	my $self           = shift;
	my $params         = shift;
	my $paramType      = shift;
	my $optComments    = shift;
	my $choiceComments = shift;
	
	print {$self->{_ptr}} "$paramType " if defined $paramType;
	print {$self->{_ptr}} "params:\n";
	
	foreach(sort keys %{$params}) {
		if(defined $params->{$_}) {
			printf {$self->{_ptr}} ("%s: %s", $_, $params->{$_}) if ref($params->{$_}) ne 'ARRAY';
			
			# If the parameter passed is an array reference, print the array as a list
			printf {$self->{_ptr}} ("%s: %s", $_,join(",", @{$params->{$_}})) if ref($params->{$_}) eq 'ARRAY';
			
			# If there is a comment describing the option, print it
			printf {$self->{_ptr}} (" #%s", $optComments->{$_}) if defined $optComments->{$_};
			
			# If there is a comment describing the choice, print it
			printf {$self->{_ptr}} (" #%s", $choiceComments->{$params->{$_}}) if defined $choiceComments->{$params->{$_}};
			
			print  {$self->{_ptr}} "\n";
		}
	}
	$self->sectionBreak();
}

sub recordError {
#prints an error report to the log file with the attempted command, the exit message and the exit code
	my $self       = shift;
	my $exitVal    = shift;
	my $attempted  = shift;
	my $errMessage = shift;
	
	$self->sectionBreak();
	printf {$self->{_ptr}} (">>Error running %s\n", $self->{_pname});
	printf {$self->{_ptr}} ("Command attempted:\n\t%s\n", $attempted);
	printf {$self->{_ptr}} ("Error message:\n\t%s\n", $errMessage) if defined $errMessage;
	printf {$self->{_ptr}} ("Exit value %d\n", $exitVal);
}

sub recordElapsed {
# same as record, but if actionTimer is not defined, clucks a warning to the user and 
# records in log with "NO TIME PROVIDED"
	my $self        = shift;
	my $action      = shift;
	my $actionTimer = shift;
	
	if(defined $actionTimer) {
		$self->record($action, $actionTimer);
	}
	else {
		$action .= " NO TIME PROVIDED";
		Carp::cluck("Warning: In Logger::recordElapsed, Timer object is undefined");
		$self->record($action);
	}
}

sub rm {
# Removes a group of files/symlinks/dirs and logs actions.  If $parentPath is specified,
# prepends $parentPath to each item in $toRemove. If $deepRm is true and rm is passed
# a symlink, it will also remove the symlink's target
	my $self       = shift;
	my $toRemove   = shift;
	my $parentPath = shift;
	my $deepRm     = shift;
	$deepRm ||= 0;
	
	#Even if only one file is passed, makes it an array reference to keep things neat
	unless(ref($toRemove)) {
		$toRemove = [ $toRemove ];
	}
	
	if(defined $parentPath) {
		for(my $i = 0; $i < scalar @{$toRemove}; $i++) {
			$toRemove->[$i] = sprintf("%s/%s", $parentPath, $toRemove->[$i]);
		}
	}
	
	foreach my $object (@{$toRemove}) {
		if(-l $object) {
			$self->rmLink($object, undef, $deepRm);
		}
		elsif(-f $object) {
			unlink $object or $self->confess(sprintf("In Logger::rm, failed to remove %s: %s", $object, $!));
			$self->record(sprintf("Removed %s", $object));
		}
		elsif(-d $object) {
			$self->rmDir($object, undef, $deepRm);
		}
	}
}

sub rmDir {
# Removes a directory and all its contents and logs its actions.  If deepRm is specified
# and any of the files in the directory is a symbolic link, remove the links target as well
	my $self       = shift;
	my $dir        = shift;
	my $parentPath = shift;
	my $deepRm     = shift;
	$deepRm ||= 0;

	my $target = undef;
	
	$dir = sprintf("%s/%s", $parentPath, $dir) if defined $parentPath;
	
	opendir DDIR, $dir or $self->confess("In Logger::rmDir, failed to open $dir: $!");
	DEL:foreach my $file (readdir DDIR) {
	# before deleting the directory, anything in the directory must first be deleted
		next DEL if $file =~ m/^\.{1,2}$/;
		$target = sprintf("%s/%s", $dir, $file);
		
		$self->rmLink($target, undef, $deepRm) if -l $target;
		$self->rm($target, undef, $deepRm)     if -f $target;
		$self->rmDir($target, undef, $deepRm)  if -d $target;
	}
	closedir DDIR;
	rmdir $dir or $self->confess("In Logger::rmDir, failed to remove $dir: $!");
	$self->record(sprintf("Removed directory %s", $dir));
}

sub rmLink {
# Removes a symbolic link and logs actions.  If $deepRm is true, then removes the link's
# target as well
	my $self       = shift;
	my $link       = shift;
	my $parentPath = shift;
	my $deepRm     = shift;
	$deepRm ||= 0;
	
	my $target = undef;
	
	$link = sprintf("%s/%s", $parentPath, $link) if defined $parentPath;
	$target = readlink($link);
	unlink $link or $self->confess("In Logger::rmLink, failed to remove symbolic link $link: $!");
	$self->record(sprintf("Removed symbolic link %s", $link));

	if($deepRm) {
		$self->rmLink($target, undef, $deepRm) if -l $target;
		$self->rm($target, undef, $deepRm)     if -f $target;
		$self->rmDir($target, undef, $deepRm)  if -d $target;
	}
	else {
		$self->record(sprintf("PLEASE NOTE: Removed symbolic link %s, but not its target %s", $link,
													$target));
	}
}

sub symlink {
	my $self   = shift;
	my @params = @_;
	
	my $overwrite = 0;
	$overwrite = 1 if(ref($params[0]) eq "HASH" and defined $params[0]->{overwrite});
	
	my %linkFiles = %{$self->resolveSourceAndTarget("symlink", @params)};
	if(-l $linkFiles{target} and $overwrite > 0) {
		$self->rmLink($linkFiles{target}, undef, 0);
	}
	elsif($overwrite == 0) {
		$self->confess("In Logger::symlink, failed to link $linkFiles{target} to $linkFiles{source}: $linkFiles{target} already exists and overwrite not allowed");
	}
	symlink($linkFiles{source}, $linkFiles{target}) or $self->confess("In Logger::symlink, failed to link $linkFiles{target} to $linkFiles{source}: $!");
	$self->record(sprintf("Linked %s to %s", $linkFiles{target}, $linkFiles{source}));
}

sub mkDir {
	my $self   = shift;
	my $newDir = shift;
	
	if(! -d $newDir) {
		my $error;
		File::Path::mkpath($newDir) or $error = $!;
		unless(-d $newDir) {
			$self->confess("In Logger::mkDir, failed to create directory $newDir: $!");
		}
		$self->record(sprintf("Created directory %s", $newDir));
	}
	else {
		$self->record(sprintf("Directory %s already exists", $newDir))
	}
}

sub zip {
# Zips input file using IO::Compress::gzip and logs actions (can add other compression types), can pass
# a list of files, or an array reference to a list of files to be zipped
	my $self   = shift;
	my @params = @_;
	my @oldFiles;
	
	if(ref($params[0]) eq "ARRAY") {
		@oldFiles = @{$params[0]}
	}
	else {
		@oldFiles = @params;
	}
	
	foreach my $oldFile (@oldFiles) {
		my $newFile = sprintf("%s.gz", $oldFile);
		gzip($oldFile => $newFile, AutoClose => 1) or $self->confess("Gzip failed in Logger::zip: $GzipError");
		$self->record(sprintf("Gzipped %s to %s", $oldFile, $newFile));
		$self->rm($oldFile);
	}
}

sub mv {
# Moves file from source to target using perl's move command, logs actions.  Input can be two scalars,
# an array reference, or a hash reference
	my $self   = shift;
	my @params = @_;
	
	my %linkFiles = %{$self->resolveSourceAndTarget("move",@params)};
	move($linkFiles{source}, $linkFiles{target}) or $self->confess(sprintf("In Logger::move, failed to move %s to %s: %s",
																																					$linkFiles{source}, $linkFiles{target},$!));
	$self->record(sprintf("Moved %s to %s", $linkFiles{source}, $linkFiles{target}));
}

sub scp {
# Copies file from source to target using scp, logs actions.  Input can be two scalars,
# an array reference, or a hash reference
	my $self   = shift;
	my @params = @_;
		
	my %linkFiles = %{$self->resolveSourceAndTarget("scp", @params)};
	if(!defined $linkFiles{targetHost}) {
		$self->confess(sprintf("In Logger::scp, can't copy file over network with undefined hostname"));
	}
	
	my $command = sprintf("scp %s %s:%s", $linkFiles{source}, $linkFiles{targetHost}, $linkFiles{target});
	my $output  = $self->runCommand($command);
	$self->record(sprintf("Copied %s to %s on %s", $linkFiles{source}, $linkFiles{target}, $linkFiles{targetHost}));
}

sub smv {
# Moves file from source to target using scp, logs actions.  Input can be two scalars, an array reference,
# or a hash reference
	my $self   = shift;
	my @params = @_;

	my %linkFiles = %{$self->resolveSourceAndTarget("scp", @params)};
	if(!defined $linkFiles{targetHost}) {
		$self->confess(sprintf("In Logger::smv, can't copy file over network with undefined hostname"));
	}
	
	my $command = sprintf("scp %s %s:%s", $linkFiles{source}, $linkFiles{targetHost}, $linkFiles{target});
	my $output  = $self->runCommand($command);
	$self->rm($linkFiles{source});
	$self->record(sprintf("Moved %s to %s on %s", $linkFiles{source}, $linkFiles{target}, $linkFiles{targetHost}));
}

sub resolveSourceAndTarget {
# Resolve what source and target directories are for symlink, move, scp and smv.  Allows for two scalars (source and target)
# to be passed, an array reference, or a hash reference. In the future, may also allow for hash reference of 2 arrays and array of 2-member hashes
	my $self      = shift;
	my $calledBy  = shift;
	my @params    = @_;
	my $linkFiles = undef;
	

	if(ref($params[0]) eq "ARRAY")  {
	# If link files are passed as an array ref, first file is source file, second file is target file.
	# If num passed files != 2, confess
		if(scalar @{$params[0]} == 2) {
			$linkFiles->{source} = $params[0]->[0];
			$linkFiles->{target} = $params[0]->[1];
		}
		else {
			$self->confess(sprintf("Must specify a source and target file for Logger::%s", $calledBy));
		}
	}
	elsif(ref($params[0]) eq "HASH") {
	# If link files are passed as a hash ref, files are specified with 'source' and 'target' keys.
	# If either key is missing, confess.  NB: allows for additional keys, e.g. 'targetHost' in scp
	# and smv
		if(defined $params[0]->{source} and defined $params[0]->{target}) {
			$linkFiles = $params[0];
		}
		else {
			$self->confess(sprintf("In Logger::%s, hash provided with keys %s, please use target/source as hash keys",
															$calledBy,(join "/", keys(%{$params[0]}))));
		}
	}
	elsif(!ref($params[0]) and !ref($params[1])) {
	# If link files are passed as scalars, first is source file, second is target file.
	# If either file is undefined, confess
		if(defined $params[0] and defined $params[1]) {
			$linkFiles->{source} = $params[0];
			$linkFiles->{target} = $params[1];
		}
		else {
			$self->confess(sprintf("Argument to Logger:%s not passed as a hash or array, undefined values for source or target", $calledBy));
		}
	}
	else {
		$self->confess(sprintf("Did not recognize format of file parameters passed to Logger::%s", $calledBy));
	}
	return $linkFiles;
}

sub fname {
	my $self = shift;
	return $self->{_fname};
}

sub runCommand {
	# runs a command and records the output.  If an error is encountered,
	# it will be printed to the log file before ending the program in an
	# error state.
	# NB: sometimes, prog_out and prog_err do not seem to be deleted. Needs
	# to be resolved
	my $self           = shift;
	my $command        = shift;
	my $allowed_errors = shift;
	my $print_command  = shift;
	
	$print_command ||= 1;
	my $prog_name = @{[split / /, $command]}[0];
	$prog_name    = $1 if $prog_name =~ m/\/(.+)$/;
	# Removes any parent directories
	
	my $read_out = 0;
	my $read_err = 0;
	
	my $output = undef;
	my $timer  = Timer->new();
	
	# create filenames to redirect STDOUT and STDERR if necessary
	my $random   = int(rand(10000));
	my $prog_out = sprintf(".%d_%d.out",$$,$random);
	my $prog_err = sprintf(".%d_%d.err",$$,$random);
	
	if(($command !~ m/(1|12|21)>/ or $command =~ m/ *[^12]>/) and $command !~ m/>{1,2}/) {
	# checks to see if the command already redirects STDOUT
		$command .= " 1>>$prog_out";
		$read_out++;
	}
	if($command !~ m/(2|12|21)>/) {
	# checks to see if the command already redirects STDERR
		$command .= " 2>>$prog_err";
		$read_err++;
	}
	
	my $action = sprintf("Running command %s\n", $command);
	$self->record($action) if $print_command;
	
	system($command);
	
	if($read_out) {
		if(-s $prog_out) {
			open OUT, '<', $prog_out;
			$output->{stdout} = [<OUT>];
			close OUT;
		}
		else {
			$output->{stdout} = undef;
		}
		unlink $prog_out or $self->cluck("Could not remove $prog_out: $!");
	}
	
	if($read_err) {
		if(-s $prog_err) {
			open ERR, '<', $prog_err;
			$output->{stderr} = [<ERR>];
			close ERR;
		}
		unlink $prog_err or $self->cluck("Could not remove $prog_err: $!");
	}

	if($? != 0 and !isAllowedError($?, $allowed_errors)) {
		my $err_out = join "", @{$output->{stderr}} if defined $output->{stderr};

		Carp::carp(sprintf(">>Error running %s", $prog_name));
		printf STDERR ("Command attempted:\n\t%s\n", $command);
		printf STDERR ("Error message:\n\t%s\n", $err_out) if defined $err_out;
		printf STDERR ("Exit value %d\n",$?);
		$self->recordError($?,$command,$err_out);
		exit 1;
	}
		
	$self->recordElapsed($command, $timer);
	
	return $output;
}

sub isAllowedError {
	my $self           = shift;
	my $error          = shift;
	my $allowed_errors = shift;
	
	return 0 if !defined $allowed_errors->[0];
	foreach(@{$allowed_errors}) {
		return 1 if $_ == $error;
	}
	return 0;
}
