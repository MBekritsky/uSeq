package Timer;

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Basename;
use Carp qw(cluck);

my $VERSION = 1.05;

sub new {
	my $class = shift;
	my $initTime = shift;
	
	my $self = {
		_time => undef,
		_timeHash => undef,
	};
	
	bless $self, $class;
	if(defined $initTime) {
		$self->set($initTime);
	}
	else {
		$self->start();
	}
	return $self;
}

sub start {
	my $self = shift;
	$self->{_time} = time;
	$self->{_timeHash} = $self->timeHash(1);
}

sub restart {
	my $self = shift;
	$self->start();
}

sub set {
	my $self = shift;
	my $time = shift;
	$self->{_time} = shift;
	$self->{_timeHash} = $self->timeHash(1);
}

sub timeHash {
	my $self = shift;
	my $rehash = shift;
	$rehash ||= 0;
	#Tried to maintain unix date conventions for the hash keys
	#S = seconds (0-60), M = minutes (0-59), H = hours (0-23)
	#d = date of month, m = month (0-11), Y = Year (from 1900)
	#a = day of week (0-6), j = day of year, i = is DST
	if(!defined $self->{_timeHash} or $rehash == 1) {
		my $localtime_keys = ["S","M","H","d","m","Y","a","j","i"];
		my $time_hash;
		@{$time_hash}{@{$localtime_keys}} = @{[localtime($self->{_time})]};
		$self->{_timeHash} = $time_hash;
	}
	else {
		return $self->{_timeHash};
	}
}

sub time {
	my $self = shift;
	if(defined $self->{_time}) {
		return $self->{_time};
	}
	else {
		cluck("Cannot return undefined time");
	}
}

#returns a timestamp in the format [DD/MM/YY HH:MM:SS]
sub stamp {
	my $self = shift;
	my $current = shift;
	
	$self->start() if(defined $current and $current == 1);
	my $timeHash = $self->{_timeHash};
	
	my $stamp = sprintf("[%02d/%02d/%02d %02d:%02d:%02d] ", $timeHash->{d}, $timeHash->{m} + 1, $timeHash->{Y} + 1900, 
																																$timeHash->{H}, $timeHash->{M},$timeHash->{S});
	return $stamp;
}

#returns a string compatible with file naming (no spaces or punctuation) in the format DDMMYYYY_HHMMSS
#current will return the file string for the current time, otherwise, fstring will return the string for the time when the timer was started
sub fstring {
	my $self = shift;
	my $current = shift;
	$self->start() if(defined $current and $current == 1);
	
	my $timeHash = $self->{_timeHash};
	my $string = sprintf("%02d%02d%d_%02d%02d%02d", $timeHash->{d}, $timeHash->{m} + 1, $timeHash->{Y} + 1900, 
																										$timeHash->{H}, $timeHash->{M},$timeHash->{S});
	return $string;
}

#returns a readable time string in the format WWW DD/MM/YYYY HH:MM:SS
sub string {
	my $self = shift;
	my $current = shift;
	$self->start() if defined $current and $current == 1;

	my %wmap = (0 => "Sun", 1 => "Mon", 2 => "Tue", 3 => "Wed", 4 => "Thu", 5 => "Fri", 6 => "Sat");
	
	my $timeHash = $self->{_timeHash};
	
	my $string = sprintf("%s %02d/%02d/%d %02d:%02d:%02d",$wmap{$timeHash->{a}},$timeHash->{d}, $timeHash->{m} + 1, $timeHash->{Y} + 1900,
																													$timeHash->{H}, $timeHash->{M},$timeHash->{S});
	return $string;
}
				 

sub elapsed {
	my $self = shift;
	my $stop = Timer->new();
	my $elapsed = $stop->time() - $self->time();
	
	my $eHash;
	my $tElapsed = $elapsed;
	$eHash->{hour} = int($tElapsed / 3600);
	$tElapsed = $tElapsed - ($eHash->{hour} * 3600);
	$eHash->{minute} = int($tElapsed / 60);
	$tElapsed = $tElapsed - ($eHash->{minute} * 60);
	$eHash->{second} = $tElapsed;
	
	my $eString = sprintf("%02d:%02d:%02d",$eHash->{hour},$eHash->{minute},$eHash->{second});
	return $eString;
}

=head1 NAME

Timer - a module for keeping and reporting time, designed for use with the uSeq pipeline

=head1 SYNOPSIS

 use Timer.pm;

 # create a new Timer object
 my $timer = Timer->new();

 # restart the Timer object
 $timer->start();
 $timer->restart();
	
 # set a new time for the timer object
 $timer->set(time);
	
 # get the time of the timer object in epoch format
 my $time      = $timer->time();
	
 # get the time of the timer object as a hash with keys corresponding to Unix date conventions
 my %timeHash  = $timer->timeHash();
	
 # get a timestamp in the format [DD/MM/YY HH:MM:SS]
 my $timeStamp = $timer->stamp();
	
 # get a file identifier using the current time in the format DDMMYYYY_HHMMSS
 my $fstring   = $timer->fstring();
	
 # get a human readable time string in the format WWW DD/MM/YYYY HH:MM:SS
 my $string    = $timer->string();
	
 # get the time elapsed since the timer began in the format HH:MM:SS
 my $elapsed   = $timer->elapsed();

=head1 DESCRIPTION

Creates a Timer object that can be used to create timestamps for log files;
create file identifiers based on the current date and time; print strings
reporting the time; and report the time elapsed while executing a command

The time kept internally by the object uses Perl's time specification, which is
the number of non-leap seconds since the epoch, which can be system-dependent

=head1 METHODS

=over 4

=item new

=item new($time)

Creates and returns a new Timer object.  If provided with $time, the timer will
start with the specified time, otherwise Timer will use the Perl time command to
get the current time at initialization

=item start

=item restart

Resets the timer's start time to the time when the function is called, using the Perl
time command

=item set time

Sets the timer's start time to the time provided by time, must be in epoch format

=item time

Returns the timer's current start time

=item timeHash

=item timeHash($rehash)

Returns a hash for the timer's current start time.  The hash has keys corresponding
to Unix date conventions: S (second), M (minute), H (hour), d (day of month), m (month),
Y (year from 1900), a (day of week), j (day of year), i (is DST).  Please note that
the month and day of week are zero-indexed.

The optional value $rehash is a boolean that tells the Timer class to rehash the time,
which is especially important if the timer has been restarted or otherwise reset.  It
is primarily for use by other timer functions, which will call it as needed.

=item stamp

=item stamp($current)

Returns a timestamp string in the form [DD/MM/YY HH:MM:SS].  $current is an optional value described
below.

=item fstring

=item fstring($current)

Returns a file identifier string in the format DDMMYYYY_HHMMSS.  
$current is an optional value described below.

=item string

=item string($current)

Returns a human-readable time string with in the format WWW DD/MM/YYYY HH:MM:SS.  $current is an 
optional value described below.

=item $current

A boolean value provided to stamp, fstring, or string.  When $current is true (1) $current
will return a string giving the current time.  If current is not specified or set to false
(0), it will return a string using the time stored in the Timer object.

=item elapsed

Returns the time elapsed since the timer was last started (either with new, start, restart, or set)

=back

=head1 REQUIRES

Cwd, File::Basename, Carp

=head1 AUTHOR

Mitchell Bekritsky (mitchell.bekritsky@gmail.com)
