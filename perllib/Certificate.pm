package Certificate;

use strict;
use warnings;
use Carp qw(confess cluck);
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
my $VERSION = 1.25;

sub new {
	my $class     = shift;
	my $targetDir = shift;
	
	my $certDir = sprintf("%s/certificates", $targetDir);
	
	my $self = {
		_certDir => $certDir,
		_timer	 => Timer->new(),
	};
	
	if(! -d $certDir) {
		# create a new certificate directory if it doesn't already exist
		if(!(mkdir $certDir)) {
			if(-d $certDir) {
				# sometimes two programs will simultaneously try to create a 
				# certificate directory simultaneously if one creates the 
				# directory before the other, then this will just print a 
				# message saying that the directory has already been created
				# ALTERNATIVE: have a front-end program that creates all 
				# directories that will be used by the pipeline
				Carp::cluck("Certificate directory was created by another program")	
			}
			else {
				Carp::confess("Failed to create ceritificate directory $certDir: $!");
			}
		}
	}
	
	bless $self, $class;
	return $self;
}

sub completed {
	# assign a completion certificate
	my $self   = shift;
	my $step   = shift;
	my $suffix = shift;
	$self->certify($step, $suffix, "complete");
}

sub created {
	# assign a creation certificate 
	my $self   = shift;
	my $step   = shift;
	my $suffix = shift;
	$self->certify($step, $suffix, "create");
}

sub certify {
	# assign an arbitrary certificate
	my $self   					 = shift;
	my $step   					 = shift;
	my $suffix 					 = shift;
	my $certified_action = shift;
	
	my $timestring       = $self->{_timer}->fstring(1);
	
	my $certificate_name = sprintf("%s/%s_%s_%s", $self->{_certDir}, $step,
																	$timestring, $certified_action);
	$certificate_name .= "_$suffix" if defined $suffix;
	
	open CERT, '>', $certificate_name or die "Could not create certificate $certificate_name: $!";
	close CERT;
}

sub getCertDir {
	my $self = shift;
	return $self->{_certDir};
}

sub certificateTimeHash {
	# get the time hash for another certificate
	my $self             = shift;
	my $certificate_name = shift;
	
	my $cert_details = [split /_/, $certificate_name];
	my $cert_time_hash;
	$cert_time_hash->{d} = substr($cert_details->[1], 0, 2);
	$cert_time_hash->{m} = substr($cert_details->[1], 2, 2) - 1;
	$cert_time_hash->{Y} = substr($cert_details->[1], 4) - 1900;
	$cert_time_hash->{H} = substr($cert_details->[2], 0, 2);
	$cert_time_hash->{M} = substr($cert_details->[2], 2, 2);
	$cert_time_hash->{S} = substr($cert_details->[2], 4);
	#a,j,i (day of week, year and if DST) are not defined here

	return $cert_time_hash;
}

sub compareTimes {
	#returns true if certificate time is more recent than the time provided in time_hash
	my $self      = shift;
	my $certFile  = shift; #This is wrong, but a quick fix for now
	my $time_hash = shift;
	
	if(ref($time_hash) ne "HASH" and $time_hash =~ m/^\d+$/) {
			Carp::carp("Passed time is not in time_hash format, attempting to convert it now");
			$time_hash = getTimeHash($time_hash);
	}
	my $cert_time_hash = $self->certificateTimeHash($certFile);
	
	foreach(qw(Y m d H M S)) {
			return 1 if $cert_time_hash->{$_} > $time_hash->{$_};
			return 0 if $cert_time_hash->{$_} < $time_hash->{$_};
	}
	return 1;
}

=head1 NAME

Certificate - a module for certifying pipeline steps, designed for use with the uSeq pipeline

=head1 SYNOPSIS

 use Certificate.pm;

 # create a new Certificate object
 my $certifier = Certificate->new($targetDir);

 # certify the creation of a file
 $certifier->created($step);
 $certifier->created($step, $suffix);
	
 # certify the completion of a pipeline step
 $certifier->completed($step);
 $certifier->completed($step, $suffix);
	
 # certify an arbitrary action
 $certifier->certify($step, undef, $certifiedAction);
 $certifier->certify($step, $suffix, $certifiedAction);

 # get the directory where certificates are written to
 my $certDir = $certifier->getCertDir();

 # extract a time hash from a certificate
 my $certificateTimeHash = $certifier->certificateTimeHash($certificateName);

 # determine whether a certificate was created before or after a certain time
 my $isMoreRecent = $certifier->compareTimes($certificateFile, $timeHash);

=head1 DESCRIPTION

Creates a Certificate object that can be used to certify when a file has been
created or when a stage of a pipeline has been completed.  This is particularly
useful when determining whether a pipeline component has successfully completed,
or when trying to assess whether a file has been created before or after the
initiation of the current pipeline instance.

Certificates are empty files with a name that describes when an event occurred.
The typical certificate filename is structured as {pipeline component calling certifier}_{date/time}_{event being certified}

Certificate uses Timer.pm to get the times for its timestamps.

=head1 METHODS

=over 4

=item new($targetDir)

Creates and returns a new Certificate instance.  $targetDir specifies the desired
directory to write certificates to. The certificate directory will be $targetDir/certificates

=item getCertDir

Return a string specifying the certificate directory.

=item created($step)

=item created($step, $suffix)

Certify the creation of $step.  $suffix is an optional value
described below.

=item completed($step)

=item completed($step, $suffix)

Certify the completion of $step.  $suffix is an
optional value described below.

=item certify($step, undef, $certifiedAction)

=item certify($step, $suffix, $certifiedAction)

Certify an arbitrary event for $step.  $certifiedAction describes the action being certified.
$suffix is an optional value described below.

=item certificateTimeHash($certificateName)

Extract the time hash from a certificate

=item compareTimes($certificateFile, $timeHash)

Compare the time specified in the name of $certificateFile to the time provided in $timeHash.
$timeHash has the structure outlined in the documentation of Timer.pm.  This function is actually
designed to use a time hash provided by Timer.pm.

=back

=head1 REQUIRES

Cwd, File::Basename, Carp, Timer

=head1 SEE ALSO

Timer.pm

=head1 AUTHOR

Mitchell Bekritsky (mitchell.bekritsky@gmail.com)
