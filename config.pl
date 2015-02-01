#!/usr/bin/perl

# Sets up a config file for uSeq
# the config file should contain paths
# to uSeq's perllib, BWA, SamTools, and
# Picard
# This script can be run to auto-detect
# paths, or it can be supplied paths
# If it is being run on auto-detect, it
# must be run from the same directory as
# the uSeq source code

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Basename;
use File::Which;
use Sys::Hostname;
use Getopt::Long;

my $paths;
my $bwaPath;
my $picardPath;
my $samtoolsPath;
my $help;

GetOptions(
	"useq:s"			=> \$paths->{uSeq},
	"perllib:s"		=> \$paths->{uSeqPerllib},
	"bwa:s"				=> \$paths->{bwa},
	"picard:s"		=> \$paths->{picard},
	"samtools:s"	=> \$paths->{samtools},
	"help|?"			=> \$help
);
printHelp() if $help;

my $currPath = dirname(abs_path($0));
my $configFile  = sprintf("%s/config.txt", $currPath);

my $currHost = [split(/\./, hostname)];
$currHost = $currHost->[0];


if(-e $configFile) {
	print "$configFile already exists.  Overwrite(y/n)? ";
	
	while(defined(my $response = <>)){
		last if(lc(substr($response,0,1)) eq "y");
		exit if(lc(substr($response,0,1)) eq "n");
		print "Unrecognized response, please enter \"y\" or \"n\: ";
	}
	
}

if(!defined $paths->{uSeqPerllib}) {
	$paths->{uSeqPerllib} = sprintf("/mnt/%s/%s/perllib/", $currHost,
													$currPath);
}
	
if(!defined $paths->{uSeq}) {
	$paths->{uSeq} = sprintf("/mnt/%s/%s/", $currHost,
													$currPath);
}
	
	
foreach my $p (keys %{$paths}){
	next if $p eq "uSeqPerllib" or $p eq "uSeq";
	$paths->{$p} = findProgram($p) if !defined $paths->{$p};
}

foreach my $p (keys %{$paths}){
	testPathExists($paths->{$p});
}

$paths->{config} = sprintf("/mnt/%s/%s", $currHost, $configFile);

open(my $cfh, ">", $configFile) or die "Could not create $configFile: $!\n";
print $cfh "{\n";
foreach my $p (keys %{$paths}){
	print $cfh "$p\t=>\t\"",$paths->{$p},"\",\n" ;
}
print $cfh "}\n";
close($cfh);


sub findProgram {
	my $programName = shift;
	
	my $programPath = which($programName);
	if(!defined $programPath){
		print STDERR "$programName could not be found in the system path. ";
		print STDERR "Please ensure $programName is installed and in your\n";
		print STDERR "system's path.  If $programName is installed and not ";
		print STDERR "in your system's path, you may specify its\n";
		print STDERR "location using the --$programName flag\n";
		exit;
	}
	return $programPath;
}

sub testPathExists {
	my $path = shift;

	if(! -d $path){
		die "The path $path could not be found\n";
	}
}

sub printHelp {
	print STDERR "Sets up a config file for uSeq with paths to uSeq and its custom";
	print STDERR " perlllib directory, as well as paths to BWA, SamTools, and\n";
	print STDERR "Picard.  The config file will try to find BWA, SamTools, and";
	print STDERR " Picard in your system path.  Make sure these programs are\n";
	print STDERR "installed before running config.pl.  If any of these programs";
	print STDERR " are installed on your cluster, but are not in the system path,\n";
	print STDERR "you may specify them with --bwa, --samtools, or --picard.\n\n";
	print STDERR "PLEASE NOTE: you must specify the paths to these executables, not";
	print STDERR " the files themselves.  For instance, if bwa was found at /path/to/bwa/bwa\n";
	print STDERR "the proper argument for the --bwa flag would be /path/to/bwa/.\n";
	exit;
}
