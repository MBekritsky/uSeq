package DirOps;

use strict;
use warnings;
use Exporter;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION = 1.05;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(checkDirExist checkTrailingSlash);
%EXPORT_TAGS = ( All 		 => [@EXPORT_OK],
								 all 		 => [@EXPORT_OK],
								 DEFAULT => [qw(checkDirExist)],
							 );

sub checkDirExist {
	my $dir     = shift;
	my $dirDesc = shift;
	
	my $dieString = "Error!  The specified ";
	$dieString .= "$dirDesc " if defined $dirDesc;
	$dieString .= "directory, $dir, does not exist\n";
		
	die "Error!  The specified $dirDesc directory, $dir, does not exist\n" if ! -d $dir;
}

sub checkTrailingSlash {
	my $dir = shift;
	
	$dir =~ m/\/$/ ? return $dir : return $dir.'/';
}

=head1 NAME

DirOps - a module for handling directory operations for the uSeq pipeline

=head1 SYNOPSIS

 use DirOps.pm;
 
 # check if a directory exists
 checkDirExists($dir);
 checkDirExists($dir, $dirDesc);
 
 # add a trailing slash to a directory name
 # if it's not already there
 checkTrailingSlash($dir);
 
=head1 DESCRIPTION

A collection of functions for managing directories.  By default, the only
function that is exported is checkDirExist.

=head1 METHODS

=over 4

=item checkDirExists($dir)

=item checkDirExists($dir, $dirDesc)

Checks to see if a directory exists.  If it does not, prints a message to STDERR
and calls 'die'.  If provided with the optional $dirDesc, a directory description, 
it will print a very short description of the directory that could not be found.

=item checkTrailingSlash($dir)

Checks to see if a directory name has a trailing slash.  If it does not have a
trailing slash, then it adds one.  In either case, it will return the directory name.

=head1 REQUIRES

None

=head1 AUTHOR

Mitchell Bekritsky (mitchell.bekritsky@gmail.com)
