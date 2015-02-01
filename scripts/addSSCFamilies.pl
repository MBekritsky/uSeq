#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my $targetDir;
my $inReportFile;

GetOptions(
 "report=s"  => \$inReportFile,
 "dir=s"     => \$targetDir
);
 
$targetDir ||= ".";
 
open(my $RFH, "<", $inReportFile) or die "Can't open $inReportFile: $!";
 
my $outReportFile = sprintf("%s/addReport_%s.txt",$targetDir,
 							  strftime("%d%m%Y_%H%M%S", localtime));

open(my $OFH, ">", $outReportFile) or die "Could not create $outReportFile: $!";

my $numLines = 0;
my @fields;
my ($sampleID, $familyID, $relation);
my (@seqInfo, $flowcell, $lane, $barcode);
my ($filepath, $linkPrefix, $linkname, $localdir, $relFile);
my $familydir;
my $numFiles;
my @seqFiles;
my $sourceFile;
my $targetFile;
my $oldSourceFile;

my %personInfo;

LINE:while(my $line = <$RFH>) {
	$numLines++;
	next LINE if $numLines == 1;
	
	chomp $line;
	@fields = split("\t",$line);
	
	$familyID = $fields[25];
	$relation = $fields[27];
	$sampleID = $fields[30];
	$barcode  = $fields[2];
	$filepath = $fields[12];

	@seqInfo  = split("-",$fields[1]);
	$flowcell = $seqInfo[0];
	$lane     = $seqInfo[1];
	
	if($sampleID !~ m/^SSC\d+$/ or $familyID !~ m/^auSSC\d+$/) {
		print $OFH "Skipping $sampleID in $familyID\n";
		next LINE;
	} 
	
	$linkPrefix = sprintf("%s%s",$flowcell,$lane);
	
	
	$localdir = sprintf("%s/%s/%s",$targetDir,$familyID,$sampleID);
	$linkname = sprintf("%s/%s.bam",$localdir,$linkPrefix);
	
	$personInfo{$sampleID}->{numFiles}++;
	$personInfo{$sampleID}->{newFile} = 0 if !defined $personInfo{$sampleID}->{newFile};
	$personInfo{$sampleID}->{dir} = $localdir if !defined $personInfo{$sampleID}->{localdir};

	if (-d $localdir) {
		$personInfo{$sampleID}->{newDir} = 0;
	}
	else {
		$personInfo{$sampleID}->{newDir} = 1;
		$familydir = sprintf("%s/%s",$targetDir,$familyID);
		if(! -d $familydir){
			mkdir($familydir) or die "Could not create $familydir: $!";
			print $OFH "Created $familydir\n";
		}
		mkdir($localdir) or die "Could not create $localdir: $!";
		print $OFH "Created $localdir\n"
	}
	
	$relFile = sprintf("%s/relationship.txt",$localdir);
	if(! -f $relFile) {
		open(my $REL, ">", $relFile) or die "Could not create $relFile: $!";
		print $REL $relation,"\n";
		close $REL;
	}
	
	$filepath = sprintf("%s/bc%s",$filepath,$barcode);
	if($filepath !~ m/^\/mnt.*/) {
		print $OFH "Skipping invalid file path $filepath for $sampleID in $familyID\n";
		next LINE;
	}
	
	my $DIR;
	if(!opendir($DIR, $filepath))
	{
	 print $OFH "Skipping...Could not open $filepath: $!\n";
	 next LINE;
	}
	@seqFiles = grep {/\.fastq\.gz/ || /read\.bam/} readdir $DIR;
	$numFiles = length(@seqFiles);
	closedir($DIR);
	
	if($numFiles > 1) {
		print $OFH "More than one sequence data file in $filepath, skipping...\n";
		next LINE;
	}
	
	foreach my $file (@seqFiles) {
		$sourceFile = sprintf("%s/%s",$filepath, $file);
		$targetFile = $linkname;
		if(! -l $linkname) {
			symlink($sourceFile,$targetFile) or die "Could not link $sourceFile to $targetFile: $!";
			print $OFH "Linked $sourceFile to $targetFile\n";
			$personInfo{$sampleID}->{newFile}++;
		}
		else {
			$oldSourceFile = readlink($linkname) or die "Could not read link to $linkname: $!";
			if(-l $targetFile and $oldSourceFile ne $sourceFile)
			{
				unlink $targetFile or die "Could not remove out of date link $targetFile: $!";
				symlink($sourceFile,$targetFile) or die "Could not link $sourceFile to $targetFile: $!";
				print $OFH "Updated $sourceFile from $oldSourceFile for $targetFile\n";
				$personInfo{$sampleID}->{newFile}++;
			}
			else
			{
				print $OFH "$targetFile is up to date\n";
			}
		}
	}
}
close($OFH);
close($RFH);

my $runFile = sprintf("%s/newPeopleToRun_%s.txt",$targetDir,
 							  strftime("%d%m%Y_%H%M%S", localtime));
open(my $FH, ">", $runFile) or die "Could not open $runFile: $!";

foreach my $key(keys %personInfo) {
	if($personInfo{$key}->{newDir} > 0 or $personInfo{$key}->{newFile} > 0) {
		print $FH $personInfo{$key}->{dir},"\t","scan\n";
	}
}
close $FH;

print "Finished adding new links and directories for SSC families from $inReportFile.\n";
print "Report can be found at $outReportFile\n";
print "A list of people to be processed by uSeq can be found at $runFile\n";
 
