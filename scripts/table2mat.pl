#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $tableFile;
my $locusFile;
my $minCount;
my $minPopPct;
my $meanPopCov;

GetOptions(
	"counts=s" => \$tableFile,
	"index=s"  => \$locusFile,
	"minCount:i" => \$minCount,
	"minPopPct:f" => \$minPopPct,
	"meanPopCov:f" => \$meanPopCov,
);

$minCount ||= 10;
$minPopPct ||= 0.5;
$meanPopCov ||= 0;

open LOC,'<',$locusFile or die "Could not open $locusFile: $!";
my $loci = [<LOC>];
close LOC;

my $locusInfo;
foreach(@{$loci})
{
	chomp;
	push @{$locusInfo},[split '\t', $_];
}

my ($popInd,$topInd,$sumInd) = (4,5,6);

my $splitPath = [split /\//, $tableFile];
my $numDirs = scalar @{$splitPath};
my $parentDir = join "/", @{$splitPath}[0..($numDirs - 2)];

if($numDirs == 1)
{
	$parentDir = ".";
}

my $alleleFile = $parentDir."/allele_matrix_info.txt";
open ALLELE,'>', $alleleFile or warn "Could not open $alleleFile: $!\n";
my $matrixFile = $parentDir."/allele_matrix.txt";
open MATRIX,'>', $matrixFile or warn "Could not open $matrixFile: $!\n";

print "Writing count matrix to $matrixFile and matrix row index to $alleleFile\n";

my $lineNum = 0;

open TAB,'<', $tableFile or die "Could not open $tableFile: $!";
while(my $line = <TAB>)
{
	chomp $line;
	my $info = [split '\t', $line];
	my $numPeople = scalar @{$info};
	my $minPop = $minPopPct * $numPeople;
	
	if($locusInfo->[$lineNum]->[$topInd] >= $minCount and $locusInfo->[$lineNum]->[$popInd] >= $minPop and 
		$locusInfo->[$lineNum]->[$sumInd]/$locusInfo->[$lineNum]->[$popInd] >= $meanPopCov)
	{
	  my $locusCounts;
	  for(my $h = 0; $h < $numPeople; $h++)
	  {
		  $locusCounts->{0}->[$h] = 0;
	  }

	  my $alleleCounts;
      my $numAlleles = 0;

	  PERSON:for(my $i = 0; $i < scalar @{$info}; $i++)
	  {
		  if($info->[$i] eq '0;;')
		  {
			  $locusCounts->{0}->[$i] = 0;
			  next PERSON;
		  }
		  my $personInfo = [split ';', $info->[$i]];
		  my $numModes = $personInfo->[0];
		  my $modes = [split ',', $personInfo->[1]];
		  my $counts = [split ',', $personInfo->[2]];
		  my $totalCount = 0;
		  for(my $j = 0; $j < $numModes; $j++)
		  {
			  $totalCount += $counts->[$j];
			  if(!defined $locusCounts->{$modes->[$j]})
			  {
                  $numAlleles++;
				  for(my $h = 0; $h < $numPeople; $h++)
				  {
					  $locusCounts->{$modes->[$j]}->[$h] = 0;
					  $alleleCounts->{$modes->[$j]}->{'top'} = $counts->[$j];
				  }
			  }
			  $locusCounts->{$modes->[$j]}->[$i] = $counts->[$j];
			  $alleleCounts->{$modes->[$j]}->{'pop'}++;
			  $alleleCounts->{$modes->[$j]}->{'sum'}+= $counts->[$j];
			  $alleleCounts->{$modes->[$j]}->{'top'} = $counts->[$j] if $counts->[$j] > $alleleCounts->{$modes->[$j]}->{'top'};

		  }
		  $locusCounts->{0}->[$i] = $totalCount;
	  }
		
	  foreach(sort {$a <=> $b} keys %{$locusCounts})
	  {
		  if($_ == 0)
		  {
			  print ALLELE join "\t", @{$locusInfo->[$lineNum]};
			  print ALLELE "\t$_\t$numAlleles\n";
		  }
		  else
		  {
			  print ALLELE join "\t", @{$locusInfo->[$lineNum]}[0..3];
			  print ALLELE "\t",$alleleCounts->{$_}->{'pop'};
			  print ALLELE "\t",$alleleCounts->{$_}->{'top'};
			  print ALLELE "\t",$alleleCounts->{$_}->{'sum'};
			  print ALLELE "\t$_\t$numAlleles\n";
		  }
		  print MATRIX join ",", @{$locusCounts->{$_}};
		  print MATRIX "\n";
	  }
	}
	$lineNum++;
}

close TAB;
close ALLELE;
close MATRIX;
