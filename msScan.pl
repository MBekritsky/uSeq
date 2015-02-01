#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
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

	print $config->{uSeqPerllib}, "\n";
	die "Could not parse config file: $@" if $@;
	die "Could not read config file: $!" 	unless defined $config;
	die "Could not run config file"				unless $config;
	unshift @INC, $config->{uSeqPerllib};
}

use MicroSeqFunc qw(createLocalMicroSeqDir getProgramPaths); #change to config file var
use Logger;
use Certificate;

sub findMicrosatellites;
sub alignReads;
sub getSamFile;
sub reindexReads;
sub splitReadsByChr;
sub changeFormatStrings;
sub convertToBam;

sub tetrascanParamComments;
sub filesToString;
sub makeTetrascanCmd;
sub makeAlnCmd;
sub makeSamCmd;
sub makeReindexCmd;
sub makeSplitCmd;
sub makeBamConvertCmd;
sub generateFileNames;

my $prog_paths = getProgramPaths();
if(defined $ENV{'JOB_ID'}) {
	$prog_paths = getProgramPaths($ENV{'USEQ_CONFIG'});
}

my $prog_name = "msScan.pl";
my $step_name = "msScan";

my $tetrascan_params;
my $gen_params;
my $files;

my $created;
my $action;

GetOptions
(
	"pair_format=s"				=> \$files->{pair_f},
	"solo_format=s"				=> \$files->{solo_f},
	"fastq=s"							=> \$files->{fastq_base},
	"dir=s"								=> \$files->{dir},
	"target_dir:s"				=> \$files->{target_dir},
	"paired"							=> \$files->{isPaired},
	"suffix=s"						=> \$files->{suffix},
	"zipped"							=> \$files->{zipped},
	
	"clean"								=> \$gen_params->{clean},
	"zip_output"					=> \$gen_params->{zip_output},
	"keep_all"						=> \$gen_params->{keep_all},
	"decentralized"				=> \$gen_params->{decentralized},
	"project_id:s"				=> \$gen_params->{project_id},
	
	"print_all_MS=i"			=> \$tetrascan_params->{M},
	"max_n=i"							=> \$tetrascan_params->{N},
	"min_unit=i"					=> \$tetrascan_params->{u},
	"min_length=i"				=> \$tetrascan_params->{l},
	"min_usize=i"					=> \$tetrascan_params->{m},
	"max_usize=i"					=> \$tetrascan_params->{x},
	"min_qual=i"					=> \$tetrascan_params->{Q},
	"barcodeLength=i"			=> \$tetrascan_params->{t},
	"qual_format=s"				=> \$tetrascan_params->{qual_format},
	"trim_type=s"					=> \$tetrascan_params->{trim_type},
);

my $hostname = `hostname`;
my $nodename = [split /\./,$hostname];
$nodename = $nodename->[0];

print "Running msScan on ",$hostname,"\n";

changeFormatStrings($files);
my $log = Logger->new($files->{target_dir},$prog_name,sprintf($files->{pair_f}, $files->{fastq_base})); #need explicit name here in case job is submitted to SGE
my $certifier = Certificate->new($files->{target_dir});

$files->{local_dir} = $files->{target_dir};
if($gen_params->{decentralized})
	{
		$files->{local_dir} = createLocalMicroSeqDir($gen_params->{project_id},$nodename,$log);
		print "Local directory: ",$files->{local_dir},"\n";
	}
	
#get filenames
generateFileNames($files,$gen_params,$log);
print "Log file for this run can be found at ", $log->fname(),"\n";
my ($tOptExplanations,$tChoiceExplanations) = tetrascanParamComments();
$log->recordParams($files,"File");
$log->recordParams($gen_params,"General");
$log->recordParams($tetrascan_params,"Tetrascan",$tOptExplanations,$tChoiceExplanations);

my $lane_name = sprintf($files->{pair_f}, $files->{fastq_base});
printf("msScan Processing %s\n", $lane_name);

findMicrosatellites($prog_paths->{uSeq},$tetrascan_params,$files,$log);
print STDERR "Completed microsatellite scan\n";
$log->rm($files->{ms_fastq});

$files->{target_counts_dir} = sprintf("%s/counts",$files->{target_dir});
$files->{local_counts_dir}	= sprintf("%s/counts",$files->{local_dir});

foreach (qw(1 2))
{
	my $fastq_stub=sprintf($files->{solo_f},$files->{fastq_base},$_);
	my $fullSourcePath = sprintf("%s/%s.scan.count",$files->{local_counts_dir},$fastq_stub);
	my $fullTargetPath = sprintf("%s/%s.scan.count",$files->{target_counts_dir},$fastq_stub);
	$log->mv($fullSourcePath,$fullTargetPath);
}

$certifier->completed($step_name,$lane_name);

sub changeFormatStrings
{
	my $params = shift;
	
	$params->{pair_f} =~ s/\@/%/g;
	$params->{solo_f} =~ s/\@/%/g;
}

sub generateFileNames
{
	my $files = shift;
	my $gen_params = shift;
	my $log = shift;
	my $file_prefix;
	my $file_stub;

	print "local dir is now ",$files->{local_dir},"\n";

	#Create names for output files
	if(! $files->{isPaired})
	{
		if(!$files->{zipped})
			{
				push @{$files->{fastq}}, sprintf("%s%s%s",					$files->{dir},$files->{fastq_base},$files->{suffix});
			}
		else
			{
				push @{$files->{fastq}}, sprintf("%s%s%s.gz",				$files->{dir},$files->{fastq_base},$files->{suffix});
			}
		push @{$files->{mod_fastq}}, sprintf("%s%s.mod.ms.txt", $files->{local_dir},$files->{fastq_base});
		push @{$files->{ms_fastq}},	 sprintf("%s%s.ms.txt",			$files->{local_dir},$files->{fastq_base});
		push @{$files->{hdr}}, 			 sprintf("%s%s.hdr.txt",		$files->{local_dir},$files->{fastq_base});
		push @{$files->{sai}}, 			 sprintf("%s%s.mod.sai",		$files->{local_dir},$files->{fastq_base});
	}
	else
	{
		for (my $i = 1; $i < 3; $i++)
		{
			$file_prefix = sprintf($files->{solo_f}, $files->{fastq_base},$i);
			if($files->{suffix} !~ m/\.bam$/)
				{
					if(!$files->{zipped})
						{
							push @{$files->{fastq}}, sprintf("%s%s%s",					$files->{dir},$file_prefix,$files->{suffix});
						}
					else
						{
							push @{$files->{fastq}}, sprintf("%s%s%s.gz",				$files->{dir},$file_prefix,$files->{suffix});
						}
				}
			else
				{
					$files->{isBam} = 1;
					$file_stub = sprintf($files->{pair_f}, $files->{fastq_base});
					$files->{fastq}->[0] = sprintf("%s%s%s",$files->{dir},$file_stub,$files->{suffix});
				}
			push @{$files->{mod_fastq}}, sprintf("%s%s.mod.ms.txt", $files->{local_dir},$file_prefix);
			push @{$files->{ms_fastq}},	 sprintf("%s%s.ms.txt",			$files->{local_dir},$file_prefix);
			push @{$files->{hdr}}, 			 sprintf("%s%s.hdr.txt",	  $files->{local_dir},$file_prefix);
			push @{$files->{sai}}, 			 sprintf("%s%s.mod.sai",	  $files->{local_dir},$file_prefix);
		}
	}
	
	for(my $i = 0; $i < scalar @{$files->{fastq}}; $i++)
		{
			if(-l $files->{fastq}->[$i])
				{
					my $file_path = [split /\//, $files->{fastq}->[$i]];
					my $new_file_name = sprintf("%s/%s",$files->{local_dir},$file_path->[$#{$file_path}]);
					if(! (-e $new_file_name and readlink($new_file_name) eq readlink($files->{fastq}->[$i])))
					{
						print "linking ", readlink($files->{fastq}->[$i])," to $new_file_name\n";
						$log->symlink({source => readlink($files->{fastq}->[$i]),target => $new_file_name, overwrite => 1});
					}
					$files->{fastq}->[$i] = $new_file_name;
				}
			#if the file is a symbolic link in original dir, link the original file to the local directory
		}

	$file_prefix = sprintf($files->{pair_f}, $files->{fastq_base});
	$files->{sam} = sprintf("%s%s.mod.sam",$files->{local_dir},$file_prefix);
	$files->{msi} = sprintf("%s%s.msi.sam",$files->{local_dir},$file_prefix);
}

sub  tetrascanParamComments
{
	#Explanations for different option names
	my $optExplanations = {
		"M" => "print only reads with microsatellites",
		"N" => "max number of Ns in read",
		"u" => "min repeats of microsatellite motif",
		"l" => "min microsatellite tract length",
		"m" => "min microsatellite motif length",
		"x" => "max microsatellite motif length",
		"Q" => "min basecall qual (independent of Phred score)",
		"t" => "length of barcode sequence at the beginning of every read"
	};

	#Explanations for different option choices (e.g. which qual format was specified)
	my $choiceExplanations = {
		"a" => "Sanger format (Phred + 33)",
		"s" => "Illumina format (Phred + 64)",
		"B" => "BWA-style trimming",
		"g" => "trim to last base pair greater than Q",
		"b" => "trim to first base pair less than Q",
		"n" => "do not trim",
	};

	return ($optExplanations,$choiceExplanations);
}

sub findMicrosatellites
{
	my $uSeq_path = shift;
	my $tetrascan_params = shift;
	my $files = shift;
	my $log = shift;
	
	my $tetrascan_cmd = makeTetrascanCmd($uSeq_path,$tetrascan_params,$files);
	my $tetrascan_out = $log->runCommand($tetrascan_cmd);
	printf("%s tetrascan\n",$files->{fastq_base});
	printf("\tstdout:\n%s\n",join "\t", @{$tetrascan_out->{stdout}}) if defined $tetrascan_out->{stdout};
	printf("\tstderr:\n%s\n",join "\t", @{$tetrascan_out->{stderr}}) if defined $tetrascan_out->{stderr};
}

sub makeTetrascanCmd
{
	my $path = shift;
	my $params = shift;
	my $files = shift;

	my $command = undef;


	$command = "bash -c \'"	if $files->{zipped};

	$command .= sprintf("%s/tetrascan ", $path);

	foreach(keys %{$params})
	{
		if($_ !~ m/(M|qual_format|trim_type)/)
		{
			$command .= sprintf(" -%s %d", $_, $params->{$_});
		}
	}

	$command .= " -M" if $params->{M} == 0;
	$command .= sprintf(" -%s", $params->{$_}) foreach(qw(qual_format trim_type));
	if(!$files->{isBam})
		{
			if(!$files->{zipped})
				{
					$command .= sprintf(" -q %s", $files->{fastq}->[0]);
					$command .= sprintf(" -p %s", $files->{fastq}->[1]) if $files->{isPaired};
				}
			else
				{
					$command .= sprintf(" -q <(gunzip -c %s)", $files->{fastq}->[0]);
					if($files->{isPaired})
						{
							#this will not work for both file naming conventions
							my $output_prefix = sprintf($files->{pair_f},$files->{fastq_base});
							$output_prefix = $files->{local_dir}.$output_prefix;
							$command .= sprintf(" -p <(gunzip -c %s) -O %s", $files->{fastq}->[1],$output_prefix);
						}
					#Need to come up with alternative if files are unpaired
					$command .= "\'";
				}
		}
	else
		{	
			my $output_prefix = sprintf($files->{pair_f},$files->{fastq_base});
			$output_prefix = $files->{local_dir}.$output_prefix;
			$command .= sprintf(" -b %s -O %s",$files->{fastq}->[0], $output_prefix);
		}
	return $command;
}
