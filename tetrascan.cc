/*
tetrascan.cc

ABSTRACT:
A program for scanning through sequence data
in fasta, fastq, or BAM formats and 
extracting microsatellite sequences.

CREATION DATE:
20.07.2010

LAST REVISION:
25.10.2011

AUTHOR:
Mitchell Bekritsky
*/

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <errno.h>

#include "api/BamReader.h"

#include "Timer.h"
#include "FileOps.h"
#include "Microsatellite.h"
#include "MicrosatelliteDetector.h"
#include "FastaFile.h"
#include "FastqFile.h"
#include "Counter.h"

#include "uScanFE.h"

using namespace std;
using namespace BamTools;

//Constants
const int BUFFER_SIZE = 100;
//Output file suffixes 
const string FASTA_OUTPUT_SUFFIX = ".ms.fa", MODIFIED_FASTA_SUFFIX = ".mod.fa", OFFSET_FILE_SUFFIX = ".osi";
const string FASTQ_OUTPUT_SUFFIX = ".ms.txt",MODIFIED_FASTQ_SUFFIX = ".mod.ms.txt", HDR_SUFFIX = ".hdr.txt";
const string COUNTER_SUFFIX = ".scan.count";
unsigned int barcodeLength = 0;

int main (int argc, char * argv [])
{
  	
  stringstream helptext;
  helptext << 
    "\n"
    "tetrascan: Find perfect microsatellite sequences\n"
    "\n"
    "Synopsis: tetrascan [-f <fasta file> -q <fastq file PE1> -p <fastq file PE2, optional> -b <bam file>]...\n"
    "\n"
    "Options:\n"
    "    General Options:\n"
    "        -f <fasta>\t\tFASTA file to scan (can contain multiple headers)\n"
    "        -q <fastq>\t\tFASTQ file to scan\n"
		"        -p <fastq>\t\tmate pair FASTQ file to scan (use when you want to print out reads with microsatellites and their mate pairs\n"
		"        -b <BAM file>\t\tBAM file to scan\n"
    "        -k <nt>\t\t\tflanking nucleotides to print with reported microsatellites from FASTA files (default: " << flankLength << ")\n"
    "        -o \t\t\tredirect output to stdout\n"
		"        -T <pathname>\t\tspecifies an alternate path to write tetrascan output files to (default is to write to same path as the input file)\n"
		"        -O <file prefix>\twith output flags, overrides default naming (can specify path in prefix or with -T, this option does not work with -p)\n"
		"        -M \t\t\tonly print reads with microsatellites in them (default is to print all reads)\n"
		"        -N \t\t\tset a maximum undefined nucleotide threshold (i.e. Ns) for sequence reads (default: " << nThres << ")\n"
		"        -d \t\t\t do not print a modified sequence file with microsatellites removed\n"
    "        -h/-H\t\t\tDisplay this page\n"
    "    Detector Options:\n"
    "        -u <units>\t\tminimum number of units in reported microsatellites (default: " << minUnits << ")\n"
    "        -l <length>\t\tminimum nt length in reported microsatellites (default: " << minLength << ")\n"
    "        -m <len>\t\tminimum repeat unit length of reported microsatellites (default: " << minUnitLen << ")\n"
    "        -x <len>\t\tmaximum repeat unit length of reported microsatellites (default: " << maxUnitLen << ")\n"
    "    Sequencing Read Scan Options:\n"
    "      Quality score flags:\n"
    "        -a \t\t\tflag the quality scores as being in Sanger format (i.e. phred + 33) (default)\n"
    "        -s \t\t\tflag the quality scores as being log-odds scores, as was done in Solexa up to v1.3\n"
    "      Trimming options:\n"
		"        -t \t\t\ttrim first n reads of a read before detection, independent of their quality (default: " << barcodeLength << ") N.B. suitable for removing barcode sequences\n"
    "        -Q <quality score>\tminimum quality value to use when trimming read (default: " << minQual << ")\n"
    "        -B \t\t\ttrim reads using BWA style trimming before scanning (default)\n"
    "        -g \t\t\ttrim reads to last good quality nucleotide before scanning\n"
    "        -F \t\t\ttrim reads to first bad quality nucleotide before scanning\n"
    "        -n \t\t\tdo not trim reads before scanning\n"
    "\n"
    "Output format:\n"
    "\tFASTA input:\n"
    "\t\t<repeat number> <coordinate> <start> <end> <length> <unit> <5' flank> <3' flank> <snapshot>\n"
    "\tFASTQ input:\n"
    "\t\t@<header>\t<microsatellite info>...\n"
    "\t\t<sequence>\n"
    "\t\t+<header>\n"
    "\t\t<quality scores>\n"
    "\tmicrosatellite info is in format <repeat unit>:<start>:<length>\n"
    "\n";
  
  bool errflg = false;
  int ch;
  optarg = NULL;
	
  while (!errflg && ((ch = getopt (argc, argv, "f:q:p:l:k:u:m:x:b:codasnBLFN:hHMQ:T:O:t:")) != EOF))
    {
      switch (ch)
				{
					case 'f' : fastaName	 					= optarg; 	        						break;
					case 'q' : fastqName						= optarg;	        							break;
					case 'p' : fastqName2						= optarg;												break;
					case 'b' : bamName							= optarg;												break;
					case 'T' : targetDir						= optarg;												break;
					case 'O' : outputPrefix					= optarg;												break;
					case 'l' : minLength	 					= atoi (optarg);	        			break;
					case 'u' : minUnits	  					= atoi (optarg);	        			break;
					case 'x' : maxUnitLen		 				= atoi (optarg);	        			break;
					case 'm' : minUnitLen		 				= atoi (optarg); 								break;
					case 'k' : flankLength	 				= atoi (optarg);	        			break;
					case 'Q' : minQual              = atoi (optarg);                break;
					case 'o' : outputToFile					= false;     										break;
					case 'd' : modifySeq						= false; 												break;
					case 'c' : printCountsToFile		= false;												break;
					case 'M' : printUsatOnly				= true;													break;
					case 'N' : nThres 							= atoi(optarg);									break;
					case 'a' : quality_t            = SolexaRead::QUAL_SANGER;      break;
					case 's' : quality_t            = SolexaRead::QUAL_SOLEXA13;    break;
					case 'B' : trim_t               = SolexaRead::TRIM_BWA_STYLE;   break;
					case 'L' : trim_t               = SolexaRead::TRIM_LAST_GOOD;   break;
					case 'F' : trim_t               = SolexaRead::TRIM_FIRST_BAD;   break;
					case 'n' : trim_t               = SolexaRead::TRIM_NONE;        break;
					case 't' : barcodeLength				= atoi (optarg);								break;
					case 'h' : errflg               = true;                         break;
					case 'H' : errflg               = true;                         break;
					  
					case '?':
		  			fprintf (stderr, "Unrecognized option -%c\n", optopt);
					  
					default:
				  	errflg = true;
				}
      
			if (errflg)
				{
	  			cout << helptext.str();
				  exit (EXIT_FAILURE);
				}
		}
  	
  if (fastaName . empty() && fastqName . empty() && bamName . empty())
    {
      cerr << "You must specify a fasta, fastq, or BAM file with the -f, -q or -b flag!" << endl;
      return 1;
    }
  
  if (maxUnitLen > BUFFER_SIZE)
    {
      cerr << "Maximum repeat unit length must be less than " << BUFFER_SIZE << "!" << endl;
      return 1;
    }
	
	MicrosatelliteDetector detector (minLength, minUnits, minUnitLen, maxUnitLen, modifySeq);
	detector . checkSensitivity();
  unsigned int tUsat = 0, seqUsat;
	Timer_t kmerTimer;
	Counter counter,counter1,counter2;
  FILE * usatO = stdout, * modO = NULL, * offsetO = NULL, * hdrO = NULL, * counter0  = NULL;
  FILE * usatO1 = stdout, * modO1 = NULL, * hdrO1 = NULL, * counter01 = NULL;
  FILE * usatO2 = stdout, * modO2 = NULL, * hdrO2 = NULL, * counter02 = NULL;

	string counterDir;
		
	if (fastaName . size())
		{
     	FastaFile fasta (fastaName);
      cerr << "Processing sequence from " << fastaName << endl;
      
			if(barcodeLength > 0)
				{
					cerr << "Barcode trimming specified for FASTA sequences, ignoring..." << endl;
				}
			
			if(printCountsToFile)
				counterDir = createCounterDir(fastaName,outputPrefix,targetDir);
			
			if(outputPrefix . empty())
				{
					usatName 		= replaceSuffix (fastaName, FASTA_OUTPUT_SUFFIX);
					modName		 	= replaceSuffix (fastaName, MODIFIED_FASTA_SUFFIX);
					offsetName 	= replaceSuffix (fastaName, OFFSET_FILE_SUFFIX);
					counterName = replaceSuffix (fastaName, COUNTER_SUFFIX);
				}
			else
				{
					usatName 		= addSuffix (outputPrefix, FASTA_OUTPUT_SUFFIX);
					modName 	 	= addSuffix (outputPrefix, MODIFIED_FASTA_SUFFIX);
					offsetName 	= addSuffix (outputPrefix, OFFSET_FILE_SUFFIX);
					counterName = addSuffix (outputPrefix, COUNTER_SUFFIX);
				}
			if(targetDir . size())
				{
					usatName 		= replacePath(usatName, 	 targetDir);
					modName 	 	= replacePath(modName, 	 	 targetDir);
					offsetName 	= replacePath(offsetName,  targetDir);
				}
      counterName = replacePath(counterName,counterDir);
			
			if(outputToFile)
			  usatO = openAndTestFile (usatName, "w");
			if(printCountsToFile)
				counter0 = openAndTestFile (counterName, "w");
				
      if(modifySeq)
				{
					modO 		= openAndTestFile (modName, "w");
					offsetO = openAndTestFile (offsetName, "w");
				}
			
      kmerTimer . start();
			string hdr;
      while (fasta . getNextSequence ())
				{
					fasta . getHdr(hdr);
					counter.increment("in");
					cerr << "Processing \"" << hdr << "\"... ";
					seqUsat = scanFasta (fasta, detector, usatO, modO, offsetO);
					counter.increment("Microsatellites found",seqUsat);
					tUsat  += seqUsat;
					cerr << kmerTimer . elapsed() << "; found " << seqUsat << " microsatellites" << endl;
				}
    }
  else if (fastqName . size() && fastqName2 . empty())
    {
			if(barcodeLength > 0)
				{
					cerr << "Trimming first " << barcodeLength << " bases from each sequence" << endl;
				}

      FastqFile fastq (fastqName, quality_t, trim_t, minQual, barcodeLength);
      cerr << "Processing sequence reads from " << fastqName << endl;
			cerr << "Trimming selected: "  << fastq . displayTrimType() << endl;

			if(printCountsToFile)
				counterDir = createCounterDir(fastqName,outputPrefix,targetDir);
      
			if(outputPrefix . empty())
				{
				  usatName 		= replaceSuffix (fastqName, FASTQ_OUTPUT_SUFFIX);
					modName 		= replaceSuffix (fastqName, MODIFIED_FASTQ_SUFFIX);
					hdrName 		= replaceSuffix (fastqName, HDR_SUFFIX);
					counterName = replaceSuffix (fastqName, COUNTER_SUFFIX);
				}
			else
				{
					usatName 		= addSuffix (outputPrefix, FASTQ_OUTPUT_SUFFIX);
					modName 		= addSuffix (outputPrefix, MODIFIED_FASTQ_SUFFIX);
					hdrName 		= addSuffix (outputPrefix, HDR_SUFFIX);
					counterName = addSuffix (outputPrefix, COUNTER_SUFFIX);
				}

			if(targetDir . size())
				{
					usatName = replacePath(usatName, targetDir);
					modName = replacePath(modName, targetDir);
					hdrName = replacePath(hdrName, targetDir);
				}
			counterName = replacePath(counterName, counterDir);

      if (outputToFile)
			  usatO = openAndTestFile (usatName, "w");
			if (printCountsToFile)
				counter0 = openAndTestFile (counterName, "w");
			
		 if(modifySeq)
			 {
					modO = openAndTestFile (modName, "w");
					hdrO = openAndTestFile (hdrName, "w");
			 }
      kmerTimer . start();
      while (fastq . getNextSequence())
			{
				counter.increment("in");
				seqUsat = scanFastq (fastq, detector, usatO, modO, hdrO,counter);
				counter.increment("Microsatellites found",seqUsat);
				tUsat += seqUsat;
				nRead++;
				if((nRead % 2000000) == 0)
				{
					cerr << "Processed read number " << nRead << "... " << kmerTimer . elapsed() << endl;
					fprintf(stderr,"\t%d reads submitted for detection, %d  total bases trimmed\n",nGoodLength,tTrimmed);
					fprintf(stderr,"\t%d total reads passed filter\n",nPassed);
					fprintf(stderr,"\tFound %d microsatellites\n",tUsat);
				}
			}
    }
	else if (fastqName . size() && fastqName2 . size())
		{
			if(barcodeLength > 0)
				{
					cerr << "Trimming first " << barcodeLength << " bases from each sequence" << endl;
				}

			FastqFile fastq1 (fastqName,  quality_t, trim_t, minQual, barcodeLength);
			FastqFile fastq2 (fastqName2, quality_t, trim_t, minQual, barcodeLength);
						
      cerr << "Processing sequence reads from " << fastqName << " and " << fastqName2 << endl;
			cerr << "Trimming selected: "  << fastq1 . displayTrimType() << endl;

			if(printCountsToFile)
				counterDir = createCounterDir(fastqName,outputPrefix,targetDir);
      
			if(outputPrefix . size())
				{
					outputPrefix1 = outputPrefix + "_1";
					outputPrefix2 = outputPrefix + "_2";
				}
				
			if(outputPrefix . empty())
				{
				  usatName1 	 = replaceSuffix (fastqName,  FASTQ_OUTPUT_SUFFIX);
					usatName2 	 = replaceSuffix (fastqName2, FASTQ_OUTPUT_SUFFIX);
					modName1 		 = replaceSuffix (fastqName,  MODIFIED_FASTQ_SUFFIX);
					modName2 		 = replaceSuffix (fastqName2, MODIFIED_FASTQ_SUFFIX);
					hdrName1 		 = replaceSuffix (fastqName,  HDR_SUFFIX);
					hdrName2 		 = replaceSuffix (fastqName2, HDR_SUFFIX);
					counterName1 = replaceSuffix (fastqName,  COUNTER_SUFFIX);
					counterName2 = replaceSuffix (fastqName2, COUNTER_SUFFIX);
				}
			else
				{
				//Need to fix output prefix to handle different formats (e.g. s_1_1_sequence or sequence_1)?
					usatName1 	 = addSuffix (outputPrefix1, FASTQ_OUTPUT_SUFFIX);
					usatName2 	 = addSuffix (outputPrefix2, FASTQ_OUTPUT_SUFFIX);						
					modName1 		 = addSuffix (outputPrefix1, MODIFIED_FASTQ_SUFFIX);
					modName2 		 = addSuffix (outputPrefix2, MODIFIED_FASTQ_SUFFIX);
					hdrName1 		 = addSuffix (outputPrefix1, HDR_SUFFIX);
					hdrName2 		 = addSuffix (outputPrefix2, HDR_SUFFIX);
					counterName1 = addSuffix (outputPrefix1, COUNTER_SUFFIX);
					counterName2 = addSuffix (outputPrefix2, COUNTER_SUFFIX);
				}
			if(targetDir . size())
				{
					usatName1 = replacePath(usatName1, targetDir);
					usatName2 = replacePath(usatName2, targetDir);
					modName1 = replacePath(modName1, targetDir);
					modName2 = replacePath(modName2, targetDir);
					hdrName1 = replacePath(hdrName1, targetDir);
					hdrName2 = replacePath(hdrName2, targetDir);
				}
			counterName1 = replacePath(counterName1, counterDir);
			counterName2 = replacePath(counterName2, counterDir);

      if (outputToFile)
				{
				  usatO1 = openAndTestFile (usatName1, "w");
					usatO2 = openAndTestFile (usatName2, "w");
				}
			
			if(printCountsToFile)
				{
					counter01	= openAndTestFile(counterName1, "w");
					counter02	= openAndTestFile(counterName2, "w");
				}
			
		 if(modifySeq)
			 {
					modO1 = openAndTestFile (modName1, "w");
					modO2 = openAndTestFile (modName2, "w");
					hdrO1 = openAndTestFile (hdrName1, "w");
					hdrO2 = openAndTestFile (hdrName2, "w");
			 }

      kmerTimer . start();
      while (fastq1 . getNextSequence())
			{
				nRead1++;
				counter1 . increment("in");
				fastq2 . getNextSequence();
				nRead2++;
				counter2 . increment("in");
				nRead += 2;
				tUsat += scanPairedFastq (fastq1, fastq2, detector, usatO1, usatO2, modO1, modO2, hdrO1, hdrO2, counter1, counter2);
				if((nRead % 2000000) == 0)
				{
					cerr << "Processed read number " << nRead << "... " << kmerTimer . elapsed() << endl;
					fprintf(stderr,"\t%d reads submitted for detection, %d total bases trimmed\n",(nGoodLength1 + nGoodLength2),(tTrimmed1 + tTrimmed2));
					fprintf(stderr,"\t%d total read pairs passed filter\n",nPassed);
					fprintf(stderr,"\tFound %d microsatellites\n",tUsat);
				}
			}
		}
  
	else if (bamName . size())
		{
			BamReader reader;
			reader . Open(bamName);
						
      cerr << "Processing sequence reads from " << bamName << endl;
//			cerr << "Trimming selected: "  << fastqFile1 . displayTrimType() << endl;
      
			if(barcodeLength > 0)
				{
					cerr << "Trimming first " << barcodeLength << " bases from each sequence" << endl;
				}

			if(printCountsToFile)
				counterDir = createCounterDir(bamName,outputPrefix,targetDir);

			if(outputPrefix . size())
				{
					outputPrefix1 = outputPrefix + "_1";
					outputPrefix2 = outputPrefix + "_2";
				}
			else
				{
					outputPrefix1 = getBaseWithPath(bamName) + "_1";
					outputPrefix2 = getBaseWithPath(bamName) + "_2";
				}

			usatName1 	 = addSuffix (outputPrefix1, FASTQ_OUTPUT_SUFFIX);
			usatName2 	 = addSuffix (outputPrefix2, FASTQ_OUTPUT_SUFFIX);
			modName1 		 = addSuffix (outputPrefix1, MODIFIED_FASTQ_SUFFIX);
			modName2 		 = addSuffix (outputPrefix2, MODIFIED_FASTQ_SUFFIX);
			hdrName1 		 = addSuffix (outputPrefix1, HDR_SUFFIX);
			hdrName2 		 = addSuffix (outputPrefix2, HDR_SUFFIX);
			counterName1 = addSuffix (outputPrefix1, COUNTER_SUFFIX);
			counterName2 = addSuffix (outputPrefix2, COUNTER_SUFFIX);
			if(targetDir . size())
				{
					usatName1 = replacePath(usatName1, targetDir);
					usatName2 = replacePath(usatName2, targetDir);
					modName1 = replacePath(modName1, targetDir);
					modName2 = replacePath(modName2, targetDir);
					hdrName1 = replacePath(hdrName1, targetDir);
					hdrName2 = replacePath(hdrName2, targetDir);
				}
			counterName1 = replacePath(counterName1, counterDir);
			counterName2 = replacePath(counterName2, counterDir);

      if (outputToFile)
				{
				  usatO1 = openAndTestFile (usatName1, "w");
					usatO2 = openAndTestFile (usatName2, "w");
				}
			
			if(printCountsToFile)
				{
					counter01	= openAndTestFile(counterName1, "w");
					counter02	= openAndTestFile(counterName2, "w");
				}
			
		 if(modifySeq)
			 {
					modO1 = openAndTestFile (modName1, "w");
					modO2 = openAndTestFile (modName2, "w");
					hdrO1 = openAndTestFile (hdrName1, "w");
					hdrO2 = openAndTestFile (hdrName2, "w");
			 }

      kmerTimer . start();
			BamAlignment al1, al2;
			
      while (reader . GetNextAlignment(al1))
			{
				nRead1++;
				counter1 . increment ("in");
				reader . GetNextAlignment(al2);
				nRead2++;
				counter2 . increment ("in");
				nRead += 2;
				
				tUsat += scanBam (al1, al2, detector, usatO1, usatO2, modO1, modO2, hdrO1, hdrO2, quality_t, trim_t, counter1, counter2, barcodeLength);
				if((nRead % 2000000) == 0)
				{
					cerr << "Processed read number " << nRead << "... " << kmerTimer . elapsed() << endl;
					cerr << "\t"  << (nGoodLength1 + nGoodLength2) << " reads submitted for detection, " << (tTrimmed1 + tTrimmed2) << " total bases trimmed" << endl;
					fprintf(stderr,"\t%d total bases trimmed from first mate-pair reads (%.02f bases per read)\n", tTrimmed1, (float)(tTrimmed1 / nRead1));
					fprintf(stderr,"\t%d total bases trimmed from second mate-pair reads (%.02f bases per read)\n", tTrimmed2, (float)(tTrimmed2 / nRead2));
					cerr << "\tFound " << tUsat << " microsatellites" << endl;
				}
			}
		}
	
  cerr << "Found a total of " << tUsat << " microsatellites" << endl;
	if(fastaName . size())
		{
			if(outputToFile)
				{
					cerr << "File with microsatellite information can be found in " << usatName <<endl;
					fclose(usatO);
				}
			if(modifySeq)
				{
					cerr << "File with offset microsatellite coordinate for " << fastaName << " can be found in " << offsetName << endl;
					cerr << "A file with a modified genome can be found in " << modName << endl;
					fclose(offsetO);
					fclose(modO);
				}
			if(printCountsToFile)
				{
					counter.print(counter0);
					cerr << "File with counts for " << fastaName << " can be found in " << counterName << endl;
					fclose(counter0);
				}
		}
		
	else if(fastqName . size()  && fastqName2 . empty())
		{
			fprintf(stderr,"Processed %d reads, %d (%.02f%%) passed filter\n",nRead,nPassed,(float)((nPassed*100)/nRead));
			if(nPassed > 0)
				{
					fprintf(stderr,"\t%d total bases trimmed from passed reads (%.02f bases per read), %d submitted for detection\n",tTrimmed,(float)(tTrimmed / nPassed),nGoodLength);
					fprintf(stderr,"\t%d  reads were too short after trimming\n",nPassed - nGoodLength);
				}
			if(outputToFile)
				{
					cerr << "File with microsatellite information can be found in " << usatName <<endl;
					fclose(usatO);
				}
			if(modifySeq)
				{
					cerr << "File with read headers and microsatellite information can be found in " << hdrName << endl;
					cerr << "A file with modified reads can be found in " << modName << endl;
					fclose(hdrO);
					fclose(modO);
				}
			if(printCountsToFile)
				{
					counter.print(counter0);
					cerr << "File with counts for " << fastqName << " can be found in " << counterName << endl;
					fclose(counter0);
				}
		}
		
	else if ((fastqName . size()  && fastqName2 . size()) || bamName . size())
		{
			nGoodLength = nGoodLength1 + nGoodLength2;
			nRead = nRead1 + nRead2;

			if(fastqName . size()  && fastqName2 . size())
				{
					cerr << "In " << fastqName << ":" << endl;
					fprintf(stderr,"Processed %d reads, %d (%.02f%%) passed filter\n",nRead1,nPassed1, (float)((nPassed*100)/nRead));
					if(nPassed1 > 0)
						{
							fprintf(stderr,"\t%d total bases trimmed from passed reads (%.02f bases per read), %d submitted for detection\n",tTrimmed1,(float)(tTrimmed1 / nPassed1),nGoodLength1);
							fprintf(stderr,"\t%d  reads were too short after trimming\n",nPassed1 - nGoodLength1);
						}
					cerr << "In " << fastqName2 << endl;
					fprintf(stderr,"Processed %d reads, %d (%.02f%%) passed filter\n",nRead2,nPassed2, (float)((nPassed*100)/nRead));
					if(nPassed2 > 0)
						{
							fprintf(stderr,"\t%d total bases trimmed from passed reads (%.02f bases per read), %d submitted for detection\n",tTrimmed2,(float)(tTrimmed2 / nPassed2),nGoodLength2);
							fprintf(stderr,"\t%d  reads were too short after trimming\n",nPassed2 - nGoodLength2);
						}
				}
			else if(bamName . size())
				{
					cerr << "In " << bamName << endl;
					cerr << "\tProcessed " << nRead << " reads, " << nGoodLength << " submitted for detection" << endl;
					if(nRead1 > 0 && nRead2 > 0)
					{
						fprintf(stderr,"\t%d total bases trimmed from first mate-pair reads (%.02f bases per read)\n", tTrimmed1, (float)(tTrimmed1 / nRead1));
						fprintf(stderr,"\t%d total bases trimmed from second mate-pair reads (%.02f bases per read)\n", tTrimmed2, (float)(tTrimmed2 / nRead2));
						cerr << "\t" << nRead - nGoodLength << " reads were too short after trimming" << endl;
					}
					else
					{
						fprintf(stderr,"%s is empty!\n",bamName.c_str());
						return 1; //returns an error status of 1, which can be used to kill downstream processes (not implemented yet)
					}
				}

			if(outputToFile)
				{
					fprintf(stderr,"Files with microsatellite information can be found at %s and %s\n", usatName1 . c_str(), usatName2 . c_str());
					fclose(usatO1);
					fclose(usatO2);
				}
			if(modifySeq)
				{
					fprintf(stderr,"Files with read headers can be found at %s and %s\n",hdrName1 . c_str(),hdrName2 . c_str());
					fprintf(stderr,"Files with modified reads can be found at %s and %s\n",modName1 . c_str(),modName2 . c_str());
					fclose(hdrO1);
					fclose(hdrO2);
					fclose(modO1);
					fclose(modO2);
				}
			if(printCountsToFile)
				{
					counter1.print(counter01);
					counter2.print(counter02);
					fprintf(stderr,"Files with counts can be found at %s and %s\n",counterName1 . c_str(), counterName2 . c_str());
					fclose(counter01);
					fclose(counter02);

				}
	}
 cerr << "Tetrascan completed in " << kmerTimer.elapsed() << endl;

  return 0;
}
