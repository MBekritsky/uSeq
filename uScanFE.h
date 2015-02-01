#ifndef __USCANFE_H__
#define __USCANFE_H__

#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>

#include "api/BamReader.h"

#include "Microsatellite.h"
#include "MicrosatelliteDetector.h"
#include "FastaFile.h"
#include "FastqFile.h"
#include "Counter.h"

using namespace std;
using namespace BamTools;

//Constants
extern const unsigned int MIN_READ_LENGTH;

 //File variables
extern string fastaName, fastqName, fastqName2, bamName;
extern string usatName, modName, offsetName, hdrName, counterName;
extern string usatName1, usatName2, modName1, modName2, hdrName1, hdrName2, counterName1, counterName2;
//Path variables
extern string targetDir;
extern string outputPrefix, outputPrefix1, outputPrefix2;

//parameters and defaults
//quality score and read trimming parameters with defaults
extern SolexaRead::Trim_t trim_t;
extern SolexaRead::Quality_t quality_t;
extern unsigned int nThres, minQual;

 //Microsatellite detector parameters with defaults
extern int minLength, minUnits, minUnitLen, maxUnitLen;

//other defaults for tetrascan
extern unsigned int flankLength;
extern bool printUsatOnly, modifySeq, outputToFile, printCountsToFile;

//counters
extern unsigned int nGoodLength, nRead, tTrimmed, nPassed;
extern unsigned int nGoodLength1, nRead1, tTrimmed1, nPassed1;
extern unsigned int nGoodLength2, nRead2, tTrimmed2, nPassed2;

//detector frontend for different input file types
int scanFasta (FastaFile & fasta, MicrosatelliteDetector & detector, FILE * usatO, FILE * modO, FILE * offsetO);
int scanFastq (FastqFile & fastq, MicrosatelliteDetector & detector, FILE * usatO, FILE * modO, FILE * headerO, 
							 Counter & counter);
int scanPairedFastq (FastqFile & fastq1, FastqFile & fastq2, MicrosatelliteDetector & detector, FILE * usatO1, 
										 FILE * usatO2, FILE * modO1, FILE * modO2, FILE * hdrO1, FILE * hdrO2, Counter & counter1,
										 Counter & counter2);
int scanBam(const BamAlignment & al1, const BamAlignment & al2, MicrosatelliteDetector & detector, FILE * usatO1, 
						FILE * usatO2, FILE * modO1, FILE * modO2, FILE * hdrO1, FILE * hdrO2, SolexaRead::Quality_t qualityType, 
						SolexaRead::Trim_t trimType, Counter & counter1, Counter & counter2, unsigned int barcodeTrim);

//Functions for printing modified genome and offset file (may be rolled into fasta file class)
void getModSeq(string & modSeq, const string & seq, const vector< Microsatellite > & usatInSeq);
void createOffsetFile(const vector< Microsatellite > & microsatellitesInSequence, FILE * offsetO);
void printCoordOffsets(FILE * offsetO, const string & coord, const vector< int > & positions, const vector< int > & offsets);

//auxiliary functions
void readInfo(const SolexaRead & read, string & hdr, string & seq, string & qual);


#endif /* uScanFE.h */
