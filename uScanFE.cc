#include "uScanFE.h"

const unsigned int MIN_READ_LENGTH = 10; //shortest possible length for detection

 // File variables
string fastaName, fastqName, fastqName2, bamName;
string usatName, modName, offsetName, hdrName, counterName;
string usatName1, usatName2, modName1, modName2, hdrName1, hdrName2, counterName1, counterName2;

// Output path variables
string targetDir;
string outputPrefix, outputPrefix1, outputPrefix2;

// Quality score and read trimming parameters with defaults
SolexaRead::Trim_t trim_t = SolexaRead::TRIM_BWA_STYLE;
SolexaRead::Quality_t quality_t = SolexaRead::QUAL_SANGER;
unsigned int nThres = 3, minQual = 20;

 // Microsatellite detection parameters with defaults
int minLength = 8, minUnits = 3, minUnitLen = 1, maxUnitLen = 6;

// Other tetrascan defaults
unsigned int flankLength = 30;
bool printUsatOnly = false, modifySeq = true, outputToFile = true, printCountsToFile = true;

// Counters
unsigned int nGoodLength  = 0, nRead  = 0, tTrimmed  = 0, nPassed  = 0;
unsigned int nGoodLength1 = 0, nRead1 = 0, tTrimmed1 = 0, nPassed1 = 0;
unsigned int nGoodLength2 = 0, nRead2 = 0, tTrimmed2 = 0, nPassed2 = 0;

int scanFasta (FastaFile & fasta, MicrosatelliteDetector & detector, FILE * usatO, FILE * modO,
							 FILE * offsetO)
{
	string seq, hdr;
	fasta . getSeq (seq);
	fasta . getHdr (hdr);

	vector< Microsatellite > usatInSeq = detector . detectPerfectMicrosatellites(seq, hdr, flankLength);

	for(unsigned int i = 0; i < usatInSeq . size(); ++i)
		usatInSeq[i] . printFastaStyle (usatO, seq);
	
	if(modifySeq)
	{
		string modSeq;
		getModSeq(modSeq, seq, usatInSeq);
		
		unsigned int ntPerLine      = fasta . getNtPerLine();
		unsigned int tSeqNt         = modSeq . length();
		unsigned int seqNumFullLine = tSeqNt / ntPerLine;

		fprintf(modO,">%s\n",hdr . c_str());
		for(unsigned int i = 0; i < seqNumFullLine; ++i)
		{
			fprintf(modO, "%s\n", (modSeq . substr((ntPerLine * i), ntPerLine)) . c_str());
		}
		fprintf (modO,"%s\n",(modSeq . substr ((ntPerLine * seqNumFullLine), 
												 (tSeqNt - (seqNumFullLine * tSeqNt)))) . c_str());
		createOffsetFile(usatInSeq, offsetO);
	}

	return usatInSeq . size();
}

void getModSeq(string & modSeq, const string & seq, const vector< Microsatellite > & usatInSeq)
{
	unsigned int nUsat = usatInSeq . size();
	vector< unsigned int > pos (nUsat + 1);
	vector< unsigned int > length (nUsat + 1);
	unsigned int readEnd  = seq . length();
	unsigned int nOverlap = 0;
	
	for(unsigned int i = 0, j = 0; i < nUsat; ++i, ++j)
	{
		// if this microsatellite starts before last stop recorded (i.e. overlap the preceding
		// microsatellite)
		if(usatInSeq [i] . start() - 1 <= pos[j])
		{
			pos[j] = usatInSeq[i] . stop();
			pos . resize(pos . size() - 1);
			length . resize(length . size() - 1);
			j--;
			nOverlap++;
		}
		else
		{
			length[j]  = (usatInSeq[i] . start()) - pos[j] - 1;
			pos[j + 1] =  usatInSeq[i] . stop();
		}
		//if microsatellite is at end of sequence
		if(usatInSeq [i] . stop() == readEnd)
		{
			pos . resize(pos . size() - 1);
			length . resize (length . size() - 1);
			if(usatInSeq [i] . start() != 1)
				readEnd = usatInSeq [i - nOverlap] .start() - 1;
		}
		nOverlap = 0;
	}
		
	if(pos . size() > 0)
	{
		length[pos . size() - 1] = readEnd - pos[pos . size() - 1];
		for(unsigned int i = 0; i < pos . size(); ++i)
			{
				modSeq . append (seq . substr (pos[i], length[i]));
			}			
	}			
}

void createOffsetFile(const vector< Microsatellite > & usatInSeq, FILE * offsetO)
{
	vector< int > positions, offsets;
	string currentCoord = usatInSeq [0] . hdr();
	positions . clear();
	offsets . clear();
	int currentOffset = 0;

	for(unsigned int i = 0; i < usatInSeq . size(); ++i)
	{		
		if ((currentCoord . compare(usatInSeq [i] . hdr())) == 0)
		{
			currentOffset += usatInSeq[i] . length();
			if((i + 1) < usatInSeq . size())
			{
				// Make sure there is at least one nucleotide between the two microsatellites
				if (usatInSeq[i] . stop() < usatInSeq[i + 1] . start() - 1)
				{
					//positions are assuming that first bp in sequence is at position 1, not zero
					positions . push_back(usatInSeq[i] . stop() - currentOffset + 1);
					offsets . push_back(currentOffset);
				}
				else
				{
					currentOffset += (usatInSeq [i + 1]. start() - usatInSeq [i] . stop() - 1);
				}
			}
			//This needs to take into account terminal overlapping microsatellites
			else
			{
				positions . push_back(usatInSeq [i] . stop() - currentOffset + 1);
				offsets . push_back(currentOffset);
			}
		}
		else
		{
			printCoordOffsets(offsetO,currentCoord,positions,offsets);

			positions . clear();
			offsets . clear();

			currentCoord  = usatInSeq[i] . hdr();
			currentOffset = usatInSeq[i] . length();
			positions . push_back(usatInSeq[i] . stop() - currentOffset + 2);
			offsets . push_back(currentOffset);
		}
	}
	printCoordOffsets(offsetO,currentCoord,positions, offsets);
}

void printCoordOffsets(FILE * offsetO, const string & coord, const vector< int > & positions, const vector< int > & offsets)
{
	fprintf(offsetO, "%s\n", coord . c_str());
	for(unsigned int i = 0; i < positions . size(); ++i)
		fprintf(offsetO, "%d\t", positions [i]);
	fprintf(offsetO, "\n");
	for(unsigned int i = 0; i < offsets . size(); ++i)
		fprintf(offsetO, "%d\t", offsets[i]);
	fprintf(offsetO, "\n");
}

void readInfo(const SolexaRead & read, string & hdr, string & seq, string & qual)
{
	hdr  = read . hdr();
	seq  = read . goodSeq();
	qual = read . goodQual();
}

int scanFastq (FastqFile & fastq, MicrosatelliteDetector & detector, FILE * usatO, FILE * modO, FILE * hdrO, Counter & counter)
{
  vector< Microsatellite > usatInRead;
  SolexaRead read;
	fastq . getLatestRead(read);

	if(read . passed())
		{
			nPassed++;
			string hdr, seq, qual;
			readInfo(read, hdr, seq, qual);
			tTrimmed += (read . readLen() - seq . length());
			counter.increment("Trimmed bases",int(read . readLen() - seq . length()));

  		if (seq . length() >= MIN_READ_LENGTH && read . nCount() < nThres)
    		{
					nGoodLength++;
		  		usatInRead = detector . detectPerfectMicrosatellites (seq, hdr, flankLength);
		  		read . setMicrosatellites (usatInRead);
    		}
			else if(seq . length() < MIN_READ_LENGTH)
			{
				counter.increment("Too short");
			}
			else
			{
				counter.increment("Too many Ns");
			}

			if ((printUsatOnly && read . nUsat() > 0) || !printUsatOnly)
				{
					read . printMicrosatelliteFastqWithQual(qual,usatO);
					if(modifySeq)
						{
							//Even though reads with terminal microsatellites cannot provide microsatellite allele length, they 
							//can still be mapped by excising detected microsatellite sequences
							read . printModifiedReadWithQualNoModHeader(qual,modO);
							read . printMicrosatelliteHeader(qual, hdrO);
						}
					counter.increment("out");
				}
			return read . nUsat();
		}
	else
		{
			counter.increment("Failed sequencer QC");
		}
	return 0;		
}

int scanPairedFastq (FastqFile & fastq1, FastqFile & fastq2, MicrosatelliteDetector & detector, FILE * usatO1, FILE * usatO2, FILE * modO1, FILE * modO2, FILE * hdrO1, FILE * hdrO2, Counter & counter1, Counter & counter2)
{
  vector< Microsatellite > usatInRead;
	SolexaRead read1, read2;
	string hdr, seq, qual1, qual2;
	
	//Look for microsatellites in first read
	fastq1 . getLatestRead(read1);
	fastq2 . getLatestRead(read2);
	
	if(read1 . passed() && read2 . passed())
		{
			nPassed+=2;
			nPassed1++;
			nPassed2++;
			readInfo (read1, hdr, seq, qual1);
			tTrimmed1 += (read1 . readLen() - seq . length());	
			counter1.increment("Trimmed bases",int(read1 . readLen() - seq . length()));

  		if (seq . length() >= MIN_READ_LENGTH && read1 . nCount() < nThres)
    		{
					nGoodLength1++;
		  		usatInRead = detector . detectPerfectMicrosatellites (seq, hdr, flankLength);
		  		read1 . setMicrosatellites (usatInRead);
					counter1 . increment("Microsatellites found",read1 . nUsat());
    		}
			else if(seq . length() < MIN_READ_LENGTH)
				{
					counter1 . increment("Too short");
				}
			else
				{
					counter1 . increment("Too many Ns");
				}

			//Look for microsatellites in second read
			readInfo (read2, hdr, seq, qual2);
			tTrimmed2 += (read2 . readLen() - seq . length());
			counter2.increment("Trimmed bases",int(read2 . readLen() - seq . length()));

			if(seq . length() >= MIN_READ_LENGTH && read2 . nCount() < nThres)
				{
					nGoodLength2++;
					usatInRead = detector . detectPerfectMicrosatellites (seq, hdr, flankLength);
					read2 . setMicrosatellites(usatInRead);
					counter2 . increment("Microsatellites found", read2 . nUsat());
				}
			else if(seq . length() < MIN_READ_LENGTH)
				{
					counter2 . increment("Too short");
				}
			else
				{
					counter2 . increment("Too many Ns");
				}

			bool read1First,read2First;
			if(read1.readNumber() == 1)
				{
					read1First = true;
					read2First = false;
				}
			else
				{
					read2First = true;
					read1First = false;
				}

			if ((printUsatOnly && (read1 . nUsat() > 0 || read2 . nUsat() > 0)) || !printUsatOnly)
				{
					read1 . printMicrosatelliteFastqWithQual(qual1,usatO1);
					read2 . printMicrosatelliteFastqWithQual(qual2,usatO2);
					if(modifySeq)
						{
							read1 . printModifiedReadWithQualNoModHeader(qual1,modO1);
							read2 . printModifiedReadWithQualNoModHeader(qual2,modO2);
							read1 . printMicrosatelliteHeader(qual1, hdrO1,read1First);
							read2 . printMicrosatelliteHeader(qual2, hdrO2,read2First);
						}
					counter1 . increment("out");
					counter2 . increment("out");
				}

			return (read1 . nUsat() + read2 . nUsat());
		}
	else
		{
			counter1.increment("Failed sequencer QC");
			counter2.increment("Failed sequencer QC");
		}
	return 0;
}

int scanBam (const BamAlignment & al1, const BamAlignment & al2, MicrosatelliteDetector & detector, FILE * usatO1, FILE * usatO2, FILE * modO1, FILE * modO2, FILE * hdrO1, FILE * hdrO2, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, Counter &
counter1, Counter & counter2,unsigned int barcodeLength)
{
	vector< Microsatellite > usatInRead;
	string hdr, seq, qual1, qual2;
	
	if((al1 . Name).compare(al2 . Name) != 0)
		{
			fprintf(stderr,"Error: read names do not match: %s != %s",al1 . Name . c_str(), al2 . Name . c_str());
			exit(EXIT_FAILURE);
		}
	
	if(!al1.IsFirstMate())
		{
			fprintf(stderr,"Error: read provided as first mate is not first mate");
			exit(EXIT_FAILURE);
		}
	if(!al2.IsSecondMate())
		{
			fprintf(stderr,"Error: read provided as second mate is not second mate");
			exit(EXIT_FAILURE);
		}
	
	SolexaRead read1(al1.Name,al1.QueryBases,al1.Qualities,qualityType,trimType,barcodeLength);
	read1.trim();
	readInfo (read1, hdr, seq, qual1);
	tTrimmed1 += (read1 . getLSeq() - seq . length());	
	counter1.increment("Trimmed bases", int(read1 . getLSeq() - seq . length()));
	
  if (seq . length() >= MIN_READ_LENGTH && read1 . nCount() < nThres)
    {
			nGoodLength1++;
		  usatInRead = detector . detectPerfectMicrosatellites (seq, hdr, flankLength);
		  read1 . setMicrosatellites (usatInRead);
			counter1 . increment("Microsatellites found",read1 . nUsat());
    }
	else if (seq . length() < MIN_READ_LENGTH)
		{
			counter1 . increment("Too short");
		}
	else
		{
			counter1 . increment("Too many Ns");
		}

	SolexaRead read2(al2.Name,al2.QueryBases,al2.Qualities,qualityType,trimType,barcodeLength);
	read2.trim();
	readInfo (read2, hdr, seq, qual2);	
	tTrimmed2 += (read2 . getLSeq() - seq . length());	
	counter2.increment("Trimmed bases", int(read2 . getLSeq() - seq . length()));
  
	if (seq . length() >= MIN_READ_LENGTH && read2 . getNCount() < nThres)
    {
			nGoodLength2++;
		  usatInRead = detector . detectPerfectMicrosatellites (seq, hdr, flankLength);
		  read2 . setMicrosatellites (usatInRead);
			counter2 . increment("Microsatellites found",read2 . nUsat());
    }
	else if(seq . length() < MIN_READ_LENGTH)
		{
			counter2 . increment("Too short");
		}
	else
		{
			counter2 . increment("Too many Ns");
		}

	if ((printUsatOnly && (read1 . nUsat() > 0 || read2 . nUsat() > 0)) || !printUsatOnly)
		{
			read1 . printMicrosatelliteFastqWithQual(qual1,usatO1);
			read2 . printMicrosatelliteFastqWithQual(qual2,usatO2);
			if(modifySeq)
				{
					read1 . printModifiedReadWithQualNoModHeader(qual1,modO1);
					read2 . printModifiedReadWithQualNoModHeader(qual2,modO2);
					read1 . printMicrosatelliteHeader(qual1, hdrO1,al1.IsFirstMate());
					read2 . printMicrosatelliteHeader(qual2, hdrO2,al2.IsFirstMate());
				}
			counter1 . increment("out");
			counter2 . increment("out");
		}

	return (read1 . nUsat() + read2 . nUsat());
}
