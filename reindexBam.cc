/*
reindexBam.cc

ABSTRACT:
A program that reindexes reads and microsatellites in a SAM file that 
is obtained after alignment.  The aligned position in the SAM files 
remains unchanged, but OP and MC tags are added given the read's 
original position and microsatellite coordinates, respectively. If the
read mapped to the reverse complement of the reference sequence, the 
microsatellites in the read are re-indexed to be on the original 
strand (e.g. if a 76 bp read has A:16:10, it would now be T:61:10) 
in the MC tag only.

Finally, the read can be checked against the microsatellite database 
for the reference sequence. If it is checked against the database and
a microsatellite at that position cannot be found, the read is failed
(i.e. all microsatellites in the read will go from uppercase to lower-
case).  About 2.5% of the time, mapping errors will occur that cause 
microsatellites to fail this reference check.  The reasons for this 
can be elaborated upon elsewhere.  If you are interested in identifying
novel microsatellite loci, this check should not be performed.  I may
add a check for novel microsatellites at a later date.

CREATION DATE:
17.05.2013

LAST REVISION:
17.05.2013

AUTHOR:
Mitchell Bekritsky
*/

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>

//included from BamTools to read and write to BAM files
#include "api/BamWriter.h"
#include "api/BamReader.h"

#include "Timer.h"
#include "FileOps.h"
#include "BamMicrosatelliteInfo.h"
#include "OsiFile.h"
#include "Microsatellite.h"
#include "MsHeader.h"
#include "MsHeaderFile.h"
#include "Counter.h"
#include "OrphanFunctions.h"

using namespace std;
using namespace BamTools;

string BAM_FILE, OSI_FILE, HDR_FILE_1, HDR_FILE_2, MSI_OUTPUT_NAME, COUNTER_NAME, PREFIX, TARGET_DIR;
string MSI_SUFFIX = ".msi.bam";
string COUNT_SUFFIX = ".reindex.count";
bool toSTDOUT = false;

string reindexVersion = "1.0";
string reindexID			= "reindexBam";

int getOffset(BamAlignment & al, const string & readCoord, BamMicrosatelliteInfo & msInfo, const OsiFile & osiFile);
int findPositionIndex(int readStart, int coordIndex, const OsiFile & osiFile);
int checkOffsetIndex(int coordIndex, int vecIndex, const OsiFile & osiFile, const BamAlignment & al, const BamMicrosatelliteInfo & msInfo);
unsigned int getMissedBasePairs(const BamAlignment & al, const BamMicrosatelliteInfo & msInfo);

int getOffset(BamAlignment & al, const string & readCoord, BamMicrosatelliteInfo & msInfo, const OsiFile & osiFile)
{

	//don't re-index a read that has failed to map or is a multi-mapper 
	if(!al.IsMapped() || al.MapQuality == 0)
		return -1;

	int readStart = al.Position + 1; //Position in BamAlignment is zero-indexed, the offset index is 1-indexed
	int coordIndex = osiFile . getCoordIndex(readCoord);
	
	if(coordIndex == -2)
		{
			cerr << "Warning!!! " << readCoord << " from " << al.Name << " does not match any coordinates in " << OSI_FILE << ".  Please check to make sure you are using the correct *.osi file." << endl;
			exit(EXIT_FAILURE);
		}

	//if the current read maps to the reverse strand
	if(al.IsReverseStrand() && msInfo.nMicrosatellites() > 0)
			msInfo . reverseMicrosatellites ();

	int positionIndex = findPositionIndex(readStart, coordIndex, osiFile);
		
	if(positionIndex >= 0)
		{
			return checkOffsetIndex(coordIndex,positionIndex, osiFile, al, msInfo);
		}
	//if findPositionIndex returns -1, it means that the read's start position is before the first offset
	else if (positionIndex == -1)
		{
			return 0;
		}
	//if findPositionIndex returns a negative number greater than -1, it means that the read starts after the last offset
	else
		{
			return osiFile . getOffset (coordIndex, osiFile . getVectorSize(coordIndex) - 1);
		}

	//If the function gets to here, that means the binary search completely failed, this should not happen
	cerr << "Failed to find position vector coordinates containing " << readStart << endl;
	exit(EXIT_FAILURE);
	return -2;
}

//Modified binary search to find position of read in position vector
int findPositionIndex(int readStart, int coordIndex, const OsiFile & osiFile)
{
	unsigned int first = 0, last = osiFile . getVectorSize(coordIndex) - 1;
	unsigned int mid = 0;
		while(first < last)
		{
			mid = (last + first)/2;
			if(readStart > osiFile . getPosition (coordIndex, mid + 1))
				{
					first = mid + 1;
				}
			else if(readStart <= osiFile . getPosition (coordIndex, mid))
				{
					last = mid;
				}
			else
				{
					return mid;
				}
		}
	return -(first + 1);
}

int checkOffsetIndex(int coordIndex, int vecIndex, const OsiFile & osiFile, const BamAlignment & al, const BamMicrosatelliteInfo & msInfo)
{
	int readStart = al.Position + 1;
	
	//If the read does not start within a microsatellite, return the offset index
	if(readStart > osiFile . getPosition(coordIndex,vecIndex) && readStart < osiFile . getPosition(coordIndex, vecIndex + 1))
		{
			return osiFile . getOffset(coordIndex, vecIndex);
		}
	else if (readStart == osiFile . getPosition (coordIndex, vecIndex + 1))
		{
			//If there are no microsatellites in the read, but it hits a microsatellite coordinate exactly, it starts after the microsatellite end
			if(msInfo.nMicrosatellites() == 0)
				{
					return osiFile . getOffset(coordIndex, vecIndex + 1);
				}
			//If there is a microsatellite in the read, see if any of them are terminal and add it to the offset of the next microsatellite end
			else
				{
					return osiFile .  getOffset(coordIndex, vecIndex + 1) - getMissedBasePairs(al, msInfo);
				}
		}
	return -1;
}

unsigned int getMissedBasePairs(const BamAlignment & al, const BamMicrosatelliteInfo & msInfo)
{
	vector< Microsatellite > microsatellites;
	msInfo.getMicrosatellites(microsatellites);
	
			unsigned int terminalMicrosatelliteLengths = 0;
			int i = 0;

			while(i < (int) microsatellites . size() && microsatellites [i] . getStart() <= terminalMicrosatelliteLengths + 1)
				{
					terminalMicrosatelliteLengths += microsatellites [i] . getLength();
					if((i - 1) >= 0)
						terminalMicrosatelliteLengths += microsatellites [i] . getStart() - microsatellites [i - 1] . getStart() - microsatellites [i - 1] . getLength();
					i++;
				}
			return terminalMicrosatelliteLengths;

			cerr << "Unrecognized flag in " << al.Name << ": " << al.AlignmentFlag << endl;
			exit(EXIT_FAILURE);
}

int main(int argc, char * argv [])
{
	stringstream helptext;
	helptext << "\n";
	
	bool errflg = false;
	int ch;
	optarg = NULL;
	
	while (!errflg && ((ch = getopt (argc,argv, "B:I:1:2:P:T:ohH")) != EOF))
		{
			switch(ch)
				{
					case 'B' : BAM_FILE     = optarg; 			break;
					case 'P' : PREFIX				= optarg;				break;
					case 'T' : TARGET_DIR		= optarg;				break;
					case 'I' : OSI_FILE   	= optarg; 			break;
					case '1' : HDR_FILE_1		= optarg;				break;
					case '2' : HDR_FILE_2		= optarg;				break;
					case 'o' : toSTDOUT			= true;					break;
					case 'H' : errflg       = true;   			break;
					case 'h' : errflg       = true;   			break;
					
					case '?':
						fprintf (stderr, "Unrecognized option -%c\n", optopt);
						errflg = true;
					
					default:
						errflg = true;
				}
				
			if (errflg)
				{
					cerr << helptext . str();
					exit (EXIT_FAILURE);
				}
		}
	
	Counter counter;

	if(OSI_FILE . empty() || BAM_FILE . empty())
		{
			if(OSI_FILE . empty() && BAM_FILE . empty())
				{
					cerr << "You must specify an index file (*.osi) with the -I flag and BAM file with the -B flag!" <<	endl;
					return 1;
				}
			else if(BAM_FILE . empty())
				{
					cerr << "You must specify a BAM file with the -B flag!" << endl;
					return 1;
				}
			else
				{
					cerr << "You must specify an index file (*.osi) with the -I flag!" << endl;
					return 1;				
				}
		}
		
	string commandLine = reconstituteCommandLine(argc,argv);
 
	SamProgram reindexProgram(reindexID);
	reindexProgram.Version = reindexVersion;
	reindexProgram.Name		 = reindexID;
	reindexProgram.CommandLine = commandLine;

	Timer_t trans_timer;
	MsHeaderFile hdrFile1 (HDR_FILE_1);
	MsHeaderFile hdrFile2 (HDR_FILE_2);
	BamReader inBam;
	BamAlignment al;
	BamWriter outBam;

	if(!inBam.Open(BAM_FILE))
	{
		cerr << "Could not read from " << BAM_FILE << ": " << inBam.GetErrorString() << endl;
		exit(EXIT_FAILURE);
	}

	Timer_t loadTimer;
	OsiFile osiFile(OSI_FILE);
	cerr << "Loaded index file " << OSI_FILE << " in " << loadTimer.elapsed() << endl;
	if(toSTDOUT)
	{
		MSI_OUTPUT_NAME = "stdout";
		outBam.SetCompressionMode(BamWriter::Uncompressed);
	}
	else
	{
		MSI_OUTPUT_NAME = replaceSuffix(BAM_FILE, MSI_SUFFIX);
	}
	
	SamHeader inHeader = inBam.GetHeader();
	inHeader.Programs.Add(reindexProgram);
	if(!outBam.Open(MSI_OUTPUT_NAME,inHeader,inBam.GetReferenceData()))
	{
		cerr << "Could not write to " << MSI_OUTPUT_NAME << ": " << outBam.GetErrorString() << endl;
	}

	if(PREFIX.empty() && TARGET_DIR.empty())
	{
		COUNTER_NAME = replaceSuffix(BAM_FILE, COUNT_SUFFIX);
		COUNTER_NAME = replacePath(COUNTER_NAME, createCounterDir(BAM_FILE));
	}
	else if(!PREFIX.empty() && TARGET_DIR.empty())
	{
		COUNTER_NAME = addSuffix(PREFIX, COUNT_SUFFIX);
	}
	else if(PREFIX.empty() && !TARGET_DIR.empty() && !toSTDOUT)
	{
		COUNTER_NAME = replacePath(BAM_FILE,TARGET_DIR);
		COUNTER_NAME = replaceSuffix(COUNTER_NAME,COUNT_SUFFIX);
		COUNTER_NAME = replacePath(COUNTER_NAME, createCounterDir(COUNTER_NAME));
	}
	else if(PREFIX.empty() && !TARGET_DIR.empty() && toSTDOUT)
	{
		cerr << "Cannot create counts directory when piping to STDOUT with a target directory but no prefix" << endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		COUNTER_NAME = addSuffix(PREFIX,COUNT_SUFFIX);
		if(COUNTER_NAME.find_first_of('/') == string::npos)
		{
			COUNTER_NAME = addPath(COUNTER_NAME,TARGET_DIR);
		}
		else
		{
			COUNTER_NAME = replacePath(COUNTER_NAME,TARGET_DIR);
		}
		COUNTER_NAME = replacePath(COUNTER_NAME, createCounterDir(COUNTER_NAME));
	}

	FILE * counter_output = openAndTestFile(COUNTER_NAME, "w");

	MsHeader currentHdr;
	vector< string > tMs;
	vector< string > tMsQ;

	int offset = 0;

	cerr << "Processing " << BAM_FILE . c_str() << endl;

	BamMicrosatelliteInfo msInfo;
	RefVector refData = inBam.GetReferenceData();

	while (inBam.GetNextAlignment(al))
	{
		msInfo . plainQName(al.Name);
		counter.increment("in");
		
		hdrFile1 . getNextHdr();
		hdrFile1 . getHdr(currentHdr);

		if((currentHdr . compareQNames(al.Name)) && al.IsFirstMate() == currentHdr . isFirst())
			{
				currentHdr . getMicrosatellites(tMs);
				currentHdr . getMicrosatelliteQuals(tMsQ);
				msInfo . addMicrosatellites(tMs);
				msInfo . addMicrosatelliteQuals(tMsQ);
				msInfo . readLength(al);
			}
		else
			{
				counter.increment("Mismatched header");
				counter.print(counter_output);
				fclose(counter_output);
				cerr << "Error!! The current header in " << BAM_FILE << " does not match the current header in " << HDR_FILE_1 << endl;
				cerr << "\"" << al.Name << "\" != \"" << currentHdr . getQName() << "\"" << endl;
				cerr << "Flag: " << al.AlignmentFlag << "\t";
				cerr << "\tHDR MATE: " << currentHdr . isFirst() << "\t" << al.IsFirstMate() << endl;
				return -1;
			}
		offset = getOffset(al, refData[al.RefID].RefName, msInfo, osiFile);

		if(offset == -2)
			{
				counter.increment("No index");
				cerr << "The search failed to find any index for the read " << al.Name << endl;
			}

		//-1 is code for a read that is a multi-mapper or fails to map
		if(offset > -1)
			{
				msInfo . addIndexingTagsToRead(al,refData[al.RefID].RefName,offset);
			}
		else
			{
				counter . increment("Multimap or no map");
			}
		
		//all reads, even mapping failures and multi-mappers, are printed back to the file
		
/*		if(msInfo.nMicrosatellites() > 0)
			msInfo . updateName(al);*/
		if(!outBam.SaveAlignment(al))
			cerr << counter.value("in") << ": Error writing alignment to " << MSI_OUTPUT_NAME << ": " << outBam.GetErrorString() << endl;
		counter.increment("out");
		msInfo.flush();

		inBam.GetNextAlignment(al);
		msInfo . plainQName(al.Name);
		counter.increment("in");

		hdrFile2 . getNextHdr();
		hdrFile2 . getHdr(currentHdr);

		if((currentHdr . compareQNames(al.Name)) && al.IsFirstMate() == currentHdr . isFirst())
			{
				currentHdr . getMicrosatellites (tMs);
				currentHdr . getMicrosatelliteQuals (tMsQ);
				msInfo . addMicrosatellites (tMs);
				msInfo . addMicrosatelliteQuals (tMsQ);
				msInfo . readLength(al);
			}
		else
			{				
				counter.increment("Mismatched header");
				counter.print(counter_output);
				fclose(counter_output);

				cerr << "Error!! The current header in " << BAM_FILE << " does not match the current header in " << HDR_FILE_2 << endl;
				cerr << "\"" << al.Name << "\" != \"" << currentHdr . getQName() << "\"" << endl;
				cerr << "Flags: " << al.AlignmentFlag << "\t";
				cerr << "\tHDR MATE: " << currentHdr . isFirst() << "\t" << al.IsFirstMate() << endl;
				return -1; 
			}
		offset = getOffset (al, refData[al.RefID].RefName, msInfo, osiFile);
		
		if(offset == -2)
			{
				counter.increment("No index");
				cerr << "The search failed to find any index for the read " << al.Name << endl;
			}

		if(offset > -1)
			{
				msInfo . addIndexingTagsToRead(al,refData[al.RefID].RefName,offset);
			}
		else
			{
				counter . increment("Multimap or no map");
			}
/*		if(msInfo.nMicrosatellites() > 0)
			msInfo . updateName(al);*/
		if(!outBam.SaveAlignment(al))
			cerr << counter.value("in") << ": Error writing alignment to " << MSI_OUTPUT_NAME << ": " << outBam.GetErrorString() << endl;
		counter.increment("out");

		if(counter.value("in") % 1000000 == 0)
			{
				cerr << "Processed read " << counter.value("in") << "... " << trans_timer . elapsed() << endl;
			}
		msInfo.flush();
	}
	
	if(counter.value("in") == 0)
	{
		fprintf(stderr,"%s is empty!\n",BAM_FILE.c_str());
	}
	
	inBam.Close();
	outBam.Close();
	counter.print();
	counter.print(counter_output);
	fclose(counter_output);
	cerr << "Read " << counter.value("in") << " reads" << endl;
	cerr << "Re-indexed SAM reads can be found in " << MSI_OUTPUT_NAME  << endl;
	cerr << "Reindex counts can be found in " << COUNTER_NAME << endl;
	cerr << "Time elapsed: " << trans_timer . elapsed() << endl;
	return 0;
 }
