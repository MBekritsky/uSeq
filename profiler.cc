/*
ABSTRACT:
Once reads with microsatellites have been aligned, 
reindexed and marked for duplicates, create a profile 
for every microsatellite observed in the individual.
Also reports information on PCR duplicate concordance
and overlapping read pair concordance.  This program
is a conversion of a perl script intended to provide
more rigorous and faster code

CREATION DATE:
20.02.2013

LAST REVISION:
20.02.2013

AUTHOR:
Mitchell Bekritsky
*/

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>

#include "api/BamWriter.h"
#include "api/BamReader.h" //BamTools is a C++ API for samtools.  The samtools API
//seems to be CPU-bound, while BamTools tries to maintain a low memory overhead, 
//and is therefore IO-bound. Using BamTools here is strictly for consistency
//with the API used in tetrascan

#include "MicrosatelliteDatabase.h"
#include "Microsatellite.h"
#include "MicrosatelliteQualityFlag.h"
#include "MicrosatelliteCoordFlag.h"
#include "MicrosatelliteProfile.h"

#include "UsefulMacros.h"
#include "Timer.h"
#include "FileOps.h"
#include "StringManip.h"
#include "Counter.h"

using namespace std;
using namespace BamTools;

struct threshold_t
{
	unsigned int minimumCount;
	unsigned int minMapQ;
	unsigned int flankLength;
};

void setDefaultThresholds(threshold_t & t)
{
  t.minimumCount = 2;
  t.minMapQ	 = 30;
  t.flankLength	 = 5;
}


int main(int argc, char * argv [])
{
	//takes a set of chromosome BAM files or one BAM file, a microsatellite database
	//a prefix for the file, and optionally sets thresholds
	
	string bamFileName, msdbName, prefix;
	threshold_t thresholds;
	bool keepDuplicates = false;
	setDefaultThresholds(thresholds);
	string counterSuffix = ".profile.count", dupSuffix = ".dups.txt",profileSuffix = ".msp";
	
	const static struct option longOpts[] =
	{
		{"bam",  	 			   required_argument, 0, 'b'},
		{"msdb", 	 			   required_argument, 0, 'd'},
		{"prefix", 			   optional_argument, 0, 'p'},
		{"count",	 			   optional_argument, 0, 'c'},
		{"mapQ", 	 			   optional_argument, 0, 'q'},
		{"flankLength",    optional_argument, 0, 'l'},
		{"keepDuplicates", optional_argument, 0, 'k'},
		{NULL,0,0,0}
	};
	
	int ch;
	
	
	Timer_t pTimer;
	Counter readTracker;
	
	while((ch = getopt_long_only(argc, argv, ":b:d:p:c:q:l:k", longOpts, NULL)) != -1)
	{
		switch(ch)
		{
			case 'b':		bamFileName							 = optarg; 				break;
			case 'd':		msdbName 								 = optarg;				break;
			case 'p':		prefix 	 								 = optarg;				break;
			case 'c':		thresholds.minimumCount  = atoi(optarg);	break;
			case 'q':		thresholds.minMapQ			 = atoi(optarg);	break;
			case 'l':		thresholds.flankLength	 = atoi(optarg);	break;
			case 'k':   keepDuplicates           = true;          break;
			case ':':
				fprintf(stderr,"Option -%c requires an argument\n",optopt);
				break;
			case '?':
			default:
				fprintf(stderr,"Option -%c is invalid and is being ignored\n",optopt);
				break;
		}
	}
	
	Timer_t dbLoadTimer;
	MicrosatelliteDatabase msdb(msdbName);
	cerr << "Loaded microsatellite database in " << dbLoadTimer.elapsed() << " seconds" << endl;
	
	BamReader msiReader;
	RefVector refData;
	BamAlignment al;
	string profileFilename,dupFilename, seqNoiseFilename, countFilename,counterDir;
	FILE * profile = NULL, * dup = NULL, * counts = NULL;
	BamWriter seqNoise;
	
	Counter counter;
	vector< MicrosatelliteCoordFlag > coordFlags;
	vector< MicrosatelliteQualityFlag > qualFlags;
	string motif, coord, keyName;
	unsigned int uncondensedPosition;
	map< int, map<string, MicrosatelliteProfile, KeySorter > > profileMap; //the KeySorter functor will sort the keys of profileMap according to coordinate, genome position and unit,
	//so that keys in profileMap will be stored in genome-sorted order
	map< string, MicrosatelliteProfile, KeySorter >::iterator chrMapIt;
	map< int, map<string, MicrosatelliteProfile, KeySorter > >::iterator	profileMapIt;
	
	if(prefix.empty())
	{
		prefix = bamFileName.substr(0,bamFileName.find_first_of("."));
	}

	if(!msiReader.Open(bamFileName))
	{
		cerr << "Could not open " << bamFileName << ": " << msiReader.GetErrorString() << endl;
		exit(EXIT_FAILURE);
	}

	profileFilename 		= addSuffix(prefix,profileSuffix);
	profile 						= openAndTestFile(profileFilename,"w");
	dupFilename					= addSuffix(prefix,dupSuffix);
	dup 								= openAndTestFile(dupFilename,"w");

	counterDir 					= createCounterDir(bamFileName,prefix,"");
	countFilename				= replacePath(addSuffix(prefix,counterSuffix),counterDir);
	counts							= openAndTestFile(countFilename,"w");
	

	seqNoiseFilename		= prefix + ".seq.noise.bam";
	if(!seqNoise.Open(seqNoiseFilename,msiReader.GetHeader(),msiReader.GetReferenceData()))
	{
		cerr << "Could not open " << seqNoiseFilename << ": " << seqNoise.GetErrorString() << endl;
	}

	string currCoord, lastCoord;
	int currRefID = 0, lastRefID = 0;

	cerr << "Processing " << bamFileName << endl;		
	refData = msiReader.GetReferenceData();
	
	while(msiReader.GetNextAlignment(al) /*&& counter.value("in") < 100*/)
	{
		counter.increment("in");
		if(!al.IsMapped())
		{
			counter.increment("unmapped");
		}
		else
		{
			currRefID = al.RefID;
			currCoord = refData[currRefID].RefName;
			coordFlags.clear();
			qualFlags.clear();
			fillMicrosatelliteCoordDetails(al,coordFlags);
			fillMicrosatelliteQualityDetails(al,qualFlags);

			al.GetTag("OP",uncondensedPosition);

			if(al.MapQuality < thresholds.minMapQ)
			{
				counter.increment("lowMapQ");
			}
			else if(al.MateRefID != al.RefID)
			{
				counter.increment("disparateCoord");
			}
			else if(coordFlags.size() > 0)
			{
				counter.increment("hasMS");
				vector< MicrosatelliteCoordFlag >::iterator it;
				for(unsigned int i = 0; i < qualFlags.size(); ++i)
				{
					coordFlags[i].getMotif(motif);
					//checks to make sure that the microsatellite has no Ns (should never happen)
					//and that it has no lowercase bases (which would mean it is at the end of a read and
					//has an indeterminate length
					if(motif.find_first_of("acgtnN") == string::npos) 
					{
						//checks that there are no non-ACGT letters in the motif
						if(motif.find_first_not_of("ACGT") == string::npos)
						{
							getProfileKeyName(currRefID,coordFlags[i],keyName);
							if(profileMap[currRefID].find(keyName) == profileMap[currRefID].end())//if this is the first read covering the microsatellite, create a new microsatellite profile instance
							{
								profileMap[currRefID][keyName] = MicrosatelliteProfile(coordFlags[i]);
							}

							if(!profileMap[currRefID][keyName].overlappingReadFragments(al,coordFlags[i],counter,seqNoise)) //if this is a unique microsatellite observation (ie, from a different read pair), add it to the profile
							{
								profileMap[currRefID][keyName].addRead(al,coordFlags[i],qualFlags[i]);
							}

						}
						else
						{
							counter.increment("badUnit");
						}
					}
				}
			}

			if(lastRefID != currRefID)
			{
				for(chrMapIt = profileMap[lastRefID].begin(); chrMapIt != profileMap[lastRefID].end(); ++chrMapIt)
				{
					if((*chrMapIt).second.count() < thresholds.minimumCount)
					{
						counter.increment("lowCoverage");
					}
					else
					{
						(*chrMapIt).second.print(msdb,thresholds.flankLength,profile,dup,counter, keepDuplicates);
					}
				}
				profileMap[lastRefID].clear();
			}
		}
		lastRefID = currRefID;
		lastCoord = currCoord;
	}
	
	//print last chromosome
	for(chrMapIt = profileMap[currRefID].begin(); chrMapIt != profileMap[currRefID].end(); ++chrMapIt)
	{
		if((*chrMapIt).second.count() < thresholds.minimumCount)
		{
			counter.increment("lowCoverage");
		}
		else
		{
			(*chrMapIt).second.print(msdb,thresholds.flankLength,profile,dup,counter, keepDuplicates);
		}
	}
	cout << "Printed all profiles" << endl;
	profileMap.clear();
	
	fclose(dup);
	fclose(profile);
	msiReader.Close();
	if(seqNoise.IsOpen())
	{
		seqNoise.Close();
	}

	counter.print();
	counter.print(counts);
	fclose(counts);
	
	cerr << "Microsatellite profiles can be found in " << profileFilename << endl;
	cerr << "Duplicate information can be found in " << dupFilename << endl;
	cerr << "Read pairs with discordant microsatellite tract lengths can be found in " << seqNoiseFilename << endl;
	cerr << "Count information can be found in " << countFilename << endl;
	cerr << "Time elapsed: " << pTimer.elapsed() << endl;

	return 0;
}
