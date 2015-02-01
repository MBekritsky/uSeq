#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <dirent.h>
#include "FileOps.h"
#include "UsefulMacros.h"

#include "ProfileSipper.h"
#include "ProfileLocus.h"
#include "PopSipper.h"
#include "PopLocus.h"
#include "Timer.h"

using namespace std;

bool cstrStartsWith(char * cstr, char * prefix) {
	size_t cstrLen = strlen(cstr), preLen = strlen(prefix);
	if(cstrLen < preLen) return false;
	return strncmp(cstr,prefix,preLen) == 0;
}

//tests if cstr is either "." or ".."
bool isDots(char * cstr) {
	if(strcmp(cstr,".") == 0 || strcmp(cstr,"..") == 0) return true;
	return false;
}

string getRelationship(const string & relFilename) {
	FILE * relFileP = openAndTestFile(relFilename,"r");
	char ch = ' ';
	string rel;
	while(ch != EOF)
	{
		while((ch = fgetc(relFileP)) != EOF && ch != '\n')
			rel.push_back(ch);
		if((ch = fgetc(relFileP)) != EOF)
			ungetc(ch, relFileP);
	}
	fclose(relFileP);
	return rel;
}

struct FamilySorter{
	bool operator ()(const ProfileSipper & a, const ProfileSipper & b) const {
		return a.relCode() < b.relCode();
	}
};

int main(int argc, char * argv[]) {

	string popDirname = argv[1], famDirname,personDirname;

	string mspFilename, mspGzipFilename, relFilename;
	string familyID, sampleID, relation;
	unsigned int popSize = 0;
	
	map< string, int > chr2int;
	map< int, string > int2chr;
	int chrCounter;
	int intCounter;
	
	char key[1024];
	for(chrCounter = 1, intCounter = 0; chrCounter <= 22; ++chrCounter, ++intCounter) {
		sprintf(key,"chr%d",chrCounter);
		chr2int[key] = intCounter;
	}
	chr2int["chrX"]  = 22;
	chr2int["chrY"]  = 23;
	chr2int["chrMT"] = 24;
	
	DIR * popDirP = NULL, * famDirP = NULL;
	struct dirent * popPent = NULL, * famPent = NULL;

	vector< ProfileSipper > profiles;
	map< string,string > sample2rel;
	map< string,string > sample2fam;
	map< string,vector< string > > fam2sample;
	vector< string > familyIDs;

	vector< ProfileSipper >::iterator vpIt;

	addTrailingSlash(popDirname);

	string popIndexFilename = popDirname + "microsat.refOnly.columns.txt";
	FILE * popIndex = openAndTestFile(popIndexFilename,"w");

	popDirP = opendir(popDirname.c_str());
	if(popDirP == NULL) {
		fprintf(stderr,"Error! Failed to open \"%s\"! Exiting...\n",popDirname.c_str());
		exit(EXIT_FAILURE);
	}
	
	while((popPent = readdir(popDirP)) != NULL) {
		if(isDots(popPent->d_name)) continue;
		if(popPent == NULL) {
		//dirent error handling 
			fprintf(stderr,"Error! Failed to load dirent struct from %s! Skipping to next entry...\n",popDirname.c_str());
		}
		else {
			if(cstrStartsWith(popPent->d_name,"auSSC")) {
				vector< ProfileSipper > currFam;
				familyID = popPent->d_name;
				familyIDs.push_back(familyID);
				famDirname = popDirname + familyID + "/";
				famDirP = opendir(famDirname.c_str());
				if(famDirP == NULL) {
					fprintf(stderr,"Error!  Failed to open %s! Exiting...\n",famDirname.c_str());
					exit(EXIT_FAILURE);
				}
				while((famPent = readdir(famDirP)) != NULL) {
					if(isDots(famPent->d_name)) continue;
					if(famPent == NULL) {
						fprintf(stderr,"Error! Failed to load dirent struct from %s! Skipping to next entry...\n",famDirname.c_str());
					}
					if(cstrStartsWith(famPent->d_name,"SSC")) {
						sampleID = famPent->d_name;
						sample2fam[sampleID] = familyID;
						fam2sample[familyID].push_back(sampleID);
						
						personDirname = famDirname + sampleID + "/";
						mspGzipFilename = personDirname + sampleID + ".msp.gz";
						mspFilename = personDirname + sampleID + ".msp";
						relFilename = personDirname + "relationship.txt";
						
						if(fileExists(relFilename)) {
						  relation = getRelationship(relFilename);
							sample2rel[sampleID] = relation;
						}
						if(fileExists(mspGzipFilename)) {
							currFam.push_back(ProfileSipper(familyID,sampleID,relation,chr2int,mspGzipFilename));
						}
						else if(fileExists(mspFilename)) {
							currFam.push_back(ProfileSipper(familyID,sampleID,relation,chr2int,mspFilename));
						}
						else {
							fprintf(stderr,"No profile exists for %s in family %s\n",sampleID.c_str(),familyID.c_str());
						}
					}
				}
				
				vector< ProfileSipper > temp = currFam;
				// sorts families in the order of mother, father, self, sibling
				sort(currFam.begin(), currFam.end(), FamilySorter());
				
				tr(currFam, vpIt)
					profiles.push_back(*vpIt);					
				closedir(famDirP);
			}
		}
	}
	closedir(popDirP);
	
	popSize = profiles.size();
	
	profiles[0].getInt2Chr(int2chr);
	
	tr(profiles,vpIt)
		fprintf(popIndex,"%s\t%s\t%s\n", vpIt->sampleID().c_str(), vpIt->familyID().c_str(), vpIt->relation().c_str());
	cerr << "Wrote population information to " << popIndexFilename << endl;
	cerr << "Population size: " << popSize << endl;
	fclose(popIndex);
	
	PopSipper popSipper(profiles,chr2int);
	
	string locusIndexFilename = popDirname + "microsat.refOnly.index.txt";
	FILE * locusIndex = openAndTestFile(locusIndexFilename,"w");

	string countsFilename = popDirname + "microsat.refOnly.counts.txt";
	FILE * counts = openAndTestFile(countsFilename,"w");

	
	unsigned int total = 0, skipped = 0, added = 0; 
	PopLocus popLocus;
	
	cerr << "Creating summary table" << endl;
	Timer_t syncTimer;
	while(popSipper.nextLocus()) {
		total++;
		popSipper.getCurrPopLocus(popLocus);
		if(popLocus.refLength() > 0) {
			popLocus.printHistInfo(counts);
			popLocus.printLocusInfo(locusIndex, int2chr);
			added++;
		}
		else {
			skipped++;
		}
		
		if(total % 100000 == 0) {
			fprintf(stderr, "Processed locus %u00,000, added %u.  Time elapsed: %s.",
			        total/100000U, added, syncTimer.elapsed().c_str());
			fprintf(stderr, " Average time per 10,000 loci: %.2f sec\n",
			         (syncTimer.length() / (double) (total/10000U)) / 100000.0);
		}
	}
	
	fclose(counts);
	fclose(locusIndex);
	cerr << "Completed summary table and index.  Time elapsed: ";
	cerr << syncTimer.elapsed() << endl;
	
	cerr << "Wrote summary index to " << locusIndexFilename << endl;
	cerr << "Wrote counts to " << countsFilename << endl;
	cerr << total << " loci were processed, " << added << " were added";
	cerr << " to the summary" << endl;
	return true;
}
