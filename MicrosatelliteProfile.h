#ifndef __MICROSATELLITE_PROFILE_H__
#define __MICROSATELLITE_PROFILE_H__

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "MicrosatelliteQualityFlag.h"
#include "MicrosatelliteCoordFlag.h"
#include "MicrosatelliteDatabase.h"
#include "Microsatellite.h"
#include "UsefulMacros.h"
#include "Counter.h"
#include <vector>
#include <string>
#include <map>

using namespace std;
using namespace BamTools;

class MicrosatelliteProfile
{
	private:
		vector< MicrosatelliteCoordFlag > coordFlags_;
		vector< MicrosatelliteQualityFlag > qualFlags_;
		vector< BamAlignment > alignments_;
		vector< bool > reportable_;
		map< string, unsigned int > bareNames_;
		map< string, unsigned int >::iterator bnIt_;
		unsigned int count_;
		
		string coord_;
		string motif_;
		unsigned int genomeStart_,genomeStop_;
		string name_;
	public:
		MicrosatelliteProfile();
		MicrosatelliteProfile(const MicrosatelliteCoordFlag & coordFlag);
		MicrosatelliteProfile(const MicrosatelliteProfile & rhs);
		MicrosatelliteProfile & operator=(const MicrosatelliteProfile & rhs);
		
		void addRead(const BamAlignment & al, const MicrosatelliteCoordFlag & coordFlag, const MicrosatelliteQualityFlag & qualFlag);
		void bareName(const BamAlignment & al, string & bareName) const;
		string bareName(const BamAlignment & al) const;
		void pushBareName(const BamAlignment & al,unsigned int index);
		bool overlappingReadFragments(const BamAlignment & al, const MicrosatelliteCoordFlag & coordFlag, Counter & counter, BamWriter & seqNoise);
		
		void coordFlags(vector< MicrosatelliteCoordFlag > & coordFlags) const {coordFlags = coordFlags_;}
		void qualFlags(vector< MicrosatelliteQualityFlag > & qualFlags) const {qualFlags = qualFlags_;}
		void alignments(vector< BamAlignment > & alignments) const {alignments = alignments_;}
		void reportable(vector< bool > & reportable) const {reportable = reportable_;}
		unsigned int count() const {return count_;}
		unsigned int genomeStart() const {return genomeStart_;}
		unsigned int genomeStop()  const {return genomeStop_;}
		void motif(string & motif) const {motif = motif_;}
		void coord(string & coord) const {coord = coord_;}
		void name(string & name) const {name = name_;}
		string motif() const {return motif_;}
		string coord() const {return coord_;}
		string name() const {return name_;}

		string getDupKeyName(const BamAlignment & al) const;
		void getDupKeyName(const BamAlignment & al, string & key) const;
		void findReadsWithSameStarts(vector< unsigned int > & duplicateIndices) const;
		unsigned int getDupLengths(vector< unsigned int > & duplicateIndices, vector< unsigned int > & dupLengthHistogram) const;
		void pickRepresentativeDuplicate(const vector< unsigned int > & duplicateIndices, unsigned int numLengths);
		void getProfileHistMap(map< unsigned int, unsigned int > & histMap) const;
		void getSubprofileHistMap(const vector< unsigned int > & indices, map< unsigned int, unsigned int > & histMap) const;
		unsigned int getTractLengthMode() const;
		void printHistogram(vector< unsigned int > hist, FILE * fh, char delim) const;
		void printDuplicateInfo(const vector< unsigned int > & duplicateIndices, const vector < unsigned int > & dupLengthHistogram,
														unsigned int refLength, FILE * dups) const;
		void condenseDuplicates(unsigned int refLength, FILE * dups, Counter & counter);
		unsigned int numReportable() const;
		void getFlank(unsigned int index, string & flank, unsigned int flankLength) const;
		unsigned int populateFlankMap(unsigned int flankLength, map< string, vector< unsigned int > > & flankMap) const;
		void sortFlanks(const map< string, vector< unsigned int > > & flankMap, vector< string > & sortedFlanks) const;
		unsigned int printMicrosatelliteHist(FILE * out) const;
		void listQnames(FILE * out, char delim) const;
		unsigned int getProfileHistogram(vector< unsigned int > & profileLengthHistogram) const;
		unsigned int getSubprofileHistogram(const vector< unsigned int > & indices, vector< unsigned int > & subprofileLengthHistogram) const;
		void print(const MicrosatelliteDatabase & msdb, unsigned int flankLength, FILE * out, FILE * dups, Counter & counter, bool keepDuplicates);
		unsigned int getSubprofileHistogram(const vector< unsigned int > & indices, vector< unsigned int > & subprofileLengthHistogram, unsigned int maxLength) const;
		void fillHistVector(const map< unsigned int, unsigned int > & histMap, vector< unsigned int > & histVec, unsigned int maxLength) const;
		void barenames2Vec(vector < string > & bnVec) const;
		void listFlankQnames(const vector< unsigned int > & indices, const vector< string > & barenameVec, FILE * out, char delim) const;
		void printFlankInfo(unsigned int numFlanks, unsigned int rank, const string & flank, unsigned int flankPop, FILE * out) const;
		void printMicrosatelliteInfo(const Microsatellite & ref, FILE * out) const;
		unsigned int printMicrosatelliteSubpopHist(const vector< unsigned int > & indices, FILE * out) const;
		unsigned int printMicrosatelliteSubpopHist(const vector< unsigned int > & indices, FILE * out, unsigned int maxLength) const;

};

class KeySorter
{
	public:
		bool operator()(const string & a, const string & b) const;
};

void getProfileKeyName(int refID, const MicrosatelliteCoordFlag & coord, string & name);
string getProfileKeyName(int refID, const MicrosatelliteCoordFlag & coord);
unsigned int extractStart(const string & key);
string extractMotif(const string & key);
int extractRefID(const string & key);

unsigned int getCoordInd(const string & coord);
unsigned int getCoordInd(const string & coord, unsigned int maxAutosome);

//condensed sequence reconstruction
void reconstructSequence(const BamAlignment & al, const vector< MicrosatelliteQualityFlag > & qualFlags, string & reconstructedSequence);
void appendMicrosatellite(string & sequence,const MicrosatelliteQualityFlag & details);
void appendOverlappedMicrosatellite(string & sequence, const MicrosatelliteQualityFlag details, unsigned int lastFullStop);

void fillMicrosatelliteQualityDetails(const BamAlignment & al, map< unsigned int, MicrosatelliteQualityFlag > & msDetails);
void fillMicrosatelliteQualityDetails(const BamAlignment & al, vector< MicrosatelliteQualityFlag > & msDetails);
void qualVectorToMap(const vector< MicrosatelliteQualityFlag > & qualVector, map< unsigned int, MicrosatelliteQualityFlag > & qualMap, const BamAlignment & al);
void fillMicrosatelliteCoordDetails(const BamAlignment & al, vector< MicrosatelliteCoordFlag > & msDetails);

#endif /* MicrosatelliteProfile.h */
