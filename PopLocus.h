/*
PopLocus.h

ABSTRACT:
A class that holds locus information for a PopSipper instance

CREATION DATE:
12.25.2013

LAST REVISION:
12.25.2013

AUTHOR:
Mitchell Bekritsky
*/

#ifndef __POP_LOCUS_H__
#define __POP_LOCUS_H__

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>

#include "UsefulMacros.h"
#include "ProfileLocus.h"

using namespace std;

class PopLocus
{
  private:
    map< string, vector< unsigned int > > sipperMotif_;
	map< string, vector< bool > > sipperMotifBool_;
	vector< map< string, ProfileLocus > > sipperLoci_;
	ProfileLocus firstLocus_;
	
	int chr_;
	int pos_;
	vector< string > motifs_;
	int refLength_;
	vector< unsigned int > whoHasLocusInds_;
	map< string, vector< unsigned int > > totalCountPerPersonByMotif_;
	map< string, unsigned int > maxLength_;
	map< string, unsigned int > maxCount_;
	int popCount_;
	unsigned int popSize_;
  public:
    PopLocus();
	PopLocus(const vector< map< string, ProfileLocus > > & sipperLoci, const vector< bool > & whoHas);
	PopLocus(const PopLocus & rhs);
	PopLocus & operator=(const PopLocus & rhs);
	
	void getSipperLoci(vector< map< string, ProfileLocus > > & sipperLoci) const {sipperLoci = sipperLoci_;}
	void getSipperMotif(map< string, vector< unsigned int > > & sipperMotif) const {sipperMotif = sipperMotif_;}
	void getSipperMotifBool(map< string, vector< bool > > & sipperMotifBool) const {sipperMotifBool = sipperMotifBool_;}
	void getFirstLocus(ProfileLocus & firstLocus) const {firstLocus = firstLocus_;}
	void getMotifs(vector< string > & motifs) const {motifs = motifs_;}
	void getWhoHasLocusInds(vector< unsigned int > & whoHasLocusInds) const {whoHasLocusInds = whoHasLocusInds_;}
	void getTotalCountPerPersonByMotif(map< string, vector< unsigned int > > & totalCountPerPersonByMotif) const {totalCountPerPersonByMotif = totalCountPerPersonByMotif_;}
	void getMaxLength(map< string, unsigned int > & maxLength) const {maxLength = maxLength_;}
	void getMaxCount(map< string, unsigned int > & maxCount) const {maxCount = maxCount_;}
	
	int chr() const {return chr_;}
	int pos() const {return pos_;}
	int refLength() const {return refLength_;}
	int popCount()  const {return popCount_;}
	unsigned int popSize() const {return popSize_;}
	
	void refresh(const vector< map< string, ProfileLocus > > & sipperLoci, const vector< bool > & whoHas);
	void clear();

	void printLocusInfo(FILE * fh, const map< int, string> & int2chr) const;
	void printHistInfo(FILE * fh) const;
};

#endif /* PopLocus.h */
