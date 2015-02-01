/*
PopSipper.h

ABSTRACT:
A class that holds a LocusSipper for each person in a population
Synchronizes loci across the population

CREATION DATE:
12.25.2013

LAST REVISION:
12.25.2013

AUTHOR:
Mitchell Bekritsky
*/

#ifndef __POP_SIPPER_H__
#define __POP_SIPPER_H__

#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <numeric>

#include "UsefulMacros.h"
#include "ProfileSipper.h"
#include "ProfileLocus.h"
#include "PopLocus.h"

using namespace std;

class PopSipper
{
  private:
	vector< ProfileSipper > sippers_;
	unsigned int popSize_;
	vector< map< string, ProfileLocus > > loci_;
	PopLocus currLocus_;
	map< string, int > chr2int_;

	vector< bool > active_;
	vector< bool > toUpdate_;
	vector< int > popChr_;
	vector< int > popPos_;
	vector< string > popMotif_;
	int currChr_;
	int currPos_;
	int popAtNextLocus_;
	int numActive_;
	
  public:
    PopSipper();
	PopSipper(const map< string, int > & chr2int);
	PopSipper(const vector< ProfileSipper > & sippers, const map< string, int > & chr2int);
	PopSipper(const PopSipper & rhs);
	PopSipper & operator=(const PopSipper & rhs);
	
	void getSippers(vector< ProfileSipper > & sippers) const {sippers = sippers_;}
	void getChr2Int(map< string, int > & chr2int) const {chr2int = chr2int_;}
	void getCurrPopLocus(PopLocus & currLocus) const {currLocus = currLocus_;}

	unsigned int popSize() const {return popSize_;}
	unsigned int numActive() const {return numActive_;}

	void setSippers(const vector< ProfileSipper > & sippers);
	void setChr2Int(const map< string, int > & chr2int) {chr2int_ = chr2int;}
	
	void clear();
	void initializeVectors();
	bool updateProfileLoci();
	void advanceUpdatableSippers();
	void setLatestPopLocus();
	bool nextLocus();
};

#endif /* PopSipper.h */
