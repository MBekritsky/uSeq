/*
ProfileLocus.h

ABSTRACT:
A class that holds information about one microsatellite locus
from a microsatellite profile

CREATION DATE:
12.23.2013

LAST REVISION:
12.23.2013

AUTHOR:
Mitchell Bekritsky
*/

#ifndef __PROFILE_LOCUS_H__
#define __PROFILE_LOCUS_H__

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include "StringManip.h"
#include "UsefulMacros.h"

using namespace std;

class ProfileLocus
{
	private:
	 int chr_;
	 int pos_;
	 string motif_;
	 int motifLength_;
	 int refLength_;
	 unsigned int numFlanks_;
	 unsigned int totalCount_;
	 unsigned int maxLength_;
	 
	 vector< string > flankSeq_;
	 vector< int > flankCount_;
	 vector< vector< int > > flankHist_;
	 
	public:
	 ProfileLocus();
	 ProfileLocus(const string & record, const map< string, int > & chr2int);
/*	 ProfileLocus(unsigned int chr, unsigned int pos, const string & motif, 
	              unsigned int refLength,unsigned int numFlanks, const string & flankSeq, unsigned int totalCount,
								const string & hist, map< string, unsigned int > & chr2int);*/
	 ProfileLocus(const ProfileLocus & rhs);
	 ProfileLocus & operator=(const ProfileLocus & rhs);
	 
	 string motif()		 const {return motif_;}
	 int chr()			 const {return chr_;}
	 int pos() const {return pos_;}
	 int motifLength() const {return motifLength_;}
	 int refLength() const {return refLength_;}
	 unsigned int numFlanks() const {return numFlanks_;}
	 unsigned int totalCount() const {return totalCount_;}
	 unsigned int maxLength() const {return maxLength_;}
	 
	 void getMotif(string & motif) const {motif = motif_;}
	 void getFlankSeq(vector< string > & flankSeq) const {flankSeq = flankSeq_;}
	 void getFlankCount(vector< int > & flankCount) const {flankCount = flankCount_;}
	 void getFlankHist(vector< vector< int > > & flankHist) const {flankHist = flankHist_;}

	 void getFlankHist(vector< int > & flankHistInstance, unsigned int index) const;
	 
	 void clear();
	 void parseRecord(const string & record, const map< string, int > & chr2int);
	 void parseHist(const string & rawHist, vector< int > & intHist) const;
	 
	 void newLocus(const string & record, const map< string, int > & chr2int);
	 void addFlankToLocus(const string & record, const map< string, int > & chr2int);
	 
	 void printLocus(const map< int, string > & int2chr) const;
	 void printLocus(FILE * fh, const map< int, string > & int2chr) const;
	 
	 void histSummary(string & hist, unsigned int flankIndex) const;
	 void histSummary(FILE * fh, unsigned int flankIndex) const;
};

#endif /*ProfileLocus.h*/
