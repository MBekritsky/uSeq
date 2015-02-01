/*
ProfileSipper.h

ABSTRACT:
A class that takes a microsatellite profile as input,
and reads it line by line.  This is a class to be used
when creating summaries of microsatellite profiles
within a population

CREATION DATE:
12.23.2013

LAST REVISION:
12.23.2013

AUTHOR:
Mitchell Bekritsky
*/

#ifndef __PROFILE_SIPPER_H__
#define __PROFILE_SIPPER_H__

#include <string>
#include "zlib.h"
#include "ProfileLocus.h"

using namespace std;

class ProfileSipper
{
 private:
  string sampleID_;
  string familyID_;
	string filename_;
	string relation_;
	string currLine_;
	unsigned int relCode_;
	gzFile fh_;
	map< string, ProfileLocus > currLocus_; //separate ProfileLocus for different motifs at same locus
	map < string, int > chr2int_;
	map < int, string > int2chr_;

	unsigned int setRelCode();
	
 public:
  ProfileSipper();
	//Passing the filename is the recommended way to initialize an instance of ProfileSipper
	ProfileSipper(const string & familyID, const string & sampleID, const string & relation,
	              const map< string, int > & chr2int, const string & filename);
	//DO NOT USE THIS CONSTRUCTOR IF YOU NEED TO COPY FILES, COPY CONSTRUCTOR WITHOUT A FILENAME
	//IS NOT IMPLEMENTED
	ProfileSipper(const string & familyID, const string & sampleID, const string & relation,
	              const map< string, int > & chr2int, gzFile * fh);
								
	//the copy and operator overload constructors will restart at the beginnging of the file
	ProfileSipper(const ProfileSipper & rhs);
	ProfileSipper & operator=(const ProfileSipper & rhs);
	
	~ProfileSipper();
	
	//get functions
	string familyID() const {return familyID_;}
	string sampleID() const {return sampleID_;}
	string relation() const {return relation_;}
	string filename() const {return filename_;}
	unsigned int relCode() const {return relCode_;}
//	FILE * fh() const {return fh_;}
	void getFamilyID(string & familyID) const {familyID = familyID_;}
	void getSampleID(string & sampleID) const {sampleID = sampleID_;}
	void getRelation(string & relation) const {relation = relation_;}
	void getFilename(string & filename) const {filename = filename_;}
	void getChr2Int(map< string, int > & chr2int) const {chr2int = chr2int_;}
	void getInt2Chr(map< int, string > & int2chr) const {int2chr = int2chr_;}
	void getCurrLocus(map< string, ProfileLocus > & currLocus) const {currLocus = currLocus_;}
	void getCurrLocusByMotif(const string & motif, ProfileLocus & currLocus);
	
	//set functions
	void familyID(const string & familyID) {familyID_ = familyID;}
	void sampleID(const string & sampleID) {sampleID_ = sampleID;}
	void relation(const string & relation) {relation_ = relation;}
	void filename(const string & filename) {filename_ = filename;}
	void setChr2Int(const map< string, int > & chr2int);
//	void fh(gzGile * fh) {fh_ = fh}; //SHOULD NOT BE USED UNLESS YOU'RE ABSOLUTELY SURE YOU KNOW WHAT YOU'RE DOING!

	void complementChr2Int();
	bool nextLine(string & record);
	bool nextLocus();
	void printLocus(FILE * out) const;
	void printLocus() const;
};

#endif /*ProfileSipper.h*/

