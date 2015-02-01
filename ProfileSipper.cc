/*
ProfileSipper.cc

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
#include "ProfileSipper.h"
using namespace std;

ProfileSipper::ProfileSipper()
{
	sampleID_.clear();
	familyID_.clear();
	filename_.clear();
	relation_.clear();
	currLine_.clear();
	currLocus_.clear();
	chr2int_.clear();
	int2chr_.clear();
	relCode_ = 4;
	fh_ = NULL;
}

ProfileSipper::ProfileSipper(const string & familyID, const string & sampleID, const string & relation,
	              const map< string, int > & chr2int, const string & filename)
{
	sampleID_ = sampleID;
	familyID_ = familyID;
	relation_ = relation;
	filename_ = filename;
	relCode_ = setRelCode();
	chr2int_ = chr2int;
	complementChr2Int();
	
	fh_ = gzopen(filename_.c_str(),"rb");
	
	if(!nextLine(currLine_)) // end of file has been reached
	{
	  cerr << "File is empty!" << endl;
	}
}

ProfileSipper::ProfileSipper(const string & familyID, const string & sampleID, const string & relation,
	              const map< string, int > & chr2int, gzFile * fh)
{
	sampleID_ = sampleID;
	familyID_ = familyID;
	relation_ = relation;
	relCode_ = setRelCode();
	chr2int_ = chr2int;
	complementChr2Int();
	
	filename_.clear();
	
	fh_ = fh;
	if(!nextLine(currLine_)) // end of file has been reached
	{
	  cerr << "File is empty!" << endl;
	}
}

//copy constructor opens to beginning of file
ProfileSipper::ProfileSipper(const ProfileSipper & rhs)
{
	rhs.getSampleID(sampleID_);
	rhs.getFamilyID(familyID_);
	rhs.getRelation(relation_);
	rhs.getFilename(filename_);
	rhs.getChr2Int(chr2int_);
	rhs.getInt2Chr(int2chr_);
	relCode_ = rhs.relCode();

	currLocus_.clear();
	
	if(!filename_.empty())
	{
		fh_ = gzopen(filename_.c_str(),"rb");
	} else
	{
		fprintf(stderr,"ERROR!  Copy constructor for ProfileSipper requires a filename\n");
		exit(EXIT_FAILURE);
	}
	if(!nextLine(currLine_)) // end of file has been reached
	{
	  cerr << "File is empty!" << endl;
	}
}

ProfileSipper & ProfileSipper::operator=(const ProfileSipper & rhs)
{
	if(this != &rhs)
	{
		rhs.getSampleID(sampleID_);
		rhs.getFamilyID(familyID_);
		rhs.getRelation(relation_);
		rhs.getFilename(filename_);
		rhs.getChr2Int(chr2int_);
		rhs.getInt2Chr(int2chr_);
		relCode_ = rhs.relCode();

		currLocus_.clear();
	
		if(!filename_.empty())
		{
			fh_ = gzopen(filename_.c_str(),"rb");
		} else
		{
			fprintf(stderr,"ERROR!  Copy constructor for ProfileSipper requires a filename\n");
			exit(EXIT_FAILURE);
		}
		if(!nextLine(currLine_)) // end of file has been reached
		{
		  cerr << "File is empty!" << endl;
		}
	}
	return *this;
}

ProfileSipper::~ProfileSipper()
{
	gzclose(fh_);
}

void ProfileSipper::setChr2Int(const map< string, int > & chr2int)
{
	chr2int_ = chr2int;
	complementChr2Int();
}

void ProfileSipper::complementChr2Int()
{
	map< string, int >::iterator mIt;
	
	tr(chr2int_, mIt)
		int2chr_[mIt->second] = mIt->first;
}

unsigned int ProfileSipper::setRelCode()
{
	if(relation_ == "mother")
	{
		return 0;
	}
	if(relation_ == "father")
	{
		return 1;
	}
	if(relation_ == "proband" || relation_ == "self")
	{
		return 2;
	}
	if(relation_ == "sibling")
	{
		return 3;
	}
	return 4;
}

void ProfileSipper::getCurrLocusByMotif(const string & motif, ProfileLocus & currLocus)
{
  if(currLocus_.find(motif) != currLocus_.end())
  {
  	currLocus = currLocus_[motif];
  }
  else
  {
  	cerr << "Motif " << motif << "could not be found at current locus" << endl;
  }
}


//read in next line from profile file, set 
bool ProfileSipper::nextLine(string & record)
{
	record.clear();
	char ch = ' ';
	
	while((ch = gzgetc(fh_)) != EOF && ch != '\n')
		record.push_back(ch);
	if(ch == EOF)
		return false;

	return true;
}


bool ProfileSipper::nextLocus()
{
	currLocus_.clear();
	
	ProfileLocus temp(currLine_,chr2int_);
	string currMotif;
	temp.getMotif(currMotif);
	
	currLocus_[currMotif] = ProfileLocus(temp);		

  	int numFlanks = currLocus_[currMotif].numFlanks();
	for(int i = 0; i < numFlanks; ++i)
	{
		currLine_.clear();
		nextLine(currLine_);
		currLocus_[currMotif].addFlankToLocus(currLine_,chr2int_);
	}
	
	int currChr = currLocus_[currMotif].chr();
	int currPos = currLocus_[currMotif].pos();

	//Get next entry in msp file
	//If next entry has same chr and pos, microsatellite
	//is at same locus, but with different motif
	if(!nextLine(currLine_))
		return false;

	temp.clear();
	temp.newLocus(currLine_,chr2int_);
	int nextChr = temp.chr();
	int nextPos = temp.pos();
	string nextMotif;
	temp.getMotif(nextMotif);
	

	while(currChr == nextChr && currPos == nextPos && currMotif != nextMotif)
	{
		currLocus_[nextMotif] = ProfileLocus(temp);
	  	int numFlanks = currLocus_[nextMotif].numFlanks();
		for(int i = 0; i < numFlanks; ++i)
		{
			currLine_.clear();
			nextLine(currLine_);
			currLocus_[nextMotif].addFlankToLocus(currLine_,chr2int_);
		}
		
		if(!nextLine(currLine_))
			return false;
		
		temp.clear();
		temp.newLocus(currLine_, chr2int_);
		nextChr = temp.chr();
		nextPos = temp.pos();
		temp.getMotif(nextMotif);
	}
		
	return true;
}

void ProfileSipper::printLocus(FILE * out) const
{
	map<string, ProfileLocus >::const_iterator spIt;
	tr(currLocus_,spIt)
		spIt->second.printLocus(out,int2chr_);
}


void ProfileSipper::printLocus() const
{
	printLocus(stdout);
}
