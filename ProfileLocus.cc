/*
ProfileLocus.cc

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

#include "ProfileLocus.h"
using namespace std;

ProfileLocus::ProfileLocus()
{
	clear();
}

ProfileLocus::ProfileLocus(const string & record, const map< string, int > & chr2int)
{
	parseRecord(record, chr2int);
}

/*
ProfileLocus::ProfileLocus(unsigned int chr, unsigned int pos, const string & motif, 
	              unsigned int refLength,unsigned int numFlanks, const string & flankSeq, unsigned int totalCount,
								const string & hist, map< strin, unsigned int > & chr2int)
{
	chr_ = chr;
	pos_ = pos;
	motif_ = motif;
	motifLength_ = motif.size();
	refLength_ = refLength;
	numFlanks_ = numFlanks;
	
}*/

ProfileLocus::ProfileLocus(const ProfileLocus & rhs)
{
	rhs.getMotif(motif_);
	rhs.getFlankSeq(flankSeq_);
	rhs.getFlankCount(flankCount_);
	rhs.getFlankHist(flankHist_);
	
	chr_         = rhs.chr();
	pos_         = rhs.pos();
	motifLength_ = rhs.motifLength();
	refLength_   = rhs.refLength();
	numFlanks_   = rhs.numFlanks();
	totalCount_  = rhs.totalCount();
	maxLength_   = rhs.maxLength();
}

ProfileLocus & ProfileLocus::operator=(const ProfileLocus & rhs)
{
	if(this != &rhs)
	{
 	  rhs.getMotif(motif_);
      rhs.getFlankSeq(flankSeq_);
  	  rhs.getFlankCount(flankCount_);
  	  rhs.getFlankHist(flankHist_);
	
  	  chr_         = rhs.chr();
  	  pos_         = rhs.pos();
  	  motifLength_ = rhs.motifLength();
	  refLength_   = rhs.refLength();
	  numFlanks_   = rhs.numFlanks();
	  totalCount_  = rhs.totalCount();
	  maxLength_   = rhs.maxLength();
	}
	return *this;
}

void ProfileLocus::getFlankHist(vector< int > & flankHistInstance, unsigned int index) const
{
	if(index < flankHist_.size())
	{
		flankHistInstance = flankHist_[index];
	}
	else
	{
		cerr << "ERROR!  Attempt to access flank with incorrect index: " << index << " > ";
		cerr << flankHist_.size() << endl;
		exit(EXIT_FAILURE);
	}
}

void ProfileLocus::clear()
{
	motif_.clear();
	
	chr_ = 0;
	pos_ = 0;
	motifLength_ = 0;
	refLength_ = 0;
	numFlanks_ = 0;
	totalCount_ = 0;
	maxLength_  = 0;
	
	flankSeq_ . clear();
	flankCount_ . clear();
	flankHist_ . clear();
}

void ProfileLocus::parseRecord(const string & record, const map< string, int > & chr2int)
{
	vector< string > fields;
	splitString(record,'\t',fields);
	chr_ = chr2int.find(fields[0])->second;
	pos_ = atoi(fields[1].c_str());
	motif_ = fields[2];
	motifLength_ = motif_.size();
	refLength_ = atoi(fields[3].c_str());
	numFlanks_ = atoi(fields[4].c_str());
	totalCount_ = atoi(fields[7].c_str());
	
	flankSeq_.push_back(fields[6]);
	flankCount_.push_back(atoi(fields[7].c_str()));
	
	string rawHist = fields[8];
	vector< int > intHist;
	parseHist(rawHist,intHist);
	flankHist_.push_back(intHist);
	maxLength_ = flankHist_[0].size();
}

void ProfileLocus::parseHist(const string & rawHist, vector< int > & intHist) const
{
  	vector< string > parsedStringHist;
	splitString(rawHist,',',parsedStringHist);
	
	intHist.clear();
	vector< string >::iterator vsIt;
	tr(parsedStringHist, vsIt)
		intHist.push_back(atoi(vsIt->c_str()));
}

void ProfileLocus::newLocus(const string & record, const map< string, int > & chr2int)
{
	parseRecord(record,chr2int);
}

void ProfileLocus::addFlankToLocus(const string & record, const map< string, int > & chr2int)
{
	vector< string > fields;
	splitString(record,'\t',fields);
		
	
	if(chr2int.find(fields[0])->second == chr_ && atoi(fields[1].c_str()) == pos_ && fields[2] == motif_)
	{
	  string flank = fields[6];
	  flank = flank.substr(0,5) + " " + flank.substr(5,string::npos);
	  flankSeq_.push_back(flank);
	  flankCount_.push_back(atoi(fields[7].c_str()));

      string rawHist = fields[8];
	  vector< int > intHist;
	  parseHist(rawHist, intHist);
	  flankHist_.push_back(intHist);
	 }
	 else
	 {
	 	cerr << "Cannot add flank with different chromosome/position/motif to locus" << endl;
		exit(EXIT_FAILURE);
	 }
}

void ProfileLocus::printLocus(FILE * fh, const map< int, string > & int2chr) const
{
	char titleString[2056];
	sprintf(titleString,"%s\t%d\t%s\t%d\t%d\t%d",int2chr.find(chr_)->second.c_str(), pos_, motif_.c_str(),
													 refLength_, numFlanks_, totalCount_);

	char histString[2056];
	vector< int >::iterator vIt;
	for(unsigned int i = 0; i <= (unsigned int) numFlanks_; ++i)
	{
		histString[0] = '\0';
		tr(flankHist_[i],vIt)
		{
			sprintf(histString,"%s;%d", histString, *vIt);
		}
		fprintf(fh,"%s\t%s\t%s\n",titleString,flankSeq_[i].c_str(),histString);
	}
}

void ProfileLocus::printLocus(const map< int, string > & int2chr) const
{
	printLocus(stdout,int2chr);
}

void ProfileLocus::histSummary(string & hist, unsigned int flankIndex) const
{
  if(flankIndex < flankHist_.size())
  {
    vector< int > counts;
    vector< unsigned int > lengths;
	
	for(unsigned int i = 0; i < flankHist_[flankIndex].size(); ++i)
	{
	  if(flankHist_[flankIndex][i] > 0)
	  {
	    counts.push_back(flankHist_[flankIndex][i]);
		lengths.push_back(i + 1);
	  }
	}
	
    char countBuff[2048];
    char lengthBuff[2048];
	for(unsigned int i = 0; i < (lengths.size() - 1); ++i)
	{
	  sprintf(countBuff,"%s%d,",countBuff,counts[i]);
	  sprintf(lengthBuff,"%s%d,",lengthBuff,lengths[i]);
	}
	sprintf(countBuff,"%s%d",countBuff, (*counts.end()));
    sprintf(lengthBuff,"%s%d",lengthBuff, (*lengths.end()));
	
	char histBuff[4128];
	sprintf(histBuff,"%d;%s;%s",(int) counts.size(),lengthBuff,countBuff);
	hist = histBuff;
  }
  else
  {
    cerr << "Flank index exceeds number of flanks" << endl;
	exit(EXIT_FAILURE);
  }
}

void ProfileLocus::histSummary(FILE * fh, unsigned int flankIndex) const
{
  if(flankIndex < flankHist_.size())
  {
    vector< int > counts;
    vector< unsigned int > lengths;
	
	for(unsigned int i = 0; i < flankHist_[flankIndex].size(); ++i)
	{
	  if(flankHist_[flankIndex][i] > 0)
	  {
	    counts.push_back(flankHist_[flankIndex][i]);
		lengths.push_back(i + 1);
	  }
	}
		
    char countBuff[2048] = {0};
    char lengthBuff[2048] = {0};
	if(lengths.size() > 1) {
	  for(unsigned int i = 0; i < (lengths.size() - 1); ++i)
	  {
	    sprintf(countBuff,"%s%d,",countBuff,counts[i]);
	    sprintf(lengthBuff,"%s%d,",lengthBuff,lengths[i]);
	  }
	}
	sprintf(countBuff,"%s%d",countBuff, (*counts.rbegin()));
    sprintf(lengthBuff,"%s%d",lengthBuff, (*lengths.rbegin()));
	
	fprintf(fh,"%d;%s;%s",(int) counts.size(),lengthBuff,countBuff);
  }
  else
  {
    cerr << "Flank index exceeds number of flanks" << endl;
	exit(EXIT_FAILURE);
  }
}
