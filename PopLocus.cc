/*
PopLocus.cc

ABSTRACT:
A class that holds locus information for a PopSipper instance

CREATION DATE:
12.25.2013

LAST REVISION:
12.26.2013

AUTHOR:
Mitchell Bekritsky
*/

#include "PopLocus.h"
using namespace std;

PopLocus::PopLocus()
{
  clear();
}

PopLocus::PopLocus(const vector< map< string, ProfileLocus > > & sipperLoci, const vector< bool > & whoHas)
{
  popSize_ = sipperLoci.size();
  refresh(sipperLoci, whoHas);
}

PopLocus::PopLocus(const PopLocus & rhs)
{
  rhs.getSipperLoci(sipperLoci_);
  rhs.getSipperMotif(sipperMotif_);
  rhs.getSipperMotifBool(sipperMotifBool_);
  rhs.getFirstLocus(firstLocus_);
  rhs.getMotifs(motifs_);
  rhs.getWhoHasLocusInds(whoHasLocusInds_);
  rhs.getTotalCountPerPersonByMotif(totalCountPerPersonByMotif_);
  rhs.getMaxLength(maxLength_);
  rhs.getMaxCount(maxCount_);
  
  chr_       = rhs.chr();
  pos_       = rhs.pos();
  refLength_ = rhs.refLength();
  popCount_  = rhs.popCount();
  popSize_   = rhs.popSize();
}

PopLocus & PopLocus::operator=(const PopLocus & rhs)
{
  if(this != &rhs)
  {
	rhs.getSipperLoci(sipperLoci_);
    rhs.getSipperMotif(sipperMotif_);
    rhs.getSipperMotifBool(sipperMotifBool_);
    rhs.getFirstLocus(firstLocus_);
    rhs.getMotifs(motifs_);
    rhs.getWhoHasLocusInds(whoHasLocusInds_);
    rhs.getTotalCountPerPersonByMotif(totalCountPerPersonByMotif_);
    rhs.getMaxLength(maxLength_);
	rhs.getMaxCount(maxCount_);
  
    chr_       = rhs.chr();
    pos_       = rhs.pos();
    refLength_ = rhs.refLength();
    popCount_  = rhs.popCount();
    popSize_   = rhs.popSize();
  }
  return(*this);
}

void PopLocus::refresh(const vector< map< string, ProfileLocus > > & sipperLoci, const vector< bool > & whoHas)
{
  clear();
  sipperLoci_ = sipperLoci;
  
  if(popSize_ == 0 && sipperLoci.size() > 0)
  {
  	popSize_ = sipperLoci.size();
  }
  
  string obsMotif;
  unsigned int personWithLocus = 0;
  ProfileLocus currLocus;
  map< string, ProfileLocus >::const_iterator spIt;
  
  //for each person who has reads at this locus
  // assign them to the proper motif, report their
  // total counts for that motif at the locus, and determine
  // the longest tract length seen for the motif
  for(unsigned int i = 0; i < popSize_; ++i) {
  	if(whoHas[i]) {
	  personWithLocus = i;
      tr(sipperLoci_[i],spIt) {
		obsMotif = spIt->first;
		currLocus = spIt->second;

		sipperMotif_[obsMotif].push_back(i);
		totalCountPerPersonByMotif_[obsMotif].push_back(currLocus.totalCount());

		//see if current maxLength is maxLength for all people with this motif
		if(maxLength_.find(obsMotif) != maxLength_.end()) {
	      maxLength_[obsMotif] = max(currLocus.maxLength(),maxLength_[obsMotif]);
		}
		else {
	      maxLength_[obsMotif] = currLocus.maxLength();
		}

		//see if current totalCount is highest totalCount for all people with this motif
		if(maxCount_.find(obsMotif) != maxCount_.end()) {
	  	  maxCount_[obsMotif] = max(currLocus.totalCount(),maxCount_[obsMotif]);
		}
		else {
	  	  maxCount_[obsMotif] = currLocus.totalCount();
		}
	  }
    }
  }
  
  map< string, vector< unsigned int > >::iterator svIt;
  vector< unsigned int >::iterator vIt;
  // create a list of all the motifs seen at the locus
  // and a boolean vector of which people have the locus
  // with a particular motif
  
  tr(sipperMotif_, svIt)
  {
    motifs_.push_back(svIt->first);
	for(unsigned int i = 0; i < popSize_; ++i)
	{
		sipperMotifBool_[svIt->first].push_back(false);
	}
	tr(svIt->second,vIt)
	{
		sipperMotifBool_[svIt->first][*vIt] = true;
	}
  }
    
  // get the locus information from the last person
  // who has the current locus
  
  firstLocus_ = sipperLoci_[personWithLocus].begin()->second;
  chr_ = firstLocus_.chr();
  pos_ = firstLocus_.pos();
  refLength_ = firstLocus_.refLength();
}

//print locus info to file
void PopLocus::printLocusInfo(FILE * fh, const map< int, string > & int2chr) const
{
	string motif;
	for(unsigned int i = 0; i < motifs_.size(); ++i)
	{
		motif = motifs_[i];
		// print index information to file in following format:
		// chr	pos	motif	refLength	popCount	maxCount	sumCount
	 	fprintf(fh,"%s\t%d\t%s\t%d\t%u\t%u\t%u\n", int2chr.find(chr_)->second.c_str(), 
		         pos_,motif.c_str(), refLength_, (unsigned int) sipperMotif_.find(motif)->second.size(), 
				 maxCount_.find(motif)->second, 
				 accumulate(totalCountPerPersonByMotif_.find(motif)->second.begin(),
				            totalCountPerPersonByMotif_.find(motif)->second.end(), 0));
	}
}

void PopLocus::printHistInfo(FILE * fh) const
{
	string motif;
	ProfileLocus currLocus;
	for(unsigned int i = 0; i < motifs_.size(); ++i)
	{
		motif = motifs_[i];
		for(unsigned int j = 0; j < (popSize_ - 1); ++j)
		{
			if(sipperMotifBool_.find(motif)->second[j])
			{
				sipperLoci_[j].find(motif)->second.histSummary(fh, 0);
				fprintf(fh,"\t");
			}
			else
			{
				fprintf(fh,"0;;\t");
			}
		}
		if(sipperMotifBool_.find(motif)->second[popSize_ - 1])
		{
			sipperLoci_[popSize_ - 1].find(motif)->second.histSummary(fh, 0);
		}
		else
		{
			fprintf(fh,"0;;");
		}
		fprintf(fh,"\n");
	}
}

void PopLocus::clear()
{
  sipperMotif_.clear();
  sipperMotifBool_.clear();
  sipperLoci_.clear();
  firstLocus_.clear();
  maxLength_.clear();
  maxCount_.clear();
  
  chr_ = -1;
  pos_ = -1;
  motifs_.clear();
  refLength_ = -1;
  popCount_  = 0;
  popSize_   = 0;
  whoHasLocusInds_.clear();
  totalCountPerPersonByMotif_.clear();
}
