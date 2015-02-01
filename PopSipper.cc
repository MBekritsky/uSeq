/*
PopSipper.cc

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

#include "PopSipper.h"

using namespace std;

PopSipper::PopSipper()
{
  clear();
}

PopSipper::PopSipper(const map< string, int > & chr2int)
{
  clear();
  chr2int_ = chr2int;
}

PopSipper::PopSipper(const vector< ProfileSipper > & sippers, const map< string, int > & chr2int)
{
  clear();
  sippers_ = sippers;
  
  popSize_ = sippers_.size();
  chr2int_ = chr2int;
  initializeVectors();
}

PopSipper::PopSipper(const PopSipper & rhs)
{
  clear();
  rhs.getSippers(sippers_);
  rhs.getChr2Int(chr2int_);
  popSize_ = rhs.popSize();
  
  //when ProfileSippers are copied, they restart
  //from the beginning of the file.  Therefore,
  //all vectors need to be reinitialized
  
  initializeVectors();
}

PopSipper & PopSipper::operator=(const PopSipper & rhs)
{
  if(this != &rhs)
  {
    clear();
	rhs.getSippers(sippers_);
	rhs.getChr2Int(chr2int_);
	popSize_ = rhs.popSize();

	//when ProfileSippers are copied, they restart
	//from the beginning of the file.  Therefore,
	//all vectors need to be reinitialized

	initializeVectors();
  }
  return *this;
}

void PopSipper::setSippers(const vector< ProfileSipper > & sippers)
{
  sippers_ = sippers;
  popSize_ = sippers.size();
  initializeVectors();
}

void PopSipper::clear()
{
  sippers_.clear();
  loci_.clear();
  chr2int_.clear();
  active_.clear();
  toUpdate_.clear();
  popChr_.clear();
  popPos_.clear();
  popMotif_.clear();
  
  popSize_ = 0;
  currChr_ = -1;
  currPos_ = -1;
  popAtNextLocus_ = 0;
  numActive_ = 0;
}

void PopSipper::initializeVectors()
{
	if(popSize_ > 0)
	{
	  active_.clear();
	  active_.resize(popSize_,true);

	  toUpdate_.clear();
	  toUpdate_.resize(popSize_,true);

	  popChr_.clear();
	  popChr_.resize(popSize_,-1);

	  popPos_.clear();
	  popPos_.resize(popSize_,-1);
	  	  
	  popMotif_.clear();
	  popMotif_.resize(popSize_,"");

	  loci_.clear();
      loci_.resize(popSize_); //default constructor for ProfileLocus will be used here

	  numActive_ = popSize_;
	}
	else
	{
	  cerr << "PopSipper::initializeVectors called with no ProfilerSippers defined" << endl;
	  exit(EXIT_FAILURE);
	}
}

bool PopSipper::updateProfileLoci()
{
  map< string, ProfileLocus > tempLocusMap;
  ProfileLocus tempLocus;
  
  //if a person is active (has not reached the end of the profile)
  // and is updateable (last locus was just reported), advance sipper
  // to the next locus
  for(unsigned int i = 0; i < popSize_; ++i)
  {
    if(active_[i] && toUpdate_[i])
	{
	  if(sippers_[i].nextLocus()) //EOF not yet reached
	  {
	    tempLocusMap.clear();
		tempLocus.clear();
	    sippers_[i].getCurrLocus(tempLocusMap);
		loci_[i] = tempLocusMap;
				
		//different loci in temp locus can only differ by motif
		//not by chromosome or position
		tempLocus = tempLocusMap.begin()->second;
		popChr_[i] = tempLocus.chr();
		popPos_[i] = tempLocus.pos();
	  }
	  else //EOF reached by ProfileSipper
	  {
	    active_[i] = false;
		popChr_[i] = 10000; //nonsense chromosome number
		popPos_[i] = -1;
	  }
	}
  }
  
  numActive_ = accumulate(active_.begin(),active_.end(),0);

  if(numActive_ > 0) {
    //find the next locus to be reported
    currChr_ = *min_element(popChr_.begin(), popChr_.end());
    vector< int > tempPos;
    for(unsigned int i = 0; i < popSize_; ++i){
      if(popChr_[i] == currChr_ && active_[i]) {
	    tempPos.push_back(popPos_[i]);
	  }
    }
    currPos_ = *min_element(tempPos.begin(), tempPos.end());
  
    //find who has the next locus to be reported
    for(unsigned int i = 0; i < popSize_; ++i)
    {
      toUpdate_[i] = (popChr_[i] == currChr_) && (popPos_[i] == currPos_);
    }
  
    popAtNextLocus_ = accumulate(toUpdate_.begin(),toUpdate_.end(),0);
    numActive_ = accumulate(active_.begin(),active_.end(),0);
    return true;
  }
  return false;
}
 
void PopSipper::setLatestPopLocus()
{
  currLocus_.refresh(loci_,toUpdate_);
}

bool PopSipper::nextLocus()
{
  if(numActive_ == 0) //no more active ProfileSippers
  {
    return false;
  }
  else
  {
	bool stillActive = updateProfileLoci();
	if(stillActive) {
	  setLatestPopLocus();
	  return true;
	}
	else {
	  return false;
	}
  }
  return false;
}
