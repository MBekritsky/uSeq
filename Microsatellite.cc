/*
Microsatellite.cc

ABSTRACT:
A basic microsatellite class and functions. Microsatellites
are defined by their repeat unit, their coordinates, and 
their length.  There is also an option to specify the mean
and standard deviation for the quality of a microsatellite if it
is derived from sequence data.  There is also an option to 
retrieve flanking sequence and a microsatellite "snapshot".

CREATION DATE:
09.08.2010

LAST REVISION:
02.09.2010

AUTHOR:
Mitchell Bekritsky
*/

#include "Microsatellite.h"
#include <cstdlib>
#include <cstdio>

using namespace std;

Microsatellite::Microsatellite()
{
  start_	= 0;
  stop_		= 0;
  length_ = 0;

  setDefaultFlankLength();
	setDefaultDelims();
  blankQual();	
  unit_ . clear();
  hdr_ . clear();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & unit, unsigned int start, unsigned int length)
{
  hdr_ . clear();
  setDefaultFlankLength();
	setDefaultDelims();
  initMicrosatelliteInfo (unit,start,length);
  blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & unit, unsigned int start, unsigned int length, unsigned int flankLength)
{
  hdr_ . clear();
  setFlankLength(flankLength);
  initMicrosatelliteInfo (unit, start, length);
  blankQual();
	pass_ = true;
}

//Microsatellites that are in reads have a default flank length of 0
Microsatellite::Microsatellite (const string & unit, unsigned int start, unsigned int length, double mean, double stdDev)
{
  hdr_ . clear();
  blankFlankLength();
	setDefaultDelims();
  initMicrosatelliteInfo (unit,start,length);
  initMicrosatelliteQual (mean,stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & unit, unsigned int start, unsigned int length, double mean, double stdDev, unsigned int flankLength)
{
  hdr_ . clear();
  setFlankLength (flankLength);
	setDefaultDelims();
  initMicrosatelliteInfo (unit, start, length);
  initMicrosatelliteQual (mean, stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & unit, unsigned int start, unsigned int length)
{
  hdr_ = hdr;
  setDefaultFlankLength();
	setDefaultDelims();
  initMicrosatelliteInfo (unit,start,length);
  blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & unit, unsigned int start, unsigned int length, unsigned int flankLength)
{
  hdr_ = hdr;
  setFlankLength (flankLength);
	setDefaultDelims();
  initMicrosatelliteInfo (unit,start,length);
  blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & unit, unsigned int start, unsigned int length, double mean, double stdDev)
{
  hdr_ = hdr;
  blankFlankLength();
	setDefaultDelims();
  initMicrosatelliteInfo (unit,start,length);
  initMicrosatelliteQual (mean,stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & unit, unsigned int start, unsigned int length, double mean, double stdDev, unsigned int flankLength)
{
  hdr_ = hdr;
  setFlankLength (flankLength);
	setDefaultDelims();
  initMicrosatelliteInfo (unit, start, length);
  initMicrosatelliteQual (mean, stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & microsatelliteInfo)
{
  hdr_ . clear();
  setDefaultFlankLength();
	setDefaultDelims();
  parseMicrosatelliteString (microsatelliteInfo);
  blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & microsatelliteInfo, char IFD)
{
	hdr_ . clear();
	setDefaultFlankLength();
	setDefaultOFD();
	setIFD(IFD);
	parseMicrosatelliteString(microsatelliteInfo);
	blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & microsatelliteInfo, unsigned int flankLength)
{
  hdr_ . clear();
  setFlankLength (flankLength);
	setDefaultDelims();
  parseMicrosatelliteString (microsatelliteInfo);
  blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & microsatelliteInfo, unsigned int flankLength, char IFD)
{
	hdr_ . clear();
	setFlankLength(flankLength);
	setDefaultOFD();
	setIFD(IFD);
	parseMicrosatelliteString (microsatelliteInfo);
	blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & microsatelliteInfo, double mean, double stdDev)
{
  hdr_ . clear();
  blankFlankLength();
	setDefaultDelims();
  parseMicrosatelliteString (microsatelliteInfo);
  initMicrosatelliteQual (mean, stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & microsatelliteInfo, double mean, double stdDev, char IFD)
{
	hdr_ . clear();
	blankFlankLength();
	setDefaultOFD();
	setIFD(IFD);
	parseMicrosatelliteString (microsatelliteInfo);
	initMicrosatelliteQual (mean, stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & microsatelliteInfo, double mean, double stdDev, unsigned int flankLength)
{
  hdr_ . clear();
  setFlankLength(flankLength);
	setDefaultDelims();
  parseMicrosatelliteString (microsatelliteInfo);
  initMicrosatelliteQual (mean, stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & microsatelliteInfo, double mean, double stdDev, unsigned int flankLength, char IFD)
{
	hdr_ . clear();
	setFlankLength(flankLength);
	setDefaultOFD();
	setIFD(IFD);
	parseMicrosatelliteString (microsatelliteInfo);
	initMicrosatelliteQual (mean, stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & microsatelliteInfo)
{
  hdr_ = hdr;
  setDefaultFlankLength();
	setDefaultDelims();
  parseMicrosatelliteString (microsatelliteInfo);
  blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & microsatelliteInfo, char IFD)
{
	hdr_ = hdr;
	setDefaultFlankLength();
	setDefaultOFD();
	setIFD(IFD);
	parseMicrosatelliteString (microsatelliteInfo);
	blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & microsatelliteInfo, char IFD, char OFD)
{
	hdr_ = hdr;
	setDefaultFlankLength();
	setIFD(IFD);
	setOFD(OFD);
	parseMicrosatelliteString(microsatelliteInfo);
	blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & microsatelliteInfo, unsigned int flankLength)
{
  hdr_ = hdr;
  setFlankLength (flankLength);
	setDefaultDelims();
  parseMicrosatelliteString (microsatelliteInfo);
  blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & microsatelliteInfo, unsigned int flankLength, char IFD)
{
	hdr_ = hdr;
	setFlankLength(flankLength);
	setDefaultOFD();
	setIFD(IFD);
	parseMicrosatelliteString (microsatelliteInfo);
	blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & microsatelliteInfo, double mean, double stdDev)
{
  hdr_ = hdr;
  blankFlankLength();
	setDefaultDelims();
  parseMicrosatelliteString (microsatelliteInfo);
  initMicrosatelliteQual (mean,stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & microsatelliteInfo, double mean, double stdDev, char IFD)
{
	hdr_ = hdr;
	blankFlankLength();
	setDefaultOFD();
	setIFD(IFD);
	parseMicrosatelliteString (microsatelliteInfo);
	initMicrosatelliteQual (mean, stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & microsatelliteInfo, double mean, double stdDev, unsigned int flankLength)
{
  hdr_ = hdr;
  setFlankLength (flankLength);
	setDefaultDelims();
  parseMicrosatelliteString (microsatelliteInfo);
  initMicrosatelliteQual (mean, stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, const string & microsatelliteInfo, double mean, double stdDev, unsigned int flankLength, char IFD)
{
	hdr_ = hdr;
	setFlankLength(flankLength);
	setDefaultOFD();
	setIFD(IFD);
	parseMicrosatelliteString (microsatelliteInfo);
	initMicrosatelliteQual (mean, stdDev);
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, unsigned int index, unsigned int count, unsigned int repeatLength, const string & seq)
{
  hdr_		= hdr;
  start_	= index - count + 1;
  stop_		= index;
  length_ = count;
  unit_ = seq . substr(start_ - 1, repeatLength);
  
  blankQual();
	setDefaultDelims();
  setDefaultFlankLength();
	pass_ = true;
}

Microsatellite::Microsatellite (const string & hdr, unsigned int index, unsigned int count, unsigned int repeatLength, const string & seq, unsigned int flankLength)
{
  hdr_ = hdr;
  start_ = index - count + 1;
  stop_ = index;
  length_ = count;
  unit_ = seq . substr(start_ - 1, repeatLength);

  setFlankLength (flankLength);
	setDefaultDelims();
  blankQual();
	pass_ = true;
}

Microsatellite::Microsatellite (const Microsatellite & microsatellite)
{
	hdr_ 	       = microsatellite . getHdr();
	unit_        = microsatellite . getUnit();
	start_       = microsatellite . getStart();
	stop_        = microsatellite . getStop();
	length_      = microsatellite . getLength();
	pass_        = microsatellite . passed();
	mean_        = microsatellite . getMean();
	stdDev_      = microsatellite . getStdDev();
	flankLength_ = microsatellite . getFlankLength();
	OFD_  			 = microsatellite . getOFD();
	IFD_  			 = microsatellite . getIFD();
}

Microsatellite & Microsatellite::operator= (const Microsatellite & rhs)
{
	if(this != &rhs)
		{
			hdr_ 					= rhs . getHdr();
			unit_ 				= rhs . getUnit();
			start_ 				= rhs . getStart();
			stop_ 				= rhs . getStop();
			length_ 			= rhs . getLength();
			pass_ 				= rhs . passed();
			mean_ 				= rhs . getMean();
			stdDev_			  = rhs . getStdDev();
			flankLength_  = rhs . getFlankLength();
			OFD_					= rhs . getOFD();
			IFD_					= rhs . getIFD();
		}
	return *this;
}

void Microsatellite::initMicrosatelliteInfo(const string & unit, unsigned int start, unsigned int length)
{
  unit_		= unit;
  start_	= start;
  length_ = length;
  stop_		= start_ + length_ - 1;
}

void Microsatellite::initMicrosatelliteQual(double mean, double stdDev)
{
	mean_		= mean;
	stdDev_ = stdDev;
}

void Microsatellite::blankQual()
{
	mean_ 	= -1;
	stdDev_	= -1;
}

void Microsatellite::setDefaultDelims()
{
	setDefaultOFD();
	setDefaultIFD();
}

void Microsatellite::setBasics(const string & unit, unsigned int start, unsigned int stop, unsigned int length)
{
	unit_ = unit;
	start_ = start;
	length_ = length;
	stop_ = stop;
}

bool Microsatellite::parseMicrosatelliteString (const string & info)
{
	unsigned int strStart = 0;
	unsigned int length = 0;

	for(; info [length] != IFD_; length++)
		;
	unit_ = info . substr(strStart,length);
	
	length++;
	strStart = length;
	
	for(; info [length] != IFD_; length++)
  	;
	start_ = atoi((info . substr(strStart, length - strStart)) . c_str());
	length_ = atoi((info . substr(length + 1)) . c_str());
	stop_		= start_ + length_ - 1;
	
  return true;
}

void Microsatellite::toLower(string & word) const
{
	for(unsigned int i = 0; i < word . length(); i++)
		word [i] = tolower(word[i]);
}

void Microsatellite::printWithHdr (char msText[1024]) const
{
	string unit = unit_;
	if(!pass_)
		toLower(unit);
  sprintf (msText,"%s%c%s%c%d%c%d", hdr_ . c_str(), OFD_, unit . c_str(), OFD_, start_, OFD_, length_);
}

void Microsatellite::printNoHdr(char msText[1024]) const
{
	string unit = unit_;
	if(!pass_)
		toLower(unit);
  sprintf (msText,"%s%c%d%c%d", unit . c_str(), OFD_, start_, OFD_, length_);
}

void Microsatellite::printNoHdr(char msText[2048], const string & qual) const
{
	string unit = unit_;
	if(!pass_)
		toLower(unit);
	sprintf (msText,"%s%c%d%c%d", unit . c_str(), OFD_, start_, OFD_, length_);
	
	if(!qual . empty())
		sprintf (msText, "%s%c%s", msText, OFD_, qual . substr(start_-1, length_) . c_str());
}

string Microsatellite::getFrontFlank (const string & seq) const
{
	//check that front flank is within sequence, otherwise, return empty string
	if(start_  == 1)
		return "";
		
	unsigned int frontFlankLength = flankLength_;
	
	//check that front flank is completely contained within sequence,
	//otherwise, return what you can
	if(frontFlankLength > start_ - 1)
		frontFlankLength = start_ -1;
	
	string fFlank = seq . substr (start_ - frontFlankLength - 1, frontFlankLength);
	
	return fFlank;
}

string Microsatellite::getRearFlank(const string & seq) const
{
	unsigned int seqLength = getSeqLength (seq, "getRearFlank");
	unsigned int rearFlankLength = flankLength_;

	/*
	//Check to see that rear flank is within sequence, otherwise return empty string
	if ( (rearFlankLength + stop_ - 1) >= seqLength);
		return "";
	*/
	
	//Check to see that rear flank is completely contained within sequence,
	//otherwise, return what you can
	if (rearFlankLength > seqLength - stop_ + 1)
		rearFlankLength = seqLength - stop_ + 1;
		
	string rFlank = seq . substr(stop_, rearFlankLength);
	
	return rFlank;
}

string Microsatellite::getSnapshot(const string & seq) const
{
	string snapshot;
	unsigned int seqLength = getSeqLength (seq, "getSnapshot");
	if (start_ > 1)
		snapshot . push_back ( (char) tolower (seq [start_ - 2]));
	for(unsigned int i = start_ - 1; i < stop_; i++)
		snapshot . push_back ( (char) toupper (seq [i]));
	if (stop_ < seqLength)
		snapshot . push_back ( (char) tolower (seq [stop_]));
		
	return snapshot;
}

//checks that sequence length is >= 0, otherwise kills program
//I chose to kill the program here because looking for sequence
//outside of seq range could seg fault program
unsigned int Microsatellite::getSeqLength (const string & seq, const string & function) const
{
	unsigned int seqLength = seq . length();
	if(seqLength == 0)
		{
			cerr << function . c_str() << ": sequence provided is empty!!" << endl;
			exit (EXIT_FAILURE);
		}
	return seqLength;
}

void Microsatellite::printFastaStyle(const string & seq) const
{
  printf("%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n", hdr_ . c_str(),
	 start_, stop_, length_, unit_ . c_str(), getFrontFlank (seq) . c_str(),
	 getRearFlank (seq) . c_str(), getSnapshot (seq) . c_str());
}

void Microsatellite::printFastaStyle(FILE * output, const string & seq) const
{
  fprintf(output, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n", hdr_ . c_str(),
	 start_, stop_, length_, unit_ . c_str(), getFrontFlank (seq) . c_str(),
	 getRearFlank (seq) . c_str(), getSnapshot (seq) . c_str());

}

void Microsatellite::printFastaStyle(FILE * output, const FastaFile & fastaFile) const
{
	string seq;
	fastaFile . getSeq (seq);
	fprintf(output, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n", hdr_ . c_str(),
	start_, stop_, length_, unit_ . c_str(), getFrontFlank (seq) . c_str(),
	getRearFlank (seq) . c_str(), getSnapshot (seq) . c_str());
}

void Microsatellite::printFastaStyle()const
{
	printFastaStyle(stdout);
}

void Microsatellite::printFastaStyle(FILE * output) const
{
	fprintf(output, "%s\t%d\t%d\t%d\t%s\n", hdr_ . c_str(), start_, stop_,
	length_, unit_ . c_str());
}
