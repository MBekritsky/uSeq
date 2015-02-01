/*
Microsatellite.h

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
09.02.2010

AUTHOR:
Mitchell Bekritsky
*/

#ifndef __MICROSATELLITE_H__
#define __MICROSATELLITE_H__

#include <string>
#include <vector>
#include <iostream>
#include <ctype.h>
#include "FastaFile.h"

using namespace std;

class Microsatellite
{
 protected:
  string unit_, hdr_;
  unsigned int start_, length_, stop_;
	double mean_, stdDev_;
	unsigned int flankLength_;
	//If microsatellite is in a sequenced read and it doesn't pass based on 
	//particulatr criteria, it will be printed in lowercase
	bool pass_;
	
	char OFD_,IFD_;
	
  void initMicrosatelliteInfo (const string & unit, unsigned int start, unsigned int length);
	void initMicrosatelliteQual (double mean, double stdDev);
	void blankQual();
	void blankFlankLength() {flankLength_ = 0;}
 
 public:
  Microsatellite ();
  Microsatellite (const string & unit, unsigned int start, unsigned int length);
  Microsatellite (const string & unit, unsigned int start, unsigned int length, unsigned int flankLength);
  Microsatellite (const string & unit, unsigned int start, unsigned int length, double mean, double stdDev);
  Microsatellite (const string & unit, unsigned int start, unsigned int length, double mean, double stdDev, unsigned int flankLength);
  Microsatellite (const string & header, const string & unit, unsigned int start, unsigned int length);
  Microsatellite (const string & header, const string & unit, unsigned int start, unsigned int length, unsigned int flankLength);
  Microsatellite (const string & header, const string & unit, unsigned int start, unsigned int length, double mean, double stdDev);
  Microsatellite (const string & header, const string & unit, unsigned int start, unsigned int length, double mean, double stdDev, unsigned int flankLength);
  Microsatellite (const string & microsatelliteInfo);
	Microsatellite (const string & microsatelliteInfo, char IFD);
  Microsatellite (const string & microsatelliteInfo, unsigned int flankLength);
	Microsatellite (const string & microsatelliteInfo, unsigned int flankLength, char IFD);
  Microsatellite (const string & microsatelliteInfo, double mean, double stdDev);
	Microsatellite (const string & microsatelliteInfo, double mean, double stdDev, char IFD);
  Microsatellite (const string & microsatelliteInfo, double mean, double stdDev, unsigned int flankLength);
	Microsatellite (const string & microsatelliteInfo, double mean, double stdDev, unsigned int flankLength, char IFD);
  Microsatellite (const string & header, const string & microsatelliteInfo);
	Microsatellite (const string & header, const string & microsatelliteInfo, char IFD);
	Microsatellite (const string & header, const string & microsatelliteInfo, char IFD, char OFD);
  Microsatellite (const string & header, const string & microsatelliteInfo, unsigned int flankLength);
	Microsatellite (const string & header, const string & microsatelliteInfo, unsigned int flankLength, char IFD);
  Microsatellite (const string & header, const string & microsatelliteInfo, double mean, double stdDev);
	Microsatellite (const string & header, const string & microsatelliteInfo, double mean, double stdDev, char IFD);
  Microsatellite (const string & header, const string & microsatelliteInfo, double mean, double stdDev, unsigned int flankLength);
	Microsatellite (const string & header, const string & microsatelliteInfo, double mean, double stdDev, unsigned int flankLength, char IFD);
  Microsatellite (const string & header, unsigned int index, unsigned int count, unsigned int repeatLength, const string & seq);
  Microsatellite (const string & header, unsigned int index, unsigned int count, unsigned int repeatLength, const string & seq, unsigned int flankLength);  
	Microsatellite (const Microsatellite & microsatellite);
	
	Microsatellite & operator=(const Microsatellite & rhs);
	
	void toLower(string & word) const;

	string 				hdr()				const {return hdr_;}
	string 				unit()			const {return unit_;}
	unsigned int  length() 		const {return length_;}
	unsigned int  start()			const {return start_;}
	unsigned int  stop() 			const {return stop_;}
	unsigned int  flankLength() const {return flankLength_;}

	void					getHdr(string & hdr) const {hdr = hdr_;}
	void					getUnit(string & unit) const {unit = unit_;}

  string 				getUnit()				 const {return unit_;}
	string 				getHdr() 				 const {return hdr_;}
  unsigned int  getLength() 		 const {return length_;}
  unsigned int	getStart() 			 const {return start_;}
  unsigned int  getStop()				 const {return stop_;}
	unsigned int	getFlankLength() const {return flankLength_;}

	double mean() 	const {return mean_;}
	double stdDev()	const {return stdDev_;}

  double getMean()	 const {return mean_;}
  double getStdDev() const {return stdDev_;}

	char OFD() const {return OFD_;}
	char IFD() const {return IFD_;}

	char getOFD() const {return OFD_;}
	char getIFD() const {return IFD_;}
	
	bool passed() const {return pass_;}

  string getFrontFlank (const string & seq) const;
  string getRearFlank (const string & seq) const;
  string getSnapshot (const string & seq) const;
  unsigned int getSeqLength (const string & seq, const string & functionName) const;
  
	void fail () {pass_ = false;}
	void pass () {pass_ = true;}
	void setBasics  (const string & unit, unsigned int start, unsigned int stop, unsigned int length);
	void setHdr			(const string & hdr)		{hdr_ = hdr;}
  void setUnit		(const string & unit) 	{unit_		= unit;}
  void setLength	(unsigned int length)  	{length_	= length;}
  void setStart		(unsigned int start)   	{start_ 	= start;}
	void setStop		(unsigned int stop)			{stop_ = stop;}
  void setMean		(double mean) 	{mean_		=  mean;}
  void setStdDev	(double stdDev) {stdDev_  = stdDev;}
	
	void setOFD (char OFD) {OFD_ = OFD;}
	void setIFD (char IFD) {IFD_ = IFD;}
	
	void setDefaultOFD () {OFD_ = ':';}
	void setDefaultIFD () {IFD_ = ':';}
	void setDefaultDelims();
	
  //When the microsatellite comes from a Fasta file, it is given default mean 
  //and standard deviation values
  void setDefaultMean () {mean_ 	= 30;}
  void setDefaultStdDev(){stdDev_ = 0;}
  void setFlankLength(unsigned int flankLength) {flankLength_ = flankLength;}
  void setDefaultFlankLength() {flankLength_ = 30;}
  
  bool parseMicrosatelliteString(const string & MicrosatelliteInfo);
  
  void printWithHdr  (char buffer[]) const;
  void printNoHdr (char buffer[]) const;
	void printNoHdr (char buffer[], const string & qual) const;
  void printFastaStyle (const string & seq) const;
  void printFastaStyle (FILE * output, const string & seq) const;
  void printFastaStyle (FILE * output, const FastaFile & fastaFile) const;
	void printFastaStyle (FILE * output) const;
	void printFastaStyle () const;
};
#endif /* microsatellite.h */
