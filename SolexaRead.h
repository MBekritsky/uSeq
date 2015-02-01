/*
SolexaRead.h

ABSTRACT:
A class and functions for handling sequencing reads from
Solexa data.  Also allows for read trimming based on
quality scores and score conversion if the sequence data
was generated using Solexa v1.3 or earlier. A child class
for sequence data when looking for microsatellites that
tracks microsatellites in a read, as well as whether they
have passed quality control filters.

CREATION DATE:
09.08.2010

LAST REVISION:
26.08.2010

AUTHOR:
Mitchell Bekritsky
*/

#ifndef __SOLEXA_READ_H__
#define __SOLEXA_READ_H__

#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include "Microsatellite.h"

using namespace std;

class SolexaRead
{
 public:
  enum Quality_t {QUAL_SANGER, QUAL_SOLEXA13};
  enum Trim_t {TRIM_NONE, TRIM_BWA_STYLE, TRIM_FIRST_BAD, TRIM_LAST_GOOD};
	
 private:
  string hdr_, seq_, qual_, oqual_;
  unsigned int readThru_;
	unsigned int barcodeLength_;
  unsigned int lSeq_;
	unsigned int nCount_;
  Quality_t qualityType_;
  Trim_t trimType_;
  bool readWarning_;
	int minQual_;
	char OFD_, ORD_, IFD_, IRD_;
	unsigned int readNumber_;
	bool passFilter_;
	
	//Microsatellite variables
  unsigned int nMicrosatellite_;
  vector< Microsatellite > microsatellites_;
	
  void initRead(const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t);
  void initRead(const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual);
  void initRead(const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, unsigned int barcodeTrim);
  void initRead(const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual, unsigned int barcodeTrim);
  void blankMicrosatelliteInfo();
  
 public:
  SolexaRead ();
  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t);
  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual);
	SolexaRead (Quality_t quality_t, Trim_t trim_t, int minQual);
	SolexaRead (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD);
	SolexaRead (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD, char ORD, char OFD);
  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, const vector< string > & unparsedMicrosatellites);
	SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, const vector< string > & unparsedMicrosatellites, char IRD, char IFD);
  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual, const vector< string > & unparsedMicrosatellites);
	SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual, const vector< string > & unparsedMicrosatellites, char IRD, char IFD);
  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, const vector< Microsatellite > & parsedMicrosatellites);
  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual, const vector< Microsatellite > & parsedMicrosatellites);

  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, unsigned int barcodeLength);
  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual, unsigned int barcodeLength);
	SolexaRead (Quality_t quality_t, Trim_t trim_t, int minQual, unsigned int barcodeLength);
	SolexaRead (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD, unsigned int barcodeLength);
	SolexaRead (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD, char ORD, char OFD, unsigned int barcodeLength);
  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, const vector< string > & unparsedMicrosatellites, unsigned int barcodeLength);
	SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, const vector< string > & unparsedMicrosatellites, char IRD, char IFD, unsigned int barcodeLength);
  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual, const vector< string > & unparsedMicrosatellites, unsigned int barcodeLength);
	SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual, const vector< string > & unparsedMicrosatellites, char IRD, char IFD, unsigned int barcodeLength);
  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, const vector< Microsatellite > & parsedMicrosatellites, unsigned int barcodeLength);
  SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual, const vector< Microsatellite > & parsedMicrosatellites, unsigned int barcodeLength);
	SolexaRead (const SolexaRead & rhs);
	
	SolexaRead & operator=(const SolexaRead & rhs);

	void initContainer (Quality_t quality_t, Trim_t trim_t, int minQual);
	void initContainer (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD);
  void initContainer (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD, char ORD, char OFD);
	void initContainer (Quality_t quality_t, Trim_t trim_t, int minQual, unsigned int barcodeLength);
	void initContainer (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD, unsigned int barcodeLength);
  void initContainer (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD, char ORD, char OFD, unsigned int barcodeLength);
		
	int 				 minQual() 			const {return minQual_;}
	char 				 OFD() 					const {return OFD_;}
	char 				 ORD()        	const {return ORD_;}
	char 				 IFD() 					const {return IFD_;}
	char 				 IRD() 					const {return IRD_;}
	bool 				 readWarning()	const {return readWarning_;}
	bool 				 passed()				const {return passFilter_;}
	bool				 passedFilter() const {return passFilter_;}
	string 			 hdr() 	 				const {return hdr_;}
	string 			 seq() 	 				const {return seq_;}
	string 			 qual()  				const {return qual_;}
	string 			 oQual() 				const {return oqual_;}
	Trim_t 			 trimType() 		const {return trimType_;}
	Quality_t 	 qualityType()	const {return qualityType_;}
	unsigned int readThru() 		const {return readThru_;}
	unsigned int lSeq() 				const {return lSeq_;}
	unsigned int readLen()			const {return lSeq_;}
	unsigned int nCount() 			const {return nCount_;}
	unsigned int readNumber()		const {return readNumber_;}
	unsigned int barcodeLength() const {return barcodeLength_;}
	
	int 				 getMinQual() 		const {return minQual_;}
	char 				 getOFD() 				const {return OFD_;}
	char 				 getORD() 				const {return ORD_;}
	char 				 getIFD() 				const {return IFD_;}
	char 				 getIRD() 				const {return IRD_;}
	bool 				 getReadWarning() const {return readWarning_;}
  string 			 getHdr()   			const  {return hdr_;}
  string 			 getSeq()   			const  {return seq_;}
  string 			 getQual()  			const	 {return qual_;}
	string 			 getOQual() 			const  {return oqual_;}
	Trim_t 			 getTrimType() 		const {return trimType_;}
	Quality_t 	 getQualityType() const {return qualityType_;}
  unsigned int getReadThru() 		const {return readThru_;}
	unsigned int getLSeq() 		 		const {return lSeq_;}
	unsigned int getNCount() 	 		const {return nCount_;}
	unsigned int getBarcodeLength() const {return barcodeLength_;}
	
  void getMicrosatellites(vector< Microsatellite > & microsatellites) const {microsatellites = microsatellites_;}
  unsigned int getNumMicrosatellites() const {return nMicrosatellite_;}
	unsigned int nUsat() const {return nMicrosatellite_;}
	unsigned int getNumPassedMicrosatellites();

	void setNCount();
  void setHdr(const string & hdr);
  void setSeq(const string & seq);
  void setQual(const string & qual);
  void setTrimType(Trim_t trimType) {trimType_ = trimType;}
  void setDefaultMinQual() {minQual_ = 20;}
  void setMinQual(int minQual) {minQual_ = minQual;}
	void setOFD (char OFD) {OFD_ = OFD;}
	void setORD (char ORD) {ORD_ = ORD;}
	void setIFD (char IFD) {IFD_ = IFD;}
	void setIRD (char IRD) {IRD_ = IRD;}
	void setBarcodelength (unsigned int barcodeLength) {barcodeLength_ = barcodeLength;}
  void setQualityType(Quality_t qualityType) {qualityType_ = qualityType;}
  void setMicrosatellites(const vector< Microsatellite > & microsatellites);

	void setDefaultOFD() {OFD_ = ':';}
	void setDefaultORD() {ORD_ = '|';}
	void setDefaultIFD() {IFD_ = ':';}
	void setDefaultIRD() {IRD_ = '|';}
	void setDefaultDelims();
	
  void resetReadThru() {readThru_ = (int) seq_ . length();}
  void printReadInfo();
  
	void compactHdr();
	void parseHdr();
	
	void printFastqStyle();
  void printFastqStyle(FILE * output);
	
	void printModifiedReadSkeleton(const string & qual, FILE * output, bool printMicrosatelliteInfo);
	void printModifiedRead();
	void printModifiedRead(FILE * output);
	void printModifiedReadWithQual(const string & qual, FILE * output);
	void printModifiedReadWithQualNoModHeader(const string & qual, FILE * output);
	void printMicrosatelliteHeader(const string & qual, FILE * output);
	void printMicrosatelliteHeader(const string & qual, FILE * output, bool isFirst);
	  
	void printPassedMicrosatellitesFastq();
  void printPassedMicrosatellitesFastq(FILE* output);
  
	void printMicrosatelliteFastqSkeleton(const string & qual, FILE * output);
	void printMicrosatelliteFastq();
  void printMicrosatelliteFastq(FILE* output);
	void printMicrosatelliteFastqWithQual(const string & qual, FILE * output);
	
  string displayTrimType() const;
	string goodSeq()  const	{return getGoodSeq();}
	string goodQual() const {return getGoodQual();}
  string getGoodSeq() const;
  string getGoodQual() const;
  string getGoodPortion(const string & text) const;
	void checkMicrosatelliteCoords();
  
  char cIllumina13ToCSanger(char cIllumina);
  void adjustQualityScores();
  void matchLengths();
  
  void msMean(Microsatellite& microsatellite);
  void msStdDev(Microsatellite& microsatellite);
  
  void addMicrosatellites(const vector< string > & unparsedMicrosatellites);
  
  void trim();
  void trimBwaStyle();
  void trimToFirstBad();
  void trimToLastGood();
	void trimBarcode();

  void printReadWarning(int warningCode);
};
#endif /* SolexaRead.h */
