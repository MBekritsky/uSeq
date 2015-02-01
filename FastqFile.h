
#ifndef __FASTQ_FILE_H__
#define __FASTQ_FILE_H__

#include <vector>
#include <string>
#include "SolexaRead.h"
#include "Microsatellite.h"
#include "FileOps.h"

using namespace std;

class FastqFile
{
 private:
  vector< SolexaRead  > reads_;
	SolexaRead currentRead_;
  FILE * solexaFileP_;
  
  SolexaRead::Quality_t qualityType_;
  SolexaRead::Trim_t trimType_;
  int minQual_;
	unsigned int barcodeLength_;
	char IRD_, IFD_, ORD_, OFD_;

  int nRead_, nPass_;
	
	void blankCounters();
	
 public:
  FastqFile();
  FastqFile(FILE * fastqFile);
  FastqFile(FILE * fastqFile, SolexaRead::Quality_t quality_t);
  FastqFile(FILE * fastqFile, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t);
  FastqFile(FILE * fastqFile, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual);
	FastqFile(FILE * fastqFile, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual, char IRD, char IFD);
	FastqFile(FILE * fastqFile, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual, char IRD, char IFD, char ORD, char OFD);
  FastqFile(FILE * fastqFile, unsigned int barcodeLength);
  FastqFile(FILE * fastqFile, SolexaRead::Quality_t quality_t, unsigned int barcodeLength);
  FastqFile(FILE * fastqFile, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, unsigned int barcodeLength);
  FastqFile(FILE * fastqFile, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual, unsigned int barcodeLength);
	FastqFile(FILE * fastqFile, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual, char IRD, char IFD, unsigned int barcodeLength);
	FastqFile(FILE * fastqFile, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual, char IRD, char IFD, char ORD, char OFD, unsigned int barcodeLength);  
	FastqFile(const string & filename);
  FastqFile(const string & filename, SolexaRead::Quality_t quality_t);
  FastqFile(const string & filename, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t);
  FastqFile(const string & filename, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual);
	FastqFile(const string & filename, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual, char IRD, char IFD);
	FastqFile(const string & filename, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual, char IRD, char IFD, char ORD, char OFD);
  FastqFile(const string & filename, unsigned int barcodeLength);
  FastqFile(const string & filename, SolexaRead::Quality_t quality_t, unsigned int barcodeLength);
  FastqFile(const string & filename, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, unsigned int barcodeLength);
  FastqFile(const string & filename, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual, unsigned int barcodeLength);
	FastqFile(const string & filename, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual, char IRD, char IFD, unsigned int barcodeLength);
	FastqFile(const string & filename, SolexaRead::Quality_t quality_t, SolexaRead::Trim_t trim_t, int minQual, char IRD, char IFD, char ORD, char OFD, unsigned int barcodeLength);
  
	/*No copy copy constructor, need to add one*/
	
  ~FastqFile();
  
	void setIRD(char IRD) {IRD_ = IRD;}
	void setIFD(char IFD) {IFD_ = IFD;}
	void setORD(char ORD) {ORD_ = ORD;}
	void setOFD(char OFD) {OFD_ = OFD;}
	void setDefaultIRD () {IRD_ = '|';}
	void setDefaultIFD () {IFD_ = ':';}
	void setDefaultORD () {ORD_ = '|';}
	void setDefaultOFD () {OFD_ = ':';}
	char getIRD() {return IRD_;}
	char getIFD() {return IFD_;}
	char getORD() {return ORD_;}
	char getOFD() {return OFD_;}
	void setDefaultDelims();

	void setBarcodeLength(unsigned int barcodeLength) {barcodeLength_ = barcodeLength;}
	unsigned int getBarcodeLength() const {return barcodeLength_;}
	void barcodeLength(unsigned int barcodeLength) {barcodeLength_ = barcodeLength;}
	unsigned int barcodeLength() const {return barcodeLength_;}
	
  vector< SolexaRead > getReads () const {return reads_;}
  
  void setQualityType (SolexaRead::Quality_t qualityType) {qualityType_ = qualityType;}
  void setTrimType    (SolexaRead::Trim_t trimType)       {trimType_ = trimType;}
	
	string displayTrimType    () {return currentRead_ . displayTrimType();}
	
  void setDefaultQualityType () {qualityType_ = SolexaRead::QUAL_SANGER;}
  void setDefaultTrimType ()    {trimType_ = SolexaRead::TRIM_NONE;}
  void setDefaultMinQual ()     {minQual_ = 20;}
	
  void getLatestRead (SolexaRead & solexaRead) const {solexaRead = currentRead_;}
	void getCurrentRead (SolexaRead & solexaRead) const {solexaRead = currentRead_;}
  
  void nReadAdd (int inc){nRead_ = nRead_ + inc;}
  void nPassAdd (int inc){nPass_ = nPass_ + inc;}
  int getNPass() const {return nPass_;}
	int getNRead() const {return nRead_;}
	
  bool loadWholeFastq ();
  bool getNextSequence ();
  
  void printFastqStyle ();
  void printFastqStyle (FILE * output);
};

#endif /* FastqFile.h */
