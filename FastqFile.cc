#include "FastqFile.h"

using namespace std;
#include <cstdio>
#include <cstdlib>

FastqFile::FastqFile()
{
  solexaFileP_ = NULL;
  reads_ . clear();
  setDefaultQualityType();
  setDefaultTrimType();
	setDefaultMinQual();
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_);
	setDefaultDelims();
	blankCounters();
}

FastqFile::FastqFile(FILE * fp)
{
  solexaFileP_ = fp;
  setDefaultQualityType();
  setDefaultTrimType();
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_);
	setDefaultDelims();
	blankCounters();
}

FastqFile::FastqFile(FILE * fp, SolexaRead::Quality_t qualityType)
{
  solexaFileP_ = fp;
  qualityType_ = qualityType;
  setDefaultTrimType();
  setDefaultMinQual();
	setDefaultDelims();
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_);
	blankCounters();
}

FastqFile::FastqFile(FILE * fp, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType)
{
  solexaFileP_ = fp;
  qualityType_ = qualityType;
  trimType_ = trimType;
  setDefaultMinQual();
	setDefaultDelims();
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_);
	blankCounters();
}

FastqFile::FastqFile(FILE * fp, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual)
{
  solexaFileP_ = fp;
  qualityType_ = qualityType;
  trimType_ = trimType;
  minQual_ = minQual;
	setDefaultDelims();
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_);
	blankCounters();
}

FastqFile::FastqFile(FILE * fp, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual, char IRD, char IFD)
{
	solexaFileP_ = fp;
	qualityType_ = qualityType;
	trimType_ = trimType;
	minQual_ = minQual;
	setDefaultDelims();
	setIRD(IRD);
	setIFD(IFD);
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_, IRD, IFD);
	blankCounters();
}

FastqFile::FastqFile(FILE * fp, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual, char IRD, char IFD, char ORD, char OFD)
{
	solexaFileP_ = fp;
	qualityType_ = qualityType;
	trimType_ = trimType;
	minQual_ = minQual;
	setIRD(IRD);
	setIFD(IFD);
	setORD(ORD);
	setOFD(OFD);
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_, IRD, IFD, ORD, OFD);
	blankCounters();
}

FastqFile::FastqFile(const string & filename)
{
  solexaFileP_ = openAndTestFile(filename, "r");
  setDefaultQualityType();
  setDefaultTrimType();
  setDefaultMinQual();
	setDefaultDelims();
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_);
	blankCounters();
}

FastqFile::FastqFile(const string & filename, SolexaRead::Quality_t qualityType)
{
  solexaFileP_ = openAndTestFile(filename, "r");
  qualityType_ = qualityType;
  setDefaultTrimType();
  setDefaultMinQual();
	setDefaultDelims();
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_);
	blankCounters();
}

FastqFile::FastqFile(const string & filename, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType)
{
  solexaFileP_ = openAndTestFile (filename, "r");
  qualityType_ = qualityType;
  trimType_ = trimType;
  setDefaultMinQual();
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_);
	blankCounters();
}

FastqFile::FastqFile(const string & filename, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual)
{
  solexaFileP_ = openAndTestFile(filename, "r");
  qualityType_ = qualityType;
  trimType_ = trimType;
  minQual_ = minQual;
	setDefaultDelims();
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_);
	blankCounters();
}

FastqFile::FastqFile(const string & filename, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual, char IRD, char IFD)
{
	solexaFileP_ = openAndTestFile(filename, "r");
	qualityType_ = qualityType;
	trimType_ = trimType;
	minQual_ = minQual;
	setDefaultDelims();
	setIRD(IRD);
	setIFD(IFD);
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_, IRD, IFD);
	blankCounters();
}

FastqFile::FastqFile(const string & filename, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual, char IRD, char IFD, char ORD, char OFD)
{
	solexaFileP_ = openAndTestFile(filename, "r");
	qualityType_ = qualityType;
	trimType_ = trimType;
	minQual_ = minQual;
	setIRD(IRD);
	setIFD(IFD);
	setORD(ORD);
	setOFD(OFD);
	barcodeLength_ = 0;
	currentRead_ . initContainer (qualityType_, trimType_ , minQual_, IRD, IFD, ORD, OFD);
	blankCounters();
}

FastqFile::FastqFile(FILE * fp, unsigned int barcodeLength)
{
  solexaFileP_ = fp;
  setDefaultQualityType();
  setDefaultTrimType();
	barcodeLength_ = barcodeLength;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_,barcodeLength_);
	setDefaultDelims();
	blankCounters();
}

FastqFile::FastqFile(FILE * fp, SolexaRead::Quality_t qualityType, unsigned int barcodeLength)
{
  solexaFileP_ = fp;
  qualityType_ = qualityType;
  setDefaultTrimType();
  setDefaultMinQual();
	setDefaultDelims();
	barcodeLength_ = barcodeLength_;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_,barcodeLength_);
	blankCounters();
}

FastqFile::FastqFile(FILE * fp, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, unsigned int barcodeLength)
{
  solexaFileP_ = fp;
  qualityType_ = qualityType;
  trimType_ = trimType;
  setDefaultMinQual();
	setDefaultDelims();
	barcodeLength_ = barcodeLength_;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_,barcodeLength_);
	blankCounters();
}

FastqFile::FastqFile(FILE * fp, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual, unsigned int barcodeLength)
{
  solexaFileP_ = fp;
  qualityType_ = qualityType;
  trimType_ = trimType;
  minQual_ = minQual;
	setDefaultDelims();
	barcodeLength_ = barcodeLength_;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_,barcodeLength_);
	blankCounters();
}

FastqFile::FastqFile(FILE * fp, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual, char IRD, char IFD, unsigned int barcodeLength)
{
	solexaFileP_ = fp;
	qualityType_ = qualityType;
	trimType_ = trimType;
	minQual_ = minQual;
	setDefaultDelims();
	setIRD(IRD);
	setIFD(IFD);
	barcodeLength_ = barcodeLength_;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_, IRD, IFD,barcodeLength_);
	blankCounters();
}

FastqFile::FastqFile(FILE * fp, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual, char IRD, char IFD, char ORD, char OFD, unsigned int barcodeLength)
{
	solexaFileP_ = fp;
	qualityType_ = qualityType;
	trimType_ = trimType;
	minQual_ = minQual;
	setIRD(IRD);
	setIFD(IFD);
	setORD(ORD);
	setOFD(OFD);
	barcodeLength_ = barcodeLength_;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_, IRD, IFD, ORD, OFD,barcodeLength_);
	blankCounters();
}

FastqFile::FastqFile(const string & filename, unsigned int barcodeLength)
{
  solexaFileP_ = openAndTestFile(filename, "r");
  setDefaultQualityType();
  setDefaultTrimType();
  setDefaultMinQual();
	setDefaultDelims();
	barcodeLength_ = barcodeLength;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_,barcodeLength_);
	blankCounters();
}

FastqFile::FastqFile(const string & filename, SolexaRead::Quality_t qualityType, unsigned int barcodeLength)
{
  solexaFileP_ = openAndTestFile(filename, "r");
  qualityType_ = qualityType;
  setDefaultTrimType();
  setDefaultMinQual();
	setDefaultDelims();
	barcodeLength_ = barcodeLength;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_,barcodeLength_);
	blankCounters();
}

FastqFile::FastqFile(const string & filename, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, unsigned int barcodeLength)
{
  solexaFileP_ = openAndTestFile (filename, "r");
  qualityType_ = qualityType;
  trimType_ = trimType;
  setDefaultMinQual();
	barcodeLength_ = barcodeLength;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_,barcodeLength_);
	blankCounters();
}

FastqFile::FastqFile(const string & filename, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual, unsigned int barcodeLength)
{
  solexaFileP_ = openAndTestFile(filename, "r");
  qualityType_ = qualityType;
  trimType_ = trimType;
  minQual_ = minQual;
	setDefaultDelims();
	barcodeLength_ = barcodeLength;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_,barcodeLength_);
	blankCounters();
}

FastqFile::FastqFile(const string & filename, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual, char IRD, char IFD, unsigned int barcodeLength)
{
	solexaFileP_ = openAndTestFile(filename, "r");
	qualityType_ = qualityType;
	trimType_ = trimType;
	minQual_ = minQual;
	setDefaultDelims();
	setIRD(IRD);
	setIFD(IFD);
	barcodeLength_ = barcodeLength;
	currentRead_ . initContainer (qualityType_, trimType_, minQual_, IRD, IFD,barcodeLength_);
	blankCounters();
}

FastqFile::FastqFile(const string & filename, SolexaRead::Quality_t qualityType, SolexaRead::Trim_t trimType, int minQual, char IRD, char IFD, char ORD, char OFD, unsigned int barcodeLength)
{
	solexaFileP_ = openAndTestFile(filename, "r");
	qualityType_ = qualityType;
	trimType_ = trimType;
	minQual_ = minQual;
	setIRD(IRD);
	setIFD(IFD);
	setORD(ORD);
	setOFD(OFD);
	barcodeLength_ = barcodeLength;
	currentRead_ . initContainer (qualityType_, trimType_ , minQual_, IRD, IFD, ORD, OFD,barcodeLength_);
	blankCounters();
}

FastqFile::~FastqFile()
{
  if (solexaFileP_)
  {
    fclose (solexaFileP_);  
  }
}

void FastqFile::setDefaultDelims()
{
	setDefaultIRD();
	setDefaultIFD();
	setDefaultORD();
	setDefaultOFD();
}

void FastqFile::blankCounters()
{
	nRead_ = 0;
	nPass_ = 0;
}

bool FastqFile::loadWholeFastq()
{
  //perform get next sequence operation until it returns false
  while (getNextSequence ())
	{
  	SolexaRead temp = currentRead_;
	  reads_ . push_back(temp);
	}
 //check that there is at least one loaded read
  if(reads_ . size() >= 1)
    return true;
  return false;
}

bool FastqFile::getNextSequence()
{
  char ch;
  string hdr, seq, qual;
  vector< string > microsatelliteStrings;
  string microsatelliteString;
  
  //Skip to next header (should never be a problem, but just in case)
  while ( (ch = fgetc (solexaFileP_)) != EOF && ch != '@')
    ;	
  if(ch == EOF)
    return false;
  
  //Load header
  while ( (ch = fgetc (solexaFileP_)) != EOF && ch != '\n' && ch != IRD_)
    hdr . push_back ( (char) ch);
  ungetc (ch, solexaFileP_);
  
  //Load microsatellite info (if any)
  while( (ch = fgetc (solexaFileP_)) != EOF && ch != '\n')
    {
      while ( (ch = fgetc (solexaFileP_)) != EOF && ch != '\n' && ch != IRD_)
				microsatelliteString . push_back ( (char) ch);
      ungetc(ch, solexaFileP_);
      
      microsatelliteStrings . push_back (microsatelliteString);
      microsatelliteString . clear();
    }
  
  //Skip any spaces trailing microsatellite(s)
  while( (ch = fgetc (solexaFileP_)) != EOF && ch == ' ')
    ;
  ungetc (ch, solexaFileP_);
  
  //Load sequence
  while ( (ch = fgetc (solexaFileP_)) != EOF && ch != '\n')
    seq . push_back ( (char) ch);
  
  //Skip quality line header
  while ( (ch = fgetc (solexaFileP_)) != EOF && ch != '\n')
    ;
  
  //Load base-calling qualities
  while ( (ch = fgetc (solexaFileP_)) != EOF && ch != '\n')
    qual . push_back ( (char) ch);
  
  currentRead_ . setHdr(hdr);
	currentRead_ . setSeq(seq);
	currentRead_ . setQual(qual);
	currentRead_ . addMicrosatellites(microsatelliteStrings);
	currentRead_ . trimBarcode();
  currentRead_ . trim();
  return true;
}

void FastqFile::printFastqStyle()
{
	printFastqStyle(stdout);
}

void FastqFile::printFastqStyle(FILE* output)
{
  for(unsigned int i = 0; i < reads_ . size(); i++)
    reads_ [i] . printFastqStyle (output);
}
