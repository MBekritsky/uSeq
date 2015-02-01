/*
SolexaRead.cc

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

#include "SolexaRead.h"
#include <cstdio>
#include <cstdlib>

using namespace std;

SolexaRead::SolexaRead ()
{
	hdr_ . clear();
	seq_ . clear();
	qual_ . clear();
	oqual_ . clear();
	lSeq_ = 0;
	readThru_ = 0;
	nCount_ = 0;
	
	qualityType_ = QUAL_SANGER;
	trimType_ = TRIM_NONE;
	setDefaultMinQual();
	setDefaultDelims();
	blankMicrosatelliteInfo();	

	readWarning_ = 0;
}

SolexaRead::SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType)
{
	initRead (hdr, seq, qual, qualityType, trimType);
	setDefaultDelims();
	blankMicrosatelliteInfo();
}

SolexaRead::SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, int minQual)
{
	initRead(hdr, seq, qual, qualityType, trimType, minQual);
	setDefaultDelims();
	blankMicrosatelliteInfo();
}

SolexaRead::SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, unsigned int barcodeLength)
{
	initRead (hdr, seq, qual, qualityType, trimType,barcodeLength);
	setDefaultDelims();
	blankMicrosatelliteInfo();
}

SolexaRead::SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, int minQual, unsigned int barcodeLength)
{
	initRead(hdr, seq, qual, qualityType, trimType, minQual, barcodeLength);
	setDefaultDelims();
	blankMicrosatelliteInfo();
}

void SolexaRead::initContainer (Quality_t qualityType, Trim_t trimType, int minQual)
{
	hdr_   . clear();
	seq_   . clear();
	qual_  . clear();
	oqual_ . clear();

	lSeq_ 		= 0;
	readThru_ = 0;
	nCount_ 	= 0;

	qualityType_ = qualityType;
	trimType_ 	 = trimType;
	minQual_  	 = minQual;
	
	blankMicrosatelliteInfo();
	setDefaultDelims();

	readWarning_ = 0;
}

void SolexaRead::initContainer (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD)
{
	initContainer (quality_t, trim_t, minQual);
	setIRD(IRD);
	setIFD(IFD);
}

void SolexaRead::initContainer (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD, char ORD, char OFD)
{
	initContainer (quality_t, trim_t, minQual);
	setIRD(IRD);
	setIFD(IFD);
	setORD(ORD);
	setOFD(OFD);
}

void SolexaRead::initContainer (Quality_t qualityType, Trim_t trimType, int minQual, unsigned int barcodeLength)
{
	hdr_   . clear();
	seq_   . clear();
	qual_  . clear();
	oqual_ . clear();

	lSeq_ 		= 0;
	readThru_ = 0;
	nCount_ 	= 0;

	qualityType_ = qualityType;
	trimType_ 	 = trimType;
	minQual_  	 = minQual;
	barcodeLength_ = barcodeLength;
	
	blankMicrosatelliteInfo();
	setDefaultDelims();

	readWarning_ = 0;
}

void SolexaRead::initContainer (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD, unsigned int barcodeLength)
{
	initContainer (quality_t, trim_t, minQual,barcodeLength);
	setIRD(IRD);
	setIFD(IFD);
}

void SolexaRead::initContainer (Quality_t quality_t, Trim_t trim_t, int minQual, char IRD, char IFD, char ORD, char OFD, unsigned int barcodeLength)
{
	initContainer (quality_t, trim_t, minQual,barcodeLength);
	setIRD(IRD);
	setIFD(IFD);
	setORD(ORD);
	setOFD(OFD);
}

SolexaRead::SolexaRead(const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, const vector< string > & unparsedMicrosatellites)
{
	initRead (hdr, seq, qual, qualityType, trimType);
	setDefaultDelims();
	addMicrosatellites (unparsedMicrosatellites);
}

SolexaRead::SolexaRead(const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, const vector< string > & unparsedMicrosatellites, char IRD, char IFD)
{
	initRead (hdr, seq, qual, qualityType, trimType);
	setDefaultDelims();
	setIRD(IRD);
	setIFD(IFD);
	addMicrosatellites (unparsedMicrosatellites);
}

SolexaRead::SolexaRead(const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, int minQual, const vector< string > & unparsedMicrosatellites)
{
	initRead (hdr, seq, qual, qualityType, trimType, minQual);
	setDefaultDelims();
	addMicrosatellites (unparsedMicrosatellites);
}

SolexaRead::SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, int minQual, const vector< string > & unparsedMicrosatellites, char IRD, char IFD)
{
	initRead (hdr, seq, qual, qualityType, trimType, minQual);
	setDefaultDelims();
	setIRD(IRD);
	setIFD(IFD);
	addMicrosatellites (unparsedMicrosatellites);
}

SolexaRead::SolexaRead(const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, const vector< Microsatellite > & parsedMicrosatellites)
{
	initRead (hdr, seq, qual, qualityType, trimType);
	microsatellites_ = parsedMicrosatellites;
	nMicrosatellite_ = parsedMicrosatellites . size();
	setDefaultDelims();
}

SolexaRead::SolexaRead(const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, int minQual, const vector< Microsatellite > & parsedMicrosatellites)
{
	initRead (hdr, seq, qual, qualityType, trimType, minQual);
	microsatellites_ = parsedMicrosatellites;
	nMicrosatellite_ = parsedMicrosatellites . size();
	setDefaultDelims();
}

SolexaRead::SolexaRead(const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, const vector< string > & unparsedMicrosatellites, unsigned int barcodeLength)
{
	initRead (hdr, seq, qual, qualityType, trimType, barcodeLength);
	setDefaultDelims();
	addMicrosatellites (unparsedMicrosatellites);
}

SolexaRead::SolexaRead(const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, const vector< string > & unparsedMicrosatellites, char IRD, char IFD, unsigned int barcodeLength)
{
	initRead (hdr, seq, qual, qualityType, trimType, barcodeLength);
	setDefaultDelims();
	setIRD(IRD);
	setIFD(IFD);
	addMicrosatellites (unparsedMicrosatellites);
}

SolexaRead::SolexaRead(const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, int minQual, const vector< string > & unparsedMicrosatellites, unsigned int barcodeLength)
{
	initRead (hdr, seq, qual, qualityType, trimType, minQual, barcodeLength);
	setDefaultDelims();
	addMicrosatellites (unparsedMicrosatellites);
}

SolexaRead::SolexaRead (const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, int minQual, const vector< string > & unparsedMicrosatellites, char IRD, char IFD, unsigned int barcodeLength)
{
	initRead (hdr, seq, qual, qualityType, trimType, minQual,barcodeLength);
	setDefaultDelims();
	setIRD(IRD);
	setIFD(IFD);
	addMicrosatellites (unparsedMicrosatellites);
}

SolexaRead::SolexaRead(const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, const vector< Microsatellite > & parsedMicrosatellites, unsigned int barcodeLength)
{
	initRead (hdr, seq, qual, qualityType, trimType, barcodeLength);
	microsatellites_ = parsedMicrosatellites;
	nMicrosatellite_ = parsedMicrosatellites . size();
	setDefaultDelims();
}

SolexaRead::SolexaRead(const string & hdr, const string & seq, const string & qual, Quality_t qualityType, Trim_t trimType, int minQual, const vector< Microsatellite > & parsedMicrosatellites, unsigned int barcodeLength)
{
	initRead (hdr, seq, qual, qualityType, trimType, minQual, barcodeLength);
	microsatellites_ = parsedMicrosatellites;
	nMicrosatellite_ = parsedMicrosatellites . size();
	setDefaultDelims();
}


SolexaRead::SolexaRead(const SolexaRead & rhs)
{
	hdr_  				= rhs . hdr();
	seq_  				= rhs . seq();
	qual_ 				= rhs . qual();
	oqual_  			= rhs . oQual();
	readThru_ 		= rhs . readThru();
	lSeq_ 				= rhs . lSeq();
	qualityType_  = rhs . qualityType();
	trimType_ 		= rhs . trimType();
	readWarning_  = rhs . readWarning();
	minQual_  		= rhs . minQual();
	nCount_ 			= rhs . nCount();
	passFilter_		= rhs . passed();
	readNumber_		= rhs . readNumber();
	OFD_  				= rhs . OFD();
	ORD_					= rhs . ORD();
	IFD_					= rhs . IFD();
	IRD_					= rhs . IRD();
	barcodeLength_ = rhs.barcodeLength();

	nMicrosatellite_   = rhs . nUsat();
	rhs . getMicrosatellites(microsatellites_);
}

SolexaRead & SolexaRead::operator=(const SolexaRead & rhs)
{
	if(this != &rhs)
		{
			hdr_  				= rhs . hdr();
			seq_  				= rhs . seq();
			qual_ 				= rhs . qual();
			oqual_  			= rhs . oQual();
			readThru_ 		= rhs . readThru();
			lSeq_ 				= rhs . lSeq();
			qualityType_  = rhs . qualityType();
			trimType_ 		= rhs . trimType();
			readWarning_  = rhs . readWarning();
			minQual_  		= rhs . minQual();
			nCount_ 			= rhs . nCount();
			passFilter_ 	= rhs . passed();
			readNumber_		= rhs . readNumber();
			OFD_					= rhs . OFD();
			ORD_					= rhs . ORD();
			IFD_					= rhs . IFD();
			IRD_					= rhs . IRD();
			barcodeLength_ = rhs . barcodeLength();
			
			nMicrosatellite_   = rhs . nUsat();
			rhs . getMicrosatellites(microsatellites_);
		}
	return *this;
}

void SolexaRead::initRead(const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t)
{
	setHdr(hdr);
	setSeq(seq);
	setQual(qual);
	
	qualityType_ = quality_t;
	trimType_ = trim_t;
	setDefaultMinQual();
	
	readWarning_ = 0;
	
	matchLengths();
}

void SolexaRead::initRead(const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual)
{
	initRead(hdr, seq, qual, quality_t, trim_t);
	minQual_ = minQual;
}

void SolexaRead::initRead(const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, unsigned int barcodeLength)
{
	setHdr(hdr);
	setSeq(seq);
	setQual(qual);
	
	qualityType_ = quality_t;
	trimType_ = trim_t;
	barcodeLength_ = barcodeLength;
	setDefaultMinQual();
	
	readWarning_ = 0;
	
	
	matchLengths();
	trimBarcode();
}

void SolexaRead::initRead(const string & hdr, const string & seq, const string & qual, Quality_t quality_t, Trim_t trim_t, int minQual, unsigned int barcodeLength)
{
	initRead(hdr, seq, qual, quality_t, trim_t, barcodeLength);
	minQual_ = minQual;
}

void SolexaRead::blankMicrosatelliteInfo()
{
	nMicrosatellite_ = 0;
	microsatellites_ . clear();
}

void SolexaRead::setDefaultDelims()
{
	setDefaultORD();
	setDefaultOFD();
	setDefaultIRD();
	setDefaultIFD();
}

void SolexaRead::setQual(const string & qual)
{
	qual_ = qual;
	oqual_ = qual;
	adjustQualityScores();
}

void SolexaRead::setSeq(const string & seq)
{
	seq_ = seq;
	readThru_ = seq . length();
	lSeq_ = (int) seq . length();
	setNCount();
}

void SolexaRead::setHdr(const string & hdr)
{
	hdr_ = hdr;
	parseHdr();
	compactHdr();
}

void SolexaRead::setMicrosatellites(const vector< Microsatellite> & microsatellites)
{
  microsatellites_ = microsatellites;
	checkMicrosatelliteCoords();
  nMicrosatellite_ = microsatellites_ . size();
	checkMicrosatelliteCoords();
}

void SolexaRead::checkMicrosatelliteCoords()
{
	for(unsigned int i = 0; i < microsatellites_ . size(); i++)
	{
		if((microsatellites_ [i] . getStart()) == 1 || (microsatellites_ [i] . getStop()) == readThru_)
			microsatellites_ [i] . fail();
	}
}

void SolexaRead::compactHdr()
{
	size_t pos = 0;
	while((pos = hdr_ . find(" ", pos)) != string::npos)
		hdr_ . replace(pos, 1, "_");
}

void SolexaRead::parseHdr()
{
	size_t first_space = hdr_ . find_first_of(" ");
	if (first_space != string::npos)
		{
			readNumber_ = atoi(hdr_ . substr(first_space + 1,1) . c_str());
			if(hdr_ . substr(first_space + 3, 1) . compare("Y") == 0)
				{
					passFilter_ = false;
				}
			else
				{
					passFilter_ = true;
				}
			hdr_ = hdr_ . substr(0,first_space);
		}
	else
		{
			size_t first_slash = hdr_ . find_first_of("/");
			readNumber_ = atoi(hdr_ . substr(hdr_ . length() - 1, 1) . c_str());
			hdr_ = hdr_ . substr(0,first_slash);
			passFilter_ = true;
		}	
}

void SolexaRead::setNCount()
{
	size_t nLoc = 0;
	nCount_ = 0;
	do
		{
			nLoc = seq_ . find('N',nLoc);
			if(nLoc != string::npos)
				{
					nCount_++;	
					nLoc++;
				}
		} while(nLoc != string::npos);
}

void SolexaRead::printReadInfo()
{
  printf ("%s\n%s\n%s\n%s\n%s\n%s\n", hdr_ . c_str(), getGoodSeq() . c_str(), seq_ . c_str(), getGoodQual() . c_str(), qual_ . c_str(), oqual_ . c_str());
}

void SolexaRead::printFastqStyle()
{
	printFastqStyle(stdout);
}

void SolexaRead::printFastqStyle(FILE * output)
{
  fprintf (output,"@%s\n%s\n+%s\n%s\n", hdr_ . c_str(), seq_ . c_str(), hdr_ . c_str(), oqual_ . c_str());
}

string SolexaRead::displayTrimType() const
{
  switch(trimType_)
    {
    case TRIM_BWA_STYLE:
      return "BWA style trimming";
    case TRIM_LAST_GOOD:
      return "Trimming to last good nucleotide";
    case TRIM_FIRST_BAD:
      return "Trimming to first bad nucleotide";
    case TRIM_NONE:
    default:
      return "No trimming applied";
    }
}

string SolexaRead::getGoodSeq() const
{
	return getGoodPortion (seq_);
}

string SolexaRead::getGoodQual() const
{
	return getGoodPortion (oqual_);
}

string SolexaRead::getGoodPortion(const string & text) const
{
	return text.substr (0, readThru_);
}

char SolexaRead::cIllumina13ToCSanger(char cIllumina)
{
  return (char) int( -10 * log10 (1 / (pow (10.0, ( (cIllumina - 64.0) / 10.0) ) + 1.0)) + 33);
}

void SolexaRead::adjustQualityScores()
{
  if(qualityType_ == QUAL_SOLEXA13)
  {
    string qual;
    for(unsigned int i = 0; i < qual_ . length(); i++)
      qual . push_back (cIllumina13ToCSanger (qual_ [i]));

    qual_ = qual;
  }
}

void SolexaRead::matchLengths()
{
  if (lSeq_ != qual_ . length())
    printReadWarning(1);
}

void SolexaRead::trim()
{
  switch(trimType_)
    {
    case TRIM_NONE:
      break;
    case TRIM_BWA_STYLE:
      trimBwaStyle();
      break;
    case TRIM_FIRST_BAD:
      trimToFirstBad();
      break;
    case TRIM_LAST_GOOD:
      trimToLastGood();
      break;
    }
}

void SolexaRead::trimBwaStyle()
{
  int pos = lSeq_ - 1, maxPos = lSeq_;
  int area = 0, maxArea = 0;

  while(pos >= 0 && area >= 0)
  {
    area += minQual_ - (qual_ [pos] - 33);
    if(area > maxArea)
    {
      maxArea = area;
      maxPos = pos;
    }
    pos--;
  }
  if(pos == 0)
  {
    readThru_ = 0;
  }
  else
  {
    readThru_ = maxPos;
  }
}

void SolexaRead::trimToFirstBad()
{
  unsigned int firstBad = 0;
  while (firstBad <  lSeq_ && (qual_ [firstBad] - 33) >= minQual_)
    firstBad++;

  readThru_ = firstBad;
}

void SolexaRead::trimToLastGood()
{
  int lastGood = lSeq_;
  while (lastGood >= 0 && (qual_ [lastGood] - 33) < minQual_)
    lastGood--;
  lastGood++;

  readThru_ = lastGood;
}

void SolexaRead::trimBarcode()
{
	if(barcodeLength_ < seq_ . length())
	{
		seq_ = seq_.substr(barcodeLength_);
		qual_  = qual_ . substr(barcodeLength_);
		oqual_ = oqual_ . substr(barcodeLength_);
		readThru_ = readThru_ - barcodeLength_;
		lSeq_ = lSeq_ - barcodeLength_;
		setNCount();
	}
	else //if the sequence length is shorter than the barcode, truncate sequence to a single N with low quality, lSeq_ and readThru_ go to zero
	{
		seq_ = "N";
		qual_ = "#";
		oqual_ = "#";
		readThru_ = 0;
		lSeq_ = 0;
		setNCount();
	}

}

void SolexaRead::printReadWarning(int warnCode)
{
  readWarning_ = true;
  cerr << "Warning!!! ";
  switch(warnCode)
    {
    case 1:
      cerr << "The sequence and quality lines of read " << hdr_ . c_str() << " have different lengths!" << endl;
      break;

    default:
      cerr << "There is an unidentified problem with read " << hdr_ . c_str() << endl;
      break;
    }
}

void SolexaRead::addMicrosatellites(const vector< string > & unparsedMicrosatellites)
{
  for(unsigned int i = 0; i < unparsedMicrosatellites . size(); i++)
    microsatellites_ . push_back ( Microsatellite (getHdr(), unparsedMicrosatellites[i], IFD_, OFD_) );
	nMicrosatellite_ = unparsedMicrosatellites . size();
	
	checkMicrosatelliteCoords();
}

unsigned int SolexaRead::getNumPassedMicrosatellites()
{
	unsigned int count = 0;
	for(unsigned int i = 0; i < nMicrosatellite_; i++)
		{
			if(microsatellites_ [i] . passed())
				count++;
		}
	return count;
}

void SolexaRead::msMean(Microsatellite &microsatellite)
{
  unsigned int start = (microsatellite . getStart() - 1), stop = microsatellite . getStop(), length = microsatellite . getLength();
  double sum = 0;
  for(unsigned int i = start; i < stop; i++)
    sum += qual_ [i] - 33;

  microsatellite . setMean( (sum/length));
}

void SolexaRead::msStdDev(Microsatellite &microsatellite)
{
  msMean(microsatellite);
  unsigned int start = (microsatellite . getStart() - 1), stop = microsatellite . getStop(), length = microsatellite . getLength();

  double sum = 0;
  for(unsigned int i = start; i < stop; i++)
    sum += pow( ( (qual_ [i] - 33) - microsatellite . getMean()), 2);
  microsatellite . setStdDev ( pow((sum/length),0.5) );
}

void SolexaRead::printPassedMicrosatellitesFastq()
{
	printPassedMicrosatellitesFastq(stdout);
}

void SolexaRead::printPassedMicrosatellitesFastq(FILE *output)
{
	if(readThru_ > 0)
		{
		  fprintf(output,"@%s", hdr_ . c_str());
		  char buffer[1024];
		  for(unsigned int i = 0; i < nMicrosatellite_; i++)
		    {
					if(microsatellites_ [i] . passed())
						{
				      microsatellites_ [i] . printNoHdr(buffer);
				      fprintf(output,"%c%s",ORD_, buffer);
    				}
				}
		  fprintf(output,"\n");
		  fprintf(output,"%s\n+%s\n%s\n", seq_ . c_str(), hdr_ . c_str(), oqual_ . c_str());
		}
}

void SolexaRead::printMicrosatelliteFastq()
{
	printMicrosatelliteFastqSkeleton(NULL, stdout);
}

void SolexaRead::printMicrosatelliteFastq(FILE * output)
{
	printMicrosatelliteFastqSkeleton(NULL, output);
}

void SolexaRead::printMicrosatelliteFastqWithQual(const string & qual, FILE * output)
{
	printMicrosatelliteFastqSkeleton(qual, output);
}

void SolexaRead::printMicrosatelliteFastqSkeleton(const string & qual, FILE * output)
{
	if(readThru_ > 0)
		{
			fprintf(output, "@%s", hdr_ . c_str());
			char buffer[1024];
			for(unsigned int i = 0; i < microsatellites_ . size(); i++)
				{
					microsatellites_ [i] . printNoHdr(buffer, qual);
					fprintf(output, "%c%s", ORD_, buffer);
				}
			fprintf(output, "\n");
			fprintf (output, "%s\n+%s\n%s\n", getGoodSeq() . c_str(), hdr_ . c_str(), getGoodQual() . c_str());
		}
}

void SolexaRead::printModifiedRead()
{
	printModifiedReadSkeleton(NULL,stdout, true);
}

void SolexaRead::printModifiedRead(FILE * output)
{
	printModifiedReadSkeleton(NULL,output, true);
}

void SolexaRead::printModifiedReadWithQual(const string & qual, FILE * output)
{
	printModifiedReadSkeleton(qual, output, true);
}

void SolexaRead::printModifiedReadWithQualNoModHeader(const string & qual, FILE * output)
{
	printModifiedReadSkeleton(qual, output, false);
}

void SolexaRead::printMicrosatelliteHeader(const string & qual, FILE * output)
{
	char buffer[1024];
	fprintf (output, "%s", hdr_ . c_str());
	for(unsigned int i = 0; i < nMicrosatellite_; i++)
		{
			microsatellites_ [i] . printNoHdr(buffer, qual);
			fprintf(output, "%c%s", ORD_, buffer);
		}
	fprintf(output,"\n");
}

void SolexaRead::printMicrosatelliteHeader(const string & qual, FILE * output, bool isFirst)
{
	char buffer[1024];
	int pairID;
	if(isFirst)
		{
			pairID = 1;
		}
	else
		{
			pairID = 2;
		}
		
	fprintf (output, "%s/%d", hdr_ . c_str(),pairID);
	for(unsigned int i = 0; i < nMicrosatellite_; i++)
		{
			microsatellites_ [i] . printNoHdr(buffer, qual);
			fprintf(output, "%c%s", ORD_, buffer);
		}
	fprintf(output,"\n");
}

void SolexaRead::printModifiedReadSkeleton(const string & qual, FILE * output, bool printMicrosatelliteInfo)
{
//A modified read will note all passing microsatellites in the read header, but it will 
//excise all microsatellite sequences and quality scores for any microsatellite in the 
//read, including terminal microsatellites

	char buffer[2048];
	string modQual, modSeq;
	modSeq . clear();
	modQual . clear();

	vector< unsigned int > pos (nMicrosatellite_ + 1);
	vector< unsigned int > length (nMicrosatellite_ + 1);

	unsigned int readEnd = readThru_;
	unsigned int nOverlap = 0;
	
	for(unsigned int i = 0, j = 0; i < nMicrosatellite_; i++, j++)
		{
			//if MS overlap or MS starts at beginning
			if(microsatellites_ [i] . getStart() - 1 <= pos [j])
				{
					pos[j] = microsatellites_ [i] . getStop();
					pos . resize(pos . size() - 1);
					length . resize(length . size() - 1);
					nOverlap++;
					j--;
				}
			else
				{
					length [j] = (microsatellites_ [i] . getStart()) - pos [j] - 1;
					pos [j + 1] = microsatellites_ [i] . getStop();
				}
			if(microsatellites_ [i] . getStop() == readEnd)
				{
					pos . resize (pos . size() - 1);
					length . resize (length . size() - 1);
					//if the 3' terminal microsatellite is also 5' terminal, then there is no microsatellites_ [i - nOverlap]
					if(microsatellites_ [i] . getStart() != 1)
						readEnd = microsatellites_ [i - nOverlap] . getStart() - 1;
				}
			nOverlap = 0;
		}

	if(pos . size() > 0)
		{
			length [pos . size() - 1] = (readEnd - pos [pos . size() - 1]);
				
			//print header with all microsatellite information and excise MS sequences and quality scores
			for(unsigned int i = 0; i < pos . size(); i++)
				{
					modQual . append (oqual_ . substr (pos[i], length[i]));
					modSeq . append(seq_ . substr (pos[i], length[i]));
				}
				
			if(modSeq . length() > 0)
				{
			   fprintf(output, "@%s", hdr_ . c_str());
				 if(printMicrosatelliteInfo)
					{
				   for(unsigned int i = 0; i < nMicrosatellite_; i++)
  			 		 {
							microsatellites_ [i] . printNoHdr (buffer, qual);
							fprintf (output, "%c%s", ORD_, buffer);
  		 		 		}
					}
			   fprintf(output,"\n");
  			 fprintf(output,"%s\n+%s\n%s\n", modSeq . c_str(), hdr_ . c_str(), modQual . c_str());
				}
			else
				{
					fprintf(output, "@%s\nN\n+%s\nB\n", hdr_ . c_str(), hdr_ . c_str());
				}
  	}
	else
		{
			fprintf(output, "@%s\nN\n+%s\nB\n", hdr_ . c_str(), hdr_ . c_str());
		}
}

