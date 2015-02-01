#include "MicrosatelliteDetector.h"

MicrosatelliteDetector::MicrosatelliteDetector()
{
  setDefaultMinLength();
  setDefaultMinUnits();
  setDefaultMinUnitLength();
  setDefaultMaxUnitLength();
	modifySeq_ = false;
	seq_ . clear();
	modSeq_ . clear();
}

MicrosatelliteDetector::MicrosatelliteDetector(unsigned int minLength, unsigned int minUnits, unsigned int minUnitLength, unsigned int maxUnitLength)
{
  setMinLength(minLength);
  setMinUnits(minUnits);
  setMinUnitLength(minUnitLength);
  setMaxUnitLength(maxUnitLength);
	modifySeq_ = false;
	seq_ . clear();
	modSeq_ . clear();
}

MicrosatelliteDetector::MicrosatelliteDetector(unsigned int minLength, unsigned int minUnits, unsigned int minUnitLength, unsigned int maxUnitLength, bool modifySeq)
{
	setMinLength(minLength);
	setMinUnits(minUnits);
	setMinUnitLength(minUnitLength);
	setMaxUnitLength(maxUnitLength);
	modifySeq_ = modifySeq;
	seq_ . clear();
	modSeq_ . clear();
	initModSeqIndices();
}

void MicrosatelliteDetector::initModSeqIndices()
{
	previousStop_ = 0;
	thisStart_ = 0;
}

void MicrosatelliteDetector::resetModSeqVars()
{
	modSeq_ . clear();
	previousStop_ = 0;
	thisStart_ = 0;
}

void MicrosatelliteDetector::setMinLength(unsigned int minLength)
{
  if(minLength >= 4)
    {
      minLength_ = minLength;
    }
  else
    {
      setDefaultMinLength();
      cerr << "The minimum microsatellite length must be >= 4!! Setting minimum length to default (" << defaultMinLength_ << " nt)" << endl;
    }
}

void MicrosatelliteDetector::setMinUnits(unsigned int minUnits)
{
  if(minUnits >= 2)
    {
      minUnits_ = minUnits;
    }
  else
    {
      setDefaultMinUnits();
      cerr << "The number of units in a microsatellite must be >= 2!! Setting minimum number of units to default (" << defaultMinUnits_ << " units)" << endl;
    }
}

void MicrosatelliteDetector::setMinUnitLength(unsigned int minUnitLength)
{
  if(minUnitLength >= 1 && minUnitLength <= 100)
    {
      minUnitLength_ = minUnitLength;
    }
  else
    {
      setDefaultMinUnitLength();
      cerr << "The minimum repeat unit length must be > 0 and <= 100!! Setting minimum repeat unit length to default (" << defaultMinUnitLength_ << " nt)" << endl;
    }
}

void MicrosatelliteDetector::setMaxUnitLength(unsigned int maxUnitLength)
{
  if(maxUnitLength >= 1 && maxUnitLength <= 100)
    {
      maxUnitLength_ = maxUnitLength;
    }
  else
    {
      setDefaultMaxUnitLength();
      cerr << "The maximum repeat unit length must be > 0 and <= 100!! Setting maximum repeat unit length to default (" << defaultMaxUnitLength_ << " nt)" << endl;
    }
}

void MicrosatelliteDetector::checkSensitivity()
{
	unsigned int minPossibleUnits, actualMinUnits, actualMinLength, unitLengthSensitivity;
	
	for(unsigned int i = minUnitLength_; i <= maxUnitLength_; i++)
		{
			minPossibleUnits = (unsigned int) ceil (minLength_ / i);
			actualMinUnits = max (minPossibleUnits, minUnits_);
			actualMinLength = max ( (actualMinUnits * i),minLength_);
			unitLengthSensitivity = max (actualMinLength, (maxUnitLength_ + i));
			if(actualMinLength < unitLengthSensitivity && unitLengthSensitivity > minLength_)
				{
					cerr << "Warning!! With the specified parameters, repeats with unit length " << i << " will only be detected if they are " << unitLengthSensitivity << " nt or longer" << endl;
//					cerr << ". The smallest max repeat unit length that can detect repeats with unit length " << i << " with length " << minLength_ << " is " << (minLength_ - i - 1) << "." << endl;
				}
		}
}

unsigned int MicrosatelliteDetector::intPow(unsigned int base, unsigned int power)
{
  unsigned int product = 1;

  for(unsigned int i = 1; i <= power; i++)
    product *= base;

  return product;
}

bool MicrosatelliteDetector::isNuc (char nuc)
{
  switch (nuc)
    {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
    case 'a':
    case 'c':
    case 'g':
    case 't':
      return 1;
      
    default:
      return 0;
    }
}

int MicrosatelliteDetector::convertNucToNumber (char nuc)
{
  switch (nuc)
    {
    case 'A':
    case 'a':
      return 0;

    case 'C':
    case 'c':
      return 1;

    case 'G':
    case 'g':
      return 2;

    case 'T':
    case 't':
      return 3;

    default:
      return -1;
    }
}

int MicrosatelliteDetector::getKmerCode (const unsigned int stop)
{
	int code = 0;
		
	//addition of 1 to i and 2 to stop-maxUnitLength_ is just bookkeeping to keep the numbers above -1 so I
	//can keep i as an unsigned int
	
	shiftNucNum();
	nucNum_ [maxUnitLength_ - 1] = convertNucToNumber( (char) seq_ [stop]);
	
	for(unsigned int i = 0, j = maxUnitLength_ - 1; i < maxUnitLength_; i++, j--)
		{
			if(nucNum_[i] == -1)
				return -2;
			code += intPow(4, j) * nucNum_ [i];
		}
/* 		
	for (unsigned int i = stop, j = (maxUnitLength_ - 1); (i + 1) >= (stop - maxUnitLength_ + 2); i--, j--)
		{
			if (!isNuc (seq_ [i]))
				return -2;
			code += intPow(4, j) * convertNucToNumber ( (char) seq_ [i]);
		}
 */ return code;
}

void MicrosatelliteDetector::shiftNucNum()
{
	for(unsigned int i = 0; i < maxUnitLength_ - 1; i++)
		nucNum_ [i] = nucNum_ [i + 1];
}

void MicrosatelliteDetector::shiftKmerIndex()
{
  for(unsigned int i = 0; i < maxUnitLength_ - 1; i++)
    kmerIndex_ [i] = kmerIndex_ [i + 1];
}

void MicrosatelliteDetector::initKmerIndex()
{
  for(unsigned int i = 0; i < maxUnitLength_; i++)
    kmerIndex_ [i] = -1;
}

void MicrosatelliteDetector::initNucNum(unsigned int pos)
{
	nucNum_[0]  = -1;
	for(unsigned int i = pos, j = 1; j < maxUnitLength_; i++, j++)
		{
			nucNum_ [j] = convertNucToNumber(seq_ [i]);
		}
}

vector< Microsatellite > MicrosatelliteDetector::detectPerfectMicrosatellites(const string & seq, const string & tag, unsigned int flankLength)
{
  vector< Microsatellite > foundMicrosatellites;
	seq_ = seq;

  seqLength_ = seq_ . length();
  indexToCheck_ = -1, count_ = 0, repeatLength_ = 0, nPeat_ = 0;
	previousStop_ = 0, thisStart_ = 0, currCode_ = 0;
  unsigned int scanStart = 0;

  while (scanStart < seqLength_ && !isNuc (seq [scanStart]))
    scanStart++;
  initKmerIndex();
	initNucNum(scanStart);

  for(unsigned int i = scanStart + maxUnitLength_ - 1; i < seqLength_; i++ )
    {
      //Skip to next called nucleotide
      while(!isNuc (seq [i]))
				{
					//If there was a microsatellite ending at the beginning of the N-run, check if it
					//passes the thresholds
				  if(indexToCheck_ >= 0)
				    detectMicrosatelliteEnd (foundMicrosatellites, tag, seq, i, flankLength);
				  i++;
				  while ( (i + maxUnitLength_ - 1) < seqLength_ && ! isNuc (seq [i]))
	    			i++;
					initNucNum(i);
				  if ( (i + maxUnitLength_ - 1) >= seqLength_)
				    break;

				  //Reset kmer index
				  initKmerIndex();
				  i += maxUnitLength_ - 1;
				}
/*			cout << i << "\t";
			cout << seq_ [i] << "\t";
*/     
  		currCode_ = getKmerCode(i);
			
/*			for(unsigned int j = 0; j < maxUnitLength_; j++)
				cout << kmerIndex_ [j]<< ";";
			cout << "\t";
			cout << "\t" << currCode_ << "\t" << indexToCheck_ << endl;
*/			
      if (indexToCheck_ == -1)
				{
				  detectMicrosatelliteStart();
				  if (indexToCheck_ >= 0)
				    {
	    			  initMicrosatelliteTrackingVariables();
				    }
				}
      else
				{
				  if (kmerIndex_ [indexToCheck_] == currCode_)
	    			{
				      count_++;
				    }
				  else
				    {
				      detectMicrosatelliteEnd (foundMicrosatellites, tag, seq, i, flankLength);
				    }
				}
      shiftKmerIndex();
      kmerIndex_[maxUnitLength_ - 1] = currCode_;
    }
	//This is in case there is a microsatellite at the very end of the sequence
	if(indexToCheck_ != -1)
		detectMicrosatelliteEnd(foundMicrosatellites, tag, seq, seqLength_, flankLength);
  return foundMicrosatellites;
}

void MicrosatelliteDetector::detectMicrosatelliteEnd(vector< Microsatellite > & microsatellites, const string tag, const string seq, unsigned int & index, unsigned int flankLength)
{
  repeatLength_ = maxUnitLength_ - indexToCheck_;
  //Check that repeat has the proper minimum unit length
  if(repeatLength_ >= minUnitLength_)
    {
      count_ += maxUnitLength_ + repeatLength_ - 1;
      //Check that microsatellite meets threshold for length and repeat units
      if(count_ >= minLength_ && (count_/repeatLength_) >= minUnits_)
				{
				  Microsatellite microsatellite (tag, index, count_, repeatLength_, seq, flankLength);
				  microsatellites . push_back (microsatellite);
				}
    }
  resetMicrosatelliteTrackingVariables();
}

void MicrosatelliteDetector::resetMicrosatelliteTrackingVariables()
{
  indexToCheck_ = -1;
  count_ = 0;
  repeatLength_ = 0;
}

void MicrosatelliteDetector::initMicrosatelliteTrackingVariables()
{
  count_ = 1;
  repeatLength_ = maxUnitLength_ - indexToCheck_;
}

void MicrosatelliteDetector::detectMicrosatelliteStart()
{
  for(unsigned int i = 0; i < maxUnitLength_; i++)
    {
      if(kmerIndex_ [i] == currCode_ && currCode_ >= 0)
					indexToCheck_ = i;
	  }
}
