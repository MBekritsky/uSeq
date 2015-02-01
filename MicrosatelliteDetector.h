#ifndef __MICROSATELLITE_DETECTOR_H__
#define __MICROSATELLITE_DETECTOR_H__

#include "Microsatellite.h"
#include <cmath>
#include <algorithm>

class MicrosatelliteDetector
{
 private:
  static const int BUFFER_SIZE_ = 100;
  unsigned int minLength_, minUnits_, minUnitLength_, maxUnitLength_;
	string seq_;
  unsigned int seqLength_;
  int kmerIndex_ [BUFFER_SIZE_];
	short int nucNum_ [BUFFER_SIZE_];
	int indexToCheck_, currCode_;
  unsigned int count_, repeatLength_, nPeat_;
  static const unsigned int defaultMinLength_ = 8, defaultMinUnits_ = 3, defaultMinUnitLength_ = 1, defaultMaxUnitLength_ = 4;
	bool modifySeq_;
	unsigned int previousStop_, thisStart_;
	string modSeq_;

	void initNucNum(unsigned int pos);
	void shiftNucNum();
  void shiftKmerIndex ();
  void initKmerIndex ();
  bool isNuc (char nuc);
  int convertNucToNumber (char nuc);
	int getKmerCode (const unsigned int stop);
  void detectMicrosatelliteStart();
  void detectMicrosatelliteEnd(vector< Microsatellite > & microsatellites, const string tag, const string seq, unsigned int & index, unsigned int flankLength);
  void initMicrosatelliteTrackingVariables();
  void resetMicrosatelliteTrackingVariables();
	void initModSeqIndices();
	void resetModSeqVars();
  unsigned int intPow(unsigned int base, unsigned int power);

 public:
  MicrosatelliteDetector ();
  MicrosatelliteDetector (unsigned int minLength, unsigned int minUnits, unsigned int minUnitLength, unsigned int maxUnitLength);
	MicrosatelliteDetector (unsigned int minLenght, unsigned int minUnits, unsigned int minUnitLength, unsigned int maxUnitLength, bool modifySeq);

  void displayParameters ();
  void displayParameters (FILE * output);
	void checkSensitivity();

  void setMinLength (unsigned int minLength);
  void setMinUnits (unsigned int minUnits);
  void setMinUnitLength (unsigned int minUnitLength);
  void setMaxUnitLength (unsigned int maxUnitLength);

  void setDefaultMinLength() {minLength_ = 8;}
  void setDefaultMinUnits() {minUnits_ = 3;}
  void setDefaultMinUnitLength() {minUnitLength_ = 1;}
  void setDefaultMaxUnitLength() {maxUnitLength_ = 4;}

  unsigned int getMinLength () const {return minLength_;}
  unsigned int getMinUnits () const {return minUnits_;}
  unsigned int getMinUnitLength () const {return minUnitLength_;}
  unsigned int getMaxUnitLength () const {return maxUnitLength_;}

  vector< Microsatellite > detectPerfectMicrosatellites (const string & seq, const string & tag, unsigned int flankLength);
};

#endif /* MicrosatelliteDetector.h */
