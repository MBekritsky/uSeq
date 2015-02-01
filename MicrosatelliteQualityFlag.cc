#include "MicrosatelliteQualityFlag.h"

using namespace std;

MicrosatelliteQualityFlag::MicrosatelliteQualityFlag()
{
	readStart_ 	 = 0;
	readStop_	 	 = 0;
	length_		 	 = 0;
	motifLength_ = 0;
	
	motif_ = "";
	qualityString_ = "";
}

MicrosatelliteQualityFlag::MicrosatelliteQualityFlag(const string & microsatelliteQualityStr)
{
	vector< string > splitQualityString;
	splitString(microsatelliteQualityStr,":",splitQualityString,3);
	
	motif_ = splitQualityString[0];
	readStart_ = atoi(splitQualityString[1].c_str());
	qualityString_ = splitQualityString[2];
	
	motifLength_ = motif_.length();
	length_ = qualityString_.length();
	readStop_ = readStart_ + length_ - 1;
}

MicrosatelliteQualityFlag::MicrosatelliteQualityFlag(const MicrosatelliteQualityFlag & rhs)
{
	readStart_   = rhs.readStart();
	readStop_	   = rhs.readStop();
	motifLength_ = rhs.motifLength();
	length_			 = rhs.length();
	rhs.getMotif(motif_);
	rhs.getQualityString(qualityString_);
}

MicrosatelliteQualityFlag & MicrosatelliteQualityFlag::operator=(const MicrosatelliteQualityFlag & rhs)
{
	if(this != &rhs)
	{
		readStart_   = rhs.readStart();
		readStop_	   = rhs.readStop();
		motifLength_ = rhs.motifLength();
		length_			 = rhs.length();
		rhs.getMotif(motif_);
		rhs.getQualityString(qualityString_);
	}
	return *this;
}
