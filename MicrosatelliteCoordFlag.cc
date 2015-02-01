#include "MicrosatelliteCoordFlag.h"

using namespace std;

MicrosatelliteCoordFlag::MicrosatelliteCoordFlag()
{
	genomeStart_ 	 = 0;
	genomeStop_	 	 = 0;
	length_		 	 = 0;
	motifLength_ = 0;
	
	motif_ = "";
	coord_ = "";
//	name_ = "";
}

MicrosatelliteCoordFlag::MicrosatelliteCoordFlag(const string & microsatelliteCoordStr)
{
	vector< string > splitCoordString;
	splitString(microsatelliteCoordStr,":",splitCoordString,4);
	
	coord_ = splitCoordString[0];
	motif_ = splitCoordString[1];
	genomeStart_ = atoi(splitCoordString[2].c_str());
	length_ = atoi(splitCoordString[3].c_str());
	
	motifLength_ = motif_.length();
	genomeStop_ = genomeStart_ + length_ - 1;

/*	ostringstream nameSStream;
	nameSStream << coord_ << motif_ << genomeStart_;
	name_ = nameSStream.str();
	*/
}

MicrosatelliteCoordFlag::MicrosatelliteCoordFlag(const MicrosatelliteCoordFlag & rhs)
{
	genomeStart_   = rhs.genomeStart();
	genomeStop_	   = rhs.genomeStop();
	motifLength_ = rhs.motifLength();
	length_			 = rhs.length();
	rhs.getMotif(motif_);
	rhs.getCoord(coord_);
//	rhs.name(name_);
}

MicrosatelliteCoordFlag & MicrosatelliteCoordFlag::operator=(const MicrosatelliteCoordFlag & rhs)
{
	if(this != &rhs)
	{
		genomeStart_   = rhs.genomeStart();
		genomeStop_	   = rhs.genomeStop();
		motifLength_ = rhs.motifLength();
		length_			 = rhs.length();
		rhs.getMotif(motif_);
		rhs.getCoord(coord_);
//		rhs.name(name_);
	}
	return *this;
}
