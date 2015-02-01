#ifndef __MICROSATELLITE_COORD_FLAG_H__
#define __MICROSATELLITE_COORD_FLAG_H__

#include <string>
#include <sstream>
#include "StringManip.h"

using namespace std;

class MicrosatelliteCoordFlag
{
	private:
		unsigned int genomeStart_, genomeStop_, length_, motifLength_;
		string motif_, coord_/*, name_*/;
	public:
		MicrosatelliteCoordFlag();
		MicrosatelliteCoordFlag(const string & MicrosatelliteCoordStr);
		MicrosatelliteCoordFlag(const MicrosatelliteCoordFlag & rhs);
		
		MicrosatelliteCoordFlag & operator=(const MicrosatelliteCoordFlag & rhs);
		
		void genomeStart(unsigned int genomeStart) 	 {genomeStart_ = genomeStart;}
		void genomeStop(unsigned int genomeStop) 	 	 {genomeStop_ = genomeStop;}
		void motifLength(unsigned int motifLength) {motifLength_ = motifLength;}
		void length(unsigned int length) 		 		 	 {length_ = length;}
		void motif(const string & motif) 		 		 	 {motif_ = motif;}
		void coord(const string & coord) 					 {coord_ = coord;}
		
		unsigned int genomeStart() const {return genomeStart_;}
		unsigned int genomeStop()  const {return genomeStop_;}
		unsigned int motifLength() const {return motifLength_;}
		unsigned int length() 		 const {return length_;}
		
		string motif() const {return motif_;}
		string coord() const {return coord_;}
		
		void getMotif(string & motif) const {motif = motif_;}
		void getCoord(string & coord) const {coord = coord_;}
};

#endif /* MicrosatelliteQualityFlag.h */
