#ifndef __MICROSATELLITE_QUALITY_FLAG_H__
#define __MICROSATELLITE_QUALITY_FLAG_H__

#include <string>
#include "StringManip.h"

using namespace std;

class MicrosatelliteQualityFlag
{
	private:
		unsigned int readStart_, readStop_, motifLength_, length_;
		string motif_, qualityString_;
	public:
		MicrosatelliteQualityFlag();
		MicrosatelliteQualityFlag(const string & MicrosatelliteQualityStr);
		MicrosatelliteQualityFlag(const MicrosatelliteQualityFlag & rhs);
		
		MicrosatelliteQualityFlag & operator=(const MicrosatelliteQualityFlag & rhs);
		
		void readStart(unsigned int readStart) 	 		 {readStart_ = readStart;}
		void readStop(unsigned int readStop) 	 		 {readStop_ = readStop;}
		void motifLength(unsigned int motifLength) 		 {motifLength_ = motifLength;}
		void length(unsigned int length) 		 		 {length_ = length;}
		void motif(const string & motif) 		 		 {motif_ = motif;}
		void qualityString(const string & qualityString) {qualityString_ = qualityString;}
		
		unsigned int readStart() 	 const {return readStart_;}
		unsigned int readStop() 	 const {return readStop_;}
		unsigned int motifLength() const {return motifLength_;}
		unsigned int length() 		 const {return length_;}
		
		string motif() const {return motif_;}
		string qualityString() const {return qualityString_;}
		
		void getMotif(string & motif) const {motif = motif_;}
		void getQualityString(string & qualityString) const {qualityString = qualityString_;}
};

#endif /* MicrosatelliteQualityFlag.h */
