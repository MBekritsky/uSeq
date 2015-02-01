#ifndef __MS_HEADER_H__
#define __MS_HEADER_H__

#include <string>
#include <vector>
#include <iostream>

using namespace std;

class MsHeader
{
	private:
		string qname_;
		vector< string > unparsedMicrosatellites_;
		vector< string > microsatelliteQuals_;
		bool isFirst_;
		char IRD_;
		char IFD_;
				
	public:
		MsHeader();
		MsHeader(const string & qname, bool isFirst, const vector< string > & unparsedMicrosatellites, const vector< string > & microsatelliteQuals);
		MsHeader(const string & hdrLine);
		MsHeader(const MsHeader & hdr);
		
		MsHeader & operator=(const MsHeader & rhs);
		
		void flushHdr();
		
		void setQName(const string & qname) {qname_ = qname;}
		void setIsFirst(bool isFirst) {isFirst_ = isFirst;}
		void setMicrosatellites(const vector< string > & microsatellites) {unparsedMicrosatellites_ = microsatellites;}
		void setMicrosatelliteQuals(const vector< string > & microsatelliteQuals) {microsatelliteQuals_ = microsatelliteQuals;}
		void setIRD(char IRD) {IRD_ = IRD;}
		void setIFD(char IFD) {IFD_ = IFD;}
		
		void setDefaultIRD() {IRD_ = '|';}
		void setDefaultIFD() {IFD_ = ':';}
		
		string getQName() const {return qname_;}
		void getQName(string & qname) const {qname = qname_;}
		bool getIsFirst() const {return isFirst_;}
		vector< string > getMicrosatellites() const {return unparsedMicrosatellites_;}
		void getMicrosatellites(vector< string > & microsatellites) const {microsatellites = unparsedMicrosatellites_;}
		vector< string > getMicrosatelliteQuals() const {return microsatelliteQuals_;}
		void getMicrosatelliteQuals(vector< string > & microsatelliteQuals) const {microsatelliteQuals = microsatelliteQuals_;}
		bool isFirst() const {return isFirst_;}
		char getIRD() const {return IRD_;}
		char getIFD() const {return IFD_;}
		
		void parseHdrLine(const string & hdrLine);
		void parseHdrLine(const string & hdrLine, bool first);
		bool determinePair(char num);
		bool compareQNames(const string & newQName);
		void printHdr(FILE * output);
		void printHdr();
};

#endif /* MsHeader.h */
