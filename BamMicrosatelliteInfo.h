#ifndef __BAM_MICROSATELLITE_INFO__
#define __BAM_MICROSATELLITE_INFO__

#include <string>
#include <vector>
#include <iostream>
#include "Microsatellite.h"
#include "StringManip.h"

#include "api/BamAlignment.h"

using namespace std;
using namespace BamTools;

class BamMicrosatelliteInfo
{
	private:
		unsigned int nMicrosatellites_;
		vector< Microsatellite > microsatellites_;
		vector< string > microsatelliteQuals_;
		string plainQName_;
		unsigned int readLength_;
		char OFD_,IFD_,ORD_,IRD_;
		
	public:
		BamMicrosatelliteInfo();
		/*BamMicrosatelliteInfo(const vector< string > & microsatellites, 
													const vector< string > & microsatelliteQuals);*/
		BamMicrosatelliteInfo(const BamMicrosatelliteInfo & rhs);
		BamMicrosatelliteInfo & operator=(const BamMicrosatelliteInfo & rhs);
		
		void flush();
		void addIndexingTagsToRead(BamAlignment & al, const string & rname, int offset) const;
		void addMUTag(BamAlignment & al) const;
		void addMCTag(BamAlignment & al, const string & rname, int offset) const;
		void addOPTag(BamAlignment & al, int offset) const;
		void readLength(const BamAlignment & al);
		void updateName(BamAlignment & al) const;
		
		void addMicrosatellites(const vector< string > & unparsedMicrosatellites);
		void addMicrosatellites(const vector< Microsatellite > & parsedMicrosatellites);
		void addMicrosatelliteQuals(const vector< string > & microsatelliteQuals) {microsatelliteQuals_ = microsatelliteQuals;}
		void plainQName(const string & plainQName) {plainQName_ = plainQName;}
		void plainQName(const BamAlignment & al) {plainQName_ = al.Name;}
		void countMicrosatellites();
		void reverseMicrosatellites();
		unsigned int adjustMicrosatelliteStart(int offset, const Microsatellite & microsatellite) const;
		unsigned int softClipAdjustment(const BamAlignment & al) const;

		void OFD (char ch) {OFD_ = ch;}
		void IFD (char ch) {IFD_ = ch;}
		void ORD (char ch) {ORD_ = ch;}
		void IRD (char ch) {IRD_ = ch;}
		
		unsigned int nMicrosatellites() const {return nMicrosatellites_;}
		void getMicrosatellites(vector< Microsatellite > & microsatellites) const {microsatellites = microsatellites_;}
		void getMicrosatelliteQuals(vector< string > & microsatelliteQuals) const {microsatelliteQuals = microsatelliteQuals_;}
		void getPlainQName(string & plainQName) const {plainQName = plainQName_;}
		char OFD() const {return OFD_;}
		char IFD() const {return IFD_;}
		char ORD() const {return ORD_;}
		char IRD() const {return IRD_;}
		unsigned int readLength() const {return readLength_;}
		
		void extractMicrosatellites(const BamAlignment & al);
		void checkMicrosatelliteQuals(const BamAlignment & al) const;
		char getPartnerNucleotide(char nt) const;
		void reverseComplement(const string & unit, string & rcUnit) const;
		unsigned int adjustMicrosatelliteStart(int offset, const Microsatellite & microsatellite, const BamAlignment & al) const;		
		void defaultOFD() {OFD_ = ':';}
		void defaultIFD() {IFD_ = ':';}
		void defaultORD() {ORD_ = '|';}
		void defaultIRD() {IRD_ = '|';}
		void defaultDelimiters();

};

#endif /* BamMicrosatelliteInfo.h */
