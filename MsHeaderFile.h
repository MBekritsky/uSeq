#ifndef __MS_HEADER_FILE_H__
#define __MS_HEADER_FILE_H__

#include <string>
#include "FileOps.h"
#include "MsHeader.h"

using namespace std;

class MsHeaderFile
{
	private:
		MsHeader currentHdr_;
		FILE * hdrFileP_;
		
	public:
		MsHeaderFile();
		MsHeaderFile(FILE * fp);
		MsHeaderFile(const string & filename);
		~MsHeaderFile();
		
		void getHdr (MsHeader & hdr) const {hdr = currentHdr_;}
		MsHeader getHdr () const {return currentHdr_;}
		void getHeader (MsHeader & hdr) const {hdr = currentHdr_;}
		MsHeader getHeader () const {return currentHdr_;}
		
		bool getNextHeader();
		bool getNextHeader(bool first);
		bool getNextHdr();
		bool getNextHdr(bool first);
};

#endif /* MsHeaderFile.h*/
