#include "MsHeader.h"

using namespace std;

MsHeader::MsHeader()
{
	flushHdr();
}

MsHeader::MsHeader(const string & qname, bool isFirst, const vector< string > & microsatellites, const vector< string > & microsatelliteQuals)
{
	qname_ = qname;
	isFirst_ = isFirst;
	setDefaultIRD();
	setDefaultIFD();
	unparsedMicrosatellites_ = microsatellites;
	microsatelliteQuals_ = microsatelliteQuals;
}

MsHeader::MsHeader(const string & hdrLine)
{
	setDefaultIRD();
	setDefaultIFD();
	parseHdrLine(hdrLine);
}

MsHeader::MsHeader(const MsHeader & hdr)
{
	IRD_ = hdr . getIRD();
	IFD_ = hdr . getIFD();
	hdr . getQName(qname_);
	isFirst_ = hdr . getIsFirst();
	hdr . getMicrosatellites(unparsedMicrosatellites_);
	hdr . getMicrosatelliteQuals(microsatelliteQuals_);
}

MsHeader & MsHeader::operator=(const MsHeader & rhs)
{
	if(this != &rhs)
		{
			IRD_ = rhs . getIRD();
			IFD_ = rhs . getIFD();
			rhs . getQName(qname_);
			isFirst_ = rhs . getIsFirst();
			rhs . getMicrosatellites(unparsedMicrosatellites_);
			rhs . getMicrosatelliteQuals(microsatelliteQuals_);
		}
	return *this;
}

void MsHeader::flushHdr()
{
	setDefaultIRD();
	setDefaultIFD();
	qname_ . clear();
	isFirst_ = false;
	unparsedMicrosatellites_ . clear();
	microsatelliteQuals_ . clear();
}

void MsHeader::parseHdrLine(const string & hdrLine)
{
	flushHdr();
	
	size_t strLength = 0;
	
	strLength = hdrLine . find_first_of("/");
	qname_ = hdrLine . substr(0, strLength);
	
	strLength++;
	isFirst_ = determinePair(hdrLine[strLength]);
	strLength++;

	if(strLength < hdrLine . size())
		{
			size_t recStart = strLength + 1;
			size_t qualStart = 0;
			size_t nextIFD = 0;
			
			while(strLength != string::npos)
				{
					strLength = hdrLine . find_first_of(IRD_, recStart);
			//Find the quality line this way since we don't know if there are ':'s in the quality line
					nextIFD = recStart;
					for(unsigned int i = 0; i < 3; i++)
						{
							nextIFD = hdrLine . find_first_of(IFD_, nextIFD);
							nextIFD++;
						}
					qualStart = nextIFD - 1;
					
					unparsedMicrosatellites_ . push_back (hdrLine . substr(recStart, qualStart - recStart));
					microsatelliteQuals_ . push_back(hdrLine . substr(qualStart + 1, strLength - qualStart - 1));
					recStart = strLength + 1;
				}
		}
}

void MsHeader::parseHdrLine(const string & hdrLine, bool first)
{
	flushHdr();
	
	size_t strLength = 0;
	
	strLength = hdrLine . find_first_of("/");
	qname_ = hdrLine . substr(0, strLength);
		
	isFirst_ = first;

	if(strLength < hdrLine . size() && strLength != string::npos)
		{
			size_t recStart = strLength + 1;
			size_t qualStart = 0;
			
			while(strLength != string::npos)
				{
					strLength = hdrLine . find_first_of(IRD_, recStart);
					qualStart = hdrLine . find_last_of(IFD_, strLength);
					
					unparsedMicrosatellites_ . push_back (hdrLine . substr(recStart, qualStart - recStart));
					microsatelliteQuals_ . push_back(hdrLine . substr(qualStart + 1, strLength - qualStart - 1));
					recStart = strLength + 1;
				}
		}
}

bool MsHeader::determinePair(char num)
{
	switch (num)
		{
			case '1':
				return true;
			case '2':
				return false;
			default:
				cerr << num << " is not a correct mate pair!" << endl;
				exit(EXIT_FAILURE);
		}
	return false;
}

bool MsHeader::compareQNames(const string & newQName)
{
	if(qname_ . compare(newQName) == 0)
		return true;
	return false;
}

void MsHeader::printHdr(FILE * output)
{
	unsigned int unparsedSize = unparsedMicrosatellites_ . size();
	
	fprintf(output,"%s",qname_ . c_str());
	if(unparsedSize > 0)
		{
			fprintf(output,":");
			for(unsigned int i = 0; i < (unparsedSize - 1); i++)
				{
					fprintf(output,"%s:%s",unparsedMicrosatellites_ [i] . c_str(),microsatelliteQuals_ [i] . c_str());
				}
			fprintf(output,"%s:%s",unparsedMicrosatellites_ [unparsedSize - 1] . c_str(),microsatelliteQuals_ [unparsedSize - 1] . c_str());
		}
	fprintf(output,"\n");
}

void MsHeader::printHdr()
{
	printHdr(stdout);
}
