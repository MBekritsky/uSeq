#include "MsHeaderFile.h"

using namespace std;

MsHeaderFile::MsHeaderFile()
{
	currentHdr_ . flushHdr();
	hdrFileP_ = NULL;
}

MsHeaderFile::MsHeaderFile(FILE * fp)
{
	currentHdr_ . flushHdr();
	hdrFileP_ = fp;
}

MsHeaderFile::MsHeaderFile(const string & filename)
{
	currentHdr_ . flushHdr();
	hdrFileP_ = openAndTestFile(filename,"r");
}

MsHeaderFile::~MsHeaderFile()
{
	if(hdrFileP_)
		fclose(hdrFileP_);
}

bool MsHeaderFile::getNextHdr()
{
	char ch;
	string hdrLine;
	
	while((ch = fgetc(hdrFileP_)) != EOF && ch != '\n')
		hdrLine . push_back ( (char) ch);
	if(ch == EOF)
		return false;
	currentHdr_ . parseHdrLine(hdrLine);
	return true;
}

bool MsHeaderFile::getNextHdr(bool first)
{
	char ch;
	string hdrLine;
	
	while((ch = fgetc(hdrFileP_)) != EOF && ch != '\n')
		hdrLine . push_back ( (char) ch);
	if(ch == EOF)
		return false;
	currentHdr_ . parseHdrLine(hdrLine, first);
	return true;
}

bool MsHeaderFile::getNextHeader()
{
	return getNextHdr();
}
