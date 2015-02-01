#include "BamMicrosatelliteInfo.h"

using namespace std;
using namespace BamTools;

BamMicrosatelliteInfo::BamMicrosatelliteInfo()
{
  flush();
	defaultDelimiters();
}

/*BamMicrosatelliteInfo::BamMicrosatelliteInfo (const vector< string > & microsatellites,
																							const vector< string > & microsatelliteQuals)
{
	addMicrosatellites(microsatellites);
	addMicrosatelliteQuals(microsatelliteQuals);
	countMicrosatellites();
	setDefaultDelimiters();
}*/

BamMicrosatelliteInfo::BamMicrosatelliteInfo(const BamMicrosatelliteInfo & rhs)
{
		rhs.getMicrosatellites(microsatellites_);
		rhs.getMicrosatelliteQuals(microsatelliteQuals_);
		rhs.getPlainQName(plainQName_);
		nMicrosatellites_ = rhs.nMicrosatellites();
		readLength_ = rhs.readLength();
		IFD_ = rhs.IFD();
		IRD_ = rhs.IRD();
		OFD_ = rhs.OFD();
		ORD_ = rhs.ORD();
}

BamMicrosatelliteInfo & BamMicrosatelliteInfo::operator=(const BamMicrosatelliteInfo & rhs)
{
	if(this != &rhs)
		{
			readLength_ = rhs.readLength();
			nMicrosatellites_ = rhs.nMicrosatellites();
			rhs . getMicrosatellites(microsatellites_);
			rhs . getMicrosatelliteQuals(microsatelliteQuals_);
			rhs.getPlainQName(plainQName_);
			IFD_ = rhs.IFD();
			IRD_ = rhs.IRD();
			OFD_ = rhs.OFD();
			ORD_ = rhs.ORD();
		}
		return *this;
}

void BamMicrosatelliteInfo::defaultDelimiters()
{
	defaultOFD();
	defaultORD();
	defaultIFD();
	defaultIRD();
}

void BamMicrosatelliteInfo::flush()
{
	readLength_ = 0;
	nMicrosatellites_ = 0;
	plainQName_ . clear();
	microsatellites_ . clear();
	microsatelliteQuals_ . clear();
}

void BamMicrosatelliteInfo::readLength(const BamAlignment & al)
{
	if(al.QueryBases.empty())
		cerr << "Warning, attempting to calculate read length for an empty read: " << al.Name << endl;
	
	readLength_ = al.QueryBases.length();
	
	if(nMicrosatellites_ > 0)
		{
			unsigned int lastStop = 0;
			
			for(unsigned int i = 0; i < nMicrosatellites_; i++)
				{
					if(microsatellites_ [i] . getStart() <= lastStop)
							readLength_ -= (lastStop - microsatellites_ [i] . getStart() + 1);
					readLength_ += microsatellites_ [i] . getLength();
					lastStop = microsatellites_ [i] . getStop();
				}
		}
}

void BamMicrosatelliteInfo::extractMicrosatellites(const BamAlignment & al)
{
	vector< string > unparsedMicrosatellites;
	string qname = al.Name;
	size_t strStart = 0, strLength = 0, qualCheckS = 0, qualCheckL = 0;
	unsigned int numDelim = 0;
	string qualTest;
	bool hasQual = false;
	
	//skip to first input record delimiter or end of qname
	strStart = qname.find_first_of(IRD_);
	//if there are microsatellites in qname, add them
	if(strStart != string::npos)
		{
			plainQName_ = qname . substr(0,strStart);
			strStart++;
			strLength = strStart;
			
			
			//Check to see if the microsatellite tags have quality scores attached to them or not
			qualCheckS = qname . find_first_of(IRD_, strStart);
			qualTest = qname . substr(strStart, qualCheckS - strStart);
			qualCheckS = 0;
			while(qualCheckL != string::npos)
				{
					qualCheckL = qualTest . find_first_of(IFD_, qualCheckS);
					if(qualCheckL != string::npos)
						{
							qualCheckL++;
							qualCheckS = qualCheckL;
							numDelim++;
						}
				}

			if(numDelim == 3)
				hasQual = true;

			splitString(qname.substr(strLength),IRD_,unparsedMicrosatellites);
			if(hasQual)
				{			
					vector< string > msDetails;
					for(unsigned int i = 0; i < unparsedMicrosatellites.size(); ++i)
					{
						splitString(unparsedMicrosatellites[i],IFD_,msDetails,4);
						//In case the field delimiter occurs
						//in the header string, split the string, limiting it to 4 total fields, the quality string,
						//even if it contains the IFD, will be at the 4th index of the array
						
						microsatelliteQuals_.push_back(msDetails[3]);
						unparsedMicrosatellites[i] = unparsedMicrosatellites[i].substr(0,
																						unparsedMicrosatellites[i].length() - (msDetails[3].length() + 1));
					}
				}
			addMicrosatellites(unparsedMicrosatellites);
		}
	else
		{
			plainQName_ = qname.substr(0);
		}
}

void BamMicrosatelliteInfo::checkMicrosatelliteQuals(const BamAlignment & al) const
{
	//if there are microsatellites in the read, either none of them have qualities
	//or all of them do
	if(nMicrosatellites_ > 0)
	{
		if(nMicrosatellites_  == microsatelliteQuals_ . size())
			{
				for(unsigned int i = 0; i < nMicrosatellites_; i++)
					{
						if(microsatelliteQuals_ [i] . length() != microsatellites_ [i] . getLength())
							{
								char buffer[1024];
								microsatellites_ [i] . printNoHdr (buffer);
								cerr << "Error!!!  Microsatellite " << buffer << " in " << al.Name << " has microsatellite quality that does not match microsatellite length\n";
								exit(EXIT_FAILURE);
							}
					}
			}
		else
			{
				cerr << cerr << "Error!!! " << al.Name << " has a discrepant number of microsatellites and microsatellite qualities" << endl;
				exit(EXIT_FAILURE);
			}
	}
	else if(nMicrosatellites_ == 0 && !microsatelliteQuals_.empty())
	{
		cerr << "Error! " << al.Name << " contains microsatellite qualities but no microsatellites" << endl;
		exit(EXIT_FAILURE);
	}
	
}

void BamMicrosatelliteInfo::addMicrosatellites(const vector < string > & unparsedMicrosatellites)
{
	for(unsigned int i = 0; i < unparsedMicrosatellites . size(); i++)
		microsatellites_ . push_back (Microsatellite (unparsedMicrosatellites [i], IFD_, OFD_));
	nMicrosatellites_ = unparsedMicrosatellites . size();
}

void BamMicrosatelliteInfo::addMicrosatellites(const vector< Microsatellite > & parsedMicrosatellites)
{
	microsatellites_ = parsedMicrosatellites;
	nMicrosatellites_ = parsedMicrosatellites . size();
}

char BamMicrosatelliteInfo::getPartnerNucleotide(char nt) const
{
	switch(nt)
		{
			case 'a':
				return 't';
				break;
			case 'A':
				return 'T';
				break;
			
			case 't':
				return 'a';
				break;
			case 'T':
				return 'A';
				break;
				
			case 'c':
				return 'g';
				break;
			case 'C':
				return 'G';
				break;
			
			case 'g':
				return 'c';
				break;
			case 'G':
				return 'C';
				break;
			
			default:
				return 'N';
				break;
		}
	return 'Q';
}

void BamMicrosatelliteInfo::reverseComplement(const string & unit, string & rcUnit) const
{
	for(int i = (unit . length() - 1); i >= 0; i--)
		{
			rcUnit . push_back(getPartnerNucleotide(unit[i]));
		}
}

void BamMicrosatelliteInfo::reverseMicrosatellites()
{
	unsigned int newStart, newStop;
	string newUnit, rcUnit;
	unsigned int unitOffset;
	for(unsigned int i = 0; i < nMicrosatellites_; i++)
		{
			newStart = readLength_ - microsatellites_ [i] . getStop() + 1;
			newStop = readLength_ - microsatellites_ [i] . getStart() + 1;
			unitOffset = microsatellites_ [i] . getLength() % (microsatellites_ [i] . getUnit()) . length();

			newUnit = (microsatellites_ [i] . getUnit()) . substr(unitOffset);
			if(unitOffset != 0)
				newUnit += (microsatellites_ [i] . getUnit()) . substr(0,unitOffset);
			reverseComplement(newUnit, rcUnit);

			microsatellites_[i] . setStart(newStart);
			microsatellites_[i] . setStop(newStop);
			microsatellites_[i] . setUnit(rcUnit);
			newUnit . clear();
			rcUnit . clear();
			
			//reverse quality scores as well
			if(microsatelliteQuals_ . size() > 0)
				reverse (microsatelliteQuals_ [i] . begin(), microsatelliteQuals_ [i] . end());
		}
}

void BamMicrosatelliteInfo::addIndexingTagsToRead(BamAlignment & al, const string & rname, int offset) const
{
	addMCTag(al,rname,offset);
	addOPTag(al,offset);
	addMUTag(al);
}

void BamMicrosatelliteInfo::addMCTag(BamAlignment & al, const string & rname, int offset) const
{
	if(nMicrosatellites_ > 0 && al.IsMapped())
		{
			static char temp[2056];

			sprintf(temp,"%s%c%s%c%d%c%d", rname . c_str(), OFD_, microsatellites_[0] . getUnit() . c_str(), OFD_, adjustMicrosatelliteStart(offset,	microsatellites_ [0],al), OFD_, microsatellites_[0] . getLength());

			if(nMicrosatellites_ > 1)
				{
					for (unsigned int i = 1; i < nMicrosatellites_; i++)
						{
							sprintf(temp,"%s%c%s%c%s%c%d%c%d", temp, ORD_, rname . c_str(), OFD_, microsatellites_[i] . getUnit() . c_str(), OFD_, adjustMicrosatelliteStart(offset, microsatellites_ [i],al), OFD_, microsatellites_[i] . getLength());
						}		
				}
			string coords (temp);
			al.AddTag("MC","Z",coords);
			coords . clear();
		}
}

//NB THIS MAY HAVE TO BE CHANGED TO ACCOUNT FOR SOFT CLIPPING ON EITHER SIDE OF THE READ
unsigned int BamMicrosatelliteInfo::adjustMicrosatelliteStart(int offset, const Microsatellite & microsatellite, const BamAlignment & al) const
{
	return al.Position + microsatellite.getStart() + offset - softClipAdjustment(al);
	//if read is mapped to forward strand and the first bases aren't clipped, then don't account for soft-clipping at beginning of read 
/* 	if(!al.IsReverseStrand() || (al.IsReverseStrand() && al.CigarData[0].Type != 'S'))
	{
		return al.Position + microsatellite . getStart() + offset;
	}
	else
	{
		unsigned int frontClip = al.CigarData[0].Length;
		return al.Position + microsatellite . getStart() + offset - frontClip;
	} */
}

unsigned int BamMicrosatelliteInfo::softClipAdjustment(const BamAlignment & al) const
{
	if(!al.IsReverseStrand() && (*al.CigarData.begin()).Type == 'S')
	{
		return (*al.CigarData.begin()).Length;
	}
	else if(al.IsReverseStrand() && (*al.CigarData.rbegin()).Type == 'S')
	{
		return (*al.CigarData.rbegin()).Length;
	}
	return 0;
}

void BamMicrosatelliteInfo::addOPTag(BamAlignment & al, int offset) const
{
	int originalPosition = al.Position + offset - softClipAdjustment(al) + 1;
	al.AddTag("OP","i",originalPosition);
}

void BamMicrosatelliteInfo::addMUTag(BamAlignment & al) const
{
	if(nMicrosatellites_ > 0)
		{
			static char temp[4112];
			
			sprintf(temp,"%s%c%d%c%s",(microsatellites_[0] . getUnit()) . c_str(), OFD_, microsatellites_[0] . getStart(), OFD_, microsatelliteQuals_ [0] . c_str());

			if(nMicrosatellites_ > 1)
				{
					for (unsigned int i = 1; i < microsatellites_ .size(); i++)
						{
							sprintf(temp,"%s%c%s%c%d%c%s",temp, ORD_, (microsatellites_[i] . getUnit()) . c_str(), OFD_, microsatellites_[i] . getStart(), OFD_, microsatelliteQuals_ [i] . c_str());
						}
				}
			string quals (temp);
			al.AddTag("MU","Z",quals);
			quals . clear();
		}

}

void BamMicrosatelliteInfo::updateName(BamAlignment & al) const
{
	string newName = plainQName_;
	if(nMicrosatellites_ > 0)
		{
			static char buffer[1024];
			for(unsigned int i = 0; i < nMicrosatellites_; ++i)
				{
					microsatellites_ [i] . printNoHdr(buffer);
					newName += ORD_;
					newName += buffer;
				}
		}
	al.Name = newName;
}
