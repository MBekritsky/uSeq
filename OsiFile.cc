#include "OsiFile.h"

using namespace std;

CoordOffsets::CoordOffsets()
{
	coordName_ . clear();
	positions_  . clear();
	offsets_ . clear();
}

CoordOffsets::CoordOffsets(string coordName, const vector< unsigned int > & positions, const vector< unsigned int > & offsets)
{
	setCoordName(coordName);
	setPositions(positions);
	setOffsets(offsets);
}

CoordOffsets::CoordOffsets(string coordName, const string & unparsedPositions, const string & unparsedOffsets)
{
	setCoordName(coordName);
	string2vec(positions_, unparsedPositions);
	string2vec(offsets_, unparsedOffsets);
}

CoordOffsets::CoordOffsets(const CoordOffsets & rhs)
{
	coordName_ = rhs . getCoordName();
	rhs . getPositionVector(positions_);
	rhs . getOffsetVector(offsets_);
}

CoordOffsets & CoordOffsets::operator=(const CoordOffsets & rhs)
{
	if(this != &rhs)
		{
			coordName_ = rhs . getCoordName();
			rhs . getPositionVector(positions_);
			rhs . getOffsetVector(offsets_);
		}
	return *this;
}

void CoordOffsets::string2vec(vector< unsigned int > & vector, const string & string)
{
	char delimiter = '\t';
	
	for(unsigned int i = 0; i < string . length(); i++)
		{
			unsigned int j = i;
			for(; string [j] != delimiter && string [j] != '\n'; j++)
				;
			vector . push_back (atoi(string . substr (i, j - i) . c_str()));
			i = j;
		}
}

void CoordOffsets::checkVecLengths()
{
	if(positions_ . size() != offsets_ . size())
		{
			cerr << "Error!!! The lengths of the offset and position vectors for " << coordName_ << " do not match!";
			exit(EXIT_FAILURE);
		}
}

int CoordOffsets::getPositionVectorIndex(unsigned int position) const
{
	return getPositionVectorIndex(position, 0);
}

int CoordOffsets::getPositionVectorIndex(unsigned int position, unsigned int start) const
{
	if(position < positions_ [0])
		return -1;
	for(unsigned int i = start; i < positions_ . size(); i++)
		{
			if(position > positions_[i] && position <= positions_ [i + 1])
				return i;
		}
	//-2 means that the position exceeds any position in the position vector (i.e., 
	//the position is after the last microsatellite
	return -2;
}

int CoordOffsets::getPositionAt(unsigned int index) const
{
	if(index > positions_ . size())
		return -1;
	return positions_ [index];
}

int CoordOffsets::getOffsetAt(unsigned int index) const
{
	if(index > offsets_ . size())
		return -1;
	return offsets_ [index];
}

OsiFile::OsiFile()
{
	osiFileP_ = NULL;
	coordOffsets_ . clear();
}

OsiFile::OsiFile(FILE * osiFileP)
{
	osiFileP_ = osiFileP;
	loadOffsets();
}

OsiFile::OsiFile(const string & osiFileName)
{
	osiFileP_ = openAndTestFile(osiFileName,"r");
	loadOffsets();
}

OsiFile::~OsiFile()
{
	fclose(osiFileP_);
}

bool OsiFile::loadOffsets()
{
	char ch = '0';
	string tCoord, tPosLine, tOffsetLine;
	while(ch != EOF)
		{
			while((ch = fgetc(osiFileP_)) != EOF && ch != '\n')
				tCoord . push_back(ch);
			while((ch = fgetc(osiFileP_)) != '\n')
				tPosLine . push_back(ch);
			while((ch = fgetc(osiFileP_)) != '\n' && ch != EOF)
				tOffsetLine  . push_back(ch);

			ch = fgetc(osiFileP_);
			if(ch != EOF)
				ungetc(ch, osiFileP_);
			
			CoordOffsets coordOffset(tCoord, tPosLine, tOffsetLine);
			coordOffsets_ . push_back (coordOffset);
			//check that the most recently added position and offset vectors have equal lengths
			coordOffsets_ . back() . checkVecLengths();
			tCoord . clear();
			tPosLine . clear();
			tOffsetLine . clear();
		}

	if(coordOffsets_ . size() == 0)
		return false;
	return true;
}

int OsiFile::getCoordIndex(string coordName) const
{
	for(unsigned int i = 0; i < coordOffsets_ . size(); i++)
		{
			if(coordName . compare(coordOffsets_ [i] . getCoordName()) == 0)
				return i;
		}
	return -2;
}

int OsiFile::getIndex(string coordName, unsigned int position) const
{
	int coordIndex = getCoordIndex(coordName);
	// -2 is code for non-existent coordinate, -1 is code for position preceding first MS
	if(coordIndex == -2)
		return coordIndex;
	return coordOffsets_ [coordIndex] . getPositionVectorIndex(position);
}

int OsiFile::getIndex(int coordIndex, unsigned int position) const
{
	if(coordIndex >= (int) coordOffsets_ . size())
		return -2;
	return coordOffsets_ [coordIndex] . getPositionVectorIndex(position);
}

int OsiFile::getIndex(string coordName, unsigned int position, unsigned int start) const
{
	int coordIndex = getCoordIndex(coordName);
	if(coordIndex == -2)
		return coordIndex;
	return coordOffsets_ [coordIndex] . getPositionVectorIndex(position, start);
}

int OsiFile::getIndex(int  coordIndex, unsigned int position, unsigned int start) const
{
	if(coordIndex >= (int) coordOffsets_ . size())
		return -2;
	return coordOffsets_ [coordIndex] . getPositionVectorIndex(position, start);
}

int OsiFile::getPosition(string coordName, unsigned int index) const
{
	int coordIndex = getCoordIndex(coordName);
	if(coordIndex == -2)
		return coordIndex;
	return coordOffsets_ [coordIndex] . getPositionAt (index);
}

int OsiFile::getPosition(int coordIndex, unsigned int index) const
{
	if(coordIndex >= (int) coordOffsets_ . size())
		return -2;
	return coordOffsets_ [coordIndex] . getPositionAt(index);
}

int OsiFile::getOffset(string coordName, unsigned int index) const
{
	int coordIndex = getCoordIndex(coordName);
	if(coordIndex == -2)
		return coordIndex;
	return coordOffsets_ [coordIndex] . getOffsetAt (index);
}

int OsiFile::getOffset(int coordIndex, unsigned int index) const
{
	if(coordIndex >= (int) coordOffsets_ . size())
		return -2;
	return coordOffsets_ [coordIndex] . getOffsetAt (index);
}

unsigned int OsiFile::getVectorSize(string coord) const
{
	int coordIndex = getCoordIndex(coord);
	if(coordIndex == -2)
		return 0;
	return coordOffsets_ [coordIndex] . getPositionVectorSize();
}

unsigned int OsiFile::getVectorSize(int coordIndex) const
{
	if(coordIndex > (int) coordOffsets_ . size())
		return 0;
	return coordOffsets_ [coordIndex] . getPositionVectorSize();
}

void OsiFile::printPositionVectors() const
{
	for(unsigned int i = 0; i < coordOffsets_ . size(); i++)
			coordOffsets_ [i] . printPositions();
}

void CoordOffsets::printPositions() const
{
	cout << "position: ";
	for(unsigned int i =  0; i < positions_ . size(); i++)
		cout << positions_ [i] << ",";
	cout << endl;
	cout << "offsets: ";
	for(unsigned int i = 0; i < offsets_ . size(); i++)
		cout << offsets_ [i] << ",";
	cout << endl;
}
