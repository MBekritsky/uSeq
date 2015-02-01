#ifndef __OSI_FILE_H__
#define __OSI_FILE_H__

#include <vector>
#include <string>
#include <stdlib.h>
#include "FileOps.h"

using namespace std;

class CoordOffsets
{
	private:
		string coordName_;
		vector< unsigned int > positions_, offsets_;
	
		void string2vec(vector< unsigned int > & vector, const string & string);
	public:
		CoordOffsets();
		CoordOffsets(string coordName, const vector< unsigned int > & positions, const vector< unsigned int > & offsets);
		CoordOffsets(string coordName, const string & unparsedPositions, const string & unparsedOffsets);
		CoordOffsets(const CoordOffsets & rhs);
		
		CoordOffsets& operator=(const CoordOffsets & rhs);
		
		void setCoordName(const string coordName) {coordName_ = coordName;}
		void setPositions(const vector< unsigned int > & positions) {positions_ = positions;}
		void setOffsets(const vector< unsigned int > & offsets) {offsets_ = offsets;}
		
		void checkVecLengths();

		string getCoordName() const {return coordName_;}
		void getPositionVector(vector< unsigned int > & positions) const {positions = positions_;}
		void getOffsetVector(vector< unsigned int > & offsets) const {offsets = offsets_;}
		unsigned int getPositionVectorSize() const {return positions_ . size();}
		//returns the proper index in the position array.  You have the option of specifying a place to start from in the array
		int getPositionVectorIndex(unsigned int position) const;
		int getPositionVectorIndex(unsigned int position, unsigned int start) const;
		int getPositionAt(unsigned int index) const;
		int getOffsetAt(unsigned int index) const;
		
		void printPositions() const;
};

class OsiFile
{
	private:
		FILE * osiFileP_;
		vector< CoordOffsets > coordOffsets_;

	public:
		OsiFile();
		OsiFile(FILE * osiFileP);
		OsiFile(const string & osiFileName);
		~OsiFile();
		
		bool loadOffsets();
		int getCoordIndex(string coord) const;
		int getIndex(string coord, unsigned int position) const;
		int getIndex(int coordIndex, unsigned int position) const;
		int getIndex(string coord, unsigned int position, unsigned int start) const;
		int getIndex(int coordIndex, unsigned int position, unsigned int start) const;
		int getPosition(string coord, unsigned int index) const;
		int getPosition(int coordIndex, unsigned int index) const;
		int getOffset(string coord, unsigned int index) const;
		int getOffset(int coordIndex, unsigned int index) const;
		unsigned int getVectorSize(string coord) const;
		unsigned int getVectorSize(int coordIndex) const;
		void printPositionVectors() const;
};


#endif /* OsiFile.h */
