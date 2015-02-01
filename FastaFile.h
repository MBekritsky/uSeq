/*
FastaFile.h

ABSTRACT:
A class and functions for loading FASTA format files
into a header and sequence string.

CREATION DATE:
20.08.2010

LAST REVISION:
26.08.2010

AUTHOR:
Mitchell Bekritsky
*/

#ifndef __LOAD_FASTA_H__
#define __LOAD_FASTA_H__

#include <string>
#include <iostream>
#include "FileOps.h"

using namespace std;

class FastaFile
{
 private:
  FILE * fastaFileP_;
  string seq_, hdr_;
	unsigned int ntPerLine_;
	unsigned int lSeq_;
 public:
  FastaFile();
  FastaFile(FILE * fastaFile);
  FastaFile(string filename);

  ~FastaFile();

  //get seq and hdr by address to speed up return of long sequences
  void getSeq(string & seq) const {seq = seq_;}
  void getHdr(string & hdr) const {hdr = hdr_;}
	unsigned int getSeqLength() {return lSeq_;}
	unsigned int getNtPerLine() {return ntPerLine_;}

  void setSeq(string ts) {seq_ = ts;}
  void setHdr(string th) {hdr_ = th;}
	void setSeqLength() {lSeq_ = seq_ . length();}

  char getSeqAtPos (int index);
  bool getNextSequence();
};

#endif /* LoadFasta.h */
