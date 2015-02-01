/*
FastaFile.cc

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

#include "FastaFile.h"
#include <cstdio>
#include <cstdlib>

using namespace std;

FastaFile::FastaFile()
{
  seq_ . clear();
  hdr_ . clear();
	ntPerLine_ = 50;
  fastaFileP_ = NULL;
}

FastaFile::FastaFile(FILE * fp)
{
  fastaFileP_ = fp;
	ntPerLine_ = 0;
}

FastaFile::FastaFile(string filename)
{
  fastaFileP_ = openAndTestFile(filename, "r");
	ntPerLine_ = 0;
}

FastaFile::~FastaFile()
{
  fclose(fastaFileP_);
}

bool FastaFile::getNextSequence()
{
  char ch;
  
  hdr_ . erase();
  seq_ . erase();
  
  //skip to ">" line of FASTA FILE
  while ((ch = fgetc(fastaFileP_)) != EOF && ch != '>')
    ;
  if(ch == EOF)
    return false;

  //Skip any leading whitespace
  while ((ch = fgetc(fastaFileP_)) != EOF && ch == ' ')
    ;
  if(ch == EOF)
    return false;

  ungetc(ch, fastaFileP_);

  while ((ch = fgetc(fastaFileP_)) != EOF && ch != '\n')
    hdr_ . push_back (ch);

  while((ch = fgetc(fastaFileP_)) != EOF && ch != '>')
    {
      if(ch != '\n')
				seq_ . push_back ((char) toupper(ch));
			if(ch == '\n' && ntPerLine_ == 0)
				ntPerLine_ = seq_ . length();
    }
  
  if(ch == '>')
    ungetc (ch, fastaFileP_);

	setSeqLength();

  return true;
}

char FastaFile::getSeqAtPos(int pos)
{
  return seq_[pos];
}
