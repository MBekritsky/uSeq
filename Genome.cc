#include "Genome.h"

using namespace std;

Genome::Genome()
{
	chromosomes_ . clear();
	name_ . clear();
	length_ = 0;
	defaultNtPerLine();
}

Genome::Genome(const string & fastaName)
{
	FastaFile fastaFile(fastaName);
	name_ = getBaseName(fastaName);
	loadFromFasta(fastaFile);
	defaultNtPerLine();
}

Genome::Genome(const string &fastaName, const string & genomeName)
{
	FastaFile fastaFile(fastaName);
	name_ = genomeName;
	loadFromFasta(fastaFile);
	defaultNtPerLine();
}

Genome::Genome(const string & fastaName, unsigned int ntPerLine)
{
	FastaFile fastaFile(fastaName);
	name_ = getBaseName(fastaName);
	loadFromFasta(fastaFile);
	ntPerLine_ = ntPerLine;
}

Genome::Genome(const string & fastaName, unsigned int ntPerLine, const string & genomeName)
{
	FastaFile fastaFile(fastaName);
	name_ = genomeName;
	loadFromFasta(fastaFile);
	ntPerLine_ = ntPerLine;
}

Genome::Genome(FastaFile & fastaFile, const string & genomeName)
{
	name_ = genomeName;
	loadFromFasta(fastaFile);
	defaultNtPerLine();
}

Genome::Genome(FastaFile & fastaFile, unsigned int ntPerLine, const string & genomeName)
{
	name_ = genomeName;
	loadFromFasta(fastaFile);
	ntPerLine_ = ntPerLine;
}

Genome::Genome(const Genome & genome)
{
	name_ = genome . name();
	genome . chromosomes (chromosomes_);
	length_ = genome . length();
	nChr_ = genome . nChr();
	ntPerLine_ = genome . ntPerLine();
}

Genome & Genome::operator=(const Genome & rhs)
{
	if(this != &rhs)
		{
			name_ = rhs . name();
			length_ = rhs . length();
			nChr_ = rhs  . nChr();
			ntPerLine_ = rhs . ntPerLine();
			rhs . chromosomes (chromosomes_);
		}
	return *this;
}

void Genome::loadFromFasta(FastaFile & fastaFile)
{
	Chromosome temp;
	string hdr, seq;
	while(fastaFile . getNextSequence())
		{
			fastaFile . getHdr(hdr);
			fastaFile . getSeq(seq);
			temp . set(hdr, seq);
			chromosomes_ . push_back(temp);
		}
	hdr . clear();
	seq . clear();
	tabLength();
	tabNChr();
}

void Genome::tabLength()
{
	for(unsigned int i = 0; i < chromosomes_ . size(); i++)
		length_ += chromosomes_ [i] . length();
}

unsigned int Genome::chrIndex(const string & chrName) const
{
	for(unsigned int i = 0; i < chromosomes_ . size(); i++)
		{
			if( chrName . compare (chromosomes_ [i] . name()) == 0)
				return i;
		}
	cerr << "Could not find " << chrName << " in " << name_ << endl;
	return 0;
}

string Genome::seq(const string & chrName, unsigned int start, unsigned int stop) const
{
	return chromosomes_ [chrIndex(chrName)] . seq(start, stop);
}

string Genome::seq(unsigned int chrIndex, unsigned int start, unsigned int stop) const
{
	return chromosomes_ [chrIndex] . seq (start, stop);
}

void Genome::seq(const string & chrName, unsigned int start, unsigned int stop, string & seq) const
{
	chromosomes_ [chrIndex(chrName)] . seq (start, stop, seq);
}

void Genome::seq(unsigned int chrIndex, unsigned int start, unsigned int stop, string & seq) const
{
	chromosomes_ [chrIndex] . seq (start, stop, seq);
}

char Genome::seq(const string & chrName, unsigned int pos) const
{
	return chromosomes_ [chrIndex(chrName)] . seq (pos);
}

char Genome::seq(unsigned int chrIndex, unsigned int pos) const
{
	return chromosomes_ [chrIndex] . seq (pos);
}

void Genome::insert(const string & chrName, unsigned int start, const string & seq)
{
	chromosomes_ [chrIndex(chrName)] . insert(start, seq);
}

void Genome::insert(unsigned int chrIndex, unsigned int start, const string & seq)
{
	chromosomes_ [chrIndex] . insert (start, seq);
}

void Genome::remove (const string & chrName, unsigned int start, unsigned int length)
{
	chromosomes_ [chrIndex(chrName)] . remove(start, length);
}

void Genome::remove(unsigned int chrIndex, unsigned int start, unsigned int length)
{
	chromosomes_ [chrIndex] . remove(start, length);
}

void Genome::substitute (const string & chrName, unsigned int position, char sub)
{
	chromosomes_ [chrIndex(chrName)] . substitute(position, sub);
}

void Genome::substitute (unsigned int chrIndex, unsigned int position, char sub)
{
	chromosomes_ [chrIndex] . substitute(position, sub);
}

void Genome::print(FILE * fp) const
{
	for(unsigned int i = 0; i < chromosomes_ . size(); i++)
		chromosomes_ [i] . print(fp, ntPerLine_);
}

void Genome::printChrNames() const
{
	for(unsigned int i = 0; i < chromosomes_ . size() - 1; i++)
		cout << chromosomes_ [i] . name() << ", ";
	cout << "and " << chromosomes_ . back() . name() << endl;
}

void Genome::printChrInfo() const
{
	for(unsigned int i = 0; i < chromosomes_ . size(); i++)
		cout << chromosomes_ [i] . name() << ": " << chromosomes_[i] . length() << " nt" << endl;
}

Chromosome::Chromosome()
{
	name_ . clear();
	sequence_ . clear();
	length_ = 0;
}

Chromosome::Chromosome(FastaFile & fastaFile)
{
	fastaFile . getNextSequence();
	fastaFile . getHdr(name_);
	fastaFile . getSeq(sequence_);
	length_ = sequence_ . length();
}

Chromosome::Chromosome(const string & fastaFileName)
{
	FastaFile fastaFile(fastaFileName);
	fastaFile . getNextSequence();
	fastaFile . getHdr(name_);
	fastaFile . getSeq(sequence_);
	length_ = sequence_ . length();
}

Chromosome::Chromosome(const string & name, const string & sequence)
{
	name_ = name;
	sequence_ = sequence;
	length_ = sequence_ . length();
}

Chromosome::Chromosome(const Chromosome & chromosome)
{
	name_ = chromosome . name();
	chromosome . fullSeq(sequence_);
	length_ = sequence_ . length();
}

Chromosome & Chromosome::operator=(const Chromosome & rhs)
{
	if(this != &rhs)
		{
			name_ = rhs . name();
			rhs . fullSeq (sequence_);
			length_ = sequence_ . length();
		}
	return *this;
}

void Chromosome::set(const string & name, const string & sequence)
{
	name_ = name;
	sequence_ = sequence;
}

string Chromosome::seq(unsigned int start, unsigned int stop) const
{
	return sequence_ . substr(start - 1, stop - start + 1);
}

void Chromosome::seq(unsigned int start, unsigned int stop, string & sequence) const
{
	sequence = sequence_ . substr(start - 1, stop - start + 1);
}

char Chromosome::seq(unsigned int position) const
{
	return sequence_ [position];
}

void Chromosome::insert(unsigned int start, const string & insertion)
{
	sequence_ . insert(start - 1, insertion);
	length_ += insertion . length();
}

void Chromosome::remove(unsigned int start, unsigned int length)
{
	sequence_ . erase(start, length);
	length_ -= length;
}

void Chromosome::substitute(unsigned int position, char sub)
{
	sequence_ [position] = sub;
}

void Chromosome::print(FILE * fp, unsigned int ntPerLine) const
{
	unsigned int lines = length_ / ntPerLine;
		
	fprintf(fp,">%s\n",name_ . c_str());
	
	for(unsigned int i = 0; i < lines; i++)
			fprintf(fp,"%s\n", sequence_ . substr((i * ntPerLine), ntPerLine) . c_str());

	if((length_ % ntPerLine) > 0)
		fprintf(fp,"%s\n", sequence_ . substr(lines * ntPerLine) . c_str());
}
