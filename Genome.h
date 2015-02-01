#ifndef __GENOME_H__
#define __GENOME_H__

#include <string>
#include <iostream>
#include <vector>
#include "FastaFile.h"
#include "FileOps.h"

using namespace std;

class Chromosome
{
	private:
		string name_;
		string sequence_;
		unsigned int length_;
	
	public:
		Chromosome();
		Chromosome(FastaFile & fastaFile);
		Chromosome(const string & fastaFileName);
		Chromosome(const string & name, const string & sequence);
		Chromosome(const Chromosome & chromosome);
		
		Chromosome & operator=(const Chromosome & rhs);
		
		unsigned int length() const {return length_;}
		string name() const {return name_;}
		
		void set(const string & name, const string & sequence);
		void fullSeq(string & sequence) const {sequence = sequence_;}
		
		string seq(unsigned int start, unsigned int stop) const;
		void seq(unsigned int start, unsigned int stop, string & sequence) const;
		char seq(unsigned int pos) const;
		
		void insert(unsigned int start, const string & sequence);
		void remove(unsigned int start, unsigned int length);
		void substitute(unsigned int position, char sub);
		
		void print(FILE * fp, unsigned int ntPerLine) const;
};

class Genome
{
	private:
		vector< Chromosome > chromosomes_;
		string name_;
		unsigned int length_;
		unsigned int nChr_;

		unsigned int ntPerLine_;

		void loadFromFasta(FastaFile & fastaFile);

	public:
		Genome();
		Genome(const string & FastaName);
		Genome(const string & FastaName, const string & genomeName);
		Genome(const string & FastaName, unsigned int ntPerLine);
		Genome(const string & FastaName, unsigned int ntPerLine, const string & genomeName);
		Genome(FastaFile & fastaFile, const string & genomeName);
		Genome(FastaFile & fastaFile, unsigned int ntPerLine, const string & genomeName);
		Genome(const Genome & genome);
		
		Genome & operator=(const Genome & rhs);
		
		void ntPerLine(unsigned int ntPerLine) {ntPerLine_ = ntPerLine;}
		void defaultNtPerLine() {ntPerLine_ = 50;}
		unsigned int ntPerLine() const {return ntPerLine_;}

		unsigned int length() const {return length_;}
		void chromosomes(vector< Chromosome > & chromosomes) const {chromosomes = chromosomes_;}
		unsigned int nChr() const {return nChr_;}
		unsigned int numChr() const {return nChr_;}
		void chromosome (const string & chrName, Chromosome & chromosome) const {chromosome = chromosomes_ [chrIndex(chrName)];}
		void chromosome (unsigned int chrIndex, Chromosome & chromosome) const {chromosome = chromosomes_ [chrIndex];}
		unsigned int lChr(const string & chrName) const {return chromosomes_ [chrIndex(chrName)] . length();}
		unsigned int lChr(unsigned int chrIndex) const {return chromosomes_ [chrIndex] . length();}
		string chrName(unsigned int chrIndex) const {return chromosomes_ [chrIndex] . name();}
		void chrName(unsigned int chrIndex, string & name) const {name = chromosomes_ [chrIndex] . name();}

		string name() const {return name_;}
		void name(const string & name) {name_ = name;}

		void tabLength();
		void tabNChr() {nChr_ = chromosomes_ . size();}

		unsigned int chrIndex(const string & chrName) const;
		unsigned int chromosomeIndex(const string & chrName) const {return chrIndex(chrName);}

		string seq  (const string & chrName, unsigned int start, unsigned int stop) const;
		string seq	(unsigned int chrIndex, unsigned int start, unsigned int stop) const;
		void seq    (const string & chrName, unsigned int start, unsigned int stop, string & sequence) const;
		void seq		(unsigned int chrIndex, unsigned int start, unsigned int stop, string & sequence) const;
		char seq    (const string & chrName, unsigned int pos) const;
		char seq		(unsigned int chrIndex, unsigned int pos) const;

		void insert (const string & chrName, unsigned int start, const string & seq);
		void insert (unsigned int chrIndex, unsigned int start, const string & seq);
		void remove (const string & chrName, unsigned int start, unsigned int length);
		void remove (unsigned int chrIndex, unsigned int start, unsigned int length);
		void substitute (const string & chrName, unsigned int position, char sub);
		void substitute (unsigned int chrIndex, unsigned int position, char sub);

		void print (FILE * fp) const;
		void printChrNames () const;
		void printChrInfo () const;
};

#endif /* Genome.h */
