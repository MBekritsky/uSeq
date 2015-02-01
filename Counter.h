#ifndef __COUNTER_H__
#define __COUNTER_H__

#include <iostream>
#include <string>
#include <map>
#include <errno.h>
#include <stdio.h>
#include <cstring>
#include <cstdlib>

#include "FileOps.h"

using namespace std;

class Counter
{
	private:
		map<string,int> counter;
		int in;
		int out;

	public:
		Counter();
		void increment(const string & key);
		void increment(const string & key,int increment);
		int value(const string & key) const;
/* 		void operator++(const string key);
		void operator+=(const string key,int increment);
 */		
 		void print(FILE * output) const;
		void print() const;
};

string createCounterDir(const string & filename,const string & outputPrefix,const string & targetDir);
string createCounterDir(const string & filename);

#endif /* Counter.h */
