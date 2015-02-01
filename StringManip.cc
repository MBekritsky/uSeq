#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include "StringManip.h"

using namespace std;

void splitString(const string & input, const string & delim, vector< string > & output, int limit)
{
	size_t entryStart = 0;
	size_t nextDelim = 0;
	int numTokens = 0;
	string token;
	while(nextDelim != string::npos && numTokens < limit)
	{
		nextDelim = input.find_first_of(delim,entryStart);
		if(numTokens < (limit - 1))
		{
			token = input.substr(entryStart,nextDelim - entryStart);
		}
		else
		{
			token = input.substr(entryStart);
		}
		output.push_back(token);
		entryStart = nextDelim + 1;
		numTokens++;
	}
}

void splitString(const string & input, char delim, vector< string > & output, int limit)
{
	size_t entryStart = 0;
	size_t nextDelim = 0;
	int numTokens = 0;
	string token;
	while(nextDelim != string::npos && numTokens < limit)
	{
		nextDelim = input.find_first_of(delim,entryStart);
		if(numTokens < (limit - 1))
		{
			token = input.substr(entryStart,nextDelim - entryStart);
		}
		else
		{
			token = input.substr(entryStart);
		}
		output.push_back(token);
		entryStart = nextDelim + 1;
		numTokens++;
	}
}

void splitString(const string & input, const string & delim, vector< string > & output)
{

	size_t entryStart = 0;
	size_t nextDelim = 0;
	string token;
	while(nextDelim != string::npos)
	{
		nextDelim = input.find_first_of(delim,entryStart);
		token = input.substr(entryStart,nextDelim - entryStart);
		output.push_back(token);
		entryStart = nextDelim + 1;
	}
}

void splitString(const string & input, char delim, vector< string > & output)
{

	size_t entryStart = 0;
	size_t nextDelim = 0;
	string token;
	while(nextDelim != string::npos)
	{
		nextDelim = input.find_first_of(delim,entryStart);
		token = input.substr(entryStart,nextDelim - entryStart);
		output.push_back(token);
		entryStart = nextDelim + 1;
	}
}

void upperCase(string & str)
{
	transform(str.begin(),str.end(),str.begin(),(int (*)(int)) std::toupper);
}

void lowerCase(string & str)
{
	transform(str.begin(),str.end(),str.begin(),(int (*)(int)) std::tolower);
}

void padString(string & str, char ch, unsigned int times)
{
	for(unsigned int i = 0; i < times; ++i)
	{
		str += ch;
	}
}

//currently leaves trailing delimiter as legacy, needs to be fixed
void stringList2string(string & str, vector< string > list, char delim)
{
	for(unsigned int i = 0; i < list.size(); ++i)
	{
			str += (list[i] + delim);
	}
}
