#ifndef __STRING_MANIP_H__
#define __STRING_MANIP_H__

#include <string>
#include <vector>

using namespace std;

void splitString(const string & input, const string & delim, vector< string > & output);
void splitString(const string & input, const string & delim, vector< string > & output, int limit);
void splitString(const string & input, char delim, vector< string > & output);
void splitString(const string & input, char delim, vector< string > & output, int limit);
void upperCase(string & str); //avoids collision with std::ios_base&::uppercase
void lowerCase(string & str);
void padString(string & str,char ch,unsigned int times);
void stringList2string(string & str, vector< string > list, char delim);

#endif /*StringManip.h*/
