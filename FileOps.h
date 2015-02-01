/*
FileOps.h

ABSTRACT: Some basic file operations for managing
frequently used file operations.

CREATION DATE:
23.08.2010

LAST REVISION:
27.08.2010

AUTHOR:
Mitchell Bekritsky
*/

#ifndef __FILE_OPS_H__
#define __FILE_OPS_H__

#include <string>
#include <iostream>
#include <sys/stat.h>

using namespace std;

string replaceSuffix(const string & filename,const string & suffix);
string addSuffix(const string & prefix, const string & suffix);
string replacePath(const string & filename, const string & newPath);
string addPath(const string & filename, const string & newPath);
FILE* openAndTestFile(const string & filename, const string & operation);
string getParentDirectory(const string & filename);
string getBaseName(const string & filename);
string getBaseWithPath(const string & filename);
bool fileExists(const string & filename);
bool dirExists(const string & filename);
string createFileName(const string & parentDir, const string & baseName, const string & ext);
void addTrailingSlash(string & dir);

#endif /* FileOps.h */
