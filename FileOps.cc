/*
FileOps.cc

ABSTRACT: Some basic file operations for managing
frequently used file operations.

CREATION DATE:
23.08.2010

LAST REVISION:
27.08.2010

AUTHOR:
Mitchell Bekritsky
*/

#include "FileOps.h"
#include <cstdio>
#include <cstdlib>

using namespace std;

FILE * openAndTestFile(const string & filename, const string & operation)
{
	FILE *fp = fopen(filename . c_str(), operation . c_str());
  if(!fp)
    {
      cerr << "Could not open " << filename . c_str() << endl;
      exit(EXIT_FAILURE);
    }
	return fp;
}

string replaceSuffix(const string & filename, const string & suffix)
{
  size_t lastSlash = filename . find_last_of("/");
	size_t firstDotAfterSlash;
	
	if(lastSlash != string::npos)
		{
		  firstDotAfterSlash = filename . find_first_of(".",lastSlash);
		}
	else
		{
			firstDotAfterSlash = filename . find_first_of(".");
		}
	
  if(firstDotAfterSlash == string::npos)
    {
      cerr << "Error! Cannot process filename!! No extension in filename " << filename << endl;
      exit(EXIT_FAILURE);
    }
  
  string outname = filename;
  outname . erase(firstDotAfterSlash);
  outname . append(suffix);
  
  return outname;
}

string addSuffix(const string & prefix, const string & suffix)
{
	string outname = prefix;
	outname . append(suffix);
	
	return outname;
}

string replacePath(const string & filename, const string & newPath)
{
	size_t lastSlash = filename . find_last_of("/");
	
	if(lastSlash == filename . length() - 1)
		{
			cerr << "Error!  You have entered what appears to be a directory, not a filename!!" << endl;
			exit (EXIT_FAILURE);
		}
	else
		{
			string newFile = newPath;
			if (newPath . find_last_of("/") != newPath . length() - 1)
				newFile += "/";
			newFile += filename . substr(lastSlash + 1);
			return newFile;
		}
	return "";
}

string addPath(const string & filename, const string & newPath)
{
	string outname = newPath;
	if((*newPath.rbegin()) != '/')
	{
		outname.append("/");
	}
	outname.append(filename);
	return outname;
}

string createFileName(const string & parentDir, const string & baseName, const string & ext)
{
	string fileName;
	
	fileName = parentDir;
	if(parentDir . find_last_of("/") != parentDir . length() - 1)
		{
			fileName += "/";
		}
	fileName += baseName;
	if(baseName . find_last_of(".") != parentDir . length() && ext . find_first_of(".") != 0)
		{
			fileName += ".";
		}
	fileName += ext;
	return fileName;
}

string getParentDirectory(const string & filename)
{
	size_t lastSlash = filename . find_last_of("/");
	//if entire string is filename, return entire file string
	if(lastSlash == filename . length() - 1)
		{
			return filename;
		}
	//if last slash is internal, check if trailing identifier is a file or directory
	return filename . substr(0, lastSlash + 1);
}

string getBaseName(const string & filename)
{
	size_t lastSlash = filename . find_last_of("/");
	size_t firstDotAfterSlash;
	if(lastSlash == filename . length() - 1)
		{
			return "";
		}
	else if(lastSlash == string::npos)
		{
			firstDotAfterSlash = filename . find_first_of(".");
		}
	else
		{
			firstDotAfterSlash = filename . find_first_of(".",lastSlash);
		}

	if(firstDotAfterSlash != string::npos)
		{
			return filename . substr(lastSlash + 1, firstDotAfterSlash - lastSlash - 1);
		}
	else
		{
			return filename . substr(lastSlash + 1);
		}
	return "";
		
}

string getBaseWithPath(const string & filename)
{
	size_t lastSlash = filename . find_last_of("/");
	size_t firstDotAfterSlash;
	
	if(lastSlash == filename . length() - 1)
		{
			return "";
		}
	else if (lastSlash == string::npos)
		{
			firstDotAfterSlash = filename . find_first_of(".");
		}
	else
		{
			firstDotAfterSlash = filename . find_first_of(".",lastSlash);
		}
		
	if(firstDotAfterSlash != string::npos)
		{
			return filename . substr(0, firstDotAfterSlash);
		}
	else
		{
			return filename;
		}
	return "";
}

bool fileExists(const string & filename)
{
	struct stat stFileInfo;
	int intStat;
	
	intStat = stat(filename . c_str(), &stFileInfo);
	if(intStat == 0)
		{
			if(S_ISREG(stFileInfo . st_mode))
				return true;
			return false;
		}
	return false;
}

bool dirExists(const string & filename)
{
	struct stat stFileInfo;
	int intStat;
	
	intStat = stat(filename . c_str(), &stFileInfo);
	if(intStat == 0)
		{
			if(S_ISDIR(stFileInfo . st_mode))
				return true;
			return false;
		}
	return false;
}

void addTrailingSlash(string & dir)
{
	if(*dir.rbegin() != '/')
	{
		dir = dir + "/";
	}
}
