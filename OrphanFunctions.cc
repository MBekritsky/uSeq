#include "OrphanFunctions.h"

using namespace std;

void reconstituteCommandLine(int argc, char * argv[], string & s)
{
	s = argv[0];
	for(int i = 1; i < argc; ++i)
	{
		s += " ";
		s += argv[i];
	}
}

string reconstituteCommandLine(int argc, char * argv[])
{
	string s;
	reconstituteCommandLine(argc,argv,s);
	return s;
}
