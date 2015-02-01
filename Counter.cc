#include "Counter.h"

using namespace std;

Counter::Counter()
{
	in  = 0;
	out = 0;
}

void Counter::increment(const string & key)
{
	if(key.compare("in") == 0)
	{
		in++;
	}
	else if(key.compare("out") == 0)
	{
		out++;
	}
	else if(counter.find(key) == counter.end())
	{
		counter[key] = 1;
	}
	else
	{
		counter[key]++;
	}
}

void Counter::increment(const string & key,int increment)
{
	if(key.compare("in") == 0)
	{
		in += increment;
	}
	else if(key.compare("out") == 0)
	{
		out += increment;
	}
	else if(counter.find(key) == counter.end())
	{
		counter[key] = increment;
	}
	else
	{
		counter[key] += increment;
	}
}

int Counter::value(const string & key) const
{
	map<string, int>::const_iterator it;
	it = counter.find(key);

	if(key.compare("in") == 0)
	{
		return in;
	}
	else if(key.compare("out") == 0)
	{
		return out;
	}
	else if(it != counter.end())
	{
		return it->second;
	}
	return -1;
}

/* void Counter::operator++(const string key)
{
	counter.increment(key);
}

void Counter::operator+=(const string key,int increment)
{
	counter.increment(key,increment);
} */

void Counter::print(FILE * output) const
{
	map<string,int>::const_iterator it;
	fprintf(output,"reads in: %d\n",in);
	if(counter.size() > 0)
	{
		for(it = counter.begin(); it != counter.end(); it++)
		{
			fprintf(output,"%s: %d\n",it->first.c_str(),it->second);
		}
	}
	fprintf(output,"reads out: %d\n",out);	
}

void Counter::print() const
{
	print(stderr);
}

string createCounterDir(const string & filename,const string & outputPrefix,const string & targetDir)
{
	string counterDir;
	if(outputPrefix . size())
	{
	 counterDir = addSuffix(getParentDirectory(outputPrefix),"counts");
	 if(targetDir . size())
	 {
		 counterDir = replacePath(counterDir, targetDir);
	 }
	}

	else if(targetDir.size() & !outputPrefix . size())
	{
	 counterDir = addSuffix(targetDir,"counts");
	}

	else if (filename . size())
	{
	 counterDir = addSuffix(getParentDirectory(filename),"counts");
	}
	else
	{
		cerr << "Failed to name counter directory, no parameters specified" << endl;
		exit(EXIT_FAILURE);
	}

	if(!dirExists(counterDir))
	{
	 bool countDirStatus = mkdir(counterDir . c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	 if(countDirStatus == 0)
	 {
		 cout << "Created " << counterDir << endl;
	 }
	 else
	 {
		 cout << "Failed to create " << counterDir << ": " << strerror(errno) << endl;
	 }
	}
	return(counterDir);
}

string createCounterDir(const string & filename)
{
	return(createCounterDir(filename,"",""));
}

