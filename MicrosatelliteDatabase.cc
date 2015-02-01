#include "MicrosatelliteDatabase.h"

using namespace std;

MicrosatelliteDatabase::MicrosatelliteDatabase()
{
	databaseP_ = NULL;
	databaseName_.clear();
	database_ . clear();
	totalRecords_ = 0;
}

MicrosatelliteDatabase::MicrosatelliteDatabase(FILE * databaseP)
{
	databaseP_ = databaseP;
	databaseName_.clear();
	load();
}

MicrosatelliteDatabase::MicrosatelliteDatabase(const string & databaseName)
{
	databaseP_ = openAndTestFile(databaseName, "r");
	load();
}

MicrosatelliteDatabase::MicrosatelliteDatabase(const MicrosatelliteDatabase & rhs)
{
	rhs.databaseName(databaseName_);
	if(databaseName_.size() > 0) //copying the database will not copy the file pointer, it will create a new pointer to the same file
	{
		databaseP_ = openAndTestFile(databaseName_, "r");
	}
	rhs.database(database_);
	totalRecords_ = rhs.totalRecords();
}

MicrosatelliteDatabase & MicrosatelliteDatabase::operator=(const MicrosatelliteDatabase & rhs)
{
	if(this != &rhs)
	{
		rhs.databaseName(databaseName_);
		if(databaseName_.size() > 0) //copying the database will not copy the file pointer, it will create a new pointer to the same file
		{
			databaseP_ = openAndTestFile(databaseName_, "r");
		}
		rhs.database(database_);
		totalRecords_ = rhs.totalRecords();
	}
	return *this;
}

MicrosatelliteDatabase::~MicrosatelliteDatabase()
{
	fclose(databaseP_);
}

void MicrosatelliteDatabase::load()
{
	char ch = ' ';
	string record;
	string coord;
	
	Microsatellite microsatellite;

	while(ch != EOF)
		{
			while((ch = fgetc(databaseP_)) != EOF && ch != '\n')
				record . push_back(ch);
			
			parseRecord (record, microsatellite);

			microsatellite.getHdr(coord);

			database_[coord] . push_back(microsatellite);
			record . clear();

			if((ch = fgetc(databaseP_)) != EOF)
				ungetc(ch, databaseP_);
		}
	countTotalRecords();
}

void MicrosatelliteDatabase::countTotalRecords()
{
	totalRecords_ = 0;
	map< string, vector< Microsatellite > >::const_iterator it;
	tr(database_,it)
		{
			totalRecords_ += it->second.size();
		}
}

/*void MicrosatelliteDatabase::printRecords() const
{
	for(unsigned int i = 0; i < database_ . size(); i++)
		database_ [i] . printFastaStyle();
}
*/
void MicrosatelliteDatabase::parseRecord(const string & record, Microsatellite & target) const
{
	static string coord, unit;
	static unsigned int start, stop, length;
	vector< string > splitRecord;
	
	splitString(record,"\t",splitRecord);
	coord  = splitRecord[0];
	start  = atoi(splitRecord[1].c_str());
	stop   = atoi(splitRecord[2].c_str());	
	length = atoi(splitRecord[3].c_str());
	unit   = splitRecord[4];
	
	target . setBasics(unit, start, stop, length);
	target . setHdr(coord);
}

bool MicrosatelliteDatabase::coordExists(const string & coord) const
{
	if(database_.find(coord) == database_.end())
	{
		return false;
	}
	return true;
}

/*
void MicrosatelliteDatabase::printCoords() const
{
	for(unsigned int i = 0; i < coords_ . size(); i++)
		cout << coords_ [i] << "\t" << microsatellitesPerCoord_ [i] << endl;
}
*/
unsigned int MicrosatelliteDatabase::nMicrosatellites(const string & coord) const
{
	map< string, vector< Microsatellite > >::const_iterator it;
	it = database_.find(coord);
	if(it != database_.end())
	{
		return it->second.size();
	}
	else
	{
		cerr << coord << " cannot be found in the database" << endl;
	}
	return 0;
}

bool MicrosatelliteDatabase::findRecordByStart(const string & coord, unsigned int recordStart, Microsatellite & record) const
{
	if(coordExists(coord))
	{
		int index;
		map< string, vector< Microsatellite> >::const_iterator it = database_.find(coord); //look for the coordinate (e.g. chr1) in the database
		int high = (*it).second.size() - 1;
		int low  = 0;
		index = binarySearch((*it).second,recordStart,low,high);
		if(index > -1)
		{
			record = (*it).second[index];
			return true;
		}
	}
	else
	{
		cerr << coord << " cannot be found in the database" << endl;
	}
	return false;
}

int MicrosatelliteDatabase::interpolationSearch(const vector< Microsatellite > & chrdb, unsigned int key, int low, int high) const
{
	if(high < low)
	{
		return -1;
	}
	else
	{
		unsigned int lowStart = chrdb[low].start();
		unsigned int highStart = chrdb[high].start();

		double frac = (key - lowStart)/(highStart - lowStart);
		if(frac < 0 || frac > 1)
		{
			return -1;
		}

		unsigned int mid = low + (unsigned int) (frac * (high - low));
		if(chrdb[mid].start() > key)
		{
			return interpolationSearch(chrdb,key,low,mid-1);
		}
		else if(chrdb[mid].start() < key)
		{
			return interpolationSearch(chrdb,key,mid+1,high);
		}
		else
		{
			return mid;
		}
	}
}

int MicrosatelliteDatabase::binarySearch(const vector< Microsatellite > & chrdb, unsigned int key, int low, int high) const
{
	if(high < low)
	{
		return -1;
	}
	else
	{
		int mid = ((high + low)/2);
		if(chrdb[mid].start() > key)
		{
			return binarySearch(chrdb,key,low,mid - 1);
		}
		else if(chrdb[mid].start() < key)
		{
			return binarySearch(chrdb,key,mid + 1,high);
		}
		else
		{
			return mid;
		}
	}
}
