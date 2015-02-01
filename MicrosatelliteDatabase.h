#ifndef __MICROSATELLITE_DATABASE_H__
#define __MICROSATELLITE_DATABASE_H__

#include <string>
#include <vector>
#include <map>
#include "FileOps.h"
#include "Microsatellite.h"
#include "StringManip.h"
#include "UsefulMacros.h"

using namespace std;

//Since microsatellite database files produced by uSeq are already in sorted order, the database will also be sorted

class MicrosatelliteDatabase
{
	private:
		FILE * databaseP_; //should never get passed between functions, will create conflicts between classes
		string databaseName_;
		map< string, vector< Microsatellite > > database_;
		unsigned int totalRecords_;
		
	public:
		MicrosatelliteDatabase();
		MicrosatelliteDatabase(FILE * databaseP);
		MicrosatelliteDatabase(const string & databaseName);
		MicrosatelliteDatabase(const MicrosatelliteDatabase & rhs);
		~MicrosatelliteDatabase();

		MicrosatelliteDatabase & operator=(const MicrosatelliteDatabase & rhs);
		
		void databaseName(string & databaseName) const {databaseName = databaseName_;}
		string databaseName() const {return databaseName_;}
		void database(map< string, vector< Microsatellite > > & database) const {database = database_;}
		unsigned int totalRecords() const {return totalRecords_;}
		unsigned int length() const {return totalRecords_;}
		
		void load();
		void parseRecord(const string & record, Microsatellite & target) const;
		bool coordExists(const string & coord) const;
		void countTotalRecords();		
		
		unsigned int nCoords() const {return database_ . size();}
		unsigned int nMicrosatellites(const string & chrName) const;
		
		bool findRecordByStart(const string & coord, unsigned int start, Microsatellite & record) const;
		int binarySearch(const vector< Microsatellite > & chrdb, unsigned int key, int low, int high) const;
		int interpolationSearch(const vector< Microsatellite > & chrdb, unsigned int key, int low, int high) const;
		
//		Microsatellite & getRecord(unsigned int index) {return database_[index];}
//		void getRecord(unsigned int index, Microsatellite & temp) const {temp = database_ [index];}

//		void printRecords() const;
//		void printCoords() const;
};

#endif /* MicrosatelliteDatabase.h */
