#include "MicrosatelliteProfile.h"

using namespace std;
using namespace BamTools;

MicrosatelliteProfile::MicrosatelliteProfile()
{
	coordFlags_ . clear();
	qualFlags_  . clear();
	alignments_ . clear();
	reportable_ . clear();
	count_ 			 = 0;
	genomeStart_ = 0;
	genomeStop_  = 0;
	name_ = "";
	coord_ = "";
	motif_  = "";
}

//A microsatellite profile needs to know which microsatellite from a read it is profiling
//if there are multiple microsatellites in a read, that read can be in several profiles
MicrosatelliteProfile::MicrosatelliteProfile(const MicrosatelliteCoordFlag & coordFlag)
{
	coordFlags_ . clear();
	qualFlags_  . clear();
	alignments_ . clear();
	reportable_ . clear();
	
	count_ = 0;
	genomeStart_ = coordFlag.genomeStart();
	genomeStop_  = coordFlag.genomeStop();
	coordFlag.getMotif(motif_);
	coordFlag.getCoord(coord_);
	
	char buffer[1024];
	sprintf(buffer,"%d%s%d%s%d",(int) coord_.length(),coord_.c_str(),(int) motif_.length(),motif_.c_str(),genomeStart_);
	name_ = buffer;
}

MicrosatelliteProfile::MicrosatelliteProfile(const MicrosatelliteProfile & rhs)
{
	rhs.coordFlags(coordFlags_);
	rhs.qualFlags(qualFlags_);
	rhs.alignments(alignments_);
	rhs.reportable(reportable_);
	rhs.motif(motif_);
	rhs.coord(coord_);
	rhs.name(name_);
	
	count_ = rhs.count();
	genomeStart_ = rhs.genomeStart();
	genomeStop_  = rhs.genomeStop();
}

MicrosatelliteProfile & MicrosatelliteProfile::operator=(const MicrosatelliteProfile & rhs)
{
	if(this != &rhs)
	{
		rhs.coordFlags(coordFlags_);
		rhs.qualFlags(qualFlags_);
		rhs.alignments(alignments_);
		rhs.reportable(reportable_);
		rhs.motif(motif_);
		rhs.coord(coord_);
		rhs.name(name_);

		count_ = rhs.count();
		genomeStart_ = rhs.genomeStart();
		genomeStop_  = rhs.genomeStop();
	}
	return *this;
}

void MicrosatelliteProfile::addRead(const BamAlignment & al, const MicrosatelliteCoordFlag & coordFlag, const MicrosatelliteQualityFlag & qualFlag)
{
	coordFlags_ . push_back(coordFlag);
	qualFlags_ . push_back(qualFlag);
	alignments_ . push_back(al);
	reportable_ . push_back(true);
	pushBareName(al,count_); //count will be the same as the index in the coordFlags_,qualFlags_ and alignments_ vectors
	count_++;
}

//return the sequences query name without any microsatellite info appended to it
void MicrosatelliteProfile::pushBareName(const BamAlignment & al,unsigned int index)
{
	vector< string > splitName;
	splitString(al.Name,"|",splitName);
	bareNames_[splitName[0]] = index;
}

//return the sequences query name without any microsatellite info appended to it
void MicrosatelliteProfile::bareName(const BamAlignment & al, string & bareName) const
{
	vector< string > splitName;
	splitString(al.Name,"|",splitName);
	bareName = splitName[0];
}

string MicrosatelliteProfile::bareName(const BamAlignment & al) const
{
	vector< string > splitName;
	splitString(al.Name,"|",splitName);
	return splitName[0];
}

bool MicrosatelliteProfile::overlappingReadFragments(const BamAlignment & al, const MicrosatelliteCoordFlag & coordFlag, Counter & counter, BamWriter & seqNoise)
{
	string readName;
	bareName(al,readName);
	
	bnIt_ = bareNames_.find(readName);
	if(bnIt_ == bareNames_.end()) //if the read's mate pair is not already in the list of bare names for reads with this microsatellite
	{
		return false;
	}
	else
	{
		unsigned int index = bareNames_[readName];
		//if both reads report the same information about the microsatellite, don't add the second read to the list
		if(coordFlags_[index].motif() == coordFlag.motif() && coordFlags_[index].length() == coordFlag.length() && coordFlags_[index].genomeStart() == coordFlag.genomeStart())
		{
			counter.increment("overlappingReadFragment");
		}
		else //if the reads report discordant information, print them to the discordant read file, and mark the read already in the profile as unreportable
		{
			counter.increment("discordantOverlappingReadFragment");
			unsigned int discordantReadPairInd = (*bnIt_).second;
			if(!seqNoise.SaveAlignment(alignments_[discordantReadPairInd]))
			{
				cerr << "Error writing discordant alignment with read ID " << al.Name << " to file: " << seqNoise.GetErrorString() << endl;
				exit(EXIT_FAILURE);
			}
			if(!seqNoise.SaveAlignment(al))
			{
				cerr << "Error writing discordant alignment with read ID " << al.Name << " to file: " << seqNoise.GetErrorString() << endl;
				exit(EXIT_FAILURE);
			}
			
			reportable_[discordantReadPairInd] = false; //in the perl version of the profiler, the first read in the
																									//pair would get printed
		}
	}
	return true;
}

void MicrosatelliteProfile::print(const MicrosatelliteDatabase & msdb, unsigned int flankLength, FILE * out, FILE * dups, Counter & counter, bool keepDuplicates)
{
	Microsatellite refMicrosatellite;
	bool print = true;
	if(coord_.empty())
	{
		print = false;
		counter.increment("noCoord");
	}

	if(print)
	{
		if(!msdb.findRecordByStart(coord_,genomeStart_,refMicrosatellite))
		{
			counter.increment("nonReference");
		}
		else
		{
			if(refMicrosatellite.unit() != motif_)
			{
				counter.increment("wrongUnit");
				print = false;
			}
		}
		
		if(!keepDuplicates) condenseDuplicates(refMicrosatellite.length(),dups,counter);
		
		unsigned int numReported = numReportable();
		unsigned int maxLength = 0;

		if(numReported > 1 && print) //i.e. if the entire locus isn't just covered by one duplicate set (in profiler.cc, we guaranteed that there are at least minimumCount reads
		{
			counter.increment("reported");
			counter.increment("out",numReported);
			map< string, vector< unsigned int > > flankMap;
			unsigned int numFlanks = populateFlankMap(flankLength,flankMap);

			printMicrosatelliteInfo(refMicrosatellite,out);
			printFlankInfo(numFlanks,0,string("all"),numReported,out);
			maxLength = printMicrosatelliteHist(out);
			listQnames(out,','); //this prints the names of all the reads to the file, not just the ones that are reported...not sure this is the right thing to do
			fprintf(out,"\n");

			vector< string > barenameVec;
			barenames2Vec(barenameVec);

			vector< string > flanksByCounts;
			sortFlanks(flankMap,flanksByCounts);
			string cFlank;

			for(unsigned int i = 0; i < flanksByCounts.size(); ++i)
			{
				cFlank = flanksByCounts[i];
				printMicrosatelliteInfo(refMicrosatellite,out);
				printFlankInfo(numFlanks,i+1,cFlank,flankMap[cFlank].size(),out);
				printMicrosatelliteSubpopHist(flankMap[cFlank],out,maxLength);
				listFlankQnames(flankMap[cFlank],barenameVec,out,',');
				fprintf(out,"\n");
			}
		}
	}
}

void MicrosatelliteProfile::barenames2Vec(vector < string > & bnVec) const
{
	bnVec.resize(bareNames_.size());
	tr(bareNames_,bnIt)
	{
		bnVec[(*bnIt).second] = (*bnIt).first;
	}
}

//prints general microsatellite information to file
void MicrosatelliteProfile::printMicrosatelliteInfo(const Microsatellite & ref, FILE * out) const
{
	fprintf(out,"%s\t%d\t%s\t%d\t",coord_.c_str(),genomeStart_,motif_.c_str(),ref.length()); //if no reference for microsatellite is found, refMicrosatellite.start() = 0
}

//prints flank information to profile
void MicrosatelliteProfile::printFlankInfo(unsigned int numFlanks, unsigned int rank, const string & flank, unsigned int flankPop, FILE * out) const
{
	fprintf(out,"%d\t%d\t%s\t%d\t",numFlanks,rank,flank.c_str(),flankPop);
}


void MicrosatelliteProfile::listQnames(FILE * out, char delim) const
{
	fprintf(out,"\t");
/*	tr(bareNames_,bnIt)
	{
		fprintf(out,"%s%c",(*bnIt).first.c_str(),delim);
	}*/
	for(unsigned int i = 0; i < alignments_.size(); ++i)
	{
		if(reportable_[i])
			fprintf(out,"%s%c",alignments_[i].Name.c_str(),delim);
	}
}

void MicrosatelliteProfile::listFlankQnames(const vector< unsigned int > & indices, const vector< string > & barenameVec, FILE * out, char delim) const
{
	vector< unsigned int >::const_iterator it;
	fprintf(out,"\t");
/*	tr(indices,it)
	{
		if(reportable_[(*it)])
		{
			fprintf(out,"%s%c",barenameVec[(*it)].c_str(),delim);
		}
	}*/
	for(unsigned int i = 0; i < alignments_.size(); ++i)
	{
		if(reportable_[i])
		{
			fprintf(out,"%s%c",alignments_[i].Name.c_str(),delim);
		}
	}

}

void MicrosatelliteProfile::sortFlanks(const map< string, vector< unsigned int > > & flankMap, vector< string > & sortedFlanks) const
{
	multimap< unsigned int, string > counts2flanks; //using a multimap since it will sort by number of reads with flank and allows for multiple flanks to have same number of reads
	map< string, vector< unsigned int > >::const_iterator it;
	tr(flankMap, it)
	{
		counts2flanks.insert(std::pair<unsigned int, string> ((*it).second.size(), (*it).first));
	}
	
	multimap< unsigned int, string >::const_reverse_iterator cfIt;
	for(cfIt = counts2flanks.rbegin(); cfIt != counts2flanks.rend(); ++cfIt) //since multimaps store keys in order from least to greatest, traverse multimap in reverse
	{
		sortedFlanks.push_back((*cfIt).second);
	}
}

unsigned int MicrosatelliteProfile::printMicrosatelliteHist(FILE * out) const
{
	vector< unsigned int > profileHistogram;
	unsigned int maxLength = getProfileHistogram(profileHistogram);
	printHistogram(profileHistogram,out,',');
	return maxLength;
}

unsigned int MicrosatelliteProfile::printMicrosatelliteSubpopHist(const vector< unsigned int > & indices, FILE * out) const
{
	vector< unsigned int > subprofileHistogram;
	unsigned int maxLength = getSubprofileHistogram(indices,subprofileHistogram);
	printHistogram(subprofileHistogram,out,',');
	return maxLength;
}

unsigned int MicrosatelliteProfile::printMicrosatelliteSubpopHist(const vector< unsigned int > & indices, FILE * out, unsigned int maxLength) const
{
	vector< unsigned int > subprofileHistogram;
	unsigned int oMaxLength = getSubprofileHistogram(indices,subprofileHistogram,maxLength);
	printHistogram(subprofileHistogram,out,',');
	return oMaxLength;
}

unsigned int MicrosatelliteProfile::populateFlankMap(unsigned int flankLength, map< string, vector< unsigned int > > & flankMap) const
{
	string flankingSequence;
	for(unsigned int i = 0; i < reportable_ .size(); ++i)
	{
		flankingSequence.clear();
		if(reportable_[i])
		{
			getFlank(i,flankingSequence,flankLength);
			flankMap[flankingSequence].push_back(i);
		}
	}
	return flankMap.size();
}

void MicrosatelliteProfile::getFlank(unsigned int index, string & flank, unsigned int flankLength) const
{
	string originalSequence;
	vector< MicrosatelliteQualityFlag > alQualFlags;
	fillMicrosatelliteQualityDetails(alignments_[index],alQualFlags);
	reconstructSequence(alignments_[index],alQualFlags,originalSequence);
	
	if(qualFlags_[index].readStart() > flankLength) //if the microsatellite has sufficient flanking sequence to get a complete flank
	{
		flank += originalSequence.substr(qualFlags_[index].readStart() - flankLength - 1,flankLength);
	}
	else //if the microsatellite does not have sufficient flanking sequence in the read (e.g. flank is length 5 and microsatellite starts at base 4, only 3 reads of flanking sequence)
	{
		padString(flank,'#',flankLength - qualFlags_[index].readStart() + 1);
		flank += originalSequence.substr(0,qualFlags_[index].readStart() - 1);
	}
	
	//the 3' flank always starts at the same spot, and string::substr will take as many letters as possible from flankLength
	//if there is not flankLength sequence after the microsatellite in the read, we can pad it afterwards
	flank += originalSequence.substr(qualFlags_[index].readStop(),flankLength);
	if( (qualFlags_[index].readStop() + flankLength - 1) >= originalSequence.length())
	{
		padString(flank,'#',flankLength - (originalSequence.length() - qualFlags_[index].readStop()));
	}
}

unsigned int MicrosatelliteProfile::numReportable() const
{
	unsigned int numReportable = 0;
	vector< bool >::const_iterator it;
	tr(reportable_, it)
	{
		numReportable += (*it);
	}
	return numReportable;
}

void MicrosatelliteProfile::condenseDuplicates(unsigned int refLength, FILE * dups, Counter & counter) //this could be checking OP tags instead, difference is probably negligible
{
	map< string, bool > readPairEncountered;
	string dupKey;
	for(unsigned int i = 0; i < alignments_.size(); ++i)
	{
		dupKey.clear();
		getDupKeyName(alignments_[i],dupKey);
		if(alignments_[i].IsDuplicate() && readPairEncountered.find(dupKey) == readPairEncountered.end()) //if Picard MarkDups thinks it's found PCR duplicates, then find all reads at the microsatellite
		//locus that meet the criteria for being a PCR duplicate, i.e. have the same start for both the read and its mate pair, provided we haven't already looked for the duplicates
		{
			readPairEncountered[dupKey] = true;
			vector< unsigned int > duplicateIndices(1,i);
			findReadsWithSameStarts(duplicateIndices);
			vector< unsigned int > dupLengthHistogram;
			unsigned int numLengths = getDupLengths(duplicateIndices,dupLengthHistogram); //lengths in histogram are zero-indexed--count for length 8 is found at index 7
			unsigned int numReads = duplicateIndices.size();

			if(numReads > 1) //if more than one read exists in the duplicate pair (not sure this conditional is really necessary)
			{
				counter.increment("msInPCRDuplicates",numReads);
				pickRepresentativeDuplicate(duplicateIndices,numLengths);
				printDuplicateInfo(duplicateIndices,dupLengthHistogram,refLength,dups);
			}
		}
	}
}

void MicrosatelliteProfile::printDuplicateInfo(const vector< unsigned int > & duplicateIndices, const vector < unsigned int > & dupLengthHistogram, 
																							 unsigned int refLength, FILE * dups) const
{
	unsigned int profileMode = getTractLengthMode();
	unsigned int dupMode = 0;
	unsigned int dupModeCov = 0;
	for(unsigned int i = 0; i < dupLengthHistogram.size(); ++i) //get the length in the duplicate set with the highest coverage
	{
		if(dupLengthHistogram[i] > dupModeCov)
		{
			dupModeCov = dupLengthHistogram[i];
			dupMode = i + 1;
		}
	}
	
	unsigned int dupNotModeCov = duplicateIndices.size() - dupModeCov; //get the number of reads not at the duplicate set mode
	
	fprintf(dups,"%s\t%d\t%s\t%lu\t%d\t%d",coord_.c_str(),genomeStart_,motif_.c_str(),
																			 (unsigned long) motif_.length(),refLength,alignments_[duplicateIndices[0]].Position);
	fprintf(dups," %d\t%.2f\t%lu ",profileMode,(float) dupNotModeCov/duplicateIndices.size(),(unsigned long) duplicateIndices.size());
	printHistogram(dupLengthHistogram,dups,',');
	fprintf(dups,"\n");
}

void MicrosatelliteProfile::printHistogram(vector< unsigned int > histogram, FILE * fh, char delim) const //generic histogram printer for printing 
//duplicate length and profile length histograms
{
	for(unsigned int i = 0; i < histogram.size()-1; ++i)
	{
		fprintf(fh,"%d%c",histogram[i],delim);
	}
	fprintf(fh,"%d",histogram[histogram.size() - 1]);
}

unsigned int MicrosatelliteProfile::getTractLengthMode() const
{
	map< unsigned int, unsigned int > histMap;
	map< unsigned int, unsigned int >::const_iterator it;
	getProfileHistMap(histMap);
	unsigned int maxCov = 0;
	unsigned int mode = 0;
	tr(histMap,it)
	{
		if((*it).second > maxCov)
		{
			mode = (*it).first;
			maxCov = (*it).second;
		}
	}
	return mode;
}

void MicrosatelliteProfile::getProfileHistMap(map< unsigned int, unsigned int > & histMap) const
{
	for(unsigned int i = 0; i < coordFlags_.size(); ++i)
	{
		if(reportable_[i])
		{
			histMap[coordFlags_[i].length()]++;
		}
	}
}

void MicrosatelliteProfile::getSubprofileHistMap(const vector< unsigned int > & indices, map< unsigned int, unsigned int > & histMap) const
{
	unsigned int cInd;
	for(unsigned int i = 0; i < indices.size(); ++i)
	{
		cInd = indices[i];
		if(reportable_[cInd])
		{
			histMap[coordFlags_[cInd].length()]++;
		}
	}
}

void MicrosatelliteProfile::pickRepresentativeDuplicate(const vector< unsigned int > & duplicateIndices, unsigned int numLengths)
{
	unsigned int currIndex;
	if(numLengths == 1) //if the reads in the PCR duplicate all report the same length for the microsatellite
	{
		unsigned int maxQual = 0, maxQualIndex = 0;
		for(unsigned int i = 0; i < duplicateIndices.size(); ++i) //find the highest quality read in the duplicate set that does not have a duplicate flag
		{
			currIndex = duplicateIndices[i];
			if(!alignments_[currIndex].IsDuplicate())
			{
				if(alignments_[currIndex].MapQuality > maxQual)
				{
					maxQual = alignments_[currIndex].MapQuality;
					maxQualIndex = i;
				}
			}
		}
		for(unsigned int i = 0; i < duplicateIndices.size(); ++i) //set reportable to false for all but the single best read in the duplicate set, mark all reads in this duplicate set
		{
			if(i != maxQualIndex)
			{
				currIndex = duplicateIndices[i];
				reportable_[currIndex] = false;
			}
		}
	}
	else //if all reads in the PCR duplicate set do not have the same length, do not report any of them in the profile, since we don't know which of the tract lengths is the true length
	//of the originating fragment
	{
		for(unsigned int i = 0; i < duplicateIndices.size(); ++i)
		{
			currIndex = duplicateIndices[i];
			reportable_[currIndex] = false;
		}
	}
}

unsigned int MicrosatelliteProfile::getProfileHistogram(vector< unsigned int > & profileLengthHistogram) const
{
	map< unsigned int, unsigned int > profileHistMap;
	getProfileHistMap(profileHistMap);
	
	unsigned int maxLength = profileHistMap.rbegin()->first;
	fillHistVector(profileHistMap,profileLengthHistogram,maxLength);
	return maxLength; //the size of the histogram corresponds to the the longest tract length observed in the profile
}

unsigned int MicrosatelliteProfile::getSubprofileHistogram(const vector< unsigned int > & indices, vector< unsigned int > & subprofileLengthHistogram) const
{
	map< unsigned int, unsigned int > subprofileHistMap;
	getSubprofileHistMap(indices,subprofileHistMap);
	
	unsigned int maxLength = subprofileHistMap.rbegin()->first;
	fillHistVector(subprofileHistMap,subprofileLengthHistogram,maxLength);
	return maxLength;
}

//if you want to print a subprofile histogram that has a fixed length (e.g. the maximum size in the complete profile), use this function
unsigned int MicrosatelliteProfile::getSubprofileHistogram(const vector< unsigned int > & indices, vector< unsigned int > & subprofileLengthHistogram, unsigned int maxLength) const
{
	map< unsigned int, unsigned int > subprofileHistMap;
	getSubprofileHistMap(indices,subprofileHistMap);
				
	fillHistVector(subprofileHistMap,subprofileLengthHistogram,maxLength);
	return maxLength;
}

void MicrosatelliteProfile::fillHistVector(const map< unsigned int, unsigned int > & histMap, vector< unsigned int > & histVec, unsigned int maxLength) const
{
	histVec.resize(maxLength, 0);
	
	map< unsigned int, unsigned int >::const_iterator it;
	tr(histMap,it)
	{
		histVec[(*it).first - 1] = (*it).second;
	}
}

unsigned int MicrosatelliteProfile::getDupLengths(vector< unsigned int > & duplicateIndices, vector< unsigned int > & dupLengthHistogram) const
{
	map< unsigned int, unsigned int > histMap;
	unsigned int currIndex;
	for(unsigned int i = 0; i < duplicateIndices.size(); ++i)
	{
		currIndex = duplicateIndices[i];
		histMap[coordFlags_[currIndex].length()]++;
	}

	//longest duplicate in map will be found at histMap.rbegin()
	unsigned int maxLength = histMap.rbegin()->first;
	dupLengthHistogram.resize(maxLength,0);

	map< unsigned int, unsigned int >::const_iterator it;
	tr(histMap,it)
	{
		dupLengthHistogram[(*it).first - 1] = (*it).second;
	}
	return histMap.size();
}

void MicrosatelliteProfile::findReadsWithSameStarts(vector< unsigned int > & duplicateIndices) const
{
	unsigned int firstIndex = duplicateIndices[0];
	for(unsigned int i = 0; i < alignments_.size(); ++i)
	{
		if(i != firstIndex)
		{
			if(alignments_[i].Position == alignments_[firstIndex].Position && 
				 alignments_[i].MatePosition == alignments_[firstIndex].MatePosition)
			{
				duplicateIndices.push_back(i);	 	
			}
			else if(alignments_[i].Position == alignments_[firstIndex].MatePosition && 
				 alignments_[i].MatePosition == alignments_[firstIndex].Position)
			{
				duplicateIndices.push_back(i);
			}
		}
	}
}

void MicrosatelliteProfile::getDupKeyName(const BamAlignment & al, string & key) const
{
	char buffer[1024];
	sprintf(buffer,"%d_%d",al.Position + 1,al.MatePosition + 1); //BamTools Position and MatePosition are 0-indexed
	key = buffer;
}

string MicrosatelliteProfile::getDupKeyName(const BamAlignment & al) const
{
	string key;
	getDupKeyName(al,key);
	return key;
}

//key structure is (motif length)+(motif)(genome_position)
void getProfileKeyName(int refID, const MicrosatelliteCoordFlag & coordFlag, string & name)
{
	char buffer[1024];
	sprintf(buffer,"%d_%s_%d",refID,coordFlag.motif().c_str(),coordFlag.genomeStart());
	name = buffer;
}

string getProfileKeyName(int refID, const MicrosatelliteCoordFlag & coordFlag)
{
	string name;
	getProfileKeyName(refID, coordFlag,name);
	return name;
}

int extractRefID(const string & key)
{
	return atoi(key.substr(0,key.find_first_of("_")).c_str()); 
}

unsigned int extractStart(const string & key)
{
	return atoi(key.substr(key.find_last_of("_") + 1).c_str()); 
}

string extractMotif(const string & key)
{
	return key.substr(key.find_first_of("_") + 1, key.find_last_of("_") - key.find_first_of("_") - 1);
}

bool KeySorter::operator()(const string & a, const string & b) const
{
	if(extractRefID(a) < extractRefID(b)) // A is on an earlier chromosome than B, A goes first
	{
		return true;
	}
	if(extractRefID(a)  > extractRefID(b)) // B is on an earlier chromosome than A, B goes first
	{
		return false;
	}
	if(extractStart(a) < extractStart(b)) // if A is earlier in the chromsome than B, A goes first
	{
		return true;
	}
	else if(extractStart(a) > extractStart(b)) // if B is earlier in the chromsome than A, B goes first
	{
		return false;
	}
	else if(extractMotif(a) != extractMotif(b)) //if A and B are at the same locus, but have differing units
	{
		if(extractMotif(a).length() < extractMotif(b).length()) //if A is shorter than B (e.g. A is C, B is AC), A goes first
		{
			return true;
		}
		else if(extractMotif(a).length() > extractMotif(b).length()) //if B is shorter than A, B goes first
		{
			return false;
		}
		else
		{
			for(unsigned int i = 0; i < extractMotif(a).size(); ++i) //if B and A are the same length, the lexicographically lower motif goes first (e.g. AAC precedes CAG)
			{
				if(extractMotif(a).at(i) < extractMotif(b).at(i))
				{
					return true;
				}
				if(extractMotif(a).at(i) > extractMotif(b).at(i))
				{
					return false;
				}
			}		
		}
		return true;
	}
	return false; // At this point, A and B have identical units and start positions, so they are the same microsatellite
}
unsigned int getCoordInd(const string & coord,unsigned int maxAutosome)
{
	string noChr = coord;
	if(coord.find("chr") != string::npos)
	{
		noChr = coord.substr(3);
	}
	
	if(noChr == "X")
	{
		return maxAutosome + 1;
	}
	else if(noChr == "Y")
	{
		return maxAutosome + 2;
	}
	else if(noChr == "MT")
	{
		return maxAutosome + 3;
	}

	return atoi(noChr.c_str());
}

unsigned int getCoordInd(const string & coord)
{
	string noChr = coord;
	if(coord.find("chr") != string::npos)
	{
		noChr = coord.substr(3);
	}
	
	if(noChr == "X")
	{
		return 23;
	}
	else if(noChr == "Y")
	{
		return 24;
	}
	else if(noChr == "MT")
	{
		return 25;
	}

	return atoi(noChr.c_str());
}

void reconstructSequence(const BamAlignment & al, const vector< MicrosatelliteQualityFlag > & qualFlags, string & reconstructedSequence)
{
	reconstructedSequence.clear();
	map< unsigned int, MicrosatelliteQualityFlag > msDetails; //by making the microsatellites a map where the read start of the microsatellite
	//is the key, we guarantee that microsatellites in the read are sorted by their start positions
	qualVectorToMap(qualFlags,msDetails,al);

	unsigned int lastTruncatedStop = 0, lastFullStop = 0;
	unsigned int querySequenceIntervalStop;
	string intervalSequence;

	MicrosatelliteQualityFlag instance;
	map<unsigned int, MicrosatelliteQualityFlag>::iterator it;
	tr(msDetails,it)
	{
		instance = (*it).second;
		if(instance.readStart() > lastFullStop)
		{
			querySequenceIntervalStop = instance.readStart() - lastFullStop - 1;
			intervalSequence = al.QueryBases.substr(lastTruncatedStop,querySequenceIntervalStop);
			upperCase(intervalSequence);
			reconstructedSequence.append(intervalSequence);
			
			appendMicrosatellite(reconstructedSequence,instance);
			lastTruncatedStop += instance.readStart() - lastFullStop - 1;
		}
		else
		{
			appendOverlappedMicrosatellite(reconstructedSequence,instance,lastFullStop);
		}
		lastFullStop = instance.readStop();
	}
	string remainingSequence = al.QueryBases.substr(lastTruncatedStop);
	upperCase(remainingSequence);
	reconstructedSequence.append(remainingSequence);
}

void qualVectorToMap(const vector< MicrosatelliteQualityFlag > & qualVector, map< unsigned int, MicrosatelliteQualityFlag > & qualMap, const BamAlignment & al)
{
	unsigned int readStart;
	vector< MicrosatelliteQualityFlag >::const_iterator it;
	tr(qualVector, it)
	{
		readStart = (*it).readStart();
		if(qualMap.find(readStart) != qualMap.end()) //if genomeStart is already in the map, two microsatellites have the same start position in the read
		{
			cerr << "Error!  Two reads in " << al.Name << " have the same start position in a read" << endl;
		}
		else
		{
			qualMap[readStart] = (*it);
		}
	}
}

void fillMicrosatelliteQualityDetails(const BamAlignment & al, vector< MicrosatelliteQualityFlag > & msDetails)
{
	string microsatelliteQualityTag;
	al.GetTag("MU",microsatelliteQualityTag);

	if(microsatelliteQualityTag.length() > 0)
	{
		vector< string > microsatelliteQualities;
		splitString(microsatelliteQualityTag,"|",microsatelliteQualities);

		vector< string >::iterator it;
		tr(microsatelliteQualities,it)
		{
			msDetails.push_back(MicrosatelliteQualityFlag(*it));
		}
	}
}

void appendMicrosatellite(string & sequence,const MicrosatelliteQualityFlag & details)
{
	string motif;
	details.getMotif(motif);
	upperCase(motif);
	
	for(unsigned int i = 0; i < details.length()/details.motifLength(); ++i) //appends complete motifs to the sequence
	{
		sequence.append(motif);
	}

	//if the microsatellite terminates in an incomplete microsatellite (e.g. ACACACACACACA), add the trailing microsatellite to the sequence
	string partialMotifSequence = motif.substr(0,details.length() % details.motifLength());
	sequence.append(partialMotifSequence);
}

void appendOverlappedMicrosatellite(string & sequence, const MicrosatelliteQualityFlag details, unsigned int lastFullStop)
{
	unsigned int overlappingMicrosatelliteBases = lastFullStop - details.readStart() + 1;
	unsigned int printedUnits = (int) ((lastFullStop - details.readStart() + details.motifLength())/details.motifLength());
	unsigned int incompleteUnits = overlappingMicrosatelliteBases % details.motifLength(); //will always return a positive number or 0

	string motif;
	details.getMotif(motif);
	upperCase(motif);

	//if the first few units of a microsatellite were overlapped, but not completed, by a preceding microsatellite, fill in the remainder of any shared units
	//(e.g. AACAACAACATATATATATAT)
	if(incompleteUnits != 0) 
	{
		sequence.append(motif.substr(incompleteUnits));
	}
	
	//append remaining complete units
	for(unsigned int i = 0; i < (details.length() - (details.motifLength() * printedUnits)) / details.motifLength(); ++i)
	{
		sequence.append(motif);
	}
	string partialMotifSequence = motif.substr(0,details.length() % details.motifLength());
	sequence.append(partialMotifSequence);
}

void fillMicrosatelliteQualityDetails(const BamAlignment & al, map< unsigned int, MicrosatelliteQualityFlag > & msDetails)
{
	string microsatelliteQualityTag;
	al.GetTag("MU",microsatelliteQualityTag);

	if(microsatelliteQualityTag.length() > 0)
	{
		vector< string > microsatelliteQualities;
		splitString(microsatelliteQualityTag,"|",microsatelliteQualities);

		vector< string >::iterator it;
		tr(microsatelliteQualities,it)
		{
			MicrosatelliteQualityFlag dummy(*it);

			if(msDetails.find(dummy.readStart()) != msDetails.end())
			{
				cerr << "Sequence reconstruction error: two microsatellites start at the same point in the read" << endl;
			}

			msDetails[dummy.readStart()] = dummy;
		}
	}
}


void fillMicrosatelliteCoordDetails(const BamAlignment & al, vector< MicrosatelliteCoordFlag > & msDetails)
{
	string microsatelliteCoordTag;
	al.GetTag("MC",microsatelliteCoordTag);

	if(microsatelliteCoordTag.length() > 0)
	{
		vector< string > microsatelliteCoords;
		splitString(microsatelliteCoordTag,"|",microsatelliteCoords);
		vector< string >::iterator it;
		tr(microsatelliteCoords,it)
		{
			msDetails.push_back(MicrosatelliteCoordFlag(*it));
		}
	}
}
