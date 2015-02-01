CXX = g++
CXXFLAGS = -Wall -pedantic -O3
CXXFLAGS_BAMTOOLS = -Wall -pedantic -Wno-long-long -O3
PROFFLAGS = -Wall -pedantic -pg

BAMTOOLS_INCLUDE = -I ./bamtools/include/
BAMTOOLS_LIB = -L ./bamtools/lib -lbamtools -lz
LZ_LIB = -lz

DEFAULT_LIBS = FileOps.cc Timer.cc Microsatellite.cc HelperFunc.cc Counter.cc StringManip.cc OrphanFunctions.cc
TETRASCAN_LIBS = FastqFile.cc FastaFile.cc SolexaRead.cc MicrosatelliteDetector.cc uScanFE.cc
REINDEX_BAM_LIBS = BamMicrosatelliteInfo.cc OsiFile.cc MsHeader.cc MsHeaderFile.cc
PROFILER_LIBS = MicrosatelliteDatabase.cc MicrosatelliteQualityFlag.cc MicrosatelliteCoordFlag.cc MicrosatelliteProfile.cc
SUMMARIZER_LIBS = ProfileSipper.cc ProfileLocus.cc PopLocus.cc PopSipper.cc

DB_TEST_LIBS = MicrosatelliteDatabase.cc
SM_TEST_LIBS = StringManip.cc

all: tetrascan profiler reindexBam summarizer

tetrascan: tetrascan.cc
	$(CXX) $(CXXFLAGS_NOPEDANT) $(DEFAULT_LIBS) $(TETRASCAN_LIBS) $(BAMTOOLS_INCLUDE) tetrascan.cc $(BAMTOOLS_LIB) -o tetrascan
reindexBam: reindexBam.cc
	$(CXX) $(CXXFLAGS_BAMTOOLS) $(DEFAULT_LIBS) $(REINDEX_BAM_LIBS) $(BAMTOOLS_INCLUDE) reindexBam.cc $(BAMTOOLS_LIB) -o reindexBam
profiler: profiler.cc
	$(CXX) $(CXXFLAGS_BAMTOOLS) $(DEFAULT_LIBS) $(PROFILER_LIBS) $(BAMTOOLS_INCLUDE) profiler.cc $(BAMTOOLS_LIB) -o profiler
summarizer: uSeqSummarizePop.cc
	$(CXX) $(CXXFLAGS) $(DEFAULT_LIBS) $(SUMMARIZER_LIBS) uSeqSummarizePop.cc $(LZ_LIB) -o $@

dbTest: microsatelliteDatabaseTest.cc
	$(CXX) $(CXXFLAGS) $(DEFAULT_LIBS) $(DB_TEST_LIBS) microsatelliteDatabaseTest.cc -o $@
stringManipTest: stringManipTest.cc
	$(CXX) $(CXXFLAGS) $(SM_TEST_LIBS) stringManipTest.cc -o $@

clean:
	rm tetrascan reindexBam profiler

summarizer: ProfileSipper.h ProfileSipper.cc ProfileLocus.cc ProfileLocus.h PopSipper.cc PopSipper.h PopLocus.cc PopLocus.h
profiler: MicrosatelliteDatabase.cc MicrosatelliteDatabase.h Counter.cc Counter.h MicrosatelliteQualityFlag.cc MicrosatelliteQualityFlag.h MicrosatelliteCoordFlag.cc MicrosatelliteCoordFlag.h UsefulMacros.h
profiler: MicrosatelliteProfile.cc MicrosatelliteProfile.h profiler.cc
tetrascan: uScanFE.cc uScanFE.h
tetrascan reindexBam summarizer: Timer.cc Timer.h
tetrascan reindexBam summarizer: FileOps.cc FileOps.h
tetrascan reindexBam summarizer: Microsatellite.cc Microsatellite.h
tetrascan reindexBam summarizer: OrphanFunctions.cc OrphanFunctions.h
tetrascan: FastqFile.cc FastqFile.h
tetrascan: FastaFile.cc FastaFile.h
tetrascan: MicrosatelliteDetector.cc MicrosatelliteDetector.h
tetrascan: SolexaRead.cc SolexaRead.h
summarizer: HelperFunc.h HelperFunc.cc
reindexBam: MsHeader.cc MsHeader.h
reindexBam: MsHeaderFile.cc MsHeaderFile.h
reindexBam: OsiFile.cc OsiFile.h
reindexBam: BamMicrosatelliteInfo.cc BamMicrosatelliteInfo.h
tetrascan: Counter.cc Counter.h
summarizer: StringManip.cc StringManip.h
summarizer: OrphanFunctions.h OrphanFunctions.cc


# addMergedBamFileReadGroups: AddMergedBamFileReadGroups.cc
#	$(CXX) $(CXXFLAGS_BAMTOOLS) $(DEFAULT_LIBS) $(BAMTOOLS_INCLUDE) AddMergedBamFileReadGroups.cc $(BAMTOOLS_LIB) -o $@
# split: splitByChr.cc
#	$(CXX) $(CXXFLAGS) $(DEFAULT_LIBS) $(SPLIT_LIBS) splitByChr.cc -o splitByChr
# simulate: MicroSim.cc
#	$(CXX) $(CXXFLAGS) $(DEFAULT_LIBS) $(SIMULATOR_LIBS) $(RNG_LIBS) MicroSim.cc -o MicroSim
# checkSim: checkSim.cc
#	$(CXX) $(CXXFLAGS) $(DEFAULT_LIBS) $(REINDEX_LIBS) checkSim.cc -o checkSim
# sort: SortSam.cc
#	$(CXX) $(CXXFLAGS) $(DEFAULT_LIBS) $(SORT_LIBS) SortSam.cc -o SortSam
# switch: switchDelims.cc
#	$(CXX) $(CXXFLAGS) $(DEFAULT_LIBS) $(SWITCH_LIBS) switchDelims.cc -o switch
# uscan: tetrascan.cc
#	$(CXX) $(CXXFLAGS_NOPEDANT) $(DEFAULT_LIBS) $(TETRASCAN_LIBS) $(BAMTOOLS_INCLUDE) tetrascan.cc $(BAMTOOLS_LIB) -o uScan
# profilerT: profiler.cc
#	$(CXX) $(CXXFLAGS_BAMTOOLS) $(DEFAULT_LIBS) $(PROFILER_LIBS) $(BAMTOOLS_INCLUDE) profiler.cc $(BAMTOOLS_LIB) -o profilerT
# reindex: reindexSam.cc
#	$(CXX) $(CXXFLAGS) $(DEFAULT_LIBS) $(REINDEX_LIBS) reindexSam.cc -o reindexSam
# tetrascan_prof: tetrascan.cc
#	$(CXX) $(PROFFLAGS) $(DEFAULT_LIBS) $(TETRASCAN_LIBS) $(BAMTOOLS_INCLUDE) tetrascan.cc $(BAMTOOLS_LIB) -o tetrascan_prof
# translate checkSim reindex split sort switch: SamRead.h SamRead.cc
# translate checkSim reindex split sort switch: SamFile.h SamFile.cc
# filter: MicrosatelliteFilter.cc MicrosatelliteFilter.h
# dbTest: MicrosatelliteDatabase.cc MicrosatelliteDatabase.h
# SPLIT_LIBS = SamFile.cc SamRead.cc
# SIMULATOR_LIBS = MicrosatelliteDatabase.cc Mutator.cc Genome.cc SolexaRead.cc FastaFile.cc
# FILTER_LIBS = SolexaRead.cc MicrosatelliteFilter.cc FastqFile.cc

RNG_LIBS = RandomNumber.cc RandomDistribution.cc
SORT_LIBS = SamFile.cc SamRead.cc
SWITCH_LIBS = SamFile.cc SamRead.cc
REINDEX_LIBS = SamRead.cc SamFile.cc OsiFile.cc MsHeaderFile.cc MsHeader.cc
