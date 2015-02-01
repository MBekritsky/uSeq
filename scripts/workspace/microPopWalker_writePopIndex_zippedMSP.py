'''
Created on Jun 20, 2011

@author: levy
'''

import numpy as np
import os
from collections import defaultdict
from collections import Counter
import time
import subprocess
import signal
import gzip

#A function to go out to the browser
def goldenFox(chrom, start, end):
    posString = "%s:%i-%i" % (chrom, start, end)
    url = """http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=%s""" % posString
    command = "firefox \"%s\"" % url
    os.system(command)

# Class to store a single locus from a single individual
class microLocus:
    # constructor builds from the first string
    # additional flanks are added in the "add entry" method
    def __init__(self, firstString,chr2int,sampleID = None):
        self.sampleID = sampleID
        entry = firstString.split("\t")
        self.chrom      = chr2int[entry[0]]
        self.pos        = int(entry[1])
        self.unit       = entry[2]
        self.unitLength = len(self.unit)
        self.refLength  = int(entry[3])
        self.numFlanks  = int(entry[4])
        self.totalCount = int(entry[7])

        self.flankSeq   = []
        self.flankHist  = []
        self.flankCount  = []

        self.flankSeq.append(entry[6])
        self.flankCount.append(int(entry[7]))
        self.flankHist.append(map(int, entry[8].split(",")))
        self.maxLength = len(self.flankHist[0])
            
    # adds the next entry
    # loads the flank information and histogram
    def addEntry(self, nextString):
        entry = nextString.split("\t")
        flank = entry[6]
        flank = "%s  %s" %(flank[:5], flank[5:])
        self.flankSeq.append(flank)
        self.flankCount.append(int(entry[7]))
        self.flankHist.append(map(int, entry[8].split(",")))
        
    # returns the histogram for a given flank index
    def getArray(self, index):
        return np.array(self.hist[index])      
        
    # print locus as a string
    def __str__(self):
        return self.toString()
    
    # prints the string depending on the "population size" of the flank -- 
    # i.e. how many reads are in the flank
    def toString(self, minPopSize=0):
        lines = []
        titleString = "\t".join(map(str, [self.chrom, self.pos, self.unit, self.refLength, 
                                          self.numFlanks, self.totalCount]))
        
        if self.sampleID != None:
            titleString = "\t".join([self.sampleID, titleString])
        
        lines.append(titleString)
        for i in xrange(len(self.flankSeq)):
            if self.flankCount[i] >= minPopSize: 
                histString = "\t".join([self.flankSeq[i]] + map(str, self.flankHist[i]))
                if i == 0:
                    histString = "\t" + histString
                lines.append(histString)        
        return "\n".join(lines)
        

# a walker for a single individual
# meaning it takes as its input a gzipped msp file
# and creates an iterator that returns one locus at a time
class microWalker:
    def __init__(self, filename, sampleID = None):
        self.infile = gzip.open(filename,'rb')
#        self.infile = subprocess.Popen(("gunzip","-c",filename),stdout=subprocess.PIPE,
#                                       preexec_fn=lambda: signal.signal(signal.SIGPIPE,signal.SIG_DFL))
        self.sampleID = sampleID
        try:
            self.nextLine = self.infile.readline().strip()
        except:
            print 'empty file %s' % filename
        
    def loadNext(self,chr2int):
        try:
            nextLoci = []
            nextLocus = microLocus(self.nextLine, chr2int, self.sampleID)
            for i in xrange(nextLocus.numFlanks):
                nextLocus.addEntry(self.infile.readline().strip())
            nextLoci.append(nextLocus)
            
            locusChrom  = nextLocus.chrom
            locusPos    = nextLocus.pos
            
            self.nextLine = self.infile.readline().strip()
            entry = self.nextLine.split("\t")
            nextChrom = chr2int[entry[0]]
            nextPos   = int(entry[1])
            while (nextChrom == locusChrom) and (locusPos == nextPos):
                nextLocus = microLocus(self.nextLine, self.sampleID)
                for i in xrange(nextLocus.numFlanks):
                    nextLocus.addEntry(self.infile.readline().strip())
                nextLoci.append(nextLocus)    
                self.nextLine = self.infile.readline().strip()
                entry = self.nextLine.split("\t")
                nextChrom = chr2int[entry[0]]
                nextPos   = int(entry[1])            
            return nextLoci
            
        except:
            if len(nextLoci) > 0:
                return nextLoci
            else:
                return None
    
    def __iter__(self):
        return self
    
    def next(self,chr2int):
        nextEntry = self.loadNext(chr2int)
        if nextEntry == None:
            self.infile.close()
            raise StopIteration
        else:
            return nextEntry
 

# a walker that iterates across the whole population
# it makes use of microWalkers and keeps them in step so that only a single locus
# is reported at a time.
# presumption is that all input msp files are sorted by chrom and pos
class microPopWalker:
    def __init__(self, microWalkers,chr2int):
        self.microWalkers = microWalkers
        self.length = len(microWalkers)
        self.loci     = []
        self.chr2int = chr2int
        
        for i in xrange(self.length):
            self.loci.append(None)
            
        self.active      = np.ones(self.length, dtype = bool)
        self.chrom       = np.zeros(self.length, dtype=int)
        self.pos         = np.zeros(self.length, dtype=int)
        self.nextEntry   = np.ones(self.length, dtype=bool)
        self.update(self.nextEntry)
        

    # called after getNext(), this function moves up the individual walkers
    # that contributed to the last locus.
    def update(self, toUpdate):
        index = 0
        for walker, toUpdate, act in zip(self.microWalkers, toUpdate, self.active):
            if act and toUpdate:
                try:
                    locus = walker.next(self.chr2int)
                    self.loci[index] = locus
                    self.chrom[index] = locus[0].chrom
                    self.pos[index] = locus[0].pos
                except:
                    self.active[index] = False
                    self.loci[index] = None
                    self.chrom[index] = 0
                    self.pos[index] = 0
            index += 1
        self.nextChrom = min(self.chrom[self.active])
        chromMatch    = self.chrom == self.nextChrom
        self.nextPos  = min(self.pos[chromMatch])
        chromPosMatch = np.logical_and(chromMatch, self.pos == self.nextPos)

        self.nextEntry = chromPosMatch
        self.sizeOfNext = np.sum(self.nextEntry)


    ## returns the next locus
    def getNext(self):
        ans, whoHas = [], []
        for locus, hasNext in zip(self.loci, self.nextEntry):
            if hasNext:
                ans += locus
                whoHas.append(True)
            else:
                whoHas.append(False)
        return ans, whoHas
    
    def __iter__(self):
        return self
  
    ## return next locus and advance the relevant walkers  
    def next(self):
        if sum(self.nextEntry) == 0:
            raise StopIteration
        else:
            loci, whoHas = self.getNext()
            info = [self.nextChrom, self.nextPos, self.sizeOfNext]
            self.update(self.nextEntry)
            return popLocus(loci)
            

# keeps all the relevant information about a particular locus across a population
# with methods for splitting by unit/flank/etc.
class popLocus:
    def __init__(self, loci):
        self.unitIndex = defaultdict(list)
        
        for index, locus in zip(range(len(loci)), loci):
            self.unitIndex[locus.unit].append(index)
        
        self.loci = loci
        firstLocus  = loci[0]
#        print firstLocus.toString()
        self.chrom  = firstLocus.chrom
        self.pos    = firstLocus.pos
        self.unit   = self.unitIndex.keys()
        self.refLength = firstLocus.refLength
        self.whoHas = np.unique([x.sampleID for x in loci])
        self.count  = [x.totalCount for x in loci]
        self.maxLength = max([x.maxLength for x in loci])

        
    ## returns a list of popLoci where each locus has a single unit.
    def splitByUnit(self):
        ans = []
        for unit in self.unit:
            ans.append(popLocus([self.loci[index] for index in self.unitIndex[unit]]))
        return ans
    
    ## returns a dictionary (of lists).
    ## the keys are the flanks
    ## the values are a list of pairs: (sampleID, flankHistogram)
    def splitByFlank(self):
        ans = defaultdict(list)
        for locus in self.loci:
            for flank, histogram in zip(locus.flankSeq, locus.flankHist):
                ans[flank].append((locus.sampleID, histogram))      
        return ans
    
    ##
    ## returns the histogram for "all" for the whole population
    ## in an array with shape (popSize, maxLengthObserved)
    ##
    def getAllHistogram(self, popIndex, popSize):
        allHist = np.zeros(shape=(popSize, self.maxLength), dtype=int)
        for locus in self.loci:
            sample = locus.sampleID
            hist = locus.flankHist[0]
            allHist[popIndex[sample], :len(hist)] = hist
        return allHist
        
    
    ## returns a dictionary where the keys are flanks
    ## and the values are the histogram array over the population
    ## if completeOnly=True, returns only flanks with no "#"s
    def flankHistograms(self, popIndex, popSize, completeOnly=False):
        allFlanks = self.splitByFlank()
        ans = {}
        for flank in allFlanks:
            if (not completeOnly) or (flank.find("#") < 0):
                thisHist = np.zeros(shape=(popSize, self.maxLength), dtype=int)
                for sample, hist in allFlanks[flank]:
                    thisHist[popIndex[sample], :len(hist)] = hist
                ans[flank] = thisHist
        return ans
                    
    def __str__(self):
        return "\t".join(map(str, [self.chrom, self.pos, self.unit, self.refLength]))
        
    ## makes the histogram plot
##    def plotHistograms(self, flankHistograms, minCount):
#        ## the first half of the code
#        ## goes through the flank histograms and determines number of reads and maximum data
#        ## so that histograms can be presented in "sorted" order (from most populous to least)
#        ## and to set the scale of the data
#        plotIndex = 1
#        flankSize = []
#        flankSeq  = []
#        dataMax   = []
#        for flank in flankHistograms:
#            hist = flankHistograms[flank]
#            flankSeq.append(flank)
#            flankSize.append(np.max(np.sum(hist, 1)))
#            dataMax.append(np.max(hist))
## loads position, referenece info, flank, etc
#            
#        flankSize = np.array(flankSize)
#        dataMax = np.array(dataMax)
#        
#        minCount = min(minCount, max(flankSize))
#        orderedList = np.argsort(-flankSize)
#        numPlots = sum(flankSize >= minCount)
#        dataMax = max(dataMax[flankSize >= minCount])
#        
#        fig = plt.figure(figsize=(25,12*numPlots))
#        plt.subplots_adjust(left=0.05, right=.90, bottom=0.02, top=0.98, wspace=0.05, hspace=0.05)
#        
#        ## for each non-trivial flank
#        for index in orderedList[:numPlots]:
#            flank = flankSeq[index]
#            popArray = flankHistograms[flank]
#            
#            ## this determines a dictionary for finding the base at a particular
#            ## length in the repeat (for indexing the x-axis by base)
#            unit = self.unit[0]
#            unitLength = len(unit)
#            unitDic = {}
#            for i in xrange(unitLength):
#                unitDic[i] = unit[i]
#
#            plotLength = max(self.refLength, self.maxLength)
#            xtickLabels = []
#            xtickLocation = []
#            for i in xrange(plotLength):
#                xtickLabels.append(unitDic[i % unitLength])
#                xtickLocation.append(i+1)
#            
#            ## make the plot
#            ax = fig.add_subplot(numPlots, 1, plotIndex)
#            data = np.hstack([np.zeros((len(popArray),1), dtype=int), popArray])
#            cax = ax.imshow(data, cmap="Blues", aspect='auto',interpolation='nearest', origin='upper', vmin=0, vmax=dataMax)
#            ax.plot([self.refLength-0.5, self.refLength-0.5], [-0.5, popSize-0.5], '--r', lw=3)
#            ax.plot([self.refLength+0.5, self.refLength+0.5], [-0.5, popSize-0.5], '--r', lw=3)
#    
#            ax.set_xlim(-0.5, plotLength+0.5)
#            ax.set_ylim(-0.5, popSize-0.5)
#            
#            ## set the x-axis ticks and the mark for reference length
#            plt.xticks(xtickLocation, xtickLabels, fontsize=15)
#            ax.text(self.refLength, -0.65, str((self.refLength)), color='r', fontsize=15)
#            ## for each rtp class, set the color
#            for i in xrange(popSize):
#                color = ['red', 'blue', 'green', 'orange'][i % 4]
#                for j in xrange(self.maxLength+1):
#                    value = data[i,j]
#                    alpha = min(1, float(value) / 10.)
#                    if value != 0:                        
#                        ax.text(j-0.1, i, str(data[i,j]), color=color, fontsize=10, alpha=alpha)
#            ax.set_title(flank)
#            plotIndex += 1
#
#        figtitle = "%s:%s %s %i" % (self.chrom, self.pos, unit, self.refLength)
#        fig.suptitle(figtitle, fontsize=18)
#        return fig


## some functions for acting on a population wide histogram
class popHistogram:    
    def __init__(self, histogram):
        self.histogram = histogram
        self.totalCount = np.sum(self.histogram, 1)
        
    def mainModes(self):
        argData = np.argsort(-self.histogram, 1)
        mode1  = argData[:, 0]
        mode2  = argData[:, 1]
        count1 = self.histogram[np.arange(len(self.histogram)), mode1 ]
        count2 = self.histogram[np.arange(len(self.histogram)), mode2 ]
        return mode1, mode2, count1, count2
    
    def getModes(self, percentCutoff, minReads=0):
        hist2 = np.array(self.histogram, dtype=float) / np.clip(self.totalCount[:, np.newaxis], 1, np.inf)
        locs = np.where(hist2 >= percentCutoff)

        modes = defaultdict(list)
        for x, y in zip(locs[0], locs[1]):
            if self.histogram[x,y] >= minReads:
                modes[y].append(x)
        return modes



def histSummary(hist):
    lengths = []
    counts = []
    for i in xrange(len(hist)):
        if hist[i] > 0:
            lengths.append(str(i+1))
            counts.append(str(hist[i]))
    return ";".join([str(len(lengths)), ",".join(lengths), ",".join(counts)])
    

## the main code for walking through the data.        

parentDir = "/mnt/wigclust10/data/safe/bekritsk/simons/Exome/workingSet/06182013_completeWiglerSSCQuads"
indexFilename  = "%s/microSat.refOnly.index.txt" % parentDir
countsFilename = "%s/microSat.refOnly.counts.txt" % parentDir
columnFilename = "%s/microSat.refOnly.columns.txt" % parentDir

indexFile  = file(indexFilename, 'w')
columnFile = file(columnFilename, 'w')
countsFile = file(countsFilename, 'w')
t0 = time.time()
t1 = time.time()

famIDs = []
for x in os.listdir(parentDir):
    if x[:5] == "auSSC":
        path = os.path.join(parentDir, x)
        if os.path.isdir(path):
            famIDs.append(x)

#famIDs = [famID]
#print famIDs

chr2int = {"chr1": 1,"chr2": 2,"chr3": 3,"chr4": 4,"chr5": 5,"chr6": 6,"chr7": 7,"chr8": 8,"chr9": 9,"chr10": 10,
           "chr11": 11,"chr12": 12,"chr13": 13,"chr14": 14,"chr15": 15,"chr16": 16,"chr17": 17,"chr18": 18,
           "chr19": 19,"chr20": 20,"chr21": 21,"chr22": 22,"chrX": 23,"chrY": 24,"chrMT": 25}
int2chr ={v:k for k,v in chr2int.items()}

sampleWalkers = []

numPeople = 0
sample2fam = {}
fam2sample = defaultdict(list)
count = 0
for famID in famIDs:
    famPath = os.path.join(parentDir, famID)
    for x in os.listdir(famPath):
        if x[:3] == "SSC":
            sampleID = x
            microGzipFilename = "%s/%s/%s/%s.msp.gz" % (parentDir, famID, sampleID, sampleID)
            relFilename = "%s/%s/%s/relationship.txt" % (parentDir, famID, sampleID)
            refFile = file(relFilename, 'r')
            
            rtp = refFile.next().strip()
            refFile.close()
            sample2fam[sampleID] = (famID, rtp)
            fam2sample[famID].append(sampleID)
            
            if os.path.exists(microGzipFilename):
                numPeople += 1
                sampleWalkers.append(microWalker(microGzipFilename, sampleID))
            else:
                print "%s/%s does not have a gzipped msp file" % (famID,sampleID)
         
popIndex = {}
popList = []
index = 0
for famID in fam2sample.keys():
    samples  = np.array(fam2sample[famID])
    relation = np.array([sample2fam[x][1] for x in samples])
    ordered  = np.array([samples[relation == 'mother'][0], 
                         samples[relation == 'father'][0], 
                         samples[relation == 'self'][0],
                         samples[relation == 'sibling'][0]])

    relation2 = np.array([sample2fam[x][1] for x in ordered])
    for x in ordered:
        popIndex[x] = index  
        popList.append(x)
        index += 1

for x in popList:
    columnFile.write("\t".join([x, sample2fam[x][0], sample2fam[x][1]]) + "\n")

columnFile.close()

popSize = len(popList)
walker = microPopWalker(sampleWalkers,chr2int)
print popSize
print 'walking'
counter = 0
plotCounter = 1
totalCounter = 0
skipped = 0
added = 0
no_ref = []
for pop_locus_raw in walker:
#    if totalCounter == 10000:
#        break
    totalCounter += 1
    if len(pop_locus_raw.unit) == 1:
        unitLoci = [pop_locus_raw]
    else:
        unitLoci = pop_locus_raw.splitByUnit()
    
    for pop_locus in unitLoci:
#        print pop_locus
#        print 'number of people', len(pop_locus.whoHas)
#        print 'maximum size', max(pop_locus.count)
#        print pop_locus.maxLength
        if pop_locus.refLength == 0:
            skipped += 1
            no_ref.append("{}\t{}".format(pop_locus.chrom, pop_locus.pos))
        if pop_locus.refLength != 0:
            added += 1
            indexLine = "\t".join(map(str, [int2chr[pop_locus.chrom], pop_locus.pos, pop_locus.unit[0], pop_locus.refLength,len(pop_locus.whoHas), max(pop_locus.count), np.sum(pop_locus.count)] ))
            ans = defaultdict(list)
            for locus in pop_locus.loci:
                ans[popIndex[locus.sampleID]] = histSummary(locus.flankHist[0])
                
            words = [ans[i] if i in ans else '0;;' for i in range(popSize)]
            countsLine = "\t".join(words)
            indexFile.write(indexLine + "\n")
            countsFile.write(countsLine + "\n")
            counter += 1
            
            if counter % 10000 == 0:
                print pop_locus.chrom
                print "\t".join(map(str, [counter, indexLine, "%0.2f min/10K" % ((time.time() - t1)/60), "%0.2f min total" % ((time.time() - t0)/60)]))
                print "Scanned %d total loci" % totalCounter
                t1 = time.time()
indexFile.close()
countsFile.close()
skippedFilename = "%s/skippedLociFromTable.txt" % parentDir
skippedFile = file(skippedFilename,"w")
skippedFile.write("\n".join(no_ref))

print totalCounter
print skipped
print added
