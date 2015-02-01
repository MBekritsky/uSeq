#!/usr/bin/Rscript

get.invoked.dir<-function()
{
  invoked.file<-""
  #when running this script on SGE, file will have some uninformative script name, so invoked.file must be overridden with script_file
  if(sum(grepl("script.file",commandArgs(),fixed=TRUE)) > 0)
  {
    cat("Script file specified\n",commandArgs(),"\n")
    cat(grep("script.file=",commandArgs(),value=TRUE),"\n")
    invoked.file<-unlist(strsplit(grep("script.file=",commandArgs(),value=TRUE),"=",fixed=TRUE))[2]
    cat(invoked.file,"\n")
  } else
  {
    invoked.file<-unlist(strsplit(grep("file=",commandArgs(),value=TRUE),"=",fixed=TRUE))[2]
  }
  cat(invoked.file,"\n")
  invoked.parsed.dir<-unlist(strsplit(invoked.file,"/",fixed=TRUE))
  invoked.dir<-sprintf("%s/",paste(invoked.parsed.dir[1:(length(invoked.parsed.dir) - 1)],sep="/",collapse="/"))
  return(invoked.dir)
}

source.genotyper.functions<-function(src.dir,trace=TRUE,...)
{
  for(src.file in list.files(src.dir,full.names=TRUE,pattern=".R$"))
  {
    if(trace) cat(src.file, ":")
    source(src.file)
    if(trace) cat("\n")
  }
}

check.for.local.R.dir<-function()
{
  local.R.dir<-sprintf("~/R/%s-library/%d.%d",version$platform,as.integer(version$major),as.integer(version$minor) %/% 1)
  if(!file.exists(local.R.dir)) dir.create(local.R.dir,recursive=TRUE)
  cat("Created local R directory",local.R.dir,"\n")
  
  if(!(local.R.dir %in% .libPaths())) .libPaths(local.R.dir)
  cat("Added",local.R.dir,"to this R session's library tree\nLibrary paths are now ",paste(.libPaths(),collapse=", "),"\n")
}

invoked.dir<-get.invoked.dir()
config.file<-sprintf("%sconfig.txt",invoked.dir)
config.vars<-read.table(config.file,colClasses="character")

src.dir<-config.vars[1,2]
anno.dir<-config.vars[2,2]
cran.repos<-config.vars[3,2]
source.genotyper.functions(src.dir,trace=FALSE)
cat("Loaded all source files from",src.dir,"\n")

check.for.local.R.dir()

check.package.installation("getopt",cran.repos=cran.repos,quietly=TRUE)
check.package.installation("fields",cran.repos=cran.repos,quietly=TRUE)
check.package.installation("colorRamps",cran.repos=cran.repos,quietly=TRUE)

genotyper.specs<-matrix(c(
  'directory'       ,'d',1,"character",
  "script.file"     ,'s',2,"character",
  'help'            ,'h',0,"logical",
  'verbose'         ,'v',2,"integer"
),byrow=TRUE,ncol=4)

opt<-getopt(genotyper.specs)

if(!is.null(opt$help))
{
  cat(getopt(genotyper.specs,usage=TRUE))
  q(status=1)
}

if(is.null(opt$directory))
{
  cat("You must specify a directory with the --directory flag\n",getopt(genotyper.specs,usage=TRUE))
  q(status=1)
}

opt$fig.directory<-sprintf("%s/SummaryFigures/",opt$directory)
dir.create(opt$fig.directory,showWarnings=FALSE) #if directory exists, dir.create does nothing

opt$data.directory<-sprintf("%s/SummaryData/",opt$directory)

if(is.null(opt$verbose)) opt$verbose<-1

opt$start.locus<-1
opt$num.loci<-1e4

locus.info.file=sprintf("%supdatedLocusInfo.RData",opt$directory)
allele.info.file=sprintf("%sallele_matrix.txt",opt$directory)
family.info.file=sprintf("%sperson_info.txt",opt$directory)
report.file=sprintf("%sreport_sampleSet_20130630.txt",opt$directory)
coverage.estimator.file=sprintf("%scoverageEstimator.RData",opt$directory)

#load locus info
load(locus.info.file)


#load person info
person.info<-load.person.info(family.info.file)
person.info<-add.genders(person.info,report.file)
if(opt$verbose >= 1) cat("Loaded person information\n")

#find families where anyone has DNA sequence not derived from whole blood (specific to the SSC project)
ss.report<-load.sample.set.report(report.file)
non.whole.blood.indices<-get.non.whole.blood.family.members(ss.report,person.info,verbose=opt$verbose)

load(coverage.estimator.file)
if(opt$verbose >= 1) cat("Loaded coverage estimator\n")

#exclude people with low coverage, poor correlation, or with DNA not derived from whole blood
if(opt$verbose >= 1) cat("The following people have very low coverage and their families are being excluded from analysis:",
                         paste(as.character(person.info$individual.id[coverage.estimator$low.cov]),collapse=", "),"\n")
if(opt$verbose >= 1) cat("The following people have poor correlation to their expected coverage and their families are being excluded from analysis:",
                         paste(as.character(person.info$individual.id[coverage.estimator$bad.cor]),collapse=", "),"\n")
excluded.by.cm<-unique(c(as.character(person.info$family.id[coverage.estimator$low.cov]),as.character(person.info$family.id[coverage.estimator$bad.cor])))
inds.excluded.by.cm<-which(person.info$family.id %in% excluded.by.cm)
all.excluded.people<-unique(c(inds.excluded.by.cm,non.whole.blood.indices$inds))

cat("Excluding the following families from analysis:",
    paste(unique(as.character(person.info$family.id[all.excluded.people])),collapse=", ",sep=", "),"\n")
person.info<-exclude.people.from.person.df(person.info,all.excluded.people)

pop<-nrow(person.info)
n.fams<-pop/4

#annotate exome loci
exon.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.exons.merged.bed",anno.dir)
intron.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.introns.merged.bed",anno.dir)
mirna.anno.file<-sprintf("%s/hg19.ms.miRNA.merged.bed",anno.dir)
utr.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.utr.merged.bed",anno.dir)
gene.id.file<-sprintf("%s/hg19GeneName.ccds.rg.geneID.txt",anno.dir)

gene.ids<-load.gene.ids(gene.id.file)
annotations<-get.annotations(exon.file=exon.anno.file,intron.file=intron.anno.file,mirna.file=mirna.anno.file,utr.file=utr.anno.file,
                        	 locus.info=updated.locus.info,verbose=opt$verbose)
load(sprintf("%slocusCoverage.RData",opt$data.directory))

y.step<-0.025

cat("Plotting coverage info\n")

#png(sprintf("%scoverage.png",opt$fig.directory),height=1000,width=1000)
#par(cex=2)
pdf(sprintf("%scoverage.pdf",opt$fig.directory),height=8,width=8)
h.break.step<-1
h.break.max<-max(coverage)
h.disp.step<-1
h.disp.max<-ifelse(h.break.max > 100, 100, h.break.max)

h<-hist(coverage,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of all genotypes",xlab="Coverage",main=sprintf("Coverage histogram\nmedian: %d",median(coverage)),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,y.max,0.01),lwd=4,col="white")
abline(v=seq(0,h.disp.max + 10,h.break.step),lwd=2,col="white")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,length.out=6),labels=seq(0,y.max*100,length.out=6),las=2)
garbage<-dev.off()

h.break.step<-1
exon.break.max<-max(coverage[which(!is.na(annotations$exon)),])
exon.hist<-hist(coverage[which(!is.na(annotations$exon)),],breaks=seq(0,exon.break.max,h.break.step),plot=FALSE)
exon.hist$counts<-exon.hist$counts/sum(exon.hist$counts)

non.exon.break.max<-max(coverage[which(is.na(annotations$exon)),])
non.exon.hist<-hist(coverage[which(is.na(annotations$exon)),],breaks=seq(0,non.exon.break.max,h.break.step),plot=FALSE)
non.exon.hist$counts<-non.exon.hist$counts/sum(non.exon.hist$counts)

h.disp.step<-1
h.disp.max<-ifelse(max(c(exon.break.max,non.exon.break.max)) > 100, 100, h.break.max)

y.max<-((max(c(exon.hist$counts,non.exon.hist$counts)) %/% y.step) * y.step) + y.step

# png(sprintf("%sexonCoverage.png",opt$fig.directory),height=1800,width=1800,res=300)
pdf(sprintf("%sexonCoverage.pdf",opt$fig.directory),height=8,width=8)
plot(exon.hist,border=NA,col="grey50",ylab="% of exon loci in SSC individuals",xlab="Coverage",main=sprintf("Exon coverage histogram\nmedian: %d",median(coverage[which(!is.na(annotations$exon)),])),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,y.max,0.01),lwd=4,col="white")
abline(v=seq(0,h.disp.max + 10,h.break.step),lwd=2,col="white")
axis(side=1,at=seq(0,h.disp.max,length.out=6))
axis(side=2,at=seq(0,y.max,length.out=6),labels=seq(0,y.max*100,length.out=6),las=2)
lines(x=non.exon.hist$mids,y=non.exon.hist$counts,col="red",lwd=5)
legend(x="topright",col=c(rgb(0,0,0,0.5),"red"),lwd=10,legend=c("Exon","Non-exon"))
garbage<-dev.off()

# png(sprintf("%snonExonCoverage.png",opt$fig.directory),height=1800,width=1800,res=300)
pdf(sprintf("%snonExonCoverage.pdf",opt$fig.directory),height=8,width=8)

plot(non.exon.hist,border=NA,col="grey50",ylab="% of non-exon loci in SSC individuals",xlab="Coverage",main=sprintf("Non-exon coverage histogram\nmedian: %d",median(coverage[which(is.na(annotations$exon)),])),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,y.max,0.01),lwd=4,col="white")
abline(v=seq(0,h.disp.max + 10,h.break.step),lwd=2,col="white")
axis(side=1,at=seq(0,h.disp.max,length.out=6))
axis(side=2,at=seq(0,y.max,length.out=6),labels=seq(0,y.max,length.out=6),las=2)
lines(x=exon.hist$mids,y=exon.hist$counts,col="red",lwd=5)
legend(x="topright",col=c(rgb(0,0,0,0.5),"red"),lwd=10,legend=c("Non-exon","Exon"))
garbage<-dev.off()

png(sprintf("%spersonMedianCoverage.png",opt$fig.directory),height=1000,width=1000)
par(cex=2)
y.step<-0.05
h.break.step<-1
h.break.max<-max(apply(coverage,2,median))
h.disp.step<-1
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(apply(coverage,2,median),breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of SSC individuals",xlab="Coverage",main=sprintf("Median coverage per person histogram\nmedian: %d",median(apply(coverage,2,median))),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,length.out=6),labels=seq(0,y.max*100,length.out=6),las=2)
garbage<-dev.off()

png(sprintf("%spersonMedianExonCoverage.png",opt$fig.directory),height=1000,width=1000)
par(cex=2)
y.step<-0.05
h.break.step<-1
temp<-apply(coverage[which(!is.na(annotations$exon)),],2,median)
h.break.max<-max(temp)
h.disp.step<-1
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(temp,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of SSC individuals",xlab="Coverage",main=sprintf("Median coverage per person histogram in exons\nmedian: %d",median(temp)),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,length.out=6))
garbage<-dev.off()

png(sprintf("%spersonMedianNonExonCoverage.png",opt$fig.directory),height=1000,width=1000)
par(cex=2)
y.step<-0.05
h.break.step<-1
temp<-apply(coverage[which(is.na(annotations$exon)),],2,median)
h.break.max<-max(temp)
h.disp.step<-1
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(temp,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of SSC individuals",xlab="Coverage",main=sprintf("Median coverage per person histogram in non-exons\nmedian: %d",median(temp)),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,length.out=6),labels=seq(0,y.max*100,length.out=6),las=2)
garbage<-dev.off()

png(sprintf("%slocusMedianCoverage.png",opt$fig.directory),height=1000,width=1000)
par(cex=2)

h.break.step<-1
h.break.max<-ceiling(max(apply(coverage,1,median)))
h.disp.step<-1
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(apply(coverage,1,median),breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of well-covered loci",xlab="Coverage",main=sprintf("Median coverage per locus histogram\nmedian: %d",median(apply(coverage,1,median))),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step),labels = seq(0,y.max*100,y.step*100),las=2)
garbage<-dev.off()

png(sprintf("%slocusExonMedianCoverage.png",opt$fig.directory),height=1000,width=1000)
par(cex=2)

h.break.step<-1
temp<-apply(coverage[which(!is.na(annotations$exon)),],1,median)
h.break.max<-ceiling(max(temp))
h.disp.step<-1
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(temp,breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of well-covered loci",xlab="Coverage",main=sprintf("Median coverage per locus histogram in exons\nmedian: %d",median(temp)),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max*100,y.step*100),las=2)
garbage<-dev.off()

png(sprintf("%slocusNonExonMedianCoverage.png",opt$fig.directory),height=1000,width=1000)
par(cex=2)

h.break.step<-1
temp<-apply(coverage[which(is.na(annotations$exon)),],1,median)
h.break.max<-ceiling(max(temp))
h.disp.step<-1
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(temp,breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of well-covered loci",xlab="Coverage",main=sprintf("Median coverage per locus histogram in non-exons\nmedian: %d",median(temp)),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max*100,y.step*100),las=2)
garbage<-dev.off()
