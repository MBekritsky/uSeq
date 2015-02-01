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
  'mendel.threshold','m',2,"double",
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

opt$geno.directory<-sprintf("%s/GenotypeInformation/",opt$directory)
opt$fig.directory<-sprintf("%s/SummaryFigures/",opt$directory)
dir.create(opt$fig.directory,showWarnings=FALSE) #if directory exists, dir.create does nothing

if(is.null(opt$verbose)) opt$verbose<-1
if(is.null(opt$mendel.threshold)) opt$mendel.threshold<-5e-3

opt$start.locus<-1
opt$num.loci<-1e4

locus.info.file=sprintf("%supdatedLocusInfo.RData",opt$directory)
allele.info.file=sprintf("%sallele_matrix.txt",opt$directory)
family.info.file=sprintf("%sperson_info.txt",opt$directory)
report.file=sprintf("%sreport_sampleSet_20130630.txt",opt$directory)
coverage.estimator.file=sprintf("%scoverageEstimator.RData",opt$directory)

#load locus information
load(locus.info.file)
updated.locus.info<-remove.chromosomes(updated.locus.info,c("chrX","chrY","chrMT"))
updated.locus.info<-update.locus.info.locators(updated.locus.info)

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
coverage.estimator <- exclude.people.from.exp.coverage.matrix(coverage.estimator=coverage.estimator,inds=all.excluded.people)

#annotate exome loci
exon.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.exons.merged.bed",anno.dir)
intron.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.introns.merged.bed",anno.dir)
mirna.anno.file<-sprintf("%s/hg19.ms.miRNA.merged.bed",anno.dir)
utr.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.utr.merged.bed",anno.dir)
gene.id.file<-sprintf("%s/hg19GeneName.ccds.rg.geneID.txt",anno.dir)

gene.ids<-load.gene.ids(gene.id.file)
annotations<-get.annotations(exon.file=exon.anno.file,intron.file=intron.anno.file,mirna.file=mirna.anno.file,utr.file=utr.anno.file,
                        	 locus.info=updated.locus.info,verbose=opt$verbose)

pop<-nrow(person.info)
n.fams<-pop/4
by.fams<-get.fam.mat(person.info)
mom.and.pop<-which(person.info$relation %in% c("mother","father"))


pdf(sprintf("%sexpVsObsCor.pdf",opt$fig.directory),height=8,width=8)
h<-hist(coverage.estimator$cor,breaks=seq(0,1,0.005),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% 0.1) * 0.1) + 0.1
plot(h,border=NA,col="grey50",ylab="% of SSC individuals",xlab="Correlation",xaxt="n",yaxt="n",
     main="Expected vs observed coverage correlation",ylim=c(0,y.max),xlim=c(0.9,1))
abline(h=seq(0,0.5,0.05),col="white",lwd=4)
abline(v=seq(0.9,1,0.005),col="white",lwd=2)
axis(side=1,at=seq(0.9,1,length.out=6))
axis(side=2,at=seq(0,0.5,length.out=6),labels=seq(0,50,length.out=6),las=2)
garbage<-dev.off()

opt$data.directory<-sprintf("%s/SummaryData/",opt$directory)

load(sprintf("%sbiasRates.RData",opt$data.directory))
load(sprintf("%serrorRates.RData",opt$data.directory))
load(sprintf("%ssingleNullFreq.RData",opt$data.directory))
load(sprintf("%sdoubleNullFreq.RData",opt$data.directory))
load(sprintf("%sgenoConfidence.RData",opt$data.directory))
load(sprintf("%sgenoAlleleFit.RData",opt$data.directory))
load(sprintf("%sgenoNoiseFit.RData",opt$data.directory))
load(sprintf("%sgenoNullProbs.RData",opt$data.directory))
load(sprintf("%slocusCoverage.RData",opt$data.directory))
load(sprintf("%sprobandDeNovoScores.RData",opt$data.directory))
load(sprintf("%ssiblingDeNovoScores.RData",opt$data.directory))
load(sprintf("%sallelesPerLocus.RData",opt$data.directory))
load(sprintf("%spctHetPerLocus.RData",opt$data.directory))
load(sprintf("%snumEMItersPerLocus.RData",opt$data.directory))

y.step<-0.1

cat("Plotting bias parameters\n")

pdf(sprintf("%stotalBiasRates.pdf",opt$fig.directory),height=8,width=8)
# png(sprintf("%stotalBiasRates.png",opt$fig.directory),height=1800,width=1800,res=300)
bias.temp<-bias.rates[which(bias.rates[,2] > 0 & bias.rates[,3] > 0 & !is.na(bias.rates[,4])),4]

h.break.step<-0.05
h.break.max<-((max(bias.temp) %/% h.break.step) * h.break.step) + h.break.step
h.disp.step<-0.25
h.disp.max<-ifelse(h.break.max > 3, 3, ((h.break.max %/% h.disp.step) * h.disp.step) + h.disp.step)

h<-hist(bias.temp,breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of all alleles",xlab="Bias parameter",main="Bias parameter\nall alleles",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max * 100,y.step* 100),las=2)
title(sub=sprintf("Mean: %.2f; median: %.2f", mean(bias.temp),median(bias.temp)))
garbage<-dev.off()

pdf(sprintf("%srefBiasRates.pdf",opt$fig.directory),height=8,width=8)
# png(sprintf("%srefBiasRates.png",opt$fig.directory),height=1800,width=1800,res=300)
bias.temp<-bias.rates[bias.rates[,2] > 0 & bias.rates[,3] > 0 & !is.na(bias.rates[,4]) & bias.rates[,1] == bias.rates[,2],4]

h.break.step<-0.05
h.break.max<-((max(bias.temp) %/% h.break.step) * h.break.step) + h.break.step
h.disp.step<-0.25
h.disp.max<-ifelse(h.break.max > 2, 2, ((h.break.max %/% h.disp.step) * h.disp.step) + h.disp.step)

h<-hist(bias.temp,breaks=seq(0,h.break.max,h.break.step), plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of all reference alleles",xlab="Bias parameter",
     main="EM-estimated allele bias parameter\nreference alleles only",ylim=c(0,y.max),xlim=c(0,h.disp.max),
     xaxt="n",yaxt="n")
abline(h=seq(0,y.max,0.025),col="white",lwd=4)
abline(v=seq(0,h.disp.max,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,h.disp.step),labels=NA)
axis(side=1,at=seq(0,h.disp.max,2 * h.disp.step),lwd=0,lwd.tick=0)
axis(side=2,at=seq(0,y.max,length.out=6),labels=seq(0,y.max * 100,length.out=6),las=2)
title(sub=sprintf("Mean: %.2f; median: %.2f", mean(bias.temp),median(bias.temp)))
garbage<-dev.off()

pdf(sprintf("%saltBiasRates.pdf",opt$fig.directory),height=8,width=8)
# png(sprintf("%saltBiasRates.png",opt$fig.directory),height=1800,width=1800,res=300)
bias.temp<-bias.rates[which(bias.rates[,2] > 0 & bias.rates[,3] > 0 & !is.na(bias.rates[,4]) & bias.rates[,1] != bias.rates[,2]),4]

h<-hist(bias.temp,breaks=seq(0,h.break.max,h.break.step), plot=FALSE)
h$counts<-h$counts/sum(h$counts)
plot(h,border=NA,col="grey50",ylab="% of all non-reference alleles",xlab="Bias parameter",
     main="Bias parameters\nnon-reference alleles only",ylim=c(0,y.max),xlim=c(0,h.disp.max),
     xaxt="n",yaxt="n")
abline(h=seq(0,y.max,0.025),col="white",lwd=4)
abline(v=seq(0,h.disp.max,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,h.disp.step),labels=NA)
axis(side=1,at=seq(0,h.disp.max,2 * h.disp.step),lwd=0,lwd.tick=0)
axis(side=2,at=seq(0,y.max,length.out=6),labels=seq(0,y.max * 100,length.out=6),las=2)
title(sub=sprintf("Mean: %.2f; median: %.2f", mean(bias.temp),median(bias.temp)))
garbage<-dev.off()

cat("Plotting null frequencies\n")

h.break.step<-0.005
h.break.max<-1
h.disp.step<-0.2
h.disp.max<-1

h<-hist(double.null.freq,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
x.max <- h$breaks[min(which(h$counts < 0.005))]
h.disp.max <- ((x.max %/% 0.1)*0.1) + 0.1

pdf(sprintf("%sdoubleNullFreq.pdf",opt$fig.directory),height=8,width=8)
plot(h,border=NA,col="grey50",ylab="% of well-covered loci",xlab="Double null call frequency",
     main="Double null call frequency by locus",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,h.disp.max,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,length.out=6))
axis(side=2,at=seq(0,y.max,0.1),labels=seq(0, y.max * 100, 10),las=2)
garbage<-dev.off()

h<-hist(single.null.freq,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
pdf(sprintf("%ssingleNullFreq.pdf",opt$fig.directory),height=8,width=8)
plot(h,border=NA,col="grey50",ylab="% of well-covered loci",xlab="Single null genotype frequency",
     main="Single null call frequency by locus",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,h.disp.max,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,length.out=6))
axis(side=2,at=seq(0,y.max,0.1),labels=seq(0, y.max * 100, 10),las=2)
garbage<-dev.off()

cat("Plotting genotype confidence scores\n")

pdf(sprintf("%sgenotypeConfidence.pdf",opt$fig.directory),height=8,width=8)
h<-hist(confidence,breaks=seq(0,1,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
x.min <- h$breaks[max(which(h$counts < 0.005))]
h.disp.min <- ((x.min %/% 0.1)*0.1)
plot(h,border=NA,col="grey50",ylab="% of all genotypes",xlab="Genotype confidence",
     main="Genotype confidence in all SSC individuals\nat all well-covered loci",
     ylim=c(0,y.max),xlim=c(h.disp.min,1),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(h.disp.min,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(h.disp.min,1,length.out=6))
axis(side=2,at=seq(0,y.max,0.1),labels=seq(0,y.max * 100, 10),las=2)
garbage<-dev.off()

pdf(sprintf("%smedianLocusConfidence.pdf",opt$fig.directory),height=8,width=8)
h<-hist(apply(confidence,1,median),breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
x.min <- h$breaks[max(which(h$counts < 0.005))]
h.disp.min <- ((x.min %/% 0.1)*0.1)
plot(h,border=NA,col="grey50",ylab="% of all well-covered loci",xlab="Median genotype confidence",
     main="Median genotype confidence by locus",ylim=c(0,y.max),xlim=c(h.disp.min,1),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(h.disp.min,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(h.disp.min,1,length.out=6))
axis(side=2,at=seq(0,y.max,0.1),labels=seq(0,y.max * 100,10),las=2)
garbage<-dev.off()

pdf(sprintf("%smedianPersonConfidence.pdf",opt$fig.directory),height=8,width=8)
h<-hist(apply(confidence,2,median),breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
x.min <- h$breaks[max(which(h$counts < 0.005))]
h.disp.min <- ((x.min %/% 0.1)*0.1)
plot(h,border=NA,col="grey50",ylab="% of all SSC individuals",xlab="Median genotype confidence",
     main="Median genotype confidence by person",ylim=c(0,y.max),xlim=c(h.disp.min,1),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(h.disp.min,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(h.disp.min,1,length.out=6))
axis(side=2,at=seq(0,y.max,0.1),labels=seq(0,y.max * 100,10),las=2)
garbage<-dev.off()


cat("Plotting genotyping EM info\n")
png(sprintf("%snumEMIterations.png",opt$fig.directory),height=1400,width=1400)
par(cex=3.5)
freqs<-tabulate(num.em.iters)/length(num.em.iters)
y.max<-((max(freqs) %/% y.step) * y.step) + y.step

bp<-barplot(freqs,border=NA,col="grey50",ylab="% of well-covered loci",
            xlab="EM iterations until convergence",
            main="Number of EM iterations until convergence per locus",
            ylim=c(0,y.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
axis(side=1,at=bp[seq(5,max(num.em.iters),5)],labels=seq(5,max(num.em.iters),5))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max*100,y.step*100),las=2)
garbage<-dev.off()


cat("Plotting coverage info\n")

png(sprintf("%scoverage.png",opt$fig.directory),height=1400,width=1400)
par(cex=3.5)

h.break.step<-1
h.break.max<-max(coverage)
h.disp.step<-5
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(coverage,breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Coverage",
     main=sprintf("Coverage histogram\nmedian: %d",median(coverage)),
     ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

png(sprintf("%spersonMedianCoverage.png",opt$fig.directory),height=1400,width=1400)
par(cex=3.5)
y.step<-0.05
h.break.step<-1
h.break.max<-(max(apply(coverage,2,median)) %/% 1) + 1
h.disp.step<-5
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(apply(coverage,2,median),breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Coverage",
     main=sprintf("Median coverage per person histogram\nmedian: %d",median(apply(coverage,2,median))),
     ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,length.out=6))
garbage<-dev.off()

png(sprintf("%spersonMedianExonCoverage.png",opt$fig.directory),height=1200,width=1200)
par(cex=3.5)
y.step<-0.05
h.break.step<-1
temp<-apply(coverage[which(!is.na(annotations$exon)),],2,median)
h.break.max<-(max(temp) %/% 1) + 1
h.disp.step<-5
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(temp,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Coverage",
     main=sprintf("Median coverage per person histogram in exons\nmedian: %d",median(temp)),
     ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,length.out=6))
garbage<-dev.off()

png(sprintf("%slocusMedianCoverage.png",opt$fig.directory),height=1400,width=1400)
par(cex=3.5)

h.break.step<-1
h.break.max<-(max(apply(coverage,1,median)) %/% 1) + 1
h.disp.step<-5
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(apply(coverage,1,median),breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Coverage",
     main=sprintf("Median coverage per locus histogram\nmedian: %d",median(apply(coverage,1,median))),
     ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,length.out=6))
garbage<-dev.off()

cat("Plotting allele fit info\n")

# fake.allele.people <- sample(ncol(coverage),5e5,replace=TRUE)
# fake.allele.loci   <- sample(nrow(coverage),5e5,replace=TRUE)
# fake.allele.r      <- coverage.estimator$exp.cov[cbind(fake.allele.loci,fake.allele.people)]
# fake.allele.x      <- rpois(5e5,lambda=fake.allele.r)
# fake.allele.p.vals <- vectorized.poisson.p.val(x=fake.allele.x, r=fake.allele.r, 
#                                                alternative="two.sided")
# 
# png(sprintf("%srandomAllelefit.png",opt$fig.directory),height=1400,width=1400)
# 
# h.break.step<-0.01
# h.break.max<-1
# h.disp.step<-0.1
# h.disp.max<-1
# 
# h<-hist(fake.allele.p.vals,breaks=seq(0,h.break.max,h.break.step))
# h$counts<-h$counts/sum(h$counts)
# y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
# plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Allele fit",main="Allele fit for randomly generated coverage",
#      ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
# axis(side=1,at=seq(0,h.disp.max,h.disp.step))
# axis(side=2,at=seq(0,y.max,y.step))
# garbage<-dev.off()

pdf(sprintf("%sallelefit.pdf",opt$fig.directory),height=8,width=8)
h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.1
h.disp.max<-1

h<-hist(allele.fits,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of all genotypes",xlab="Allele fit",
     main="Allele goodness-of-fit",ylim=c(0,y.max),xlim=c(0,1),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,1,length.out=6))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max * 100, y.step * 100),las=2)
garbage<-dev.off()

pdf(sprintf("%spersonMedianAlleleFit.pdf",opt$fig.directory),height=8,width=8)
h<-hist(apply(allele.fits,2,median),breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of SSC individuals",xlab="Allele fit",
     main=sprintf("Median allele fit per person histogram\nmedian:%0.2f",median(apply(allele.fits,2,median))),
     ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max * 100, y.step * 100),las=2)
garbage<-dev.off()

pdf(sprintf("%slocusMedianAlleleFit.pdf",opt$fig.directory),height=8,width=8)
h<-hist(apply(allele.fits,1,median),breaks=seq(0,h.break.max,h.break.step), plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of well-covered loci",xlab="Allele fit",
     main=sprintf("Median allele fit per locus histogram\nmedian: %0.2f",median(apply(allele.fits,1,median))),
     ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max * 100, y.step * 100),las=2)
garbage<-dev.off()

cat("Plotting noise fit info\n")

h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.1
h.disp.max<-1
h<-hist(noise.fits,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
pdf(sprintf("%snoisefit.pdf",opt$fig.directory),height=8,width=8)
plot(h,border=NA,col="grey50",ylab="% of all genotypes",xlab="Noise fit",
     main="Noise goodness-of-fit",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,0.1),labels=seq(0,y.max * 100, 10),las=2)
garbage<-dev.off()

# fake.noise.people <- sample(ncol(coverage),5e5,replace=TRUE)
# fake.noise.loci   <- sample(nrow(coverage),5e5,replace=TRUE)
# fake.noise.n      <- coverage[cbind(fake.noise.loci,fake.noise.people)]
# fake.noise.n[fake.noise.n < 1] <- 1
# fake.noise.p      <- error.rates[fake.noise.loci]
# fake.noise.x      <- rbinom(5e5,size=fake.noise.n,prob=fake.noise.p)
# fake.noise.p.vals <- vectorized.binom.p.val(x=fake.noise.x, n=fake.noise.n, 
#                                             p=fake.noise.p, alternative="greater")
# 
# png(sprintf("%srandomNoisefit.png",opt$fig.directory),height=1400,width=1400)
# par(cex=3.5)
# 
# h.break.step<-0.01
# h.break.max<-1
# h.disp.step<-0.1
# h.disp.max<-1
# 
# h<-hist(fake.noise.p.vals,breaks=seq(0,h.break.max,h.break.step))
# h$counts<-h$counts/sum(h$counts)
# y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
# plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Noise fit",main="Noise fit for randomly generated noise",
#      ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
# axis(side=1,at=seq(0,h.disp.max,h.disp.step))
# axis(side=2,at=seq(0,y.max,y.step))
# garbage<-dev.off()

pdf(sprintf("%spersonMedianNoiseFit.pdf",opt$fig.directory),height=8,width=8)

h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.1
h.disp.max<-1

h<-hist(apply(noise.fits,2,median),breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of SSC individuals",xlab="Noise fit",
     main=sprintf("Median noise fit per person histogram\nmedian:%0.2f",median(apply(noise.fits,2,median))),
     ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max * 100, y.step * 100),las=2)
garbage<-dev.off()

pdf(sprintf("%slocusMedianNoiseFit.pdf",opt$fig.directory),height=8,width=8)

h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.1
h.disp.max<-1

h<-hist(apply(noise.fits,1,median),breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of well-covered loci",xlab="Noise fit",
     main=sprintf("Median noise fit per locus histogram\nmedian: %0.2f",median(apply(noise.fits,1,median))),
     ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max * 100, y.step * 100),las=2)
garbage<-dev.off()

cat("Plotting null genotype probabilities\n")

pdf(sprintf("%snullProbability.pdf",opt$fig.directory),height=8,width=8)

h.break.step<-0.005
h.break.max<-1
h.disp.step<-0.1
h.disp.max<-0.2
y.step <- 0.1

h<-hist(p.nulls,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of all genotypes",xlab="Null probability",
     main="Marginal null probability",ylim=c(0,y.max),
     xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max * 100, y.step * 100),las=2)
garbage<-dev.off()

pdf(sprintf("%spersonMedianNullProbability.pdf",opt$fig.directory),height=8,width=8)

h.break.step<-0.001
h.break.max<-1
h.disp.step<-0.01
h.disp.max<-0.05

h<-hist(apply(p.nulls,2,median),breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of SSC individuals",xlab="Null probability",
     main=sprintf("Median null probability per person histogram\nmedian:%0.2f",median(apply(p.nulls,2,median))),
     ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max * 100, y.step * 100),las=2)
garbage<-dev.off()

pdf(sprintf("%slocusMedianNullProbability.pdf",opt$fig.directory),height=8,width=8)

h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.05
h.disp.max<-0.3

h<-hist(apply(p.nulls,1,median),breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of well-covered loci",xlab="Null probability",
     main=sprintf("Median null probability per locus histogram\nmedian: %0.2f",median(apply(p.nulls,1,median))),
     ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,1,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step),labels=seq(0,y.max * 100, y.step * 100),las=2)
garbage<-dev.off()

cat("Plotting heterozygosity info\n")

pdf(sprintf("%spctHetPerLocus.pdf",opt$fig.directory),height=8,width=8)
h.break.min <- min(log10(pct.het.per.locus[pct.het.per.locus > 0])) %/% 1
h.break.step <- 0.1
h.disp.min <- h.break.min

h<-hist(log10(pct.het.per.locus),breaks=seq(h.break.min,0,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col="grey50",ylab="% of well-covered loci",xlab="Heterozygous genotype call frequency",
     main="Heterozygous genotype call frequency per locus",ylim=c(0,y.max),xlim=c(h.break.min,0),
     xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(h.disp.min,0,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(h.disp.min,0),labels=sapply(seq(h.disp.min,0),function(x) sprintf("%.0e",10^x)))
axis(side=2,at=seq(0,y.max,0.1),labels=seq(0,y.max * 100, 10),las=2)
garbage<-dev.off()

pdf(sprintf("%snumAllelesPerLocus.pdf",opt$fig.directory),height=8,width=8)
het.freqs<-tabulate(num.alleles.per.locus)/length(num.alleles.per.locus)
y.max <- 1
y.step <- 0.1

bp<-barplot(het.freqs[het.freqs >= 1e-3],border=NA,col="grey50",
            ylab="% of well-covered loci",xlab="Number of alleles",main="Number of alleles per locus",
            ylim=c(0,y.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
axis(side=1,at=bp,labels=seq(1,length(bp),1),lwd=0,lwd.tick=2)
axis(side=2,at=seq(0,y.max,length.out=6),labels=seq(0,y.max * 100, length.out=6),las=2)
garbage<-dev.off()

cat("Plotting Mendel scores\n")

pdf(sprintf("%scombinedDeNovoScores.pdf",opt$fig.directory),height=8,width=8)

t <- -10*log10(c(pro.dn.score[pro.dn.score > 0],sib.dn.score[sib.dn.score > 0]))
h.break.step <- 1
h.break.min  <- 0
h.break.max  <- ((max(t) %/% h.break.step) * h.break.step) + h.break.step

h<-hist(t,breaks=seq(h.break.min,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max <- ((max(h$counts) %/% 0.1) * 0.1) + 0.1
x.max <- ((h$breaks[max(which(h$counts > 0.001))] %/% 10) * 10) + 10
plot(h,col="grey50",ylab="% of SSC trios",xlab="Mendel obedience scores",main="Mendel obedience scores for all SSC trios",
     ylim=c(0,y.max),xlim=c(0,20),xaxt="n",yaxt="n",border=NA)
abline(h=seq(0,y.max,0.05),col="white",lwd=4)
abline(v=seq(0,20,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,20,5))
axis(side=2,at=seq(0,y.max,0.1),labels=seq(0,y.max * 100,10),las=2)
garbage<-dev.off()
