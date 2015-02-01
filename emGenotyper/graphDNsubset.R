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

opt$fig.directory<-sprintf("%s/SummaryFigures/",opt$directory)
dir.create(opt$fig.directory,showWarnings=FALSE) #if directory exists, dir.create does nothing

if(is.null(opt$verbose)) opt$verbose<-1
if(is.null(opt$mendel.threshold)) opt$mendel.threshold<-5e-3

opt$start.locus<-1
opt$num.loci<-1e4

locus.info.file=sprintf("%sallele_matrix_info.txt",opt$directory)
allele.info.file=sprintf("%sallele_matrix.txt",opt$directory)
family.info.file=sprintf("%sperson_info.txt",opt$directory)
report.file=sprintf("%sreport_sampleSet_20130630.txt",opt$directory)
coverage.estimator.file=sprintf("%scoverageEstimator.RData",opt$directory)

#load locus information
locus.info<-load.locus.info(locus.info.file)
locus.info<-remove.chromosomes(locus.info,c("chrX","chrY"))
locus.info<-update.locus.info.locators(locus.info)

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
by.fams<-get.fam.mat(person.info)
mom.and.pop<-which(person.info$relation %in% c("mother","father"))


#determine stop locus
opt$max.locus<-max(locus.info$locus.row.ind)
opt$intervals<-floor(opt$max.locus/opt$num.loci)

error.rates<-rep(0,opt$max.locus)
single.null.freq<-rep(0,opt$max.locus)
double.null.freq<-rep(0,opt$max.locus)
num.alleles.per.locus<-rep(0,opt$max.locus)
pct.het.per.locus<-rep(0,opt$max.locus)
num.em.iters<-rep(0,opt$max.locus)
confidence<-matrix(0,nrow=opt$max.locus,ncol=pop)
allele.fits<-matrix(0,nrow=opt$max.locus,ncol=pop)
coverage<-matrix(0,nrow=opt$max.locus,ncol=pop)
pro.mendel.score<-matrix(0,nrow=opt$max.locus,ncol=n.fams)
sib.mendel.score<-matrix(0,nrow=opt$max.locus,ncol=n.fams)
bias.rates<-matrix(0,nrow=nrow(locus.info),ncol=4)

for(i in 0:opt$intervals)
{
	chunk.start.locus<-(1 + (i * opt$num.loci))
	
	chunk.stop.locus<-((i+1) * opt$num.loci)
	if(chunk.stop.locus > max(locus.info$locus.row.ind))
	{
		chunk.stop.locus<-max(locus.info$locus.row.ind)
	}

	if(opt$num.loci < 0)
	{
		chunk.stop.locus<-max(locus.info$locus.row.ind)
	}	
	
	#set locus information for this chunk
	chunk.locus.info<-locus.info[which(locus.info$locus.row.ind >= chunk.start.locus & locus.info$locus.row.ind <= chunk.stop.locus),]
	chunk.start.allele<-which(locus.info$locus.row.ind == chunk.start.locus & locus.info$allele.no == 0)
	if(opt$verbose >= 1) cat("Loaded locus information\n")
	
	#set number of alleles to load for this chunk
	chunk.num.alleles<-nrow(chunk.locus.info)
	chunk.stop.allele<-chunk.start.allele + chunk.num.alleles - 1
	#load alleles
	chunk.alleles<-load.alleles(allele.info.file,chunk.locus.info,exclude.XY=TRUE,exclude.MT=FALSE,start.locus=chunk.start.allele,num.alleles=chunk.num.alleles)
	if(opt$verbose >= 1) cat("Loaded allele counts\n")

	#remove sex chromosomes from locus info
	chunk.locus.info<-remove.chromosomes(chunk.locus.info,c("chrX","chrY"))

	#exclude families
	chunk.locus.info<-exclude.people.from.locus.df(chunk.locus.info,chunk.alleles,all.excluded.people)
	chunk.alleles<-exclude.people.from.allele.df(chunk.alleles,all.excluded.people)
	chunk.locus.info<-update.locus.info.locators(chunk.locus.info)

	#annotate loci in this chunk
#	exon.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.exons.merged.bed",anno.dir)
#	intron.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.introns.merged.bed",anno.dir)
#	mirna.anno.file<-sprintf("%s/hg19.ms.miRNA.merged.bed",anno.dir)
#	utr.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.utr.merged.bed",anno.dir)
#	gene.id.file<-sprintf("%s/hg19GeneName.ccds.rg.geneID.txt",anno.dir)

#	gene.ids<-load.gene.ids(gene.id.file)
#	annotations<-get.annotations(exon.file=exon.anno.file,intron.file=intron.anno.file,mirna.file=mirna.anno.file,utr.file=utr.anno.file,
#                           	 locus.info=chunk.locus.info,verbose=opt$verbose)

	chunk.locus.info.file<-sprintf("%slocusInfo_%d-%d.RData",opt$directory,chunk.start.locus,chunk.stop.locus)
	load(chunk.locus.info.file)
	if(opt$verbose >= 1) cat("Loaded chunk locus info from",chunk.locus.info.file,"\n")

	em.geno.file<-sprintf("%semGenotypes_%d-%d.RData",opt$directory,chunk.start.locus,chunk.stop.locus)
	load(em.geno.file)
	if(opt$verbose >= 1) cat("Loaded EM genotypes from",em.geno.file,"\n")

#	new.locus.info.file<-sprintf("%slocusInfo_%d-%d.RData",opt$directory,chunk.start.locus,chunk.stop.locus)
#	load(new.locus.info.file)
#	if(opt$verbose >= 1) cat("Loaded locus info with genotype statistics from",new.locus.info.file,"\n")

#	max.alleles<-6 #number of most common alleles within family to consider when calculating the Mendel score
#	mendel.ind<-get.mendel.trio.indices(max.alleles)
#	n.ind<-nrow(mendel.ind) #total number of genotype combinations (of 4^max.alleles) that are Mendelian

	fam.genotypes.file<-sprintf("%sfamilyGenotypes_%d-%d.RData",opt$directory,chunk.start.locus,chunk.stop.locus)
	load(fam.genotypes.file)
	if(opt$verbose >= 1) cat("Loaded family genotypes from",fam.genotypes.file,"\n")

	bias.rates[chunk.start.allele:chunk.stop.allele,]<-as.matrix(chunk.locus.info[,c(4,8,12,13)])
	error.rates[chunk.start.locus:chunk.stop.locus]<-em.geno$locus.error.rates
	single.null.freq[chunk.start.locus:chunk.stop.locus]<-apply(em.geno$genotypes,1, function(x) sum(x[,1] != -1 & x[,2] == -1))/pop
	double.null.freq[chunk.start.locus:chunk.stop.locus]<-apply(em.geno$genotypes,1, function(x) sum(x[,1] == -1 & x[,2] == -1))/pop
	confidence[chunk.start.locus:chunk.stop.locus,]<-em.geno$genotypes[,,3]
	allele.fits[chunk.start.locus:chunk.stop.locus,]<-em.geno$genotypes[,,4]
	coverage[chunk.start.locus:chunk.stop.locus,]<-chunk.alleles$alleles[which(chunk.locus.info$allele.no == 0),]
	pro.mendel.score[chunk.start.locus:chunk.stop.locus,]<-fam.genotypes[,,10]
	sib.mendel.score[chunk.start.locus:chunk.stop.locus,]<-fam.genotypes[,,11]
	num.alleles.per.locus[chunk.start.locus:chunk.stop.locus]<-apply(em.geno$genotypes,1,function(x) sum(unique(unique(x[,1]),unique(x[,2])) > 0))
	pct.het.per.locus[chunk.start.locus:chunk.stop.locus]<-apply(em.geno$genotypes,1, function(x) sum(x[,1] != x[,2] & x[,2] > 0))/pop
	num.em.iters[chunk.start.locus:chunk.stop.locus]<-em.geno$n.em.iter
}
cat("Writing summary data to files\n")

opt$data.directory<-sprintf("%s/SummaryData/",opt$directory)
dir.create(opt$data.directory,showWarnings=FALSE) #if directory exists, dir.create does nothing

cat(bias.rates,file=sprintf("%sbiasRates.txt",opt$data.directory),sep="\n")
save(bias.rates,file=sprintf("%sbiasRates.RData",opt$data.directory))

cat(error.rates,file=sprintf("%serrorRates.txt",opt$data.directory),sep="\n")
save(error.rates,file=sprintf("%serrorRates.RData",opt$data.directory))

cat(single.null.freq,file=sprintf("%ssingleNullFreq.txt",opt$data.directory),sep="\n")
save(single.null.freq,file=sprintf("%ssingleNullFreq.RData",opt$data.directory))

cat(double.null.freq,file=sprintf("%sdoubleNullFreq.txt",opt$data.directory),sep="\n")
save(double.null.freq,file=sprintf("%sdoubleNullFreq.RData",opt$data.directory))

write.table(confidence,file=sprintf("%sgenoConfidence.txt",opt$data.directory),sep="\t")
save(confidence,file=sprintf("%sgenoConfidence.RData",opt$data.directory))

write.table(allele.fits,file=sprintf("%sgenoAlleleFit.txt",opt$data.directory),sep="\t")
save(allele.fits,file=sprintf("%sgenoAlleleFit.RData",opt$data.directory))

write.table(coverage,file=sprintf("%slocusCoverage.txt",opt$data.directory),sep="\t")
save(coverage,file=sprintf("%slocusCoverage.RData",opt$data.directory))

write.table(pro.mendel.score,file=sprintf("%sprobandMendelScores.txt",opt$data.directory),sep="\t")
save(pro.mendel.score,file=sprintf("%sprobandMendelScores.RData",opt$data.directory))

write.table(sib.mendel.score,file=sprintf("%ssiblingMendelScores.txt",opt$data.directory),sep="\t")
save(sib.mendel.score,file=sprintf("%ssiblingMendelScores.RData",opt$data.directory))

cat(num.alleles.per.locus,file=sprintf("%sallelesPerLocus.txt",opt$data.directory),sep="\t")
save(num.alleles.per.locus,file=sprintf("%sallelesPerLocus.RData",opt$data.directory))

cat(pct.het.per.locus,file=sprintf("%spctHetPerLocus.txt",opt$data.directory),sep="\t")
save(pct.het.per.locus,file=sprintf("%spctHetPerLocus.RData",opt$data.directory))

cat(num.em.iters,file=sprintf("%snumEMItersPerLocus.txt",opt$data.directory),sep="\t")
save(num.em.iters,file=sprintf("%snumEMItersPerLocus.RData",opt$data.directory))

y.step<-0.1

cat("Plotting error rates for",sum(error.rates > 0),"loci\n")

png(sprintf("%serrorRateHist.png",opt$fig.directory),height=1000,width=1000)
h<-hist(error.rates,breaks=seq(0,1,0.01))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Locus error rate",main="Locus error rate histogram",ylim=c(0,y.max))
garbage<-dev.off()

png(sprintf("%slogErrorRateHist.png",opt$fig.directory),height=1000,width=1000)
h.step<-0.05
h.break.min<-(min(log10(error.rates)) %/% h.step) - h.step
h.disp.step<-1
h.disp.min<- ifelse(h.break.min < -6, -6, ((h.break.min %/% h.disp.step) * h.disp.step) - h.disp.step)
h<-hist(log10(error.rates),breaks=seq(h.break.min,0,h.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Locus error rate",main="Locus error rate",xaxt="n",yaxt="n",xlim=c(h.disp.min,0),ylim=c(0,y.max))
axis(side=1,at=seq(h.disp.min,0,1),labels=10^(seq(h.disp.min,0,1)))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

cat("Plotting bias parameters\n")

png(sprintf("%stotalBiasRates.png",opt$fig.directory),height=1000,width=1000)
bias.temp<-bias.rates[which(bias.rates[,2] > 0 & bias.rates[,3] > 0 & !is.na(bias.rates[,4])),4]

h.break.step<-0.05
h.break.max<-((max(bias.temp) %/% h.break.step) * h.break.step) + h.break.step

h.disp.step<-0.25
h.disp.max<-ifelse(h.break.max > 3, 3, ((h.break.max %/% h.disp.step) * h.disp.step) + h.disp.step)

h<-hist(bias.temp,breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Bias parameter",main="Bias parameter histogram, all alleles",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

png(sprintf("%srefBiasRates.png",opt$fig.directory),height=1000,width=1000)
bias.temp<-bias.rates[which(bias.rates[,2] > 0 & bias.rates[,3] > 0 & !is.na(bias.rates[,4]) & bias.rates[,1] == bias.rates[,2]),4]

h.break.step<-0.05
h.break.max<-((max(bias.temp) %/% h.break.step) * h.break.step) + h.break.step

h.disp.step<-0.25
h.disp.max<-ifelse(h.break.max > 3, 3, ((h.break.max %/% h.disp.step) * h.disp.step) + h.disp.step)

h<-hist(bias.temp,breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Bias parameter",main="Bias parameter histogram, reference alleles only",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

png(sprintf("%saltBiasRates.png",opt$fig.directory),height=1000,width=1000)
bias.temp<-bias.rates[which(bias.rates[,2] > 0 & bias.rates[,3] > 0 & !is.na(bias.rates[,4]) & bias.rates[,1] != bias.rates[,2]),4]

h.break.step<-0.05
h.break.max<-((max(bias.temp) %/% h.break.step) * h.break.step) + h.break.step

h.disp.step<-0.25
h.disp.max<-ifelse(h.break.max > 3, 3, ((h.break.max %/% h.disp.step) * h.disp.step) + h.disp.step)

h<-hist(bias.temp,breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Bias parameter",main="Bias parameter histogram, non-reference alleles only",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

cat("Plotting null frequencies\n")

png(sprintf("%ssingleNullFreq.png",opt$fig.directory),height=1000,width=1000)

h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.2
h.disp.max<-1

h<-hist(single.null.freq,breaks=seq(0,h.break.max,h.break.step))
h$counts<-log10(h$counts)
y.max<-(max(h$counts) %/% 1) + 1
plot(h,border=NA,col=rgb(1,0,0,0.5),ylab="Frequency",xlab="Single null call frequency",main="Single null call frequency by locus histogram",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max),labels=10^seq(0,y.max))
abline(h=1:y.max,col="black",lty=2)
garbage<-dev.off()

png(sprintf("%sdoubleNullFreq.png",opt$fig.directory),height=1000,width=1000)

h<-hist(double.null.freq,breaks=seq(0,h.break.max,h.break.step))
h$counts<-log10(h$counts)
y.max<-(max(h$counts) %/% 1) + 1
plot(h,border=NA,col=rgb(1,0,0,0.5),ylab="Frequency",xlab="Double null call frequency",main="Double null call frequency by locus histogram",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max),labels=10^seq(0,y.max))
abline(h=1:y.max,col="black",lty=2)
garbage<-dev.off()

cat("Plotting genotype confidence scores\n")

png(sprintf("%sgenotypeConfidence.png",opt$fig.directory),height=1000,width=1000)

h<-hist(confidence,breaks=seq(0,h.break.max,h.break.step))
h$counts<-log10(h$counts)
y.max<-(max(h$counts) %/% 1) + 1
plot(h,border=NA,col=rgb(1,0,0,0.5),ylab="Frequency",xlab="Genotype confidence",main="Genotype confidence histogram",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max),labels=10^seq(0,y.max))
abline(h=1:y.max,col="black",lty=2)
garbage<-dev.off()

png(sprintf("%smedianLocusConfidence.png",opt$fig.directory),height=1000,width=1000)
h<-hist(apply(confidence,1,median),breaks=seq(0,h.break.max,h.break.step))
h$counts<-log10(h$counts)
y.max<-(max(h$counts) %/% 1) + 1
plot(h,border=NA,col=rgb(1,0,0,0.5),ylab="Frequency",xlab="Median genotype confidence",main="Median genotype confidence by locus histogram",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max),labels=10^seq(0,y.max))
abline(h=1:y.max,col="black",lty=2)
garbage<-dev.off()

png(sprintf("%smedianPersonConfidence.png",opt$fig.directory),height=1000,width=1000)
h<-hist(apply(confidence,2,median),breaks=seq(0,h.break.max,h.break.step))
h$counts<-log10(h$counts)
y.max<-(max(h$counts) %/% 1) + 1
plot(h,border=NA,col=rgb(1,0,0,0.5),ylab="Frequency",xlab="Median genotype confidence",main="Median genotype confidence by person histogram",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max),labels=10^seq(0,y.max))
abline(h=1:y.max,col="black",lty=2)
garbage<-dev.off()

cat("Plotting heterozygosity info\n")

png(sprintf("%spctHetPerLocus.png",opt$fig.directory),height=1000,width=1000)

h.break.step<-0.025

h<-hist(pct.het.per.locus,breaks=seq(0,h.break.max,h.break.step))
h$counts<-log10(h$counts)
y.max<-(max(h$counts) %/% 1) + 1
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Percent heterozygous genotype calls",main="Percent heterozygous genotype calls per locus histogram",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max),labels=10^seq(0,y.max))
abline(h=1:y.max,col="black",lty=2)
garbage<-dev.off()

png(sprintf("%snumAllelesPerLocus.png",opt$fig.directory),height=1000,width=1000)

het.freqs<-tabulate(num.alleles.per.locus)/opt$max.locus
y.max<-((max(het.freqs) %/% y.step) * y.step) + y.step

bp<-barplot(het.freqs,border=NA,col=rgb(228,108,10,200,max=255),ylab="Frequency",xlab="Number of alleles",main="Number of alleles per locus histogram",ylim=c(0,y.max),xaxt="n",yaxt="n")
axis(side=1,at=bp,labels=seq(1,max(num.alleles.per.locus),1))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

cat("Plotting genotyping EM info\n")
png(sprintf("%snumEMIterations.png",opt$fig.directory),height=1000,width=1000)

freqs<-tabulate(num.em.iters)/opt$max.locus
y.max<-((max(freqs) %/% y.step) * y.step) + y.step

bp<-barplot(freqs,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="EM iterations until convergence",main="Number of EM iterations per locus until convergence histogram",ylim=c(0,y.max),xaxt="n",yaxt="n")
axis(side=1,at=bp,labels=seq(1,max(num.em.iters),1))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

cat("Plotting Mendel scores\n")

png(sprintf("%sprobandMendelScores.png",opt$fig.directory),height=1000,width=1000)

h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.2
h.disp.max<-1

h<-hist(pro.mendel.score,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-log10(h$counts)
y.max<-(max(h$counts) %/% 1) + 1
plot(h,border=NA,col=rgb(1,0,0,0.5),ylab="Frequency",xlab="Mendel scores",main="Proband Mendel scores histogram",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max),labels=10^seq(0,y.max))
abline(h=1:y.max,col="black",lty=2)
garbage<-dev.off()

png(sprintf("%ssiblingMendelScores.png",opt$fig.directory),height=1000,width=1000)

h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.2
h.disp.max<-1

h<-hist(sib.mendel.score,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-log10(h$counts)
y.max<-(max(h$counts) %/% 1) + 1
plot(h,border=NA,col=rgb(1,0,0,0.5),ylab="Frequency",xlab="Mendel scores",main="Sibling Mendel scores histogram",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max),labels=10^seq(0,y.max))
abline(h=1:y.max,col="black",lty=2)
garbage<-dev.off()

png(sprintf("%scombinedMendelScores.png",opt$fig.directory),height=1000,width=1000)

h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.2
h.disp.max<-1

h<-hist(c(sib.mendel.score,pro.mendel.score),breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-log10(h$counts)
y.max<-(max(h$counts) %/% 1) + 1
plot(h,border=NA,col=rgb(1,0,0,0.5),ylab="Frequency",xlab="Mendel scores",main="Sibling Mendel scores histogram",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max),labels=10^seq(0,y.max))
abline(h=1:y.max,col="black",lty=2)
garbage<-dev.off()

cat("Plotting coverage info\n")

png(sprintf("%scoverage.png",opt$fig.directory),height=1000,width=1000)

h.break.step<-1
h.break.max<-max(coverage)
h.disp.step<-1
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(coverage,breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Coverage",main=sprintf("Coverage histogram\nmedian: %d",median(coverage)),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

png(sprintf("%spersonMedianCoverage.png",opt$fig.directory),height=1000,width=1000)

h.break.step<-1
h.break.max<-max(apply(coverage,2,median))
h.disp.step<-1
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(apply(coverage,2,median),breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Coverage",main=sprintf("Median coverage per person histogram\nmedian: %d",median(apply(coverage,2,median))),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

png(sprintf("%slocusMedianCoverage.png",opt$fig.directory),height=1000,width=1000)

h.break.step<-1
h.break.max<-max(apply(coverage,1,median))
h.disp.step<-1
h.disp.max<-ifelse(h.break.max > 50, 50, h.break.max)

h<-hist(apply(coverage,1,median),breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Coverage",main=sprintf("Median coverage per locus histogram\nmedian: %d",median(apply(coverage,1,median))),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

cat("Plotting allele fit info\n")

png(sprintf("%sallelefit.png",opt$fig.directory),height=1000,width=1000)

h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.01
h.disp.max<-1

h<-hist(allele.fits,breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Allele fit",main=sprintf("Allele fit histogram\nmedian: %0.2f",median(allele.fits)),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

png(sprintf("%spersonMedianAlleleFit.png",opt$fig.directory),height=1000,width=1000)

h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.01
h.disp.max<-1

h<-hist(apply(allele.fits,2,median),breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Allele fit",main=sprintf("Median allele fit per person histogram\nmedian:%0.2f",median(apply(allele.fits,2,median))),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()

png(sprintf("%slocusMedianAlleleFit.png",opt$fig.directory),height=1000,width=1000)

h.break.step<-0.01
h.break.max<-1
h.disp.step<-0.01
h.disp.max<-1

h<-hist(apply(allele.fits,1,median),breaks=seq(0,h.break.max,h.break.step))
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
plot(h,border=NA,col=rgb(0,0,0,0.5),ylab="Frequency",xlab="Coverage",main=sprintf("Median allele fit per locus histogram\nmedian: %0.2f",median(apply(allele.fits,1,median))),ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
axis(side=1,at=seq(0,h.disp.max,h.disp.step))
axis(side=2,at=seq(0,y.max,y.step))
garbage<-dev.off()
