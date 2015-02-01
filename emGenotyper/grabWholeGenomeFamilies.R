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
check.package.installation("RColorBrewer",cran.repos=cran.repos,quietly=TRUE)
check.package.installation("plyr",cran.repos=cran.repos,quietly=TRUE)
check.package.installation("abind",cran.repos=cran.repos,quietly=TRUE)

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

if(is.null(opt$verbose)) opt$verbose<-1
if(is.null(opt$mendel.threshold)) opt$mendel.threshold<-5e-3

opt$start.locus<-1
opt$num.loci<-1e4

locus.info.file=sprintf("%sallele_matrix_info.txt",opt$directory)

#load locus information
locus.info<-load.locus.info(locus.info.file)
locus.info<-remove.chromosomes(locus.info,c("chrX","chrY"))
locus.info<-update.locus.info.locators(locus.info)

#determine stop locus
opt$max.locus<-max(locus.info$locus.row.ind)
opt$intervals<-floor(opt$max.locus/opt$num.loci)

updated.locus.info<-locus.info
updated.locus.info$num.called<-NA
updated.locus.info$bias<-NA

start.allele<-1
stop.allele<--1

person.info.file<-sprintf("%s/personInfoQC.txt",opt$directory)
person.info<-read.table(person.info.file,header=TRUE)

wg.people.ind<-which(person.info$family.id %in% c("auSSC12596","auSSC12605"))

chunk.start.locus<-1
chunk.stop.locus<-1e4
em.geno.file<-sprintf("%semGenotypes_%d-%d.RData",opt$directory,chunk.start.locus,chunk.stop.locus)
load(em.geno.file)

wg.fam.geno<-em.geno$genotypes[,wg.people.ind,]

for(i in 1:opt$intervals)
{
	chunk.start.locus<-(1 + (i * opt$num.loci))
	
	chunk.stop.locus<-((i+1) * opt$num.loci)
	if(chunk.stop.locus > opt$max.locus)
	{
		chunk.stop.locus<-opt$max.locus
	}

	if(opt$num.loci < 0)
	{
		chunk.stop.locus<-opt$max.locus
	}	
	
	em.geno.file<-sprintf("%semGenotypes_%d-%d.RData",opt$directory,chunk.start.locus,chunk.stop.locus)
	load(em.geno.file)
	if(opt$verbose >= 1) cat("Loaded genotype info from",em.geno.file,"\n")
	
	wg.fam.geno<-abind(wg.fam.geno,em.geno$genotypes[,wg.people.ind,],along=1)
}


wg.geno.rdata<-sprintf("%s/wholeGenomePeopleExomeGenotypes.RData",opt$directory)
save(wg.fam.geno,file=wg.geno.rdata)
cat("Saved exome genotypes for people with whole genome seqeuncing data to ",wg.geno.rdata,"\n")
