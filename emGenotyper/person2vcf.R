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

genotyper.specs<-matrix(c(
  'directory'       ,'d',1,"character",
  'script.file'     ,'s',2,"character",
  'person.id'       ,'p',1,"character",
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

if(is.null(opt$person.id))
{
  cat("You must specify a the person whose genotypes you want to convert to VCF format with the --person.id flag\n",getopt(genotyper.specs,usage=TRUE))
  q(status=1)
}


opt$geno.directory<-sprintf("%s/GenotypeInformation/",opt$directory)
opt$vcf.directory<-sprintf("%s/vcf/",opt$directory)
conditional.make.dir(opt$vcf.directory)

if(is.null(opt$verbose)) opt$verbose<-1
if(is.null(opt$mendel.threshold)) opt$mendel.threshold<-5e-3

opt$start.locus<-1
opt$num.loci<-1e4

locus.info.file=sprintf("%sallele_matrix_info.txt",opt$directory)
updated.family.info.file=sprintf("%sUpdatedPersonInfo.RData",opt$directory)
report.file=sprintf("%sreport_sampleSet_20130630.txt",opt$directory)
coverage.estimator.file=sprintf("%scoverageEstimator.RData",opt$directory)

#load person info
load(updated.family.info.file)

person.ind <- which(person.info$individual.id == opt$person.id)
if(length(person.ind) == 0)
  stop("Invalid person.id ",person.id)

vcf.file <- sprintf("%s%s.uSeq.vcf", opt$vcf.directory, opt$person.id)

print.vcf.meta.information(vcf.file, person.id=opt$person.id, 
                           family.id=as.character(person.info[person.ind, 'family.id']))

#load locus information
locus.info<-load.locus.info(locus.info.file)
locus.info<-remove.chromosomes(locus.info,c("chrX","chrY","chrMT"))
locus.info<-update.locus.info.locators(locus.info)

#determine stop locus
opt$max.locus<-max(locus.info$locus.row.ind)
opt$intervals<-floor(opt$max.locus/opt$num.loci)

format.string <- get.format.string()

for(i in 0:opt$intervals)
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
	
	em.geno.file<-sprintf("%semGenotypes_%d-%d.RData",opt$geno.directory,chunk.start.locus,chunk.stop.locus)
	load(em.geno.file)
	if(opt$verbose >= 1) cat("Loaded EM genotypes from",em.geno.file,"\n")

	new.locus.info.file<-sprintf("%slocusInfo_%d-%d.RData",opt$geno.directory,chunk.start.locus,chunk.stop.locus)
	load(new.locus.info.file)
	if(opt$verbose >= 1) cat("Loaded locus info with genotype statistics from",new.locus.info.file,"\n")
  
  chunk.locus.info <- chunk.locus.info[chunk.locus.info$allele.no == 0,]
  
  vcf.table <- data.frame(chr=chunk.locus.info$chr,pos=chunk.locus.info$pos,id=".",ref=get.ref.string(chunk.locus.info),
                          alt=".", stringsAsFactors=F)

  #qual here is genotype quality, NOT marginal probability of alternate allele
  vcf.table$qual <- -10 * log10(1 - em.geno$genotypes[,person.ind,3])
  vcf.table$qual[is.infinite(vcf.table$qual)] <- 100
  
  vcf.table$filter <- "."
  
  info <- data.frame(END=(chunk.locus.info$pos + chunk.locus.info$ref.length - 1),
                     RL=chunk.locus.info$ref.length, RU=chunk.locus.info$unit,
                     POP=chunk.locus.info$pop.count, SUM=chunk.locus.info$sum.count,
                     TOP=chunk.locus.info$top.count, NR=em.geno$locus.error.rates,
                     stringsAsFactors=F)
  
  vcf.table$info <- info.table.2.string(info)
	vcf.table$format.string <- format.string
  
  useq.info <- data.frame(G1=em.geno$genotypes[,person.ind,1], G2=em.geno$genotypes[,person.ind,2],
                          GC=em.geno$genotypes[,person.ind,3], AF=em.geno$genotypes[,person.ind,4],
                          NF=em.geno$genotypes[,person.ind,5], AL1=em.geno$genotypes[,person.ind,6],
                          AL2=em.geno$genotypes[,person.ind,7], ALN=em.geno$genotypes[,person.ind,8],
                          EC=em.geno$genotypes[,person.ind,9], PNULL=-10 * log10(em.geno$genotypes[,person.ind,11]),
                          RL=chunk.locus.info$ref.length, RU=chunk.locus.info$unit, stringsAsFactors=F)

  useq.info$DP <- rowSums(useq.info[,c('AL1','AL2','ALN')])
  useq.info$SEEN <- apply(useq.info[,c('G1','G2')], 1, function(x) paste(unique(x[x > 0]), collapse=","))
  
  has.alt <- (useq.info$G1 > 0 & useq.info$G1 != useq.info$RL) | 
    (useq.info$G2 > 0 & useq.info$G2 != useq.info$RL)
  
  useq.info$ALT <- NA
  useq.info$ALT[has.alt] <- apply(useq.info[has.alt,c('G1','G2','RL')], 1, function(x) 
                                         unique(x[1:2][x[1:2] != x[3] & x[1:2] > 0]))
  
  useq.info$ALTSTR <- NA
  for(i in which(has.alt)){
    motif.length <- nchar(as.character(useq.info[i,"RU"]))
    alt.string <- lapply(unlist(useq.info$ALT[i]), function(x) sprintf("%s%s", paste(rep(useq.info[i,"RU"], 
                                                                                 as.integer(x) %/% motif.length), collapse=""), 
                                                               substr(useq.info[i,"RU"], 1, as.integer(x) %% motif.length)))
    useq.info$ALTSTR[i] <- paste(alt.string,collapse=",")
    vcf.table$alt[i] <- useq.info$ALTSTR[i]
  }
	
  useq.info$GT <- "0/0"
  useq.info$GT[has.alt] <- apply(useq.info[has.alt,c("RL","ALT","G1","G2")], 1, function(x)
    paste(if(as.integer(x["G1"]) == x["RL"]) 0 else if(x["G1"] > 0) which(unlist(x["ALT"]) == x["G1"]) else -1,
          if(as.integer(x["G2"]) == x["RL"]) 0 else if(x["G2"] > 0) which(unlist(x["ALT"]) == x["G2"]) else -1, 
          sep="/"))
  useq.info$GB <- apply(useq.info[,c('G1','G2')], 1, function(x) paste(x, collapse="/"))
  
	vcf.table$useq <- useq.table.2.string(useq.info)
  write.table(vcf.table,file=vcf.file,append=T,quote=F, sep="\t",row.names=F, col.names=F)
}