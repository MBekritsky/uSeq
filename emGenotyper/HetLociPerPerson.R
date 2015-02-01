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
check.package.installation("plyr",cran.repos=cran.repos,quietly=TRUE)

require(plyr)

genotyper.specs<-matrix(c(
  'directory'       ,'d',1,"character",
  "script.file"     ,'s',2,"character",
  'help'            ,'h',0,"logical",
  'verbose'         ,'v',2,"integer"
),byrow=TRUE,ncol=4)

opt<-getopt(genotyper.specs)

if(is.null(opt$verbose)) opt$verbose <- 1

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
opt$data.directory<-sprintf("%s/SummaryData/",opt$directory)
dir.create(opt$data.directory,showWarnings=FALSE)

opt$start.locus<-1
opt$num.loci<-1e4

locus.info.file  <- sprintf("%supdatedLocusInfo.RData",opt$directory)
person.info.file <- sprintf("%sUpdatedPersonInfo.RData",opt$directory)

#load locus information
load(locus.info.file)
updated.locus.info<-remove.chromosomes(updated.locus.info,c("chrX","chrY","chrMT"))
updated.locus.info<-update.locus.info.locators(updated.locus.info)
opt$max.locus<-max(updated.locus.info$locus.row.ind)
opt$intervals<-floor(opt$max.locus/opt$num.loci)

exon.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.exons.merged.bed",anno.dir)
intron.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.introns.merged.bed",anno.dir)
mirna.anno.file<-sprintf("%s/hg19.ms.miRNA.merged.bed",anno.dir)
utr.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.utr.merged.bed",anno.dir)
gene.id.file<-sprintf("%s/hg19GeneName.ccds.rg.geneID.txt",anno.dir)

gene.ids<-load.gene.ids(gene.id.file)

#load person information
load(person.info.file)

hwe.pvals     <- rep(0,opt$max.locus)
context       <- rep(0,opt$max.locus)
num.loci.het  <- data.frame(exon=rep(0,nrow(person.info)),utr.3=rep(0,nrow(person.info)),
                            utr.5=rep(0,nrow(person.info)),intron=rep(0,nrow(person.info)),
                            intergenic=rep(0,nrow(person.info)),mirna=rep(0,nrow(person.info)))

# Only look at parents
is.parent <- person.info[,3] %in% c("mother","father")

hwe.pval <- function(genotype.slice, pval.only=FALSE){
  # Filter all genotypes containing null alleles
  non.null <- genotype.slice[,1] > 0 & genotype.slice[,2] > 0
  
  # If more than half the population has null genotypes, skip the locus
  # and return NA
  if(sum(non.null & is.parent) < (0.5 * sum(is.parent))){
    return(NA)
  }
  
  # Get allele frequencies
  allele.counts <- data.frame(table(genotype.slice[which(is.parent & non.null),1:2]))
  
  # if there's only one allele at the locus, the HWE is 1
  if(nrow(allele.counts) == 1){
    if(pval.only){
      return(1)
    } else {
      return(list(chisq.pval=1,n.parents=sum(is.parent & non.null),table=NA))
    }
  }
    
  allele.counts$Freq <- allele.counts$Freq / (2 * sum(is.parent & non.null))
  
  possible.genotypes <- unique(t(apply(expand.grid(allele.counts$Var1,allele.counts$Var1), 1, function(x) sort(as.integer(x)))))
  
  contingency.table <- data.frame(expected=apply(possible.genotypes, 1, function(x) ifelse(x[1] == x[2], (allele.counts$Freq[allele.counts$Var1 == x[1]])^2,
                                                                                           2 * allele.counts$Freq[allele.counts$Var1 == x[1]] * allele.counts$Freq[allele.counts$Var1 == x[2]])))
  
  rownames(contingency.table) <- apply(possible.genotypes, 1, function(x) paste(x,collapse="|"))
  contingency.table$expected <- contingency.table$expected * sum(is.parent & non.null)
  
  observed <- count(as.data.frame(genotype.slice[which(is.parent & non.null),1:2]), c("allele.1", "allele.2"))
  observed <- data.frame(t(apply(observed, 1, function(x) c(sort(x[1:2]),x[3]))))
  rownames(observed) <- apply(observed[,1:2],1, function(x) paste(x, collapse="|"))
  
  contingency.table$observed <- 0
  
  for(x in rownames(observed)){
    contingency.table[x,"observed"] <- observed[x,3]
  }
  
  contingency.table$sq.diff <- ((contingency.table$observed - contingency.table$expected)^2)/contingency.table$expected
  
  score <- sum(contingency.table$sq.diff)
  df <- nrow(possible.genotypes) - nrow(allele.counts)
  
  chisq.pval <- 1-pchisq(score,df)
  
  if(pval.only){
    return(chisq.pval)
  } else {
    return(list(chisq.pval=chisq.pval, df=df, score=score, n.parents=sum(is.parent & non.null), table=contingency.table))
  }
}

#get precedence
annotation.precedence<-function(locus.annotation) {
  ifelse(sum(is.na(locus.annotation['exon']))   == 0, return("exon"),
         ifelse(sum(is.na(locus.annotation['intron'])) == 0, return("intron"),
                ifelse(sum(is.na(locus.annotation['utr'])) == 0 & grepl("utr5",locus.annotation['utr']), return("5' UTR"),
                       ifelse(sum(is.na(locus.annotation['utr'])) == 0 & grepl("utr3",locus.annotation['utr']), return("3' UTR"),
                              ifelse(sum(is.na(locus.annotation['mirna']))  == 0, return("miRNA"),
                                     return("intergenic")
                              )))))
}

for(i in 0:opt$intervals)
{
  chunk.start.locus <- (1 + (i * opt$num.loci))
  
  chunk.stop.locus  <- ((i+1) * opt$num.loci)
  if(chunk.stop.locus > max(updated.locus.info$locus.row.ind))
  {
    chunk.stop.locus <- max(updated.locus.info$locus.row.ind)
  }
  
  if(opt$num.loci < 0)
  {
    chunk.stop.locus <- max(updated.locus.info$locus.row.ind)
  }	
  
  chunk.locus.info.file<-sprintf("%slocusInfo_%d-%d.RData",opt$geno.directory,chunk.start.locus,chunk.stop.locus)
  load(chunk.locus.info.file)
  if(opt$verbose >= 1) cat("Loaded chunk locus info from",chunk.locus.info.file,"\n")
  
  annotations <- get.annotations(exon.file=exon.anno.file, intron.file=intron.anno.file,
                                 mirna.file=mirna.anno.file, utr.file=utr.anno.file,
                                 locus.info=chunk.locus.info, verbose=opt$verbose)
  
  em.geno.file<-sprintf("%semGenotypes_%d-%d.RData",opt$geno.directory,chunk.start.locus,chunk.stop.locus)
  load(em.geno.file)
  if(opt$verbose >= 1) cat("Loaded EM genotypes from",em.geno.file,"\n")
    
  # Get HWE p-values for each locus
  p.vals <- apply(em.geno$genotypes,1,function(x) hwe.pval(x,TRUE))
  # Filter loci violating HWE
  hwe.obedient <- p.vals > 0.05
  
  chunk.contexts <- apply(annotations,1,annotation.precedence)
  
  #Count heterozygous genotypes per person
  num.het.per.person <- data.frame(exon=rep(0,nrow(person.info)), utr.3=rep(0,nrow(person.info)), 
                                   utr.5=rep(0,nrow(person.info)), intron=rep(0,nrow(person.info)),
                                   intergenic=rep(0,nrow(person.info)), mirna=rep(0,nrow(person.info)))
  if(any(chunk.contexts == "exon")){
    num.het.per.person$exon       <- apply(em.geno$genotypes[which(hwe.obedient & chunk.contexts == "exon"),,1:2],2,
                                                    function(x) sum(x[,1] > 0 & x[,2] > 0 & x[,1] != x[,2]))
  }
  if(any(chunk.contexts == "3' UTR")){
    num.het.per.person$utr.3      <- apply(em.geno$genotypes[which(hwe.obedient & chunk.contexts == "3' UTR"),,1:2],2,
                                           function(x) sum(x[,1] > 0 & x[,2] > 0 & x[,1] != x[,2]))
  }
  if(any(chunk.contexts == "5' UTR")){
    num.het.per.person$utr.5      <- apply(em.geno$genotypes[which(hwe.obedient & chunk.contexts == "5' UTR"),,1:2],2,
                                        function(x) sum(x[,1] > 0 & x[,2] > 0 & x[,1] != x[,2]))
  }
  if(any(chunk.contexts == "intron")){
    num.het.per.person$intron     <- apply(em.geno$genotypes[which(hwe.obedient & chunk.contexts == "intron"),,1:2],2,
                                        function(x) sum(x[,1] > 0 & x[,2] > 0 & x[,1] != x[,2]))
  }
  if(any(chunk.contexts == "intergenic")){
    num.het.per.person$intergenic <- apply(em.geno$genotypes[which(hwe.obedient & chunk.contexts == "intergenic"),,1:2],2,
                                        function(x) sum(x[,1] > 0 & x[,2] > 0 & x[,1] != x[,2]))
  }
  if(any(chunk.contexts == "mirna")){
    num.het.per.person$mirna      <- apply(em.geno$genotypes[which(hwe.obedient & chunk.contexts == "miRNA"),,1:2],2,
                                        function(x) sum(x[,1] > 0 & x[,2] > 0 & x[,1] != x[,2]))
  }
  
  hwe.pvals[chunk.start.locus:chunk.stop.locus] <- p.vals
  num.loci.het$exon       <- num.loci.het$exon       + num.het.per.person$exon
  num.loci.het$utr.3      <- num.loci.het$utr.3      + num.het.per.person$utr.3
  num.loci.het$utr.5      <- num.loci.het$utr.5      + num.het.per.person$utr.5
  num.loci.het$intron     <- num.loci.het$intron     + num.het.per.person$intron
  num.loci.het$intergenic <- num.loci.het$intergenic + num.het.per.person$intergenic
  num.loci.het$mirna      <- num.loci.het$mirna      + num.het.per.person$mirna
  context[chunk.start.locus:chunk.stop.locus] <- chunk.contexts
}

save(num.loci.het,file=sprintf("%shetGenoPerPerson.RData",opt$data.directory))
save(hwe.pvals,file=sprintf("%shwePVals.RData",opt$data.directory))
save(context,file=sprintf("%scontexts.RData",opt$data.directory))
