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
  'swap.threshold'  ,'w',2,"double",
  'help'            ,'h',0,"logical",
  'verbose'         ,'v',2,"integer"
),byrow=TRUE,ncol=4)

opt<-getopt(genotyper.specs)

if(!is.null(opt$help)) {
  cat(getopt(genotyper.specs, usage=TRUE))
  q(status=1)
}

if(is.null(opt$directory)) {
  cat("You must specify a directory with the --directory flag\n", 
      getopt(genotyper.specs, usage=TRUE))
  q(status=1)
}

if(is.null(opt$verbose)) opt$verbose <- 1
if(is.null(opt$mendel.threshold)) opt$mendel.threshold <- 5e-3
if(is.null(opt$swap.threshold))   opt$swap.threshold   <- 0.8

opt$start.locus <- 1
opt$num.loci    <- 1e4

opt$geno.directory <- sprintf("%s/GenotypeInformation/", opt$directory)

locus.info.file=sprintf("%sallele_matrix_info.txt", opt$directory)
allele.info.file=sprintf("%sallele_matrix.txt", opt$directory)
family.info.file=sprintf("%sperson_info.txt", opt$directory)
report.file=sprintf("%sreport_sampleSet_20130630.txt", opt$directory)
coverage.estimator.file=sprintf("%scoverageEstimator.RData", opt$directory)
dn.file=sprintf("%sUnsortedMendelScores.txt", opt$directory)

if(file.exists(dn.file)) {
  cat("Removing old copy of", dn.file, "\n")
  unlink(dn.file)
}

#load locus information
locus.info <- load.locus.info(locus.info.file)
locus.info <- remove.chromosomes(locus.info, c("chrX", "chrY", "chrMT"))
locus.info <- update.locus.info.locators(locus.info)

#determine stop locus
opt$max.locus <- max(locus.info$locus.row.ind)
opt$intervals <- floor(opt$max.locus / opt$num.loci)

#load person info
person.info <- load.person.info(family.info.file)
person.info <- add.genders(person.info, report.file)
if(opt$verbose >= 1) cat("Loaded person information\n")

#find families where anyone has DNA sequence not derived from whole blood (specific to the SSC project)
ss.report <- load.sample.set.report(report.file)
non.whole.blood.indices <- get.non.whole.blood.family.members(ss.report, person.info,
                                                              verbose=opt$verbose)

#load coverage estimator
load(coverage.estimator.file)
if(opt$verbose >= 1) cat("Loaded coverage estimator\n")

#exclude people with low coverage, poor correlation, or with DNA not derived from whole blood
if(opt$verbose >= 1) cat("The following people have very low coverage and their families are being excluded from analysis:",
                         paste(as.character(person.info$individual.id[coverage.estimator$low.cov]), collapse=", "), "\n")
if(opt$verbose >= 1) cat("The following people have poor correlation to their expected coverage and their families are being excluded from analysis:",
                         paste(as.character(person.info$individual.id[coverage.estimator$bad.cor]), collapse=", "), "\n")
excluded.by.cm <- unique(c(as.character(person.info$family.id[coverage.estimator$low.cov]),
                           as.character(person.info$family.id[coverage.estimator$bad.cor])))
inds.excluded.by.cm <- which(person.info$family.id %in% excluded.by.cm)
all.excluded.people <- unique(c(inds.excluded.by.cm, non.whole.blood.indices$inds))

cat("Excluding the following families from analysis:",
    paste(unique(as.character(person.info$family.id[all.excluded.people])), collapse=", ", sep=", "), "\n")
person.info        <- exclude.people.from.person.df(person.info, all.excluded.people)
coverage.estimator <- exclude.people.from.exp.coverage.matrix(coverage.estimator, all.excluded.people)

pop     <- nrow(person.info)
n.fams  <- pop / 4
by.fams <- get.fam.mat(person.info)
mom.and.pop <- which(person.info$relation %in% c("mother", "father"))
moms        <- which(person.info$relation == "mother")
dads        <- which(person.info$relation == "father")
pro.m       <- which((person.info$relation == "proband" | person.info$relation == "self") & 
                      person.info$gender == "M")
pro.f       <- which((person.info$relation == "proband" | person.info$relation == "self") &
                      person.info$gender == "F")
sib.m       <- which(person.info$relation == "sibling" & person.info$gender == "M")
sib.f       <- which(person.info$relation == "sibling" & person.info$gender == "F")


for(i in 0:opt$intervals) {
  chunk.start.locus <- (1 + (i * opt$num.loci))
  chunk.stop.locus  <- ((i + 1) * opt$num.loci)
  if(chunk.stop.locus > max(locus.info$locus.row.ind))
  {
    chunk.stop.locus <- max(locus.info$locus.row.ind)
  }
  
  if(opt$num.loci < 0)
  {
    chunk.stop.locus <- max(locus.info$locus.row.ind)
  }	
  
  #set locus information for this chunk
  chunk.locus.info   <- locus.info[which(locus.info$locus.row.ind >= chunk.start.locus &
                                         locus.info$locus.row.ind <= chunk.stop.locus), ]
  chunk.start.allele <- which(locus.info$locus.row.ind == chunk.start.locus & locus.info$allele.no == 0)
  if(opt$verbose >= 1) cat("Loaded locus information\n")
  
  #set number of alleles to load for this chunk
  chunk.num.alleles <- nrow(chunk.locus.info)
  #load alleles
  chunk.alleles     <- load.alleles(allele.info.file, chunk.locus.info, exclude.XY=TRUE, exclude.MT=FALSE,
                                    start.locus=chunk.start.allele, num.alleles=chunk.num.alleles)
  if(opt$verbose >= 1) cat("Loaded allele counts\n")
  #remove sex chromosomes and mitochondrial chromosome from locus info
  chunk.locus.info <- remove.chromosomes(chunk.locus.info, c("chrX", "chrY", "chrMT"))
  chunk.locus.info <- exclude.people.from.locus.df(chunk.locus.info, chunk.alleles, all.excluded.people)
  chunk.alleles    <- exclude.people.from.allele.df(chunk.alleles, all.excluded.people)
  
  chunk.locus.info <- update.locus.info.locators(chunk.locus.info)
  #annotate loci in this chunk
  exon.anno.file   <- sprintf("%s/hg19.ms.ccds.rg.info.exons.merged.bed", anno.dir)
  intron.anno.file <- sprintf("%s/hg19.ms.ccds.rg.info.introns.merged.bed", anno.dir)
  mirna.anno.file  <- sprintf("%s/hg19.ms.miRNA.merged.bed", anno.dir)
  utr.anno.file    <- sprintf("%s/hg19.ms.ccds.rg.info.utr.merged.bed", anno.dir)
  gene.id.file     <- sprintf("%s/hg19GeneName.ccds.rg.geneID.txt", anno.dir)
  
  gene.ids    <- load.gene.ids(gene.id.file)
  annotations <- get.annotations(exon.file=exon.anno.file, intron.file=intron.anno.file,
                                 mirna.file=mirna.anno.file, utr.file=utr.anno.file,
                                 locus.info=chunk.locus.info, verbose=opt$verbose)
  
  em.geno.file <- sprintf("%semGenotypes_%d-%d.RData", opt$geno.directory, chunk.start.locus, 
                          chunk.stop.locus)
  load(em.geno.file)
  if(opt$verbose >= 1) cat("Loaded EM genotypes from", em.geno.file, "\n")
  
  new.locus.info.file <- sprintf("%slocusInfo_%d-%d.RData", opt$geno.directory, chunk.start.locus,
                                 chunk.stop.locus)
  load(new.locus.info.file)
  if(opt$verbose >= 1) cat("Loaded locus info with genotype statistics from", new.locus.info.file, "\n")
  
  max.alleles <- 6 #number of most common alleles within family to consider when calculating the Mendel score
  mendel.ind  <- get.mendel.trio.indices(max.alleles)
  n.ind       <- nrow(mendel.ind) #total number of genotype combinations (of 4^max.alleles) that are Mendelian
  
  fam.genotypes.file <- sprintf("%sfamilyGenotypes_%d-%d.RData", opt$geno.directory, chunk.start.locus,
                                chunk.stop.locus)
  load(fam.genotypes.file)
  if(opt$verbose >= 1) cat("Loaded family genotypes from", fam.genotypes.file, "\n")
  
  # for mendel violation, everyone in family must be not null
  mendel.violators <- which((apply(fam.genotypes[, , 11:12], 1:2, function(x) sum(x <= opt$mendel.threshold)) > 0) & 
                            (apply(fam.genotypes[, , 1:8], 1:2, function(x) sum(x < 0)) == 0),
                            arr.ind=TRUE)
  num.mendel.violators <- nrow(mendel.violators)
  
  if(opt$verbose >= 1) cat("Found", num.mendel.violators, "potential Mendel violations\n")
  
  mendel.violator.table <- data.frame(locus.num=rep(0, num.mendel.violators),
                                      allele.num=rep(0, num.mendel.violators),
                                      chunk.locus.num=rep(0, num.mendel.violators),
                                      chunk.allele.num=rep(0, num.mendel.violators),
                                      family.num=rep(0, num.mendel.violators),
                                      chr=rep(0, num.mendel.violators),
                                      pos=rep(0, num.mendel.violators),
                                      motif=rep(0, num.mendel.violators),
                                      motif.length=rep(0, num.mendel.violators),
                                      ref.length=rep(0, num.mendel.violators),
                                      em.noise.rate=rep(0, num.mendel.violators),
                                      locus.diversity=rep(0, num.mendel.violators),
                                      non.null.par.genotypes=rep(0, num.mendel.violators),
                                      context=rep(0, num.mendel.violators),
                                      context.details=rep(0, num.mendel.violators),
                                      family.id=rep(0, num.mendel.violators),
                                      mom.allele.1=rep(0, num.mendel.violators),
                                      mom.allele.2=rep(0, num.mendel.violators),
                                      mom.confidence=rep(0, num.mendel.violators),
                                      mom.allele.fit=rep(0, num.mendel.violators),
                                      mom.noise.fit=rep(0, num.mendel.violators),
                                      mom.p.null=rep(0, num.mendel.violators),
                                      mom.allele.1.cov=rep(0, num.mendel.violators),
                                      mom.allele.2.cov=rep(0, num.mendel.violators),
                                      mom.noise.cov=rep(0, num.mendel.violators),
                                      mom.exp.cov=rep(0, num.mendel.violators),
                                      dad.allele.1=rep(0, num.mendel.violators),
                                      dad.allele.2=rep(0, num.mendel.violators),
                                      dad.confidence=rep(0, num.mendel.violators),
                                      dad.allele.fit=rep(0, num.mendel.violators),
                                      dad.noise.fit=rep(0, num.mendel.violators),
                                      dad.p.null=rep(0, num.mendel.violators),
                                      dad.allele.1.cov=rep(0, num.mendel.violators),
                                      dad.allele.2.cov=rep(0, num.mendel.violators),
                                      dad.noise.cov=rep(0, num.mendel.violators),
                                      dad.exp.cov=rep(0, num.mendel.violators),
                                      pro.allele.1=rep(0, num.mendel.violators),
                                      pro.allele.2=rep(0, num.mendel.violators),
                                      pro.confidence=rep(0, num.mendel.violators),
                                      pro.allele.fit=rep(0, num.mendel.violators),
                                      pro.noise.fit=rep(0, num.mendel.violators),
                                      pro.p.null=rep(0, num.mendel.violators),
                                      pro.allele.1.cov=rep(0, num.mendel.violators),
                                      pro.allele.2.cov=rep(0, num.mendel.violators),
                                      pro.noise.cov=rep(0, num.mendel.violators),
                                      pro.exp.cov=rep(0, num.mendel.violators),
                                      sib.allele.1=rep(0, num.mendel.violators),
                                      sib.allele.2=rep(0, num.mendel.violators),
                                      sib.confidence=rep(0, num.mendel.violators),
                                      sib.allele.fit=rep(0, num.mendel.violators),
                                      sib.noise.fit=rep(0, num.mendel.violators),
                                      sib.p.null=rep(0, num.mendel.violators),
                                      sib.allele.1.cov=rep(0, num.mendel.violators),
                                      sib.allele.2.cov=rep(0, num.mendel.violators),
                                      sib.noise.cov=rep(0, num.mendel.violators),
                                      sib.exp.cov=rep(0, num.mendel.violators),
                                      pro.mendel.score=rep(0, num.mendel.violators),
                                      sib.mendel.score=rep(0, num.mendel.violators),
                                      who.violates=rep(0, num.mendel.violators),
                                      obedience.score=rep(0, num.mendel.violators),
									  swap.obedience.score=rep(0, num.mendel.violators),
                                      ins.or.del=rep(0, num.mendel.violators),
                                      om.or.com=rep(0, num.mendel.violators),
                                      dn.confidence=rep(0, num.mendel.violators),
                                      dn.noise.fit=rep(0, num.mendel.violators),
                                      dn.allele.fit=rep(0, num.mendel.violators),
                                      dn.p.null=rep(0, num.mendel.violators),
                                      dn.score=rep(0, num.mendel.violators),
                                      dn.allele=rep(0, num.mendel.violators),
                                      dn.occ.in.mom=rep(0, num.mendel.violators),
                                      dn.occ.in.dad=rep(0, num.mendel.violators),
                                      dn.inh.in.pro.m=rep(0, num.mendel.violators),
                                      dn.inh.in.pro.f=rep(0, num.mendel.violators),
                                      dn.inh.in.sib.m=rep(0, num.mendel.violators),
                                      dn.inh.in.sib.f=rep(0, num.mendel.violators),
                                      dn.allele.in.non.null.parents=rep(0, num.mendel.violators),
                                      dn.geno.in.non.null.parents=rep(0, num.mendel.violators),
                                      dn.allele.freq.in.non.null.parents=rep(0, num.mendel.violators),
                                      dn.geno.freq.in.non.null.parents=rep(0, num.mendel.violators)
  )

  mendel.violator.table$chunk.locus.num  <- mendel.violators[, 1]
  mendel.violator.table$locus.num        <- mendel.violator.table$chunk.locus.num + chunk.start.locus - 1
  mendel.violator.table$allele.num       <- sapply(mendel.violator.table$locus.num, 
                                                   function(x) which(locus.info$locus.row.ind == x & locus.info$allele.no == 0))
  mendel.violator.table$chunk.allele.num <- sapply(mendel.violator.table$chunk.locus.num, 
                                                   function(x) which(chunk.locus.info$locus.row.ind == x & chunk.locus.info$allele.no == 0))
  mendel.violator.table$family.num       <- mendel.violators[, 2]
  
  if(opt$verbose >= 1) cat("Got potential Mendel violation indices\n")
  
  #remove genotypes where the family and individual genotypes are not the same
  num.concordant.genotypes <- apply(mendel.violator.table, 1, 
                                    function(x) sum(fam.genotypes[x['chunk.locus.num'], x['family.num'], 1:8] == 
                                                    matrix(em.geno$genotypes[x["chunk.locus.num"], by.fams[x["family.num"], ], 1:2],
                                                           nrow=1, byrow=TRUE)[c(1, 5, 2, 6, 3, 7, 4, 8)]))
  mendel.violator.table    <- mendel.violator.table[-1 * (which(num.concordant.genotypes < 8)), ]
  
  if(opt$verbose >= 1) cat("Removed", length(which(num.concordant.genotypes < 8)),
                           "Mendel violations with discordances between individual and family genotypes\n")
  
  mendel.violator.table$chr             <- chunk.locus.info$chr[mendel.violator.table$chunk.allele.num]
  mendel.violator.table$pos             <- chunk.locus.info$pos[mendel.violator.table$chunk.allele.num]
  mendel.violator.table$motif           <- chunk.locus.info$unit[mendel.violator.table$chunk.allele.num]
  mendel.violator.table$motif.length    <- nchar(as.character(mendel.violator.table$motif))
  mendel.violator.table$ref.length      <- chunk.locus.info$ref.length[mendel.violator.table$chunk.allele.num]
  mendel.violator.table$em.noise.rate   <- em.geno$locus.error.rates[mendel.violator.table$chunk.locus.num]
  mendel.violator.table$locus.diversity <- sapply(mendel.violator.table$chunk.locus.num, 
                                                  function(x) sum(unique(c(unique(em.geno$genotypes[x, , 1]),
                                                                           unique(em.geno$genotypes[x, , 2]))) > 0))
  mendel.violator.table$non.null.par.genotypes <- sapply(mendel.violator.table$chunk.locus.num, 
                                                         function(x) sum(em.geno$genotypes[x, mom.and.pop, 1] > 0 &
                                                                         em.geno$genotypes[x, mom.and.pop, 2] > 0))
  mendel.violator.table$context                <- sapply(mendel.violator.table$chunk.locus.num, 
                                                         function(x) annotation.precedence(annotations[x, ]))
  mendel.violator.table$context.details        <- sapply(mendel.violator.table$chunk.locus.num,
                                                        function(x) annotation.string(get.locus.annotations(x, annotations, gene.ids)))
  mendel.violator.table$family.id              <- sapply(mendel.violator.table$family.num,
                                                         function(x) person.info$family.id[by.fams[x, 1]])
  
  if(opt$verbose >= 1) cat("Got general locus info\n")
  
  mendel.violator.table$mom.allele.1     <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 1], 1])
  mendel.violator.table$mom.allele.2     <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 1], 2])
  mendel.violator.table$mom.confidence   <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 1], 3])
  mendel.violator.table$mom.allele.fit   <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 1], 4])
  mendel.violator.table$mom.noise.fit    <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1, 
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 1], 5])
  mendel.violator.table$mom.p.null       <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1, 
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 1], 11])
  mendel.violator.table$mom.allele.1.cov <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 1], 6])
  mendel.violator.table$mom.allele.2.cov <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 1], 7])
  mendel.violator.table$mom.noise.cov    <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 1], 8])
  mendel.violator.table$mom.exp.cov      <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 1], 9])
  
  mendel.violator.table$dad.allele.1     <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 2], 1])
  mendel.violator.table$dad.allele.2     <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 2], 2])
  mendel.violator.table$dad.confidence   <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 2], 3])
  mendel.violator.table$dad.allele.fit   <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 2], 4])
  mendel.violator.table$dad.noise.fit    <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 2], 5])
  mendel.violator.table$dad.p.null       <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 2], 11])
  mendel.violator.table$dad.allele.1.cov <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 2], 6])
  mendel.violator.table$dad.allele.2.cov <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 2], 7])
  mendel.violator.table$dad.noise.cov    <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 2], 8])
  mendel.violator.table$dad.exp.cov      <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 2], 9])
  
  mendel.violator.table$pro.allele.1     <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 3], 1])
  mendel.violator.table$pro.allele.2     <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 3], 2])
  mendel.violator.table$pro.confidence   <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 3], 3])
  mendel.violator.table$pro.allele.fit   <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 3], 4])
  mendel.violator.table$pro.noise.fit    <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 3], 5])
  mendel.violator.table$pro.p.null       <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 3], 11])
  mendel.violator.table$pro.allele.1.cov <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 3], 6])
  mendel.violator.table$pro.allele.2.cov <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 3], 7])
  mendel.violator.table$pro.noise.cov    <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 3], 8])
  mendel.violator.table$pro.exp.cov      <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 3], 9])
  
  mendel.violator.table$sib.allele.1     <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 4], 1])
  mendel.violator.table$sib.allele.2     <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 4], 2])
  mendel.violator.table$sib.confidence   <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 4], 3])
  mendel.violator.table$sib.allele.fit   <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 4], 4])
  mendel.violator.table$sib.noise.fit    <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 4], 5])
  mendel.violator.table$sib.p.null       <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 4], 11])
  mendel.violator.table$sib.allele.1.cov <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 4], 6])
  mendel.violator.table$sib.allele.2.cov <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 4], 7])
  mendel.violator.table$sib.noise.cov    <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 4], 8])
  mendel.violator.table$sib.exp.cov      <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) em.geno$genotypes[x[1], by.fams[x[2], 4], 9])
  
  if(opt$verbose >= 1) cat("Got genotype info\n")
  
  mendel.violator.table$pro.mendel.score <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) fam.genotypes[x[1], x[2], 11])
  mendel.violator.table$sib.mendel.score <- apply(mendel.violator.table[, c('chunk.locus.num', 'family.num')], 1,
                                                  function(x) fam.genotypes[x[1], x[2], 12])
  
  mendel.violator.table$who.violates[which(mendel.violator.table$pro.mendel.score <= opt$mendel.threshold & 
                                             mendel.violator.table$sib.mendel.score <= opt$mendel.threshold)] <- "both"
  mendel.violator.table$who.violates[which(mendel.violator.table$pro.mendel.score <= opt$mendel.threshold & 
                                             mendel.violator.table$sib.mendel.score > opt$mendel.threshold)]  <- "proband"
  mendel.violator.table$who.violates[which(mendel.violator.table$pro.mendel.score > opt$mendel.threshold & 
                                             mendel.violator.table$sib.mendel.score <= opt$mendel.threshold)] <- "sibling"
  
  mendel.violator.table$obedience.score[which(mendel.violator.table$who.violates == "proband")] <- 
    mendel.violator.table$pro.mendel.score[which(mendel.violator.table$who.violates == "proband")]
  mendel.violator.table$obedience.score[which(mendel.violator.table$who.violates == "sibling")] <- 
    mendel.violator.table$sib.mendel.score[which(mendel.violator.table$who.violates == "sibling")]
  mendel.violator.table$obedience.score[which(mendel.violator.table$who.violates == "both")] <-
    apply(mendel.violator.table[mendel.violator.table$who.violates == "both",c('pro.mendel.score','sib.mendel.score')],1,max)
  
  mendel.violator.table$swap.obedience.score[mendel.violator.table$who.violates == "proband"] <- 
    apply(mendel.violator.table[mendel.violator.table$who.violates == "proband", c('chunk.locus.num','family.num')],1,
	      function(x) fam.genotypes[x[1], x[2], 13])
  mendel.violator.table$swap.obedience.score[mendel.violator.table$who.violates == "sibling"] <- 
    apply(mendel.violator.table[mendel.violator.table$who.violates == "sibling", c('chunk.locus.num','family.num')],1,
	      function(x) fam.genotypes[x[1], x[2], 14])
  mendel.violator.table$swap.obedience.score[mendel.violator.table$who.violates == "both"] <-
    apply(mendel.violator.table[mendel.violator.table$who.violates == "proband", c('chunk.locus.num','family.num')],1,
	      function(x) max(c(fam.genotypes[x[1], x[2], 13], fam.genotypes[x[1], x[2], 14])))
  
  
  get.violation.info <- function(v.table, dn.alleles, child.indices, genotypes) {
    parent.alleles   <- v.table[child.indices, c('mom.allele.1', 'mom.allele.2', 'dad.allele.1', 'dad.allele.2')]
    trio.alleles     <- cbind(parent.alleles, dn.alleles)
    total.in.parents <- apply(trio.alleles, 1, function(x) sum(x[c(5, 6)] %in% x[1:4]))
    violation.info   <- data.frame(om.or.com=rep(NA, length(child.indices)))
    
    violation.info$om.or.com[which(total.in.parents < 2)]  <- "commission"
    violation.info$om.or.com[which(total.in.parents == 2)] <- "omission"
    
    #get.dn.allele
    violation.info$allele <- rep(NA, length(child.indices))
    violation.info$allele[which(total.in.parents < 2)]  <- apply(trio.alleles[which(total.in.parents < 2), ], 1,
                                                                 function(x) unique(x[4 + which(!(x[5:6] %in% x[1:4]))]))
    violation.info$allele[which(total.in.parents == 2)] <- apply(trio.alleles[which(total.in.parents == 2), ], 1,
                                                                 function(x) unique(unlist(x[4 + c(which(x[5:6] %in% x[3:4] & !(x[5:6] %in% x[1:2])),
                                                                                                   which(x[5:6] %in% x[1:2] & !(x[5:6] %in% x[3:4])))])))
    violation.info$allele[which(lapply(violation.info$allele,length) == 0)] <- NA
    
    violation.info$mom.occ       <- rep(NA,length(child.indices))
    violation.info$dad.occ       <- rep(NA,length(child.indices))
    violation.info$pro.m.inh     <- rep(NA,length(child.indices))
    violation.info$pro.f.inh     <- rep(NA,length(child.indices))
    violation.info$sib.m.inh     <- rep(NA,length(child.indices))
    violation.info$sib.f.inh     <- rep(NA,length(child.indices))
    violation.info$allele.par.count <- rep(NA,length(child.indices))
    violation.info$geno.par.count   <- rep(NA,length(child.indices))
    #get dn allele and genotype frequency
    for(i in seq_along(child.indices)) {
      child.dn.locus   <- v.table$chunk.locus.num[child.indices[i]]
      child.dn.alleles <- violation.info$allele[[i]]
      a.freq           <- vector()
      mom.occ       <- vector()
      dad.occ       <- vector()
      pro.m.inh     <- vector()
      pro.f.inh     <- vector()
      sib.m.inh     <- vector()
      sib.f.inh     <- vector()
      
      for(j in child.dn.alleles) {
        a.freq       <- c(a.freq, sum(genotypes[child.dn.locus, c(moms, dads), 1:2] == j & 
                                        genotypes[child.dn.locus, c(moms, dads), 1] > 0 & 
                                        genotypes[child.dn.locus, c(moms, dads), 2] > 0))
        mom.has.allele <- genotypes[child.dn.locus, moms, 1] == j | genotypes[child.dn.locus, moms, 2] == j
        dad.has.allele <- genotypes[child.dn.locus, dads, 1] == j | genotypes[child.dn.locus, dads, 2] == j

        mom.has.allele.inds <- which(mom.has.allele)
        dad.has.allele.inds <- which(dad.has.allele)
        
        mom.occ   <- c(mom.occ, sum(mom.has.allele))
        dad.occ   <- c(dad.occ, sum(dad.has.allele))
        pro.m.inh <- c(pro.m.inh, sum((mom.has.allele[(pro.m %/% 4) + 1] | dad.has.allele[(pro.m %/% 4) + 1]) & 
                     (genotypes[child.dn.locus, pro.m, 1] == j | genotypes[child.dn.locus, pro.m, 2] == j)))
        pro.f.inh <- c(pro.f.inh, sum((mom.has.allele[(pro.f %/% 4) + 1] | dad.has.allele[(pro.f %/% 4) + 1]) & 
                                       (genotypes[child.dn.locus, pro.f, 1] == j | genotypes[child.dn.locus, pro.f, 2] == j)))
        sib.m.inh <- c(sib.m.inh, sum((mom.has.allele[(sib.m %/% 4)] | dad.has.allele[(sib.m %/% 4)]) & 
                                        (genotypes[child.dn.locus, sib.m, 1] == j | genotypes[child.dn.locus, sib.m, 2] == j)))
        sib.f.inh <- c(sib.f.inh, sum((mom.has.allele[(sib.f %/% 4)] | dad.has.allele[(sib.f %/% 4)]) & 
                                        (genotypes[child.dn.locus, sib.f, 1] == j | genotypes[child.dn.locus, sib.f, 2] == j)))
      }
      violation.info$mom.occ[[i]]       <- sum(mom.occ)
      violation.info$dad.occ[[i]]       <- sum(dad.occ)
      violation.info$pro.m.inh[[i]]     <- sum(pro.m.inh)
      violation.info$pro.f.inh[[i]]     <- sum(pro.f.inh)
      violation.info$sib.m.inh[[i]]     <- sum(sib.m.inh)
      violation.info$sib.f.inh[[i]]     <- sum(sib.f.inh)
      violation.info$allele.par.count[[i]] <- sum(a.freq)
      violation.info$geno.par.count[[i]]   <- sum(genotypes[child.dn.locus, mom.and.pop, 1] == dn.alleles[i, 1] &
                                                    genotypes[child.dn.locus, mom.and.pop, 2] == dn.alleles[i, 2])
    }
    return(violation.info)
  }
  
  #get de novo alleles
  #for proband
  dn.proband         <- which(mendel.violator.table$who.violates == "proband")
  proband.dn.alleles <- mendel.violator.table[dn.proband, c('pro.allele.1', 'pro.allele.2')]
  t                  <- get.violation.info(mendel.violator.table, proband.dn.alleles, dn.proband, em.geno$genotypes)
  mendel.violator.table[dn.proband, c('om.or.com', 'dn.allele', 'dn.occ.in.mom', 'dn.occ.in.dad', 
                                      'dn.inh.in.pro.m', 'dn.inh.in.pro.f', 'dn.inh.in.sib.m', 
                                      'dn.inh.in.sib.f', 'dn.allele.in.non.null.parents', 
                                      'dn.geno.in.non.null.parents')] <- t
  mendel.violator.table[dn.proband,'dn.confidence'] <- 
    apply(mendel.violator.table[dn.proband,c('mom.confidence','dad.confidence','pro.confidence')],1,prod)
  mendel.violator.table[dn.proband,'dn.noise.fit'] <- 
    apply(mendel.violator.table[dn.proband,c('mom.noise.fit','dad.noise.fit','pro.noise.fit')],1,prod)
  mendel.violator.table[dn.proband,'dn.allele.fit'] <- 
    apply(mendel.violator.table[dn.proband,c('mom.allele.fit','dad.allele.fit','pro.allele.fit')],1,prod)
  mendel.violator.table[dn.proband,'dn.p.null'] <- 
    apply(mendel.violator.table[dn.proband,c('mom.p.null','dad.p.null','pro.p.null')],1,function(x) 1 - prod(1 - x))
  mendel.violator.table[dn.proband,'dn.score'] <- (1 - mendel.violator.table[dn.proband,'dn.p.null']) * 
                                                  mendel.violator.table[dn.proband,'dn.allele.fit'] *
                                                  mendel.violator.table[dn.proband,'dn.noise.fit'] *
                                                  mendel.violator.table[dn.proband,'dn.confidence']
        
  #for sibling
  dn.sibling         <- which(mendel.violator.table$who.violates == "sibling")
  sibling.dn.alleles <- mendel.violator.table[dn.sibling, c('sib.allele.1', 'sib.allele.2')]
  t                  <- get.violation.info(mendel.violator.table, sibling.dn.alleles, dn.sibling, em.geno$genotypes)
  mendel.violator.table[dn.sibling, c('om.or.com', 'dn.allele', 'dn.occ.in.mom', 'dn.occ.in.dad', 
                                      'dn.inh.in.pro.m', 'dn.inh.in.pro.f', 'dn.inh.in.sib.m', 
                                      'dn.inh.in.sib.f', 'dn.allele.in.non.null.parents', 
                                      'dn.geno.in.non.null.parents')] <- t
  mendel.violator.table[dn.sibling,'dn.confidence'] <- 
    apply(mendel.violator.table[dn.sibling,c('mom.confidence','dad.confidence','sib.confidence')],1,prod)
  mendel.violator.table[dn.sibling,'dn.noise.fit'] <- 
    apply(mendel.violator.table[dn.sibling,c('mom.noise.fit','dad.noise.fit','sib.noise.fit')],1,prod)
  mendel.violator.table[dn.sibling,'dn.allele.fit'] <- 
    apply(mendel.violator.table[dn.sibling,c('mom.allele.fit','dad.allele.fit','sib.allele.fit')],1,prod)
  mendel.violator.table[dn.sibling,'dn.p.null'] <- 
    apply(mendel.violator.table[dn.sibling,c('mom.p.null','dad.p.null','sib.p.null')],1,function(x) 1 - prod(1 - x))
  mendel.violator.table[dn.sibling,'dn.score'] <- (1 - mendel.violator.table[dn.sibling,'dn.p.null']) * 
                                                  mendel.violator.table[dn.sibling,'dn.allele.fit'] *
                                                  mendel.violator.table[dn.sibling,'dn.noise.fit'] *
                                                  mendel.violator.table[dn.sibling,'dn.confidence']
                                                  
  #for both
  dn.both         <- which(mendel.violator.table$who.violates == "both")
  both.dn.alleles <- mendel.violator.table[dn.both, c('pro.allele.1', 'pro.allele.2', 'sib.allele.1', 'sib.allele.2')]
  t               <- get.violation.info(mendel.violator.table, both.dn.alleles, dn.both, em.geno$genotypes)
  mendel.violator.table[dn.both, c('om.or.com', 'dn.allele', 'dn.occ.in.mom', 'dn.occ.in.dad', 
                                      'dn.inh.in.pro.m', 'dn.inh.in.pro.f', 'dn.inh.in.sib.m', 
                                      'dn.inh.in.sib.f', 'dn.allele.in.non.null.parents', 
                                      'dn.geno.in.non.null.parents')] <- t
  mendel.violator.table[dn.both,'dn.confidence'] <- 
    apply(mendel.violator.table[dn.both,c('mom.confidence','dad.confidence','pro.confidence','sib.confidence')],1,prod)
  mendel.violator.table[dn.both,'dn.noise.fit'] <- 
    apply(mendel.violator.table[dn.both,c('mom.noise.fit','dad.noise.fit','pro.noise.fit','sib.noise.fit')],1,prod)
  mendel.violator.table[dn.both,'dn.allele.fit'] <- 
    apply(mendel.violator.table[dn.both,c('mom.allele.fit','dad.allele.fit','pro.allele.fit','sib.allele.fit')],1,prod)
  mendel.violator.table[dn.both,'dn.p.null'] <- 
    apply(mendel.violator.table[dn.both,c('mom.p.null','dad.p.null','pro.p.null','sib.p.null')],1,function(x) 1 - prod(1 - x))
  mendel.violator.table[dn.both,'dn.score'] <- (1 - mendel.violator.table[dn.both,'dn.p.null']) * 
                                               mendel.violator.table[dn.both,'dn.allele.fit'] *
                                               mendel.violator.table[dn.both,'dn.noise.fit'] *
                                               mendel.violator.table[dn.both,'dn.confidence']
                                       
  
  mendel.violator.table$dn.allele.freq.in.non.null.parents <- mendel.violator.table$dn.allele.in.non.null.parents / 
                                                              (2 * mendel.violator.table$non.null.par.genotypes)
  mendel.violator.table$dn.geno.freq.in.non.null.parents   <- mendel.violator.table$dn.geno.in.non.null.parents /
                                                              mendel.violator.table$non.null.par.genotypes
  
  mendel.violator.table[,'ins.or.del'] <- apply(mendel.violator.table[,c('dn.allele','ref.length')],1,function(x)
                                                 ifelse(max(unlist(x['dn.allele'])) > x['ref.length'], "I",
                                                        ifelse(max(unlist(x['dn.allele'])) < x['ref.length'], "D",NA)))
  
  mendel.violator.table$dn.allele <- sapply(mendel.violator.table$dn.allele, FUN = paste, collapse = ",")
  mendel.violator.table$context.details  <- sapply(mendel.violator.table$context.details, FUN = paste, collapse = ",")
  
  if(i == 0)
  {
    write.table(mendel.violator.table, file=dn.file, row.names=FALSE, quote=TRUE, sep="\t")
  } else
  {
    write.table(mendel.violator.table, file=dn.file, row.names=FALSE, quote=TRUE, sep="\t", 
                append=TRUE, col.names=FALSE)
  }
  cat("Wrote", nrow(mendel.violator.table), "de novo mutations to", dn.file, "from locus",
      chunk.start.locus, "to", chunk.stop.locus, "\n")
}

cat("Completed writing de novo mutations to", dn.file, "\n")
