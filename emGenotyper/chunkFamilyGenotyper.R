#!/usr/bin/Rscript

get.invoked.dir<-function() {
  invoked.file<-""
  if(sum(grepl("script.file", commandArgs(), fixed=TRUE)) > 0) {
    # when running this script on SGE, file will have an uninformative script name
    # so invoked.file must be overridden with script_file
    cat("Script file specified\n",commandArgs(),"\n")
    invoked.file<-unlist(strsplit(grep("script.file=", commandArgs(), value=TRUE), "=", 
                                  fixed=TRUE))[2]
  } else {
    invoked.file<-unlist(strsplit(grep("file=", commandArgs(), value=TRUE), "=",
                                  fixed=TRUE))[2]
  }
  cat(invoked.file,"\n")
  invoked.parsed.dir<-unlist(strsplit(invoked.file, "/", fixed=TRUE))
  invoked.dir<-sprintf("%s/", paste(invoked.parsed.dir[1:(length(invoked.parsed.dir) - 1)],
                                    sep="/", collapse="/"))
  return(invoked.dir)
}

source.genotyper.functions<-function(src.dir,trace=TRUE,...) {
  for(src.file in list.files(src.dir, full.names=TRUE, pattern=".R$")) {
    if(trace) cat(src.file, ":")
    source(src.file)
    if(trace) cat("\n")
  }
}

check.for.local.R.dir<-function() {
  local.R.dir<-sprintf("~/R/%s-library/%d.%d", version$platform, as.integer(version$major),
                       as.integer(version$minor) %/% 1)
  if(!file.exists(local.R.dir)) dir.create(local.R.dir, recursive=TRUE)
  cat("Created local R directory", local.R.dir, "\n")
  
  if(!(local.R.dir %in% .libPaths())) .libPaths(local.R.dir)
  cat("Added", local.R.dir, "to this R session's library tree\nLibrary paths are now ",
      paste(.libPaths(), collapse=", "), "\n")
}

invoked.dir<-get.invoked.dir()
config.file<-sprintf("%sconfig.txt", invoked.dir)
config.vars<-read.table(config.file, colClasses="character")

src.dir<-config.vars[1, 2]
anno.dir<-config.vars[2, 2]
cran.repos<-config.vars[3, 2]
source.genotyper.functions(src.dir, trace=FALSE)
cat("Loaded all source files from", src.dir, "\n")

check.for.local.R.dir()

check.package.installation("getopt",     cran.repos=cran.repos, quietly=TRUE)
check.package.installation("fields",     cran.repos=cran.repos, quietly=TRUE)
check.package.installation("colorRamps", cran.repos=cran.repos, quietly=TRUE)

genotyper.specs<-matrix(c(
    'start.locus'     , 'b', 2, "integer",
    'num.loci'        , 'e', 2, "integer",
    'directory'       , 'd', 1, "character",
    "denovo.dir"      , 'o', 2, "character",
    "script.file"     , 's', 2, "character",
    'help'            , 'h', 0, "logical",
    'verbose'         , 'v', 2, "integer"
  ), byrow=TRUE, ncol=4)

opt <- getopt(genotyper.specs)

if(!is.null(opt$help)) {
  cat(getopt(genotyper.specs, usage=TRUE))
  q(status=1)
}

if(is.null(opt$directory)) {
  cat("You must specify a directory with the --directory flag\n", 
      getopt(genotyper.specs, usage=TRUE))
  q(status=1)
}

opt$geno.directory <- sprintf("%s/GenotypeInformation/", opt$directory)
conditional.make.dir(opt$geno.directory)

if(is.null(opt$verbose))          opt$verbose          <- 1
if(is.null(opt$start.locus))      opt$start.locus      <- 1
if(is.null(opt$num.loci))         opt$num.loci         <- (-1)

locus.info.file  <- sprintf("%sallele_matrix_info.txt",         opt$directory)
allele.info.file <- sprintf("%sallele_matrix.txt",              opt$directory)
family.info.file <- sprintf("%sperson_info.txt",                opt$directory)
report.file      <- sprintf("%sreport_sampleSet_20140603.txt",  opt$directory)
coverage.estimator.file <- sprintf("%scoverageEstimator.RData", opt$directory)

#load locus information
locus.info <- load.locus.info(locus.info.file)
locus.info <- remove.chromosomes(locus.info, c("chrX", "chrY", "chrMT"))
locus.info <- update.locus.info.locators(locus.info)

#determine stop locus
opt$stop.locus <- opt$start.locus + opt$num.loci - 1
if(opt$stop.locus > max(locus.info$locus.row.ind)) {
	opt$stop.locus <- max(locus.info$locus.row.ind)
}

if(opt$num.loci < 0) {
  opt$stop.locus <- max(locus.info$locus.row.ind)
}

#set locus information for this chunk
chunk.locus.info <- locus.info[which(locus.info$locus.row.ind >= opt$start.locus & 
                                     locus.info$locus.row.ind <= opt$stop.locus), ]
opt$start.allele <- which(locus.info$locus.row.ind == opt$start.locus &
                          locus.info$allele.no == 0)
if(opt$verbose >= 1) cat("Loaded locus information\n")

#set number of alleles to load for this chunk
opt$num.alleles <- nrow(chunk.locus.info)
#load alleles
chunk.alleles <- load.alleles(allele.info.file, chunk.locus.info, exclude.XY=TRUE,
                              exclude.MT=FALSE, start.locus=opt$start.allele, 
                              num.alleles=opt$num.alleles)
if(opt$verbose >= 1) cat("Loaded allele counts\n")

#remove sex chromosomes from locus info
chunk.locus.info <- remove.chromosomes(chunk.locus.info, c("chrX", "chrY"))

#load person info
person.info <- load.person.info(family.info.file)
person.info <- add.genders(person.info,report.file)
if(opt$verbose >= 1) cat("Loaded person information\n")

#find families where anyone has DNA sequence not derived from whole blood (specific to the SSC project)
ss.report               <- load.sample.set.report(report.file)
non.whole.blood.indices <- get.non.whole.blood.family.members(ss.report, person.info, verbose=opt$verbose)

#load coverage estimator
load(coverage.estimator.file)
if(opt$verbose >= 1) cat("Loaded coverage estimator\n")

#exclude people with low coverage, poor correlation, or with DNA not derived from whole blood
if(opt$verbose >= 1) cat("The following people have very low coverage and their families are being excluded from analysis:",
                         paste(as.character(person.info$individual.id[coverage.estimator$low.cov]), collapse=", "), "\n")
if(opt$verbose >= 1) cat("The following people have poor correlation to their expected coverage and their families are being excluded from analysis:",
                         paste(as.character(person.info$individual.id[coverage.estimator$bad.cor]), collapse=", "),"\n")
excluded.by.cm <- unique(c(as.character(person.info$family.id[coverage.estimator$low.cov]),
                           as.character(person.info$family.id[coverage.estimator$bad.cor])))
inds.excluded.by.cm <- which(person.info$family.id %in% excluded.by.cm)
all.excluded.people <- unique(c(inds.excluded.by.cm, non.whole.blood.indices$inds))

cat("Excluding the following families from analysis:",
    paste(unique(as.character(person.info$family.id[all.excluded.people])), collapse=", ", sep=", "), "\n")
person.info        <- exclude.people.from.person.df(person.info, all.excluded.people)
chunk.locus.info   <- exclude.people.from.locus.df(chunk.locus.info, chunk.alleles, all.excluded.people)
chunk.alleles      <- exclude.people.from.allele.df(chunk.alleles, all.excluded.people)
coverage.estimator <- exclude.people.from.exp.coverage.matrix(coverage.estimator, all.excluded.people)
coverage.estimator$exp.cov <- coverage.estimator$exp.cov[opt$start.locus:opt$stop.locus, ]

chunk.locus.info <- update.locus.info.locators(chunk.locus.info)
#annotate loci in this chunk
exon.anno.file   <- sprintf("%s/hg19.ms.ccds.rg.info.exons.merged.bed", anno.dir)
intron.anno.file <- sprintf("%s/hg19.ms.ccds.rg.info.introns.merged.bed", anno.dir)
mirna.anno.file  <- sprintf("%s/hg19.ms.miRNA.merged.bed", anno.dir)
utr.anno.file    <- sprintf("%s/hg19.ms.ccds.rg.info.utr.merged.bed", anno.dir)
gene.id.file     <- sprintf("%s/hg19GeneName.ccds.rg.geneID.txt", anno.dir)

gene.ids <- load.gene.ids(gene.id.file)
annotations <- get.annotations(exon.file=exon.anno.file, intron.file=intron.anno.file, mirna.file=mirna.anno.file,
                               utr.file=utr.anno.file, locus.info=chunk.locus.info, verbose=opt$verbose)

pop         <- ncol(chunk.alleles$alleles)
n.fams      <- length(unique(person.info$family.id))
by.fams     <- get.fam.mat(person.info)
mom.and.pop <- which(person.info$relation %in% c("mother", "father"))

em.geno.file <- sprintf("%semGenotypes_%d-%d.RData", opt$geno.directory, opt$start.locus, opt$stop.locus)
load(em.geno.file)
if(opt$verbose >= 1) cat("Loaded EM genotypes from", em.geno.file, "\n")

new.locus.info.file <- sprintf("%slocusInfo_%d-%d.RData", opt$geno.directory, opt$start.locus, opt$stop.locus)
load(new.locus.info.file)
if(opt$verbose >= 1) cat("Loaded locus info with genotype statistics from", new.locus.info.file, "\n")

# number of most common alleles within family to consider when calculating the Mendel score
max.alleles   <- 6 
mendel.ind    <- get.mendel.trio.indices(max.alleles)
# total number of genotype combinations (of 4^max.alleles) that are Mendelian
n.ind         <- nrow(mendel.ind)
fam.genotypes <- get.mendel.scores(allele.info=chunk.alleles, locus.info=chunk.locus.info, geno.list=em.geno,
                                 n.fams=n.fams, by.fam=by.fams, pop.size=pop, coverage.model=coverage.estimator,
                                 mendel.ind=mendel.ind)
fam.genotypes.file <- sprintf("%sfamilyGenotypes_%d-%d.RData", opt$geno.directory, opt$start.locus, opt$stop.locus)
save(fam.genotypes, file=fam.genotypes.file)
if(opt$verbose >= 1) cat("Saved family genotypes to", fam.genotypes.file, "\n")