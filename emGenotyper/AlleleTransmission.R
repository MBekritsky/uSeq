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
  'min.confidence'  ,'c',2,"double",
  'min.allele.fit'  ,'a',2,"double",
  'min.noise.fit'   ,'n',2,"double",
  'max.p.null'      ,'0',2,"double",
  "script.file"     ,'s',2,"character",
  'help'            ,'h',0,"logical",
  'verbose'         ,'v',2,"integer",
  'allele.freq'     ,'q',2,"double"
),byrow=TRUE,ncol=4)

opt<-getopt(genotyper.specs)
if(is.null(opt$verbose))        opt$verbose        <- 1
if(is.null(opt$min.confidence)) opt$min.confidence <- 0.99
if(is.null(opt$min.allele.fit)) opt$min.allele.fit <- 1e-5
if(is.null(opt$min.noise.fit))  opt$min.noise.fit  <- 1e-3
if(is.null(opt$max.p.null))     opt$max.p.null     <- 0.02
if(is.null(opt$allele.freq))    opt$allele.freq    <- 0.01

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
hwe.pval.file    <- sprintf("%shwePVals.RData", opt$data.directory)

load(hwe.pval.file)

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

#load person information
load(person.info.file)

by.fams <- get.fam.mat(person.info)

is.parent <- person.info$relation %in% c("mother","father")

rare.file <- sprintf("%salleleInheritance.txt",opt$directory)
cat("SSC.family.ID","chr","pos","motif","ref.length","allele.length","allele.bias",
    "hwe.pval","context","context.info","n.good.geno.parents","good.geno.parent.alleles","family.type","mom.allele.count","dad.allele.count",
    "pro.allele.count","sib.allele.count","\n", sep="\t", file=rare.file)
counter <- list(rare.allele.loci = 0, too.many.nulls = 0, passed.filter = 0)

quad.report<-read.table(sprintf("%sreport_quad_20131006.txt",opt$directory),header=TRUE,sep="\t")
#remove CHP families and non-wholeblood families
quad.report<-quad.report[grep("auSSC.*wholeblood",quad.report$quad.quad_id),]
#remove wholeblood tag from remaining quad ids
quad.report$quad.quad_id<-sub("-wholeblood","",quad.report$quad.quad_id)
gatk.bad.families<-quad.report$quad.quad_id[grep("ok",quad.report$status,ignore.case=TRUE,invert=TRUE)]
message("The following families have high de novo counts or sample mix-ups according to the SSC GATK pipeline: ",
        paste(as.character(gatk.bad.families),collapse=", "),"\n")

for(i in 0:opt$intervals){
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

  # convert any factors in chunk.locus.info to characters
  f <- sapply(chunk.locus.info, is.factor)
  chunk.locus.info[f] <- lapply(chunk.locus.info[f], as.character)
   
  # get genome annotations
  annotations <- get.annotations(exon.file=exon.anno.file, intron.file=intron.anno.file,
                                 mirna.file=mirna.anno.file, utr.file=utr.anno.file,
                                 locus.info=chunk.locus.info, verbose=opt$verbose)
  
  em.geno.file<-sprintf("%semGenotypes_%d-%d.RData",opt$geno.directory,chunk.start.locus,chunk.stop.locus)
  load(em.geno.file)
  if(opt$verbose >= 1) cat("Loaded EM genotypes from",em.geno.file,"\n")
  
  for(l in seq.int(1,dim(em.geno$genotypes)[1])){
    genotype.slice <- em.geno$genotypes[l,,]
    # Filter all genotypes containing null alleles
    non.null <- genotype.slice[,1] > 0 & genotype.slice[,2] > 0
      
    # If more than half the population has null genotypes, skip the locus
    # and return NA
    if(sum(non.null & is.parent) < (0.5 * sum(is.parent))){
      counter$too.many.nulls <- counter$too.many.nulls + 1
      next
    }
      
    # Get allele frequencies
    allele.counts <- data.frame(table(genotype.slice[which(is.parent & non.null),1:2]))
    max.parents   <- opt$allele.freq * sum(is.parent & non.null)
    
    if(any(allele.counts$Freq <= max.parents)){
      counter$rare.allele.loci <- counter$rare.allele.loci + 1
      for(j in which(allele.counts$Freq <= max.parents)){
        rare.allele <- allele.counts[j,1]
        
        # find the parents with the rare allele
        parent <- which(genotype.slice[,1:2] == rare.allele & is.parent & non.null, arr.ind = T)[,1]
        seen.families <- vector("character")
        for(p in parent){
          parent.type <- as.character(person.info$relation[p])
  
          family           <- by.fams[which(by.fams == p, arr.ind = T)[1], ]
          family.id        <- as.character(person.info[family[1],'family.id'])
          if(!(family.id %in% gatk.bad.families) && !(family.id %in% seen.families)){
            seen.families <- c(seen.families, family.id)
            family.with.rare <- genotype.slice[family,]
          
            if(all(family.with.rare[,"confidence"] > opt$min.confidence) & 
               all(family.with.rare[,"allele.fit"] > opt$min.allele.fit) &
               all(family.with.rare[,"noise.fit"]  > opt$min.noise.fit) & 
               all(family.with.rare[,"p.null"] < opt$max.p.null)){
            
            
            counter$passed.filter <- counter$passed.filter + 1
            
            #get parents that have good genotypes at locus
            good.par.geno <- is.parent & genotype.slice[,"confidence"] > opt$min.confidence & 
                             genotype.slice[,"allele.fit"] > opt$min.allele.fit &
                             genotype.slice[,"noise.fit"]  > opt$min.noise.fit & 
                             genotype.slice[,"p.null"] < opt$max.p.null
            
            good.par.count   <- sum(good.par.geno)
            good.par.alleles <- length(unique(c(genotype.slice[good.par.geno,1],genotype.slice[good.par.geno,2])))
            
            # family types refer to the gender of the proband and sibling
            family.type <- paste(person.info$gender[family[3:4]],collapse="")
            
            rare.locus.info  <- chunk.locus.info[chunk.locus.info$allele.no == rare.allele & chunk.locus.info$locus.row.ind == l,
                                                 c("chr","pos","unit","ref.length","allele.no","bias")]
            rare.hwe.p.val <- hwe.pvals[chunk.start.locus + l - 1]
            inheritance.pattern <- apply(genotype.slice[family,1:2], 1, function(x) sum(x == rare.allele))
            context <- annotation.precedence(locus.annotation=annotations[l,])
            context.info <- NA
            if(context != "intergenic") context.info <- annotation.string(get.locus.annotations(l,annotations,gene.ids))
            
            cat(as.character(person.info$family.id[p]), unlist(rare.locus.info), rare.hwe.p.val, 
                context,context.info,good.par.count,good.par.alleles,family.type, inheritance.pattern,"\n",sep="\t",
                file=rare.file,append=T)
          }
        }
      }
    }
  }
  }
}
print(counter)