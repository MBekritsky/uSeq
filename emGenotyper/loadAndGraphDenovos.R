#!/usr/bin/Rscript

get.invoked.dir<-function()
{
  invoked.file<-""
  #when running this script on SGE, file will have some uninformative script name, so invoked.file must be overridden with script_file
  if(sum(grepl("script.file",commandArgs(),fixed=TRUE)) > 0)
  {
    message("Script file specified\n",commandArgs(),"\n")
    message(grep("script.file=",commandArgs(),value=TRUE),"\n")
    invoked.file<-unlist(strsplit(grep("script.file=",commandArgs(),value=TRUE),"=",fixed=TRUE))[2]
    message(invoked.file,"\n")
  } else
  {
    invoked.file<-unlist(strsplit(grep("file=",commandArgs(),value=TRUE),"=",fixed=TRUE))[2]
  }
  message(invoked.file,"\n")
  invoked.parsed.dir<-unlist(strsplit(invoked.file,"/",fixed=TRUE))
  invoked.dir<-sprintf("%s/",paste(invoked.parsed.dir[1:(length(invoked.parsed.dir) - 1)],sep="/",collapse="/"))
  return(invoked.dir)
}

source.genotyper.functions<-function(src.dir,trace=TRUE,...)
{
  for(src.file in list.files(src.dir,full.names=TRUE,pattern=".R$"))
  {
    if(trace) message(src.file, ":")
    source(src.file)
    if(trace) message("\n")
  }
}

check.for.local.R.dir<-function()
{
  local.R.dir<-sprintf("~/R/%s-library/%d.%d",version$platform,as.integer(version$major),as.integer(version$minor) %/% 1)
  if(!file.exists(local.R.dir)) dir.create(local.R.dir,recursive=TRUE)
  message("Created local R directory ",local.R.dir,"\n")
  
  if(!(local.R.dir %in% .libPaths())) .libPaths(local.R.dir)
  message("Added ",local.R.dir," to this R session's library tree\nLibrary paths are now ",paste(.libPaths(),collapse=", "),"\n")
}

invoked.dir<-get.invoked.dir()
config.file<-sprintf("%sconfig.txt",invoked.dir)
config.vars<-read.table(config.file,colClasses="character")

src.dir<-config.vars[1,2]
anno.dir<-config.vars[2,2]
cran.repos<-config.vars[3,2]
source.genotyper.functions(src.dir,trace=FALSE)

check.for.local.R.dir()

check.package.installation("getopt",cran.repos=cran.repos,quietly=TRUE)
check.package.installation("fields",cran.repos=cran.repos,quietly=TRUE)
check.package.installation("colorRamps",cran.repos=cran.repos,quietly=TRUE)
check.package.installation("plyr",cran.repos=cran.repos,quietly=TRUE)

genotyper.specs<-matrix(c(
  'start.locus'     ,'b',2,"integer",
  'num.loci'        ,'e',2,"integer",
  'directory'       ,'d',1,"character",
  "denovo.dir"      ,'o',2,"character",
  "script.file"     ,'s',2,"character",
  'mendel.threshold','m',2,"double",
  'swap.threshold'  ,'w',2,"double",
  'trio.confidence' ,'f',2,"double",
  'trio.p.null'     ,'p',2,"double",
  'trio.noise.fit'  ,'g',2,"double",
  'trio.allele.fit' ,'a',2,"double",
  'max.dn.noise'    ,'n',2,"double",
  'min.exp.cov'     ,'c',2,"double",
  'help'            ,'h',0,"logical",
  'verbose'         ,'v',2,"integer"
),byrow=TRUE,ncol=4)

opt<-getopt(genotyper.specs)

if(!is.null(opt$help)) {
  cat(getopt(genotyper.specs,usage=TRUE))
  q(status=1)
}

if(is.null(opt$directory))
{
  error("You must specify a directory with the --directory flag\n",getopt(genotyper.specs,usage=TRUE))
}

if(is.null(opt$verbose))          opt$verbose          <- 1
if(is.null(opt$start.locus))      opt$start.locus      <- 1
if(is.null(opt$num.loci))         opt$num.loci         <- -1
if(is.null(opt$mendel.threshold)) opt$mendel.threshold <- 10^-4
if(is.null(opt$swap.threshold))   opt$swap.threshold   <- 0.8
if(is.null(opt$trio.confidence))  opt$trio.confidence  <- 0.99
if(is.null(opt$trio.p.null))      opt$trio.p.null      <- 0.01
if(is.null(opt$trio.noise.fit))   opt$trio.noise.fit   <- 10^-3
if(is.null(opt$trio.allele.fit))  opt$trio.allele.fit  <- 10^-5
if(is.null(opt$max.dn.noise))     opt$max.dn.noise     <- 0.17
if(is.null(opt$min.exp.cov))      opt$min.exp.cov      <- 1

opt$geno.directory<-sprintf("%s/GenotypeInformation/",opt$directory)

locus.info.file=sprintf("%sallele_matrix_info.txt",opt$directory)
allele.info.file=sprintf("%sallele_matrix.txt",opt$directory)
family.info.file=sprintf("%sperson_info.txt",opt$directory)
report.file=sprintf("%sreport_sampleSet_20130630.txt",opt$directory)
coverage.estimator.file=sprintf("%scoverageEstimator.RData",opt$directory)

#load locus information
locus.info<-load.locus.info(locus.info.file)
locus.info<-remove.chromosomes(locus.info,c("chrX","chrY","chrMT"))
locus.info<-update.locus.info.locators(locus.info)

#determine stop locus
opt$stop.locus<-opt$start.locus + opt$num.loci - 1
if(opt$stop.locus > max(locus.info$locus.row.ind))
{
  opt$stop.locus<-max(locus.info$locus.row.ind)
}

if(opt$num.loci < 0)
{
  opt$stop.locus<-max(locus.info$locus.row.ind)
}

#set locus information for this chunk
chunk.locus.info<-locus.info[which(locus.info$locus.row.ind >= opt$start.locus & locus.info$locus.row.ind <= opt$stop.locus),]
opt$start.allele<-which(locus.info$locus.row.ind == opt$start.locus & locus.info$allele.no == 0)
if(opt$verbose >= 1) message("Loaded locus information\n")
#set number of alleles to load for this chunk
opt$num.alleles<-nrow(chunk.locus.info)
#load alleles
chunk.alleles<-load.alleles(allele.info.file,chunk.locus.info,exclude.XY=TRUE,exclude.MT=FALSE,start.locus=opt$start.allele,num.alleles=opt$num.alleles)
if(opt$verbose >= 1) message("Loaded allele counts\n")
#remove sex chromosomes from locus info
chunk.locus.info<-remove.chromosomes(chunk.locus.info,c("chrX","chrY"))
#load person info
person.info<-load.person.info(family.info.file)
person.info<-add.genders(person.info,report.file)
if(opt$verbose >= 1) message("Loaded person information\n")

#find families where anyone has DNA sequence not derived from whole blood (specific to the SSC project)
ss.report<-load.sample.set.report(report.file)
non.whole.blood.indices<-get.non.whole.blood.family.members(ss.report,person.info,verbose=opt$verbose)

#load coverage estimator
load(coverage.estimator.file)
if(opt$verbose >= 1) message("Loaded coverage estimator\n")

#exclude people with low coverage, poor correlation, or with DNA not derived from whole blood
if(opt$verbose >= 1) message("The following people have very low coverage and their families are being excluded from analysis: ",
                         paste(as.character(person.info$individual.id[coverage.estimator$low.cov]),collapse=", "),"\n")
if(opt$verbose >= 1) message("The following people have poor correlation to their expected coverage and their families are being excluded from analysis: ",
                         paste(as.character(person.info$individual.id[coverage.estimator$bad.cor]),collapse=", "),"\n")
excluded.by.cm<-unique(c(as.character(person.info$family.id[coverage.estimator$low.cov]),as.character(person.info$family.id[coverage.estimator$bad.cor])))
inds.excluded.by.cm<-which(person.info$family.id %in% excluded.by.cm)

quad.report<-read.table(sprintf("%sreport_quad_20131006.txt",opt$directory),header=TRUE,sep="\t")
#remove CHP families and non-wholeblood families
quad.report<-quad.report[grep("auSSC.*wholeblood",quad.report$quad.quad_id),]
#remove wholeblood tag from remaining quad ids
quad.report$quad.quad_id<-sub("-wholeblood","",quad.report$quad.quad_id)
gatk.bad.families<-quad.report$quad.quad_id[grep("ok",quad.report$status,ignore.case=TRUE,invert=TRUE)]
if(opt$verbose >= 1) message("The following families have high de novo counts or sample mix-ups according to the SSC GATK pipeline: ",
                         paste(as.character(gatk.bad.families),collapse=", "),"\n")

all.excluded.people<-unique(c(inds.excluded.by.cm,non.whole.blood.indices$inds))

if(opt$verbose >= 1) message("Excluding the following families from analysis: ",
				     paste(unique(as.character(person.info$family.id[all.excluded.people])),collapse=", ",sep=", "),"\n")
person.info<-exclude.people.from.person.df(person.info,all.excluded.people)
chunk.locus.info<-exclude.people.from.locus.df(chunk.locus.info,chunk.alleles,all.excluded.people)
chunk.alleles<-exclude.people.from.allele.df(chunk.alleles,all.excluded.people)
coverage.estimator<-exclude.people.from.exp.coverage.matrix(coverage.estimator,all.excluded.people)
coverage.estimator$exp.cov<-coverage.estimator$exp.cov[opt$start.locus:opt$stop.locus,]

chunk.locus.info<-update.locus.info.locators(chunk.locus.info)
#annotate loci in this chunk
exon.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.exons.merged.bed",anno.dir)
intron.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.introns.merged.bed",anno.dir)
mirna.anno.file<-sprintf("%s/hg19.ms.miRNA.merged.bed",anno.dir)
utr.anno.file<-sprintf("%s/hg19.ms.ccds.rg.info.utr.merged.bed",anno.dir)
gene.id.file<-sprintf("%s/hg19GeneName.ccds.rg.geneID.txt",anno.dir)

gene.ids<-load.gene.ids(gene.id.file)
annotations<-get.annotations(exon.file=exon.anno.file,intron.file=intron.anno.file,mirna.file=mirna.anno.file,utr.file=utr.anno.file,
                             locus.info=chunk.locus.info,verbose=opt$verbose)

pop<-ncol(chunk.alleles$alleles)
n.fams<-pop/4
by.fams<-get.fam.mat(person.info)
mom.and.pop<-which(person.info$relation %in% c("mother","father"))

em.geno.file<-sprintf("%semGenotypes_%d-%d.RData",opt$geno.directory,opt$start.locus,opt$stop.locus)
load(em.geno.file)
if(opt$verbose >= 1) message("Loaded EM genotypes from ",em.geno.file,"\n")

new.locus.info.file<-sprintf("%slocusInfo_%d-%d.RData",opt$geno.directory,opt$start.locus,opt$stop.locus)
load(new.locus.info.file)
if(opt$verbose >= 1) message("Loaded locus info with genotype statistics from ",new.locus.info.file,"\n")

max.alleles<-6 #number of most common alleles within family to consider when calculating the Mendel score
mendel.ind<-get.mendel.trio.indices(max.alleles)
n.ind<-nrow(mendel.ind) #total number of genotype combinations (of 4^max.alleles) that are Mendelian

fam.genotypes.file<-sprintf("%sfamilyGenotypes_%d-%d.RData",opt$geno.directory,opt$start.locus,opt$stop.locus)
load(fam.genotypes.file)
if(opt$verbose >= 1) message("Loaded family genotypes from ",fam.genotypes.file,"\n")

#find quads with low mendel scores
num.null <- apply(fam.genotypes[, , 1:8], 1:2, function(x) sum(x < 0))

if(opt$verbose >= 1) message("Got number of null genotypes per family\n")

pro.low.mendel.indices     <- as.data.frame(which(fam.genotypes[, , 9]  >= opt$trio.confidence &
												  fam.genotypes[, , 11] <= opt$mendel.threshold &
                                                  fam.genotypes[, , 13] >= opt$swap.threshold & 
                                                  num.null == 0, arr.ind=T))
pro.low.mendel             <- nrow(pro.low.mendel.indices)
pro.low.mendel.indices$who <- 1

sib.low.mendel.indices     <- as.data.frame(which(fam.genotypes[, , 10] >= opt$trio.confidence &
												  fam.genotypes[, , 12] <= opt$mendel.threshold &
                                                  fam.genotypes[, , 14] >= opt$swap.threshold &
                                                  num.null == 0, arr.ind=T))
sib.low.mendel             <- nrow(sib.low.mendel.indices)
sib.low.mendel.indices$who <- 2

low.mendel.indices     <- rbind(sib.low.mendel.indices, pro.low.mendel.indices)
low.mendel.indices     <- low.mendel.indices[order(low.mendel.indices[, 1], low.mendel.indices[, 2],
                                                   low.mendel.indices[, 3]), ]

tot.low.mendel <- ddply(low.mendel.indices, .(row,col), function(x) sum(x['who']))
names(tot.low.mendel)[3] <- "who"

denovo.info<-array(NA,dim=c(nrow(low.mendel.indices),2))

if(opt$verbose >= 1) message("Filtering de novo mutations\n")

for(i in seq_along(tot.low.mendel[,1])) {
  locus<-tot.low.mendel[i,1]
  family<-tot.low.mendel[i,2]
  who.dn<-tot.low.mendel[i,3]

  fam.idx <- by.fams[family,]
  dn.idx <- if(who.dn == 1) fam.idx[1:3] else if(who.dn == 2) fam.idx[c(1,2,4)] else fam.idx

  trio.p.null      <- 1 - prod(1 - em.geno$genotypes[locus,dn.idx, 11])
  trio.noise.fit   <- prod(em.geno$genotypes[locus,dn.idx, 5])
  trio.allele.fit  <- prod(em.geno$genotypes[locus,dn.idx, 4])

  family.id<-unique(person.info$family.id[by.fams[family,]])

  if(!(family.id %in% gatk.bad.families)) {
 	  if(em.geno$locus.error.rates[locus] <= opt$max.dn.noise & trio.p.null < opt$trio.p.null &
	     trio.noise.fit > opt$trio.noise.fit & trio.allele.fit > opt$trio.allele.fit) {
        if(who.dn == 1) {
      	  denovo.info[i, 1]<-fam.genotypes[locus, family, 11]
        } else if (who.dn == 2) {
    	    denovo.info[i,1]<-fam.genotypes[locus,family,12]
        } else if (who.dn == 3) {
    	    denovo.info[i,1]<-min(fam.genotypes[locus,family,11:12])
        }
        denovo.info[i,2]<-i
	  }
  }
}

denovo.info<-denovo.info[rowSums(is.na(denovo.info)) == 0,]
denovo.info<-denovo.info[order(denovo.info[,1]),]

# plot(log10(denovo.info[,1]),main=sprintf("Top %d de novo candidate scores",nrow(denovo.info)),pch=16,ylab="log10 denovo score")

denovo.dir<-ifelse(length(opt$denovo.dir) == 0,sprintf("%sdenovos",opt$directory),sprintf("%s%s",opt$directory,opt$denovo.dir))
fig.dirs<-create.dn.figure.dirs(denovo.dir)

total.denovo<-0
dn.df<-data.frame(who=rep(NA,nrow(denovo.info)),type=rep(NA,nrow(denovo.info)),
                  context=rep(NA,nrow(denovo.info)),gender=rep(NA,nrow(denovo.info)))

if(opt$verbose >= 1) message(nrow(denovo.info)," final de novo candidates\n")
			  
dn.count<-0
if(nrow(denovo.info) > 0) {
  print("Plotting de novo mutations")
  for(i in 1:nrow(denovo.info)) {
    denovo.ind<-tot.low.mendel[denovo.info[i, 2], ]
    fam.ind<-by.fams[denovo.ind[[2]], ]
    dn.locus<-denovo.ind[[1]]
    dn.family<-denovo.ind[[2]]
    dn.person<-denovo.ind[[3]]
    
    locus.info.index<-chunk.alleles$locus.row.indices[dn.locus]
    
    dn.type<-NULL
    ifelse(dn.person == 1, dn.indices <- 5:6,
    ifelse(dn.person == 2, dn.indices <- 7:8,
                           dn.indices <- 5:8))
    ifelse(dn.person < 3, dn.count <- dn.count + 1,
           dn.count <- dn.count + 2)
    
    dn.anno<-get.locus.annotations(dn.locus,annotations,gene.ids)
    context<-annotation.precedence(dn.anno)
    
    dn.alleles <- unique(fam.genotypes[dn.locus, dn.family, dn.indices])
    graph.dir <- NULL
    if(sum(dn.alleles %in% fam.genotypes[dn.locus,dn.family, 1:4]) == length(dn.alleles)) {
      #all alleles in Mendel violation child are in parents, omission violation
      dn.df$type[i] <- "omission"
      ifelse(dn.person == 1, graph.dir <- fig.dirs[["omissions"]][["proband"]][[context]],
      ifelse(dn.person == 2, graph.dir <- fig.dirs[["omissions"]][["sibling"]][[context]],
                             graph.dir <- fig.dirs[["omissions"]][["both"]][[context]]))
    } else {
      #new allele in Mendel violation child, commission violation
      dn.df$type[i] <- "commission"
      ifelse(dn.person == 1, graph.dir <- fig.dirs[["commissions"]][["proband"]][[context]],
      ifelse(dn.person == 2, graph.dir <- fig.dirs[["commissions"]][["sibling"]][[context]],
                             graph.dir <- fig.dirs[["commissions"]][["both"]][[context]]))
    }
    
	if(sum(fam.genotypes[dn.locus,dn.family,1:8] != matrix(t(em.geno$genotypes[dn.locus,by.fams[dn.family,],1:2]),nr=1)) == 0) {
	    total.denovo<-total.denovo + 1
	    dn.df$context[i]<-context
    	ifelse(dn.person == 1, dn.df$who[i] <- "proband",
	    ifelse(dn.person == 2, dn.df$who[i] <- "sibling",
        	                   dn.df$who[i] <- "both"))
    
    	ifelse(dn.person == 1, dn.df$gender[i] <- as.character(person.info$gender[by.fams[dn.family, 3]]),
	    ifelse(dn.person == 2, dn.df$gender[i] <- as.character(person.info$gender[by.fams[dn.family, 4]]),
    	                       dn.df$gender[i] <- sprintf("%s/%s",as.character(person.info$gender[by.fams[dn.family, 3]]),
	                                           as.character(person.info$gender[by.fams[dn.family, 4]]))))
    
    	plot.family(chunk.locus.info, chunk.alleles, em.geno$genotypes, dn.locus, person.info,
        	        family.ind=by.fams[dn.family, ], fam.genotypes=fam.genotypes[dn.locus,dn.family, ],
            	    who.dn=dn.person, err.rate=em.geno$locus.error.rates[dn.locus], print.dir=graph.dir,
                	anno=dn.anno, file.type="pdf",x.max="fam")
		plot.pop.heatmap(chunk.locus.info, chunk.alleles, pop, em.geno$genotypes,dn.locus,
						 fam.genotypes=fam.genotypes[dn.locus, dn.family, ], who.dn=dn.person,
						 person.info=person.info, family.ind=by.fams[dn.family, ], print.dir=graph.dir,
						 anno=dn.anno, display.dn.freq=T, display.ref.freq=T, file.type="pdf")
	    plot.pop.allele.count(chunk.locus.info, chunk.alleles, em.geno$genotypes,dn.locus, print.dir=graph.dir,
    	               fam.genotypes=fam.genotypes[dn.locus, dn.family, ], who.dn=dn.person, anno=dn.anno,file.type="pdf")
	    plot.dn.eo.scatter(chunk.locus.info, chunk.alleles, pop, em.geno$genotypes, dn.locus, dn.family, by.fams,
    	                   coverage.estimator, who.dn=dn.person, fam.genotypes=fam.genotypes[dn.locus, dn.family, ],
        	               print.dir=graph.dir, anno=dn.anno, file.type="pdf")
	    plot.dn.cov.scatter(chunk.locus.info, chunk.alleles, pop, em.geno$genotypes, dn.locus, dn.family, by.fams, 
    	                    who.dn=dn.person, fam.genotypes=fam.genotypes[dn.locus, dn.family, ], print.dir=graph.dir,
        	                anno=dn.anno, file.type="pdf")
	    plot.dn.cov.scatter(chunk.locus.info, chunk.alleles, pop, em.geno$genotypes, dn.locus, dn.family, by.fams, 
    	                    who.dn=dn.person, fam.genotypes=fam.genotypes[dn.locus, dn.family, ], print.dir=graph.dir,
        	                anno=dn.anno, plot.3d=T, file.type="pdf")
	}
  }
  print(dn.df)
  message("Found ",dn.count," total de novo mutations, which can be found in ",denovo.dir,"\n")
} else
{
  message("No de novo mutations found\n")
}
