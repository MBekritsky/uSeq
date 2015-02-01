#!/usr/bin/Rscript

get.invoked.dir<-function()
{
  invoked.file<-unlist(strsplit(grep("file=",commandArgs(),value=TRUE),"=",fixed=TRUE))[2]
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

check.for.local.R.dir<-function() {
  local.R.dir <- sprintf("~/R/%s-library/%d.%d", version$platform, as.integer(version$major),
                         as.integer(version$minor) %/% 1)
  if(!file.exists(local.R.dir)) {
    dir.create(local.R.dir, recursive=TRUE)
    cat("Created local R directory", local.R.dir, "\n")
  }
  
  if(!(local.R.dir %in% .libPaths())) {
    .libPaths(local.R.dir)
    cat("Added", local.R.dir, "to this R session's library tree\nLibrary paths are now ",
        paste(.libPaths(), collapse=", "), "\n")
  }
  return(local.R.dir)
}

invoked.dir<-get.invoked.dir()
config.file<-sprintf("%sconfig.txt",invoked.dir)
config.vars<-read.table(config.file,colClasses="character")

src.dir<-config.vars[1,2]
anno.dir<-config.vars[2,2]
cran.repos<-config.vars[3,2]
source.genotyper.functions(src.dir,trace=FALSE)
cat("Loaded all source files from",src.dir,"\n")

local.R.dir <- check.for.local.R.dir()

check.package.installation("getopt",cran.repos=cran.repos,quietly=TRUE)
check.package.installation("irlba",cran.repos=cran.repos,quietly=TRUE)

genotyper.specs<-matrix(c(
  'directory'       ,'d',1,"character",
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

if(is.null(opt$verbose)) opt$verbose<-1

locus.info.file=sprintf("%sallele_matrix_info.txt",opt$directory)
allele.info.file=sprintf("%sallele_matrix.txt",opt$directory)
family.info.file=sprintf("%sperson_info.txt",opt$directory)
report.file=sprintf("%sreport_sampleSet_20130630.txt",opt$directory)

#load person information
person.info<-load.person.info(family.info.file)
pop<-nrow(person.info)

#load alleles
locus.cov<-load.locus.total.cov(locus.info.file,allele.info.file,pop.size=pop,exclude.XY=TRUE,exclude.MT=TRUE)
if(opt$verbose >= 1) cat("Loaded locus counts\n")
coverage.estimator<-linear.coverage.model(locus.cov,print.dir=opt$directory,verbose=opt$verbose)

total.cov<-pop.total.coverage.hist(locus.cov,opt$directory)
mean.cov<-pop.mean.coverage.hist(locus.cov,opt$directory)

h.break.step<-0.1
h.break.max<-(log10(max(total.cov)) %/% 1) + 1
h.disp.step<-0.2
h.disp.max<-h.break.max
y.step <- 0.1

h<-hist(log10(total.cov),breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
x.min <-  h$breaks[max(which(h$mids < log10(mean(total.cov)) & h$counts < 0.001))]
h.disp.min <- x.min %/% 1
opt$fig.directory<-sprintf("%s/SummaryFigures/",opt$directory)
dir.create(opt$fig.directory,showWarnings=FALSE) #if directory exists, dir.create does nothing

pdf(sprintf("%stotalCoverage.pdf",opt$fig.directory),height=8,width=8)
plot(h,border=NA,col="lightgrey",ylab="",xlab="Total coverage",
     main="Total coverage in all SSC individuals",ylim=c(0,y.max),xlim=c(h.disp.min,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,h.disp.max,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(h.disp.min,h.disp.max),labels=sapply(seq(h.disp.min,h.disp.max), function(x) sprintf("%.0e",10^x)))
axis(side=2,at=seq(0,y.max,0.1),las=2)
garbage<-dev.off()

h.break.step<-1
h.break.max<-(max(mean.cov) %/% 1) + 1
h.disp.step<-1
y.step <- 0.01

h<-hist(mean.cov,breaks=seq(0,h.break.max,h.break.step),plot=FALSE)
h$counts<-h$counts/sum(h$counts)
y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
x.max <-  h$breaks[max(which(h$counts > 0.001))]
h.disp.max <- ((x.max %/% 10)*10) + 10

pdf(sprintf("%smeanCoverage.pdf",opt$fig.directory),height=8,width=8)
plot(h,border=NA,col="lightgrey",ylab="",xlab="Mean coverage",
     main="Mean coverage in all SSC individuals",ylim=c(0,y.max),xlim=c(0,h.disp.max),xaxt="n",yaxt="n")
abline(h=seq(0,1,y.step/2),col="white",lwd=4)
abline(v=seq(0,h.disp.max,h.break.step),col="white",lwd=2)
axis(side=1,at=seq(0,h.disp.max,10))
axis(side=2,at=seq(0,y.max,0.01),las=2)
garbage<-dev.off()

#find individuals with low correlation between expected and observed allele coverage or very low observed coverage and remove them
low.cov<-as.integer(which(total.cov < (mean(total.cov) - (2*sd(total.cov)))))
bad.cor<-which(coverage.estimator$cor < 0.8)
cat("The people with the following indices have coverage that is too low:",low.cov,"\n")
cat("The people with the following indices have poor correlation with their expected coverage estimators:",bad.cor,"\n")
coverage.estimator$low.cov<-low.cov
coverage.estimator$bad.cor<-bad.cor
coverage.estimator$problem.people<-unique(low.cov,bad.cor)

coverage.estimator.file<-sprintf("%scoverageEstimator.RData",opt$directory)
save(coverage.estimator,file=coverage.estimator.file)
cat("Saved coverage estimator to",coverage.estimator.file,"\n")
