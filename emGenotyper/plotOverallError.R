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

opt$geno.directory<-sprintf("%s/GenotypeInformation/",opt$directory)
opt$fig.directory<-sprintf("%s/SummaryFigures/",opt$directory)
dir.create(opt$fig.directory,showWarnings=FALSE) #if directory exists, dir.create does nothing

if(is.null(opt$verbose)) opt$verbose<-1
if(is.null(opt$mendel.threshold)) opt$mendel.threshold<-5e-3

opt$start.locus<-1
opt$num.loci<-1e4

locus.info.file=sprintf("%sallele_matrix_info.txt",opt$directory)

#load locus information
locus.info<-load.locus.info(locus.info.file)
locus.info<-remove.chromosomes(locus.info,c("chrX","chrY","chrMT"))
locus.info<-update.locus.info.locators(locus.info)

#determine stop locus
opt$max.locus<-max(locus.info$locus.row.ind)
opt$intervals<-floor(opt$max.locus/opt$num.loci)

error.info<-data.frame(em.noise.rate=rep(0,opt$max.locus),ref.length=rep(0,opt$max.locus),motif=rep(0,opt$max.locus),motif.length=rep(0,opt$max.locus))

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

	error.info$em.noise.rate[chunk.start.locus:chunk.stop.locus]<-em.geno$locus.error.rates
	error.info$ref.length[chunk.start.locus:chunk.stop.locus]<-chunk.locus.info$ref.length[chunk.locus.info$allele.no == 0]
	error.info$motif[chunk.start.locus:chunk.stop.locus]<-as.character(chunk.locus.info$unit[chunk.locus.info$allele.no == 0])
}
error.info$motif.length<-nchar(as.character(error.info$motif))

h.step<-0.05
error.info$mid<-round(((log10(error.info$em.noise.rate) %/% h.step) * h.step),2)

# png(sprintf("%snoiseBoxplot.png",opt$fig.directory),height=1200,width=2000)
# par(cex=5,pch=16)
# layout(matrix(1:6,nrow=2,byrow=TRUE))
# for(i in 1:6)
# {
#   boxplot(log10(em.noise.rate) ~ (ref.length %/% i),data=error.info,subset=motif.length == i,ylab="EM-deduced error rate",xlab="Reference microsatellite # repeats",
#           cex=0.5,yaxt="n",cex.axis=1.5,cex.main=2,cex.lab=2,
#           main=sprintf("EM-deduced error rate as related to reference microsatellite length\n%d bp microsatellites only",i),
#           ylim=c(-5,0),col=RColorBrewer::brewer.pal(6,"Set1")[i])
#   axis(side=2,at=-5:0,labels=10^(-5:0))
# }
# garbage<-dev.off()

# png(sprintf("%snoiseHistByMotifLength.png",opt$fig.directory),height=1200,width=2000)
# par(cex=5,pch=16)
# layout(matrix(1:6,nrow=2,byrow=TRUE))
# y.step=0.1
# for(i in 1:6)
# {
#    h.step<-0.05
#    h.break.min<-(min(log10(error.info$em.noise.rate[error.info$motif.length == i])) %/% h.step) - h.step
#    h.disp.step<-1
#    h.disp.min<- ifelse(h.break.min < -6, -6, ((h.break.min %/% h.disp.step) * h.disp.step) - h.disp.step)
#    h<-hist(log10(error.info$em.noise.rate[error.info$motif.length == i]),breaks=seq(h.break.min,0,h.step),plot=FALSE)
#    h$counts<-h$counts/sum(h$counts)
#    y.max<-((max(h$counts) %/% y.step) * y.step) + y.step
#    plot(h,border=NA,col=RColorBrewer::brewer.pal(6,"Dark2")[i],ylab="Frequency",xlab="Locus error rate",main=sprintf("Locus error rate\n%d bp motifs",i),xaxt="n",yaxt="n",xlim=c(h.disp.min,0),ylim=c(0,y.max))
#    axis(side=1,at=seq(h.disp.min,0,1),labels=10^(seq(h.disp.min,0,1)))
#    axis(side=2,at=seq(0,y.max,y.step))
# }
# garbage<-dev.off()

mids<-round(seq(min(error.info$mid),0,h.step),2)
count.by.motif.length<-t(laply(mids, function(i) laply(1:6, function(x) sum(error.info$motif.length == x & error.info$mid == i))))

min.hist<-mids[which(colSums(count.by.motif.length) > 0.0005*sum(count.by.motif.length))[1]] %/% 1
min.hist.cind<-which(mids == min.hist)

y.max<-((max(colSums(count.by.motif.length)/sum(count.by.motif.length)) %/% 0.05) * 0.05) + 0.05

ml.col<-RColorBrewer::brewer.pal(6,"Set2")
# png(sprintf("%s/logLocusErrorWithMotifLength.png",opt$fig.directory),height=1600,width=1600)
# par(cex=4)
pdf(sprintf("%s/logLocusErrorWithMotifLength.pdf",opt$fig.directory),height=8,width=8)
bp<-barplot(count.by.motif.length[,min.hist.cind:ncol(count.by.motif.length)]/sum(count.by.motif.length),col=ml.col,border=NA,ylim=c(0,y.max),yaxt="n",
			xlab="Locus error rate",ylab="% of well-covered loci",main="Locus error rate with motif length")
axis(side=2,at=seq(0,y.max,length.out=6),labels=seq(0,y.max*100,length.out=6),las=2)
axis(side=1,at=bp[seq(1,length(bp),20)],labels=10^mids[min.hist.cind+seq(0,length(bp),20)])
legend(x="topright",title="Motif length (bp)",legend=1:6,fill=ml.col,border=NA,bty="n")
garbage<-dev.off()
cat("Wrote graph",sprintf("%s/logLocusErrorWithMotifLength.png",opt$fig.directory),"\n")

rl.col<-colorRampPalette(RColorBrewer::brewer.pal(8,"Set2"),bias=1.5)

for(l in 1:6)
{
	lengths<-sort(unique(error.info$ref.length[error.info$motif.length == l]))
	count.by.bin<-t(laply(mids,  function(i) laply(lengths, function(x) sum(error.info$ref.length == x & error.info$mid == i & error.info$motif.length == l))))

	good.rows<-which(rowSums(count.by.bin) > 0.001*sum(count.by.bin))
	min.hist<-mids[which(colSums(count.by.bin) > 0.0001*sum(count.by.bin))[1]] %/% 1
	min.hist<-ifelse(min.hist < -6, -6, min.hist)
	min.hist.cind<-which(mids == min.hist)
	y.max<-((max(colSums(count.by.bin)/sum(count.by.bin)) %/% 0.05) * 0.05) + 0.05

	graph.name<-sprintf("%s/%dbpLogLocusErrorWithRefLength.pdf",opt$fig.directory,l)
	pdf(graph.name,height=8,width=8)
	bp<-barplot(count.by.bin[good.rows,min.hist.cind:ncol(count.by.bin)]/sum(count.by.bin),col=rl.col(length(good.rows)),border=NA,ylim=c(0,y.max),yaxt="n",
				xlab="Locus error rate",ylab=sprintf("%% of well-covered loci with %d bp motifs", l),main=sprintf("Locus error rate with tract length for microsatellites with %d bp motifs",l))
	axis(side=2,at=seq(0,y.max,length.out=6),labels=seq(0,y.max*100,length.out=6),las=2)
	axis(side=1,at=bp[seq(1,length(bp),20)],labels=10^mids[min.hist.cind+seq(0,length(bp),20)])
	legend(x="topright",legend=lengths[good.rows],fill=rl.col(length(good.rows)),border=NA,title="tract length (bp)",bty="n")
	garbage<-dev.off()
	cat("Wrote graph",graph.name,"\n")

# 	num.repeats<-sort(unique(error.info$ref.length[error.info$motif.length == l] %/% l))
# 	count.by.bin<-t(laply(mids,  function(i) laply(num.repeats, function(x) sum((error.info$ref.length %/% l) == x & error.info$mid == i & error.info$motif.length == l))))
# 
# 	good.rows<-which(rowSums(count.by.bin) > 0.001*sum(count.by.bin))
# 	min.hist<-mids[which(colSums(count.by.bin) > 0.0001*sum(count.by.bin))[1]] %/% 1
# 	min.hist<-ifelse(min.hist < -6, -6, min.hist)
# 	min.hist.cind<-which(mids == min.hist)
# 	y.max<-((max(colSums(count.by.bin)/sum(count.by.bin)) %/% 0.05) * 0.05) + 0.05
# 
# 	graph.name<-sprintf("%s/%dbpLogLocusErrorWithRefNumRepeats.png",opt$fig.directory,l)
# 	png(graph.name,height=1600,width=1600)
# 	par(cex=4)
# 	bp<-barplot(count.by.bin[good.rows,min.hist.cind:ncol(count.by.bin)]/sum(count.by.bin),col=rl.col(length(good.rows)),border=NA,ylim=c(0,y.max),yaxt="n",
# 				xlab="Locus error rate",ylab="Frequency",main=sprintf("Locus error rate for %d bp microsatellites by reference # repeats\n%s total loci",l,formatC(sum(count.by.bin),big.mark=",")))
# 	axis(side=2,at=seq(0,y.max,length.out=6))
# 	axis(side=1,at=bp[seq(1,length(bp),20)],labels=10^mids[min.hist.cind+seq(0,length(bp),20)])
# 	legend(x="topright",legend=num.repeats[good.rows],fill=rl.col(length(good.rows)),border=NA,title="# repeats")
# 	garbage<-dev.off()
# 	cat("Wrote graph",graph.name,"\n")
}
