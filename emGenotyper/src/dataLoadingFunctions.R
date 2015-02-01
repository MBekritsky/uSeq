#Loads information about each locus being analyzed
load.locus.info<-function(filename)
{
  info<-read.table(filename)
  colnames(info)<-c("chr","pos","unit","ref.length","pop.count","top.count","sum.count","allele.no","num.lengths")
  
  locus.inds<-which(info$allele.no == 0)
  info$allele.zero.ind<-rep(locus.inds,info$num.lengths[locus.inds]+1)
  info$locus.row.ind<-rep(1:length(locus.inds),info$num.lengths[locus.inds]+1)
  
  return(info)  
}

load.locus.total.cov <- function(locus.info.file,allele.info.file,pop.size,exclude.XY=TRUE,exclude.MT=TRUE,verbose=1)
{
  locus.info <- read.table(locus.info.file, sep="\t")
  if(!exclude.XY & !exclude.MT)
  {
    locus.lines <- which(locus.info[,8] == 0)
  } else if(exclude.XY & !exclude.MT)
  {
    locus.lines <- which(locus.info[,8] == 0 & !(locus.info[,1] %in% c("chrX","chrY")))
  } else if(!exclude.XY & exclude.MT)
  {
    locus.lines <- which(locus.info[,8] == 0 & !(locus.info[,1] %in% "chrMT"))
  } else
  {
    locus.lines <- which(locus.info[,8] == 0 & !(locus.info[,1] %in% c("chrX","chrY","chrMT")))
  }
  
  rm(locus.info)
  
  if(verbose >= 1) cat("Identified locus lines\n")
  
  locus.cov <- matrix(NA,ncol=pop.size,nrow=length(locus.lines))
  
  allele.info <- file(allele.info.file, "r")
  
  allele.line.no<-1
  locus.no<-1
  while(length(allele.line <- readLines(allele.info, n=1)) > 0)
  {
    if(allele.line.no %in% locus.lines)
    {
      allele.parsed.line<-as.integer(unlist(strsplit(allele.line,split=',',fixed=TRUE)))
      locus.cov[locus.no,] <- allele.parsed.line
      locus.no <- locus.no + 1
    }
    allele.line.no <- allele.line.no + 1
  }
  
  close(allele.info)
  return(locus.cov)
}

update.locus.info.locators<-function(info)
{
  locus.inds<-which(info$allele.no == 0)
  info$allele.zero.ind<-rep(locus.inds,info$num.lengths[locus.inds]+1)
  info$locus.row.ind<-rep(1:length(locus.inds),info$num.lengths[locus.inds]+1)
  return(info)
}

remove.chromosomes<-function(locus.info,chr.list)
{
  if(sum(chr.list %in% locus.info$chr) > 0)
  {
    locus.info<-locus.info[-which(locus.info$chr %in% chr.list),]
  }
  return(locus.info)
}

#Loads information about each person being analyzed
load.person.info<-function(filename)
{
  info<-read.table(filename)
  colnames(info)<-c("individual.id","family.id","relation")
  info$rel.code<-rep(0,nrow(info))
  info$rel.code[which(info$relation %in% "mother")]<-1
  info$rel.code[which(info$relation %in% "father")]<-2
  info$rel.code[which(info$relation %in% "self")]<-3
  info$rel.code[which(info$relation %in% "sibling")]<-4
  return(info)  
}

add.genders<-function(person.info,report.file)
{
  report.table<-read.table(report.file,header=TRUE,sep="\t", fill = TRUE)
  genders<-apply(person.info,1,function(x) as.character(report.table$sample.Gender[which(as.character(report.table$sample.sample_id) == as.character(x[1]))[1]]))
  person.info$gender<-genders
  return(person.info)
}

load.sample.set.report<-function(filename)
{
  sample.set.report<-read.table(filename,sep="\t",header=TRUE, fill = TRUE)
  return(sample.set.report)
}

get.non.whole.blood.family.members<-function(sample.set.report,person.info,verbose=1)
{
  non.whole.blood.people<-grep("^SSC",unique(as.character(sample.set.report$sample.person_id[grep("whole blood",sample.set.report$sample.sample_type,fixed=TRUE,invert=TRUE)])),
                               value=TRUE)
  if(verbose >= 1) cat("The following people have DNA obtained from saliva:",paste(non.whole.blood.people,collapse=", "),"\n")
  non.whole.blood.families<-grep("^auSSC",unique(as.character(sample.set.report$sample.FamilyId[grep("whole blood",sample.set.report$sample.sample_type,
                                                                                                     fixed=TRUE,invert=TRUE)])),value=TRUE)
  if(verbose >= 1) cat("The following families will be excluded since a member has DNA sequence that was not obtained from whole blood DNA:",paste(non.whole.blood.families,collapse=", "),"\n")
  non.whole.blood.families<-which(person.info$family.id %in% non.whole.blood.families)
  return(list(inds=non.whole.blood.families,family.ids=non.whole.blood.families))
}

exclude.families<-function(allele.info,person.info,excluded.members)
{
  excluded.families<-unique(person.info$family.id[excluded.members])
  excluded.indices<-which(person.info$family.id %in% excluded.families)
  print("Excluding the following families with overdispersed members: ")
  print(paste(as.character(excluded.families),collapse=", "))
  
  allele.info$alleles<-allele.info$alleles[,-excluded.indices]
  allele.info$locus.cov<-allele.info$locus.cov[,-excluded.indices]
  allele.info$log.locus.cov<-allele.info$log.locus.cov[,-excluded.indices]
  
  person.info<-person.info[-excluded.indices,]
  
  return(list(person.info=person.info,allele.info=allele.info))
}

#Loads coverage matrix
load.alleles<-function(filename,locus.info,exclude.XY=TRUE,exclude.MT=FALSE,start.locus=1,num.alleles=-1)
{
  alleles<-read.table(filename,sep=",",skip=(start.locus-1),nrows=num.alleles)
  alleles<-as.matrix(alleles)
  if(exclude.XY)
  {
    alleles<-alleles[which(!(locus.info$chr %in% c("chrX","chrY"))),]
    locus.info<-remove.chromosomes(locus.info,c("chrX","chrY"))
  }
  if(exclude.MT)
  {
    alleles<-alleles[which(locus.info$chr != "chrMT"),]
    locus.info<-remove.chromosomes(locus.info,"chrMT")
  }
  return.list=list(alleles=alleles,num.loci=sum(locus.info$allele.no == 0))
  return(return.list)
}

add.distance.to.exon<-function(locus.info,distance.bed)
{
  distance.info<-read.table(distance.bed)
  names(distance.info)<-c("ms.chr","ms.start","ms.stop","exon.chr","exon.start","exon.stop","exon.name","distance")
  locus.info$distance<-distance.info$distance[which(distance.info$ms.chr == locus.info$chr & distance.info$ms.start == locus.info$pos)]
  return(locus.info)
}

exclude.people.from.person.df<-function(person.info,inds)
{
  return(person.info[-inds,])
}

exclude.people.from.locus.df<-function(locus.info,allele.info,inds)
{
  total.excluded.cov<-rowSums(allele.info$alleles[,inds])
  total.excluded.with.cov<-rowSums(allele.info$alleles[,inds] > 0)
  locus.info$pop.count<-locus.info$pop.count-total.excluded.with.cov
  locus.info$sum.count<-locus.info$sum.count-total.excluded.cov
  return(locus.info)
}

exclude.people.from.allele.df<-function(allele.info,inds)
{
  allele.info$alleles<-allele.info$alleles[,-inds]
  allele.info$locus.cov<-allele.info$locus.cov[,-inds]
  allele.info$log.locus.cov<-allele.info$log.locus.cov[,-inds]
  return(allele.info)
}

exclude.people.from.exp.coverage.matrix<-function(coverage.estimator,inds)
{
  coverage.estimator$exp.cov<-coverage.estimator$exp.cov[,-inds]
  coverage.estimator$cor<-coverage.estimator$cor[-inds]
  return(coverage.estimator)
}

add.num.called<-function(f.geno,f.locus.info,start.locus=1,stop.locus=NULL)
{
  f.locus.info$num.called<-NA
  
  if(is.null(stop.locus))
  {
    stop.locus<-dim(f.geno$genotypes)[1]
  }
  for(i in start.locus:stop.locus)
  {
    counts<-tabulate(f.geno$genotypes[i,,1:2][f.geno$genotypes[i,,1:2]  > 0])
    for(j in 8:length(counts))
    {
      ind<-which(f.locus.info$locus.row.ind == i & f.locus.info$allele.no == j)
      if(length(ind) > 0)
      {
        f.locus.info$num.called[ind]<-counts[j]
      }
    }
  }
  f.locus.info$num.called[is.na(f.locus.info$num.called) & f.locus.info$allele.no > 0]<-0
  return(f.locus.info)
}

add.allele.biases<-function(f.geno,f.locus.info,start.locus=1,stop.locus=NULL)
{
  f.locus.info$bias<-NA
  
  if(is.null(stop.locus))
  {
    stop.locus<-dim(f.geno$genotypes)[1]
  }
  
  for(i in start.locus:stop.locus)
  {
    for(j in 1:nrow(f.geno$allele.biases[[i]]))
    {
      ind<-which(f.locus.info$locus.row.ind == i & f.locus.info$allele.no == f.geno$allele.biases[[i]][j,1])
      f.locus.info$bias[ind]<-f.geno$allele.biases[[i]][j,2]
    }
  }
  return(f.locus.info)
}