conditional.make.dir<-function(dirnames, verbose = 0) {
  # Check if a set of directories exist.  If any of them
  # do not, create them
  
  if(missing(dirnames)) stop("dirnames is missing")
  if(length(dirnames) == 0) warning("dirnames is empty")
  
  for(dir in dirnames) {
    if(!file.exists(dir)) {
      if(verbose > 0) cat("Creating", dir, "\n",sep=" ")
      dir.create(dir)
    }
  }
}

######### OBSOLETE CODE #########

# assign.exon.status<-function(allele.info,locus.info)
# {
#   ms.exon.info<-read.table("~/Desktop/hg19.exon.bed")
#   names(ms.exon.info)<-c("chr","start","stop")
#   in.exons<-0
#   for(i in 1:allele.info$num.loci)
#   {
#     exon.ind<-which(ms.exon.info$chr == as.character(locus.info$chr[i]) & ms.exon.info$start == locus.info$pos[i])
#     if(length(exon.ind)>0)
#     {
#       locus.info$in.exon[i]<-TRUE
#     } else
#     {
#       locus.info$in.exon[i]<-FALSE
#     }
#   }
#   rm(ms.exon.info)
#   return(locus.info) 
# }
# 
# locus.genotype.breakdown<-function(f.geno,person.info,index)
# {
#  breakdown<-data.frame(unique(em.geno$genotypes[index,,1:2]))
#  names(breakdown)<-c("allele.1","allele.2")
#  print(breakdown)
#  gender.cat<-c("F","M","F","M","F","M")
#  rel.cat<-c("mother","father","self","self","sibling","sibling")
#  
#  for(i in 1:2)
#  {
#    for(j in 1:nrow(breakdown))
#    {
#      breakdown[[rel.cat[i]]][j]<-sum(as.character(person.info$relation) == as.character(rel.cat[i]) & 
#                                     f.geno$genotypes[index,,1] == breakdown$allele.1[j] & 
#                                     f.geno$genotypes[index,,2] == breakdown$allele.2[j]) / 
#                                     sum(as.character(person.info$relation) == as.character(rel.cat[i]))
#    }
#  }
# 
#  for(i in 3:6)
#  {
#    for(j in 1:nrow(breakdown))
#    {
#      cat.name<-sprintf("%s-%s",rel.cat[i],gender.cat[i])
#      breakdown[[cat.name]][j]<-length(which(as.character(person.info$relation) == as.character(rel.cat[i]) & 
#                                    as.character(person.info$gender) == as.character(gender.cat[i]) & 
#                                        f.geno$genotypes[index,,1] == breakdown$allele.1[j] & 
#                                        f.geno$genotypes[index,,2] == breakdown$allele.2[j])) / 
#        length(which(as.character(person.info$relation) == as.character(rel.cat[i]) & as.character(person.info$gender) == as.character(gender.cat[i])))
#    }
#  }
#  
#  return(breakdown)
# }
