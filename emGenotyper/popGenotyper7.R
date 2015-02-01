source("/data/safe/bekritsk/simons/scripts/genotyper/emGenotyper/src/mendelScoreFunctions2.R")
source("/data/safe/bekritsk/simons/scripts/genotyper/emGenotyper/src/modelingFunctions3.R")
source("/data/safe/bekritsk/simons/scripts/genotyper/emGenotyper/src/genotyperFunctions2.R")
source("/data/safe/bekritsk/simons/scripts/genotyper/emGenotyper/src/figureFunctions.R")
source("/data/safe/bekritsk/simons/scripts/genotyper/emGenotyper/src/dataLoadingFunctions.R")
source("/data/safe/bekritsk/simons/scripts/genotyper/emGenotyper/src/annotationFunctions.R")
source("/data/safe/bekritsk/simons/scripts/genotyper/emGenotyper/src/helperFunctions.R")

dirname="/data/safe/bekritsk/simons/Exome/workingSet/06182013_completeWiglerSSCQuads/subsample/"

locus.info.file=paste(dirname,"allele_matrix_info.txt",sep="")
allele.info.file=paste(dirname,"allele_matrix.txt",sep="")
family.info.file=paste(dirname,"person_info.txt",sep="")
report.file="/data/safe/bekritsk/simons/Exome/workingSet/06182013_completeWiglerSSCQuads/report_sampleSet_20130611.txt"


#load families
locus.info<-load.locus.info(locus.info.file)
#Exclude sex chromosomes--X and Y can be corrected for copy number (MT copy number is unknown), but genotyper 
#also needs to accomodate single allele genotypes for X and Y, depending on gender.  MT could be single allele genotype, 
#but could have many more than 2 genotypes as well
allele.info<-load.alleles(allele.info.file,locus.info,exclude.XY=TRUE,exclude.MT=FALSE)
locus.info<-remove.chromosomes(locus.info,c("chrX","chrY"))
locus.info<-add.locus.idxs.to.locus.info(locus.info,allele.info)
person.info<-load.person.info(family.info.file)

exon.anno.file<-"/mnt/wigclust4/home/bekritsk/genomes/hg19/hg19BedFiles/hg19.ms.ccds.rg.info.exons.merged.bed"
intron.anno.file<-"/mnt/wigclust4/home/bekritsk/genomes/hg19/hg19BedFiles/hg19.ms.ccds.rg.info.introns.merged.bed"
mirna.anno.file<-"/mnt/wigclust4/home/bekritsk/genomes/hg19/hg19BedFiles/hg19.ms.miRNA.merged.bed"
utr.anno.file<-"/mnt/wigclust4/home/bekritsk/genomes/hg19/hg19BedFiles/hg19.ms.ccds.rg.info.utr.merged.bed"

gene.id.file<-"/mnt/wigclust4/home/bekritsk/genomes/hg19/hg19GeneName.ccds.rg.geneID.txt"

gene.ids<-load.gene.ids(gene.id.file)
ptm<-proc.time()
annotations<-get.annotations(exon.file=exon.anno.file,intron.file=intron.anno.file,mirna.file=mirna.anno.file,utr.file=utr.anno.file,
                            locus.info=locus.info)
print(proc.time() - ptm)
save(annotations,file=sprintf("%sannotations.RData",dirname))

ptm<-proc.time()
coverage.estimator<-linear.coverage.model(allele.info$locus.cov,print.dir=dirname,people=person.info)
print(proc.time()-ptm)
save(coverage.estimator,file=sprintf("%scoverageEstimator.RData",dirname))

#plot coverage per person
total.cov<-pop.total.coverage.hist(allele.info$locus.cov,dirname)
mean.cov<-pop.mean.coverage.hist(allele.info$locus.cov,dirname)

#find individuals with low correlation between expected and observed allele coverage or very low observed coverage and remove them
problem.people<-unique(as.integer(c(which(total.cov < (mean(total.cov) - (2*sd(total.cov)))),which(coverage.estimator$cor < 0.8))))
problem.families<-which(person.info$family.id %in% person.info$family.id[problem.people])

#exclude families from analysis where a member has low correlation between expected and observed allele coverage or has low total coverage
person.info<-exclude.people.from.person.df(person.info,problem.families)
locus.info<-exclude.people.from.locus.df(locus.info,allele.info,problem.families)
allele.info<-exclude.people.from.allele.df(allele.info,problem.families)
coverage.estimator<-exclude.people.from.exp.coverage.matrix(coverage.estimator,problem.families)
person.info<-add.genders(person.info,report.file)

pop<-ncol(allele.info$alleles)
n.fams<-pop/4
by.fams<-get.fam.mat(person.info)
mom.and.pop<-which(person.info$relation %in% c("mother","father"))

em.geno<-call.genotypes(allele.info,locus.info,coverage.estimator,pop)
save(em.geno,file=sprintf("%semGenotpyes.RData",dirname))

locus.info<-add.num.called(em.geno,locus.info,stop.locus=1000)
locus.info<-add.allele.biases(em.geno,locus.info,stop.locus=1000)
save(locus.info,file=sprintf("%slocusInfo.RData",dirname))

graphBiasInfo(dirname,locus.info)
graphGenotypeOverviewStats(dirname,em.geno)

max.alleles<-6 #number of most common alleles within family to consider when calculating the Mendel score
mendel.ind<-get.mendel.trio.indices(max.alleles)
n.ind<-nrow(mendel.ind) #total number of genotype combinations (of 4^max.alleles) that are Mendelian
fam.genotypes<-get.mendel.scores(allele.info=allele.info,locus.info=locus.info,geno.list=em.geno,
                                n.fams=n.fams,by.fam=by.fams,pop.size=pop,coverage.model=coverage.estimator,mendel.ind=mendel.ind)
save(fam.genotypes,file=sprintf("%sfamilyGenotypes.RData",dirname))

low.mendel.threshold<-5e-3
max.error.rate<-1.5e-1
min.exp.cov<-5

include.na.noise<-FALSE
include.na<-FALSE
pro.low.mendel<-length(which(fam.genotypes[,,10] <= low.mendel.threshold))
pro.low.mendel.indices<-which(fam.genotypes[,,10] <= low.mendel.threshold, arr.ind=T)
pro.low.mendel.indices<-cbind(pro.low.mendel.indices,1)
sib.low.mendel<-length(which(fam.genotypes[,,11] <= low.mendel.threshold))
sib.low.mendel.indices<-which(fam.genotypes[,,11] <= low.mendel.threshold, arr.ind=T)
sib.low.mendel.indices<-cbind(sib.low.mendel.indices,2)
low.mendel.indices<-rbind(sib.low.mendel.indices,pro.low.mendel.indices)
low.mendel.indices<-low.mendel.indices[order(low.mendel.indices[,1],low.mendel.indices[,2],low.mendel.indices[,3]),]

tot.low.mendel<-vector()
denovo.loci<-unique(low.mendel.indices[,1])
for(i in denovo.loci)
{
  loci.indices<-which(low.mendel.indices[,1] == i)
  denovo.fams.at.loci<-unique(low.mendel.indices[loci.indices,2])
  for(j in denovo.fams.at.loci)
  {
    dn.fam.at.locus<-which(low.mendel.indices[,1] == i & low.mendel.indices[,2] == j)
    if(length(dn.fam.at.locus) == 1)
    {
      tot.low.mendel<-rbind(tot.low.mendel,low.mendel.indices[dn.fam.at.locus,])
    } else
    {
      tot.low.mendel<-rbind(tot.low.mendel,c(low.mendel.indices[dn.fam.at.locus[1],1:2],3))
    }
  }
}

denovo.info<-array(NA,dim=c(nrow(low.mendel.indices),2))
for(i in 1:nrow(tot.low.mendel))
{
  locus<-tot.low.mendel[i,1]
  family<-tot.low.mendel[i,2]
  who.dn<-tot.low.mendel[i,3]

  num.na<-sum(em.geno$genotypes[rbind(cbind(locus,by.fams[family,],1),cbind(locus,by.fams[family,],2))] == -1) + 
          sum(fam.genotypes[locus,family,1:8] == -1)
  good.ex.cov<-sum(coverage.estimator$exp.cov[locus,by.fams[family,]] > min.exp.cov)
  if(num.na == 0 & em.geno$locus.error.rates[locus] <= max.error.rate & good.ex.cov == 4)
  {
    if(who.dn == 1)
    {
      denovo.info[i,1]<-fam.genotypes[locus,family,10]
    } else if (who.dn == 2)
    {
      denovo.info[i,1]<-fam.genotypes[locus,family,11]
    } else if (who.dn == 3)
    {
      denovo.info[i,1]<-mean(fam.genotypes[locus,family,10:11])
    }
    denovo.info[i,2]<-i
  }
}
denovo.info<-na.omit(denovo.info)
denovo.info<-denovo.info[order(denovo.info[,1]),]

x11();plot(-10*log10(denovo.info[,1]),main=sprintf("Top %d de novo candidate scores",nrow(denovo.info)),pch=16,ylab="log10 denovo score")

denovo.dir<-sprintf("%s/denovos/",dirname)
fig.dirs<-create.dn.figure.dirs(denovo.dir)

total.denovo<-0
dn.df<-data.frame(who=rep(NA,nrow(denovo.info)),type=rep(NA,nrow(denovo.info)),
                  context=rep(NA,nrow(denovo.info)),gender=rep(NA,nrow(denovo.info)),
				  geno.freq=rep(NA,nrow(denovo.info)),par.geno.freq=rep(NA,nrow(denovo.info)))
for(i in 1:nrow(denovo.info))
{
  total.denovo<-total.denovo+1
  denovo.ind<-tot.low.mendel[denovo.info[i,2],]
  fam.ind<-by.fams[denovo.ind[2],]
  dn.locus<-denovo.ind[1]
  dn.family<-denovo.ind[2]
  dn.person<-denovo.ind[3]
  
  locus.info.index<-allele.info$locus.row.indices[dn.locus]
  
  dn.type<-NULL
  ifelse(dn.person == 1, dn.indices<-5:6,
  ifelse(dn.person == 2, dn.indices<-7:8,
                         dn.indices<-5:8))

  dn.anno<-get.locus.annotations(dn.locus,annotations,gene.ids)
  context<-annotation.precedence(dn.anno)
  
  ifelse(dn.person == 1, dn.df$who[i]<-"proband",
  ifelse(dn.person == 2, dn.df$who[i]<-"sibling",
                         dn.df$who[i]<-"both"))
  dn.df$context[i]<-context
  ifelse(dn.person == 1, dn.df$gender[i]<-as.character(person.info$gender[by.fams[dn.family,3]]),
  ifelse(dn.person == 2, dn.df$gender[i]<-as.character(person.info$gender[by.fams[dn.family,4]]),
                         dn.df$gender[i]<-sprintf("%s/%s",as.character(person.info$gender[by.fams[dn.family,3]]),
                                                  as.character(person.info$gender[by.fams[dn.family,4]]))))
												  
  ifelse(dn.person == 1, dn.df$geno.freq[i]<-sum(em.geno$genotypes[dn.locus,,1] == em.geno$genotypes[dn.locus,by.fams[dn.family,3],1] & em.geno$genotypes[dn.locus,,2] == em.geno$genotypes[dn.locus,by.fams[dn.family,3],2]),
  ifelse(dn.person == 2, dn.df$geno.freq[i]<-sum(em.geno$genotypes[dn.locus,,1] == em.geno$genotypes[dn.locus,by.fams[dn.family,4],1] & em.geno$genotypes[dn.locus,,2] == em.geno$genotypes[dn.locus,by.fams[dn.family,4],2]),
                         dn.df$geno.freq[i]<-length(unique(c(which(em.geno$genotypes[dn.locus,,1] == em.geno$genotypes[dn.locus,by.fams[dn.family,3],1] & 
						 											em.geno$genotypes[dn.locus,,2] == em.geno$genotypes[dn.locus,by.fams[dn.family,3],2]),
						 					 				 which(em.geno$genotypes[dn.locus,,1] == em.geno$genotypes[dn.locus,by.fams[dn.family,4],1] & 
															 		em.geno$genotypes[dn.locus,,2] == em.geno$genotypes[dn.locus,by.fams[dn.family,4],2]))))))
  
  
  ifelse(dn.person == 1, dn.df$par.geno.freq[i]<-sum(em.geno$genotypes[dn.locus,mom.and.pop,1] == em.geno$genotypes[dn.locus,by.fams[dn.family,3],1] & em.geno$genotypes[dn.locus,mom.and.pop,2] == em.geno$genotypes[dn.locus,by.fams[dn.family,3],2]),
  ifelse(dn.person == 2, dn.df$par.geno.freq[i]<-sum(em.geno$genotypes[dn.locus,mom.and.pop,1] == em.geno$genotypes[dn.locus,by.fams[dn.family,4],1] & em.geno$genotypes[dn.locus,mom.and.pop,2] == em.geno$genotypes[dn.locus,by.fams[dn.family,4],2]),
                         dn.df$par.geno.freq[i]<-length(unique(c(which(em.geno$genotypes[dn.locus,mom.and.pop,1] == em.geno$genotypes[dn.locus,by.fams[dn.family,3],1] & 
						 											em.geno$genotypes[dn.locus,mom.and.pop,2] == em.geno$genotypes[dn.locus,by.fams[dn.family,3],2]),
									 					  	 	 which(em.geno$genotypes[dn.locus,mom.and.pop,1] == em.geno$genotypes[dn.locus,by.fams[dn.family,4],1] & 
																 	em.geno$genotypes[dn.locus,mom.and.pop,2] == em.geno$genotypes[dn.locus,by.fams[dn.family,4],2]))))))
  
  dn.alleles<-unique(fam.genotypes[dn.locus,dn.family,dn.indices])
  graph.dir<-NULL
  if(sum(dn.alleles %in% fam.genotypes[dn.locus,dn.family,1:4]) == length(dn.alleles)) #all alleles in Mendel violation child are in parents, omission violation
  {
    dn.df$type[i]<-"omission"
    ifelse(dn.person == 1, graph.dir<-fig.dirs[["omissions"]][["proband"]][[context]],
    ifelse(dn.person == 2, graph.dir<-fig.dirs[["omissions"]][["sibling"]][[context]],
                           graph.dir<-fig.dirs[["omissions"]][["both"]][[context]]))
  } else #new allele in Mendel violation child, commission violation
  {
    dn.df$type[i]<-"commission"
    ifelse(dn.person == 1, graph.dir<-fig.dirs[["commissions"]][["proband"]][[context]],
    ifelse(dn.person == 2, graph.dir<-fig.dirs[["commissions"]][["sibling"]][[context]],
                           graph.dir<-fig.dirs[["commissions"]][["both"]][[context]]))
  }
#  plot.family(locus.info,allele.info,em.geno$genotypes,dn.locus,person.info,family.ind=by.fams[dn.family,],fam.genotypes=fam.genotypes[dn.locus,dn.family,],
#              who.dn=dn.person,print.dir=graph.dir,anno=dn.anno)
#  plot.pop.image(locus.info,em.geno$genotypes,denovo.ind[1],print.dir=graph.dir,fam.genotypes=fam.genotypes[dn.locus,dn.family,],who.dn=denovo.ind[3],anno=dn.anno)
#  plot.dn.eo.scatter(allele.info,locus.info,dn.locus,dn.family,by.fams,em.geno$genotypes,coverage.estimator,who.dn=dn.person,fam.genotypes=fam.genotypes[dn.locus,dn.family,],
#                     print.dir=graph.dir,anno=dn.anno)
#  plot.dn.cov.scatter(allele.info,locus.info,dn.locus,dn.family,by.fams,em.geno$genotypes,who.dn=dn.person,fam.genotypes=fam.genotypes[dn.locus,dn.family,],
#                      print.dir=graph.dir,anno=dn.anno)
}

table(dn.df$gender)
