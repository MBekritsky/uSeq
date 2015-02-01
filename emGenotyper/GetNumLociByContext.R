#!/usr/bin/Rscript

source("/data/safe/bekritsk/simons/scripts/genotyper/emGenotyper/src/dataLoadingFunctions.R")
source("/data/safe/bekritsk/simons/scripts/genotyper/emGenotyper/src/annotationFunctions.R")
source("/data/safe/bekritsk/simons/scripts/genotyper/emGenotyper/src/figureFunctions.R")

dirname="/data/safe/bekritsk/simons/Exome/workingSet/06182013_completeWiglerSSCQuads/"

locus.info.file=paste(dirname,"allele_matrix_info.txt",sep="")
allele.info.file=paste(dirname,"allele_matrix.txt",sep="")
family.info.file=paste(dirname,"person_info.txt",sep="")

#load families
locus.info<-load.locus.info(locus.info.file)
#Exclude sex chromosomes--X and Y can be corrected for copy number (MT copy number is unknown), but genotyper 
#also needs to accomodate single allele genotypes for X and Y, depending on gender.  MT could be single allele genotype, 
#but could have many more than 2 genotypes as well
locus.info<-remove.chromosomes(locus.info,c("chrX","chrY"))

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

precedence<-function(anno)
{
  ifelse(sum(is.na(anno['exon'])) == 0, return("Exon"),
         ifelse(sum(is.na(anno['mirna'])) == 0, return("miRNA"),
                ifelse(sum(is.na(anno['utr'])) == 0, return("UTR"),
                       ifelse(sum(is.na(anno['intron'])) == 0, return("Intron"),
                              return("Intergenic")))))
}

main.annos<-apply(annotations,1,precedence)
by.context<-as.data.frame(table(main.annos))
by.context$pct<-by.context$Freq/sum(by.context$Freq)

# png(sprintf("%s/SummaryFigures/coveredLociByContextLogSpace.png",dirname),height=1200,width=1200)
# par(cex=2)
# bp<-barplot(log10(by.context$Freq),col=addAlphaToHex(RColorBrewer::brewer.pal(5,"Dark2"),0.7),ylim=c(0,5),yaxt="n",border=NA,
#             main="Well-covered SSC microsatellites by context",ylab="Count",font=2,font.lab=2)
# axis(side=1,at=bp,labels=by.context$main.anno,lwd=0,font=2)
# axis(side=2,at=0:5,labels=10^(0:5))
# text(x=bp,y=log10(by.context$Freq)+0.1,labels=by.context$Freq,font=2)
# dev.off()

png(sprintf("%s/SummaryFigures/coveredLociByContext.png",dirname),height=1400,width=2000)
par(cex=4)
y.max<-((max(by.context$Freq) %/% 1e4) * 1e4) + 1e4
bp<-barplot(by.context$Freq,col=addAlphaToHex(RColorBrewer::brewer.pal(5,"Dark2"),0.7),ylim=c(0,y.max),yaxt="n",border=NA,
            main="Well-covered SSC microsatellites by context",font=2,font.lab=2)
axis(side=1,at=bp,labels=by.context$main.anno,lwd=0,font=2)
# axis(side=2,at=seq(0,y.max,length.out=6))
text(x=bp,y=by.context$Freq+1.5e3,labels=sapply(by.context$Freq,function(x) prettyNum(x,big.mark=",")),font=2)
dev.off()
