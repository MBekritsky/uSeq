low.obedience.scores<-read.table("~/Google Drive/Microsatellites/UnsortedMendelScores.txt",header=TRUE)

excluded.by.cm.or.ebv<-unique(c("auSSC12762", "auSSC12905", "auSSC11101", "auSSC11497","auSSC12353", "auSSC11292", "auSSC14236", "auSSC14698", "auSSC11300", "auSSC11101", "auSSC12905", "auSSC11497"))

quad.report<-read.table("~/Google Drive/Microsatellites/report_quad_20131006.txt",header=TRUE,sep="\t")
#remove CHP families and non-wholeblood families
quad.report<-quad.report[grep("auSSC.*wholeblood",quad.report$quad.quad_id),]
#remove wholeblood tag from remaining quad ids
quad.report$quad.quad_id<-sub("-wholeblood","",quad.report$quad.quad_id)
gatk.bad.families<-quad.report$quad.quad_id[grep("ok",quad.report$status,ignore.case=TRUE,invert=TRUE)]
gatk.bad.families.ind<-grep("ok",quad.report$status,ignore.case=TRUE,invert=TRUE)

low.obedience.clean<-low.obedience.scores[which(!(low.obedience.scores$family.id %in% gatk.bad.families) &
                                                !(low.obedience.scores$family.id %in% excluded.by.cm.or.ebv)),]
low.obedience.clean$obedience.score <- -10 * log10(low.obedience.clean$obedience.score)


t<-subset(low.obedience.clean,obedience.score >= 60 & context == "exon" & who.violates != "both" & om.or.com == "commission")
write.table(t,file="~/Google Drive/Microsatellites/exonDeNovos2.txt",sep="\t",quote=F,row.names=F)

plot(log10(1-low.obedience.clean$dn.confidence),low.obedience.clean$obedience.score,pch=16,cex=0.5,col=rgb(0,0,0,0.1),xlim=c(-5,0))

plot(log10(low.obedience.clean$em.noise.rate),low.obedience.clean$obedience.score,pch=16,cex=0.5,col=rgb(0,0,0,0.1),xlim=c(-5,0),ylim=c(0,200))
t2<-subset(low.obedience.clean,em.noise.rate < 0.1 & dn.confidence > 0.95 & obedience.score >= 40 & dn.noise.fit > 1e-4 & 
           dn.allele.fit > 1e-4 & dn.p.null < 0.01)

plot(low.obedience.clean$dn.confidence,low.obedience.clean$obedience.score,pch=16,cex=0.5,col=rgb(0,0,0,0.1))
plot(log10(low.obedience.clean$dn.allele.fit),low.obedience.clean$obedience.score,pch=16,cex=0.5,col=rgb(0,0,0,0.1),xlim=c(-5,0))
plot(low.obedience.clean$dn.noise.fit,low.obedience.clean$obedience.score,pch=16,cex=0.5,col=rgb(0,0,0,0.1))
plot(log10(low.obedience.clean$dn.p.null),low.obedience.clean$obedience.score,pch=16,cex=0.5,col=rgb(0,0,0,0.1),xlim=c(-5,0))


plot(log10(low.obedience.scores$denovo.p.null),-10*log10(low.obedience.scores$obedience.score),xlim=c(-5,0))
plot(log10(low.obedience.scores$dn.p.null),log10(low.obedience.scores$dn.allele.fit),pch=16,cex=0.5,col=rgb(0,0,0,0.1),
     xlim=c(-5,0),ylim=c(-5,0))
plot(low.obedience.scores$denovo.p.null,low.obedience.scores$denovo.noise.fit)
plot(low.obedience.scores$denovo.allele.fit,low.obedience.scores$denovo.noise.fit)

plot(sort(log10(low.obedience.scores$denovo.allele.fit)),ylim=c(-10,0))
plot(sort(log10(low.obedience.scores$denovo.noise.fit)),ylim=c(-10,0))
plot(sort(log10(low.obedience.scores$denovo.p.null)),ylim=c(-10,0))