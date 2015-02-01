hexToR<-function(h) {
  return(strtoi(substring(h, 1, 2), base=16))
}
hexToG<-function(h) {
  return(strtoi(substring(h, 3, 4), base=16))
}
hexToB<-function(h) {
  return(strtoi(substring(h, 5, 6), base=16))
}
cutHex<-function(h) {
  if(substring(h, 1, 1) == "#") return(substring(h, 2, nchar(h))) else return(h)
}
addAlphaToHex <- function(hex, alpha=0.5) {
  hex <- cutHex(hex)
  r   <- hexToR(hex)
  g   <- hexToG(hex)
  b   <- hexToB(hex)
  maxVal <- 255
  alpha <- maxVal * alpha
  return(rgb(r, g, b, alpha, max=maxVal))
}

get.denovo.allele <- function(who.dn, fam.genotypes) {
  parent.alleles <- fam.genotypes[1:4]
  unique.parent.alleles <- unique(fam.genotypes[1:4])
  child.alleles  <- {
    if(who.dn == 1) {
      fam.genotypes[5:6]
    } else if(who.dn == 2) {
      fam.genotypes[7:8]
    } else {
      fam.genotypes[5:8]
    }
  }
  unique.child.alleles<-unique(child.alleles)
  
  total.dn.in.parents <- sum(unique.child.alleles %in% unique.parent.alleles)
  if(total.dn.in.parents == length(unique.child.alleles)) {
    # de novo is an omission
    return(unique(c(unique.child.alleles[unique.child.alleles %in% parent.alleles[1:2] & 
                                         !(unique.child.alleles %in% parent.alleles[3:4])],
                    unique.child.alleles[!(unique.child.alleles %in% parent.alleles[1:2]) &
                                           unique.child.alleles %in% parent.alleles[3:4]])))
  } else {
    return(unique(unique.child.alleles[!(unique.child.alleles %in% parent.alleles)]))
  }
}

plot.family <- function(locus.info, allele.info, genotypes, locus.ind, person.info, individual.ind=NULL, family.ind=NULL, 
                        individual.id=NULL,family.id=NULL, palette="Blues", new.display=FALSE, print.dir, who.dn, 
                        fam.genotypes, anno, err.rate, trio.metrics=TRUE, file.type=c("png","pdf"), x.max=c("fam","pop"),
                        simple=F) {
  require(RColorBrewer)
  
  x.max <- match.arg(x.max)
  file.type <- match.arg(file.type)
  
  fam.col <- colorRampPalette(brewer.pal(brewer.pal.info[palette,'maxcolors'], palette))
  t.ind   <- which(locus.info$locus.row.ind == locus.ind & locus.info$allele.no == 0)
  a.ind   <- (t.ind + 1):(t.ind + locus.info$num.lengths[t.ind])
  pop.lengths <- locus.info$allele.no[a.ind]
  
  # get maximum allele length in population
  max.length  <- max(pop.lengths)
  
  fam.ind <- get.fam.member.indices(individual.ind, family.ind, individual.id, family.id, person.info)
  if(length(fam.ind) == 0) {
    stop("Could not find specified family")
  }
  
  fam.mat <- matrix(0, ncol=max.length, nrow=4)
  fam.mat[1:4, locus.info$allele.no[a.ind]] <- t(allele.info$alleles[a.ind, fam.ind])
  
  if(x.max == "fam"){
  # reset max.length to maximum length observed in family
    max.length <- max(which(colSums(fam.mat) > 0))
    fam.mat    <- fam.mat[, 1:max.length]
  }
  fam.mat    <- apply(fam.mat, 2, rev)
  max.count  <- max(fam.mat)
  
  seq        <- get.seq.axis(as.character(locus.info$unit[t.ind]), max.length)
  
  label.info  <- data.frame(cbind(which(fam.mat > 0, arr.ind=T),
                                 count=fam.mat[which(fam.mat > 0)]))
  if(!simple){
    label.color <- rep("black", nrow(label.info))
    label.color[which(label.info$count > (0.5 * max.count))] <- "gold"
  } else
  {
    label.color <- rep("gold", nrow(label.info))
  }
  # get expected coverage per genotyped allele
  expected.cov <- data.frame(unique(cbind(rep(4:1, 2), matrix(genotypes[locus.ind, fam.ind, 1:2], ncol=1),
                                          round(rep(genotypes[locus.ind, fam.ind, 9], 2), 1))))
  # remove duplicated expected coverage (i.e. at homozygous genotypes)
  expected.cov <- expected.cov[which(expected.cov[, 1] %in% label.info[, 1] & expected.cov[, 2] %in% label.info[, 2]), ]
  names(expected.cov) <- c("row", "col", "exp.cov")
  # expected coverage color is the same as the color for the observed coverage
  expected.cov.color  <- unlist(sapply(1:nrow(expected.cov), function(i) 
                         label.color[which(expected.cov$col[i] == label.info$col &
                                           expected.cov$row[i] == label.info$row)]))
  
  #label format is sampleID (relation, gender)\ngenotype genotype.confidence\nallele.fit,noise.fit
  if(!simple){
  ytickLabels <- sapply(1:4, function(k) sprintf("%s (%s)\n%d|%d %0.2f\n%0.1E;%0.1E",
                                              capitalize(as.character(person.info$relation[fam.ind[k]])), 
                                              person.info$gender[fam.ind[k]], 
                                              genotypes[locus.ind,fam.ind[k], 1], genotypes[locus.ind,fam.ind[k], 2], 
                                              genotypes[locus.ind,fam.ind[k], 3], genotypes[locus.ind,fam.ind[k], 4],
                                              genotypes[locus.ind,fam.ind[k], 5]))
  } else{
    ytickLabels <- sapply(1:4, function(k) sprintf("%s (%s)\n%d|%d %0.2f",
                                                   capitalize(as.character(person.info$relation[fam.ind[k]])), 
                                                   person.info$gender[fam.ind[k]], 
                                                   genotypes[locus.ind,fam.ind[k], 1], genotypes[locus.ind,fam.ind[k], 2], 
                                                   genotypes[locus.ind,fam.ind[k], 3]))
  }

  dn.score <- NULL
  if(!missing(who.dn) && !missing(fam.genotypes)) {
    dn.score <- get.dn.score(who.dn, fam.genotypes)
  }
  
  if(!missing(print.dir) && is.character(print.dir)) {
    if(file.type == "png"){
      fname<-sprintf("%s_%s.png", get.plot.name.preamble(print.dir, locus.info, t.ind, dn.score),
                                  person.info$family.id[fam.ind[1]])
      png(fname, width=(1.618 * 2000), height=2000)
      if(trio.metrics){
        layout(matrix(1:2,nr=1),widths=c(8,2))
      }
      par(cex=4)
    } else if(file.type == "pdf"){
      fname<-sprintf("%s_%s.pdf", get.plot.name.preamble(print.dir, locus.info, t.ind, dn.score),
                     person.info$family.id[fam.ind[1]])
      pdf(fname, width=13, height=8)
      if(trio.metrics){
        layout(matrix(1:2,nr=1),widths=c(8,3))
      }
    }
  } else if(new.display) {
    get.new.display(width=(1.618 * 12.5),height=12.5)
  }
  fam.mat[fam.mat == 0] <- NA
  
  dn.alleles <- get.denovo.allele(who.dn, fam.genotypes)
  
  
  if(!simple) par(mar=c(5, 11, 8, 2) + 0.1)
  if(simple) par(mar=c(5, 11, 6, 2) + 0.1)
  image(1:ncol(fam.mat), 1:nrow(fam.mat), t(fam.mat), col=fam.col(max.count), ylab="", xlab="",
        xaxt="n", yaxt="n", bty="n")
  par(font.lab=1)
  main.title <- get.locus.title(locus.info, t.ind)
  if(!missing(err.rate) && !simple) {main.title <- sprintf("%s, error rate %.1E",main.title,err.rate)}
  if(!missing(anno))                {main.title <- append.annotation(main.title,anno)}
  if(!missing(who.dn) && !simple)   {main.title <- append.dn.info(main.title, who.dn,fam.genotypes, 
                                                       genotypes[locus.ind,fam.ind,],
                                                       person.info$family.id[fam.ind[1]])}
  title(main=main.title, line=3)
  
  title(ylab="Coverage", line=9.0, font.lab=2)
  title(xlab="Tract length(nt)", font.lab=2)
  
  # label boxes in plot with the appropriate observed and expected coverages
  text(label.info$col, label.info$row, label.info$count, col=label.color, font=2)
  text(expected.cov$col, expected.cov$row - 0.2, sapply(expected.cov$exp.cov, function(x) sprintf("(%.1f)", x)), col=expected.cov.color)

  motif.length <- nchar(as.character(locus.info$unit[t.ind]))
  axis(side=1, at=seq(motif.length,ncol(fam.mat),motif.length), lwd=0, lwd.ticks=1)
  axis(side=2, at=1:4, labels=rev(ytickLabels), las=2, lwd=0, lwd.ticks=1)
  axis(side=3, 1:ncol(fam.mat), labels=seq, lwd=0, lwd.ticks=1)
  
  # draw a red box around all de novo alleles
  box.lwd <- ifelse(file.type=="png", 8, 4)
  for(i in dn.alleles)
  {
    rect(xleft=(i - 0.5), xright=(i + 0.5), ytop=4.5, ybottom=0.5, col=NA, border="red", 
         lwd=box.lwd,lty=2)
  }
  
  ref.length <- locus.info$ref.length[t.ind]
  rect(xleft=(ref.length - 0.5), xright=(ref.length + 0.5), ytop=4.5, ybottom=0.5, col=NA, border="green", 
       lwd=box.lwd, lty=2)
  
  
  if(trio.metrics){
    add.trio.metrics(who.dn, fam.genotypes, genotypes[locus.ind, fam.ind, ],file.type=file.type)
  }
  
  if(!missing(print.dir) && is.character(print.dir)) {
    dev.off()
    return(fname)
  } else {
    return(NULL)
  }
}

plot.pop.heatmap <- function(locus.info, allele.info, pop.size, genotypes, locus.ind, person.info, individual.ind=NULL, 
                             family.ind=NULL, individual.id=NULL, family.id=NULL, palette="Blues", new.display=TRUE,
                             display.dn.freq=FALSE, display.ref.freq=FALSE, print.dir, who.dn, fam.genotypes, anno, 
                             err.rate, subgroup,file.type=c("png","pdf"),scale=T,simple=T) {
  require(RColorBrewer)
  file.type <- match.arg(file.type)
  
  pop.col <- colorRampPalette(brewer.pal(brewer.pal.info[palette,'maxcolors'], palette))
  t.ind   <- which(locus.info$locus.row.ind == locus.ind & locus.info$allele.no == 0)
  a.ind   <- (t.ind + 1):(t.ind + locus.info$num.lengths[t.ind])
  pop.lengths <- locus.info$allele.no[a.ind]
  max.length  <- max(pop.lengths)
  
  group.info <- {
    if(!missing(subgroup)){
      get.pop.subgroup(person.info, subgroup)
    }  else {
      list(ind=1:pop.size, anno=NULL, size=pop.size)
    }
  }
  
  pop.mat <- matrix(0, ncol=max.length, nrow=group.info$size)
  pop.mat[1:group.info$size, locus.info$allele.no[a.ind]] <- t(allele.info$alleles[a.ind, group.info$ind])
                                                              
  max.count <- max(pop.mat[!is.na(pop.mat)])
  max.color <- ((max.count %/% 50) * 50) + 50
  
  seq <- get.seq.axis(as.character(locus.info$unit[t.ind]), max.length)

  dn.alleles <- vector(mode="numeric",length=0)
  if(!missing(fam.genotypes)) {
    dn.alleles <- get.denovo.allele(who.dn,fam.genotypes)
  }
  
  dn.score <- NULL
  if(!missing(who.dn) && !missing(fam.genotypes)) {
    dn.score <- get.dn.score(who.dn, fam.genotypes)
  }
  
  if(!missing(print.dir) && is.character(print.dir)) {
    if(file.type == "png"){
      fname <- sprintf("%s_popHm.png", get.plot.name.preamble(print.dir, locus.info, 
                                                              t.ind, dn.score))
      png(fname, width=(1.618 * 2000), height=2000)
	  if(scale) layout(mat=matrix(1:2, nr=1), widths=c(10, 1))
      par(cex=4)
    } else if(file.type == "pdf"){
      fname <- sprintf("%s_popHm.pdf", get.plot.name.preamble(print.dir, locus.info, 
                                                              t.ind, dn.score))
      pdf(fname, width=11, height=8)
      if(scale) layout(mat=matrix(1:2, nr=1), widths=c(10, 1))
    }
  } else if(new.display) {
    get.new.display(width=(1.618 * 12.5),height=12.5)
    layout(mat=matrix(1:2, nr=1), widths=c(10, 1))
  }
  if(!simple) par(mar=c(5, 5, 8, 1) + 0.1)
  if(simple) par(mar=c(5, 5, 6, 1) + 0.1)
  pop.mat[pop.mat == 0] <- NA

  image(1:ncol(pop.mat), 1:nrow(pop.mat), t(pop.mat), col=pop.col(max.color), 
        xlab="Tract length(nt)", ylab="SSC person", font.lab=2, axes=FALSE)
  
  main.title <- get.locus.title(locus.info, t.ind)
  if(!missing(err.rate) && !simple) {main.title <- sprintf("%s, error rate %.1E", main.title, err.rate)}
  if(!missing(anno))                {main.title <- append.annotation(main.title,anno)}
  main.title <- sprintf("Population coverage for %s", main.title)
#   ifelse(as.log10,main.title <- sprintf("Population log10 coverage for %s", main.title),main.title<-sprintf("Population coverage for %s",main.title))
  if(length(who.dn) > 0 && !simple) {main.title <- append.dn.info(main.title,who.dn,fam.genotypes, genotypes[locus.ind,fam.ind,],
                                                     person.info$family.id[fam.ind[1]])}
  main.title<-append.subgroup(main.title,group.info)
  title(main=main.title,line=3)
  
  parents <- person.info$relation %in% c("mother","father")
  
  box.lwd <- ifelse(file.type=="png", 8, 4)

  for(i in dn.alleles)
  {
    rect(xleft=(i - 0.5), xright=(i + 0.5), ytop=group.info$size + 0.5, ybottom=0.5, col=NA, border="red", 
        lwd=box.lwd, lty=2)
    if(display.dn.freq && !simple) {
      text(x=(i - 0.5), y=nrow(pop.mat)/4, font=2, pos=2, cex=0.75,
           labels=sprintf("Parental freq\n%0.1E", sum(genotypes[locus.ind,parents,1:2] == i) / (2 * sum(parents))))
    }
  }

  ref.length <- locus.info$ref.length[t.ind]
  rect(xleft=(ref.length - 0.5), xright=(ref.length + 0.5), ytop=group.info$size + 0.5, ybottom=0.5, col=NA, border="green", 
       lwd=box.lwd, lty=2)
  if(display.ref.freq && !simple) {
    text(x=(ref.length - 0.5), y=nrow(pop.mat)/2, font=2, pos=2, cex=0.75,
         labels=sprintf("Parental freq\n%0.1E", sum(genotypes[locus.ind,parents,1:2] == ref.length) / (2 * sum(parents))))
  }
  
  motif.length <- nchar(as.character(locus.info$unit[t.ind]))
  axis(side=1,at=seq(1,ncol(pop.mat)),labels=NA)
  axis(side=1,at=seq(motif.length,ncol(pop.mat),motif.length),lwd=0,lwd.tick=0)
  axis(side=2,at=c(1,seq(500,nrow(pop.mat),500),nrow(pop.mat)))
  axis(side=3,1:ncol(pop.mat),labels=seq)

  if(scale){
	par(mar=c(7,1,10,3) + 0.1)
	image(1,1:max.color,matrix(1:max.color,nr=1,nc=max.color),col=pop.col(max.color),
          ylab="",xlab="",xaxt="n",yaxt="n",bty="n")
	axis(side=4,at=seq(0,max.color,25), lwd=0, lwd.tick=1)
  }
  if(!missing(print.dir) && is.character(print.dir)) {
    dev.off()
    return(fname)
  } else {
    return(NULL)
  }
}

capitalize <- function(string) {
  return(paste(toupper(substr(string, 1, 1)),
               substr(string, 2, nchar(string)),
               sep=""))
}

# can pass width and height to graphic device available on OS
get.new.display <- function(width=(1.618 * 7), height=7){
  if(capabilities("aqua")){
    quartz(width=width, height=height)
  } else if(capabilities("X11")) {
    x11(width=width, height=height)
  } else{
    message("The current R build does not support X11 or quartz, cannot launch new display")
  }
}

plot.pop.allele.count <- function(locus.info, allele.info, genotypes, locus.ind, palette="Blues",
                                  new.display=TRUE, print.dir, fam.genotypes, who.dn, anno, 
                                  file.type=c("png","pdf")) {  
  require(RColorBrewer)
  file.type <- match.arg(file.type)
  
  af.col <- colorRampPalette(brewer.pal(brewer.pal.info[palette,'maxcolors'], palette))
  t.ind  <- which(locus.info$locus.row.ind == locus.ind & locus.info$allele.no == 0)
  a.ind  <- (t.ind+1):(t.ind + locus.info$num.lengths[t.ind])
  num.observed.alleles <- length(a.ind) + 1
  pop.lengths <- locus.info$allele.no[a.ind]
    
  pop.mat<-matrix(0,ncol=num.observed.alleles,nrow=num.observed.alleles)
  geno.codes<-cbind(c(-1, pop.lengths), 1:num.observed.alleles)
  label.info<-data.frame(x=vector(),y=vector(),text=vector(),col=vector())
  
  for(i in seq_len(num.observed.alleles)) {
    for(j in seq_len(num.observed.alleles)) {
      num.with.geno <- sum(genotypes[locus.ind, , 1] == geno.codes[i, 1] &
                           genotypes[locus.ind, , 2] == geno.codes[j, 1])
      pop.mat[i, j] <- num.with.geno
    }
  }
    
  max.count   <- max(pop.mat)
  label.info  <- data.frame(cbind(which(pop.mat > 0, arr.ind=T),
                                 count=pop.mat[which(pop.mat > 0)]))
  label.color <- rep("black", sum(pop.mat > 0))
  label.color[which(label.info$count > (0.5 * max.count))] <- "gold"
  
  dn.score <- NULL
  if(!missing(who.dn) && !missing(fam.genotypes)) {
    dn.score <- get.dn.score(who.dn, fam.genotypes)
  }
  
  if(!missing(print.dir) && is.character(print.dir)) {
    if(file.type == "png"){
      fname <- sprintf("%s_genoFreq.png", get.plot.name.preamble(print.dir, locus.info, t.ind, dn.score))
      png(fname, width=(1.618 * 2000), height=2000)
      par(cex=4)
    } else if(file.type == "pdf"){
      fname <- sprintf("%s_genoFreq.pdf", get.plot.name.preamble(print.dir, locus.info, t.ind, dn.score))
      pdf(fname, width=8, height=8)
    }
  } else if(new.display==TRUE)
  {
    get.new.display()
  }

  main.title <- get.locus.title(locus.info, t.ind)
  main.title <- paste(main.title, "\n", sep="")
  if(!missing(anno)) {main.title <- append.annotation(main.title, anno)}

  pop.mat[pop.mat==0] <- NA
  image(1:ncol(pop.mat), 1:nrow(pop.mat), t(pop.mat), col=af.col(max.count), main=main.title,
        xlab="Allele B", ylab="Allele A", font.lab=2, xaxt="n", yaxt="n", bty="n")
  text(x=label.info$col, y=label.info$row, labels=label.info$count, col=label.color, font=2)
  axis(side=1, at=geno.codes[, 2], labels=geno.codes[, 1], lwd=0, lwd.tick=1)
  axis(side=2, at=geno.codes[, 2], labels=geno.codes[, 1], lwd=0, lwd.tick=1)

  if(!missing(fam.genotypes) && !missing(who.dn)) {
    if(who.dn == 1) {
      dn.al.1.ind <- geno.codes[which(geno.codes[, 1] == fam.genotypes[5]), 2]  
      dn.al.2.ind <- geno.codes[which(geno.codes[, 1] == fam.genotypes[6]), 2]
      
      rect(xleft=(dn.al.2.ind - 0.5), xright=(dn.al.2.ind + 0.5),
           ytop=(dn.al.1.ind - 0.5), ybottom=(dn.al.1.ind + 0.5),
           col=NA, border="red", lwd=6, lty=2)
    } else if(who.dn == 2) {
      dn.al.1.ind <- geno.codes[which(geno.codes[, 1] == fam.genotypes[7]), 2]  
      dn.al.2.ind <- geno.codes[which(geno.codes[, 1] == fam.genotypes[8]), 2]
      
      rect(xleft=(dn.al.2.ind - 0.5), xright=(dn.al.2.ind + 0.5),
           ytop=(dn.al.1.ind - 0.5), ybottom=(dn.al.1.ind + 0.5),
           col=NA, border="red", lwd=6, lty=2)
    } else if(who.dn == 3) {
      dn.al.1.ind <- geno.codes[which(geno.codes[, 1] == fam.genotypes[5]), 2]  
      dn.al.2.ind <- geno.codes[which(geno.codes[, 1] == fam.genotypes[6]), 2]
      
      rect(xleft=(dn.al.2.ind - 0.5), xright=(dn.al.2.ind + 0.5),
           ytop=(dn.al.1.ind - 0.5), ybottom=(dn.al.1.ind + 0.5),
           col=NA, border="red", lwd=6, lty=2)
      if(fam.genotypes[5] != fam.genotypes[7] | fam.genotypes[6] != fam.genotypes[8]) {
        # if the children have different de novo mutations, highlight both of them
        dn.al.1.ind <- geno.codes[which(geno.codes[, 1] == fam.genotypes[7]), 2]  
        dn.al.2.ind <- geno.codes[which(geno.codes[, 1] == fam.genotypes[8]), 2]
        
        rect(xleft=(dn.al.2.ind - 0.5), xright=(dn.al.2.ind + 0.5),
             ytop=(dn.al.1.ind - 0.5), ybottom=(dn.al.1.ind + 0.5),
             col=NA, border="red", lwd=6, lty=2) 
      }   
    }
  }
  
  if(!missing(print.dir) && is.character(print.dir)) {
    dev.off()
    return(fname)
  } else {
    return(NULL)
  }
}

plot.dn.eo.scatter <- function(locus.info, allele.info, pop.size, genotypes, locus.ind, fam.ind, by.fams, 
                               coverage.estimator, new.display=TRUE, who.dn, fam.genotypes, print.dir,
                               anno, file.type=c("png","pdf")) {
  file.type <- match.arg(file.type)
  
  locus.details <- get.locus.details(locus.ind, locus.info)
  t.ind         <- which(locus.info$locus.row.ind == locus.ind & locus.info$allele.no == 0)
  cov.info      <- get.coverage(allele.info, locus.details,pop.size)
  obs.exp.cov.ratio <- cov.info$allele.cov / matrix(rep(coverage.estimator$exp.cov[locus.ind, ],length(locus.details$lengths)),
                                                    nrow=length(locus.details$lengths), byrow=TRUE)
  
  dn.info <- get.dn.alleles.and.indices(who.dn, locus.info, genotypes, by.fams, fam.ind, locus.details, cov.info, locus.ind)
  
  max.oe<-6 
  # this ensures that the smaller axis is still visually meaningful, we don't expect too many people to have 6+ times their expected coverage
  
  plot.colors <- get.genotype.colors(genotypes=genotypes, locus.ind)
  
  dn.score    <- NULL
  if(!missing(who.dn) && !missing(fam.genotypes)) {
    dn.score  <- get.dn.score(who.dn, fam.genotypes)
  }
  
  if(!missing(print.dir) && is.character(print.dir)) {
    if(file.type == "png"){
      fname <- sprintf("%s_%d_%d_oe.png", get.plot.name.preamble(print.dir, locus.info, t.ind, dn.score),
                       dn.info$alleles[1], dn.info$alleles[2])
      png(fname, width=2400, height=2400)
      par(cex=4)
    } else if(file.type == "pdf"){
      fname <- sprintf("%s_%d_%d_oe.pdf", get.plot.name.preamble(print.dir, locus.info, t.ind, dn.score),
                       dn.info$alleles[1], dn.info$alleles[2])
      pdf(fname, width=8, height=8)
    }
  } else if(new.display) {
    get.new.display()
  }
  
  plot.allele.pop.scatter(data=obs.exp.cov.ratio, locus.details=locus.details, x.row=dn.info$indices[1], y.row=dn.info$indices[2],
                          axis.label.info="observed / expected coverage ratio", locus.info=locus.info, plot.colors=plot.colors,
                          max.val=max.oe, anno=anno, highlight.family=by.fams[fam.ind,], file.type=file.type)
  if(!missing(print.dir) && is.character(print.dir)) {
    dev.off()
    return(fname)
  } else {
    return(NULL)
  }
}

plot.dn.cov.scatter <- function(locus.info, allele.info, pop.size, genotypes, locus.ind, fam.ind, by.fams,
                                new.display=TRUE, who.dn, fam.genotypes, print.dir, anno, plot.3d=FALSE,
                                file.type=c("png","pdf")) {
  file.type <- match.arg(file.type)
  
  t.ind         <- which(locus.info$locus.row.ind == locus.ind & locus.info$allele.no == 0)
  locus.details <- get.locus.details(locus.ind, locus.info)
  cov.info      <- get.coverage(allele.info, locus.details, pop.size)
  
  dn.info <- get.dn.alleles.and.indices(who.dn, locus.info, genotypes, by.fams, fam.ind, 
                                        locus.details, cov.info, locus.ind)
  
  max.cov     <- 1.1 * max(max(cov.info$allele.cov[dn.info$indices[1], ]),
                           max(cov.info$allele.cov[dn.info$indices[2], ]))
  if(plot.3d && (locus.details$n.alleles - 1) > 2){
    if((locus.details$n.alleles - 1) > 3){
      max.cov <- 1.1 * max(max(cov.info$allele.cov[dn.info$indices[1], ]),
                           max(cov.info$allele.cov[dn.info$indices[2], ]),
                           max(colSums(cov.info$allele.cov[-dn.info$indices, ])))
    } else {
      max.cov <- 1.1 * max(max(cov.info$allele.cov[dn.info$indices[1], ]),
                           max(cov.info$allele.cov[dn.info$indices[2], ]),
                           max(cov.info$allele.cov[-dn.info$indices, ]))
    }
  }

  plot.colors <- get.genotype.colors(genotypes=genotypes, locus.ind)
  
  dn.score <- NULL
  if(!missing(who.dn) && !missing(fam.genotypes)) {
    dn.score <- get.dn.score(who.dn, fam.genotypes)
  }
  
  if(!missing(print.dir) && is.character(print.dir)) {
    fname <- sprintf("%s_%d_%d_cov", get.plot.name.preamble(print.dir, locus.info, t.ind, dn.score),
                     dn.info$alleles[1], dn.info$alleles[2])
    if(plot.3d){
      fname <- sprintf("%s_3d",fname)
    }
    if(file.type == "png"){
      fname <- sprintf("%s.png",fname)
      png(fname,width=2400,height=2400)
      par(cex=4)
    } else if(file.type == "pdf"){
      fname <- sprintf("%s.pdf",fname)
      pdf(fname,width=8,height=8)
    }
    if(plot.3d) {
      par(mar=(par("mar") + c(0, 0, 0, 2)))
    }
  } else if(new.display) {
    get.new.display()
  }
  plot.allele.pop.scatter(data=cov.info$allele.cov, locus.details=locus.details, x.row=dn.info$indices[1], y.row=dn.info$indices[2],
                          axis.label.info="observed coverage", locus.info=locus.info, plot.colors=plot.colors, max.val=max.cov, 
                          anno=anno,highlight.family=by.fams[fam.ind,], plot.3d=plot.3d,file.type=file.type)
  
  if(!missing(print.dir) && is.character(print.dir)) {
    dev.off()
    return(fname)
  } else {
    return(NULL)
  }
}

plot.cov.scatter <- function(locus.info, allele.info, pop.size, genotypes, locus.ind, allele.1, allele.2,
                             print.dir,highlight.family,new.display=FALSE,anno, plot.3d=FALSE) {
  locus.details <- get.locus.details(locus.ind, locus.info)
  cov.info      <- get.coverage(allele.info, locus.details, pop.size)
  
  allele.1.ind <- which(locus.details$lengths == allele.1)
  allele.2.ind <- which(locus.details$lengths == allele.2)
  
  max.cov <- 1.1 * max(max(cov.info$allele.cov[allele.1.ind, ]),
                       max(cov.info$allele.cov[allele.2.ind, ]))
  
  if(plot.3d && (locus.details$n.alleles - 1) > 2){
    if((locus.details$n.alleles - 1) > 3){
      max.cov <- 1.1 * max(max(cov.info$allele.cov[allele.1.ind, ]),
                           max(cov.info$allele.cov[allele.2.ind, ]),
                           max(colSums(cov.info$allele.cov[-c(allele.1.ind,allele.2.ind), ])))
    } else {
      max.cov <- 1.1 * max(max(cov.info$allele.cov[allele.1.ind, ]),
                           max(cov.info$allele.cov[allele.2.ind, ]),
                           max(cov.info$allele.cov[--c(allele.1.ind,allele.2.ind), ]))
    }
  }
  
  plot.colors <- get.genotype.colors(genotypes=genotypes, locus.ind)
  
  if(!missing(print.dir) && is.character(print.dir)) {
 
    fname <- sprintf("%s_%d_%d_cov", get.plot.name.preamble(print.dir, locus.info, allele.info$locus.row.indices[locus.ind]),
                     allele.1, allele.2)
    if(plot.3d){
      fname <- sprintf("%s_3d",fname)
    }
    fname <- sprintf("%s.png",fname)
    
    png(fname,width=2400,height=2400)
    par(cex=4)
  } else if(new.display) {
    get.new.display()
  }
  
  if(missing(highlight.family)){
    plot.allele.pop.scatter(data=cov.info$allele.cov, locus.details=locus.details, x.row=allele.1.ind, y.row=allele.2.ind,
                            axis.label.info="observed coverage", locus.info=locus.info, plot.colors=plot.colors, 
                            max.val=max.cov, anno=anno, plot.3d=plot.3d)
  } else {
    plot.allele.pop.scatter(data=cov.info$allele.cov, locus.details=locus.details, x.row=allele.1.ind, y.row=allele.2.ind,
                            axis.label.info="observed coverage", locus.info=locus.info, plot.colors=plot.colors, 
                            max.val=max.cov, anno=anno, plot.3d=plot.3d, highlight.family=highlight.family)
  }
  if(!missing(print.dir) && is.character(print.dir)) {
    dev.off()
  }
}

plot.eo.scatter <- function(locus.info, allele.info, pop.size, genotypes, locus.ind, coverage.estimator,
                            allele.1, allele.2, print.dir, new.display=FALSE, anno) {
  locus.details <- get.locus.details(locus.ind, locus.info)
  cov.info      <- get.coverage(allele.info, locus.details, pop.size)
  obs.exp.cov.ratio <- cov.info$allele.cov / matrix(rep(coverage.estimator$exp.cov[locus.ind, ],
                                                        length(locus.details$lengths)), 
                                                    nrow=length(locus.details$lengths), byrow=TRUE)
  
  allele.1.ind <- which(locus.details$lengths == allele.1)
  allele.2.ind <- which(locus.details$lengths == allele.2)
  
  max.oe <- 6 
  # this ensures that the smaller axis is still visually meaningful, we don't expect too many people to have 6+ times their expected coverage
  plot.colors <- get.genotype.colors(genotypes=genotypes, locus.ind)
  
  if(!missing(print.dir) && is.character(print.dir)) {
    fname <- sprintf("%s_%d_%d_oe.png", get.plot.name.preamble(print.dir, locus.info,
                                                               allele.info$locus.row.indices[locus.ind]),
                     allele.1, allele.2)
    png(fname,width=2400, height=2400)
    par(cex=4)
  } else if(new.display) {
    get.new.display()
  }
  
  plot.allele.pop.scatter(data=obs.exp.cov.ratio, locus.details=locus.details, x.row=allele.1.ind, y.row=allele.2.ind,
                          axis.label.info="observed / expected coverage ratio", locus.info=locus.info, plot.colors=plot.colors,
                          max.val=max.oe, anno=anno)
  if(!missing(print.dir) && is.character(print.dir)) {
    dev.off()
  }
}

get.dn.alleles.and.indices <- function(who.dn, locus.info, genotypes, by.fams, fam.ind, locus.details, 
                                     cov.info, locus.ind) {
  who.dn.alleles <- {
    if(who.dn == 2) {
      genotypes[locus.ind, by.fams[fam.ind, 4], 1:2]
      } else {
        genotypes[locus.ind, by.fams[fam.ind, 3], 1:2]
      }
  }
  
  dn.allele.inds <- which(locus.details$length %in% who.dn.alleles)
  if(length(unique(who.dn.alleles)) == 1) {
    # if child with the de novo allele is homozygous
    if(dn.alleles[1] != locus.info$ref.length[locus.details$locus.ind]) {
      # if the child is homozygous for a non-reference allele,
      # graph scatter against the reference length allele
      dn.alleles[2]     <- locus.info$ref.length[locus.details$locus.ind]
      dn.allele.inds[2] <- which(locus.details$lengths == dn.alleles[2])
    } else  {
      # if the child is homozygous for a reference allele,
      # graph scatter against the non-reference allele with the highest total coverage
      dn.allele.inds[2] <- which(rev(order(rowSums(cov.info$allele.cov[-dn.allele.inds[1], ]))) != dn.allele.inds)[1]
      dn.alleles[2]     <- locus.details$lengths[dn.allele.inds[2]]
    }
  }
  return(list(alleles=dn.alleles, indices=dn.allele.inds))
}
  
plot.allele.pop.scatter <- function(data, locus.details, x.row, y.row, axis.label.info, locus.info,
                                  pt.color, legend.info, plot.colors, max.val, spot.pch=16,
                                  anno, highlight.family, file.type=c("png","pdf"), plot.3d=FALSE) {
  if(!missing(plot.colors)) {
    if(missing(pt.color)) {
      # if pt.color isn't specified and plot.colors is, 
      # pt.colors = plot.colors$pt.colors, otherwise, use pt.colors
      pt.color <- plot.colors$pt.colors
    } else {
      message("Specified plot colors using plot.colors and pt.color, using pt.color specification")
    }
    if(missing(legend.info)) {
      # if legend.info isn't specified and plot.colors is
      # legend.info = plot.colors$legend.info
      legend.info<-plot.colors$legend
    } else {
      message("Specified plot legend information using plot.colors and legend.info, using legend.info specification")
    }
  } else {
    #if no specification was given for plot colors,
    # make them all black with 0.6 alpha
    pt.color<-rgb(0,0,0,0.6)
  }
  
  main.title <- get.locus.title(locus.info, locus.details$locus.ind)
  main.title <- paste(main.title, "\n", sep="")
  if(!missing(anno)) {main.title <- append.annotation(main.title, anno)}
  
  if(!plot.3d) {
    plot(data[x.row, ],data[y.row, ], xlab=sprintf("%d bp allele %s",locus.details$lengths[x.row], axis.label.info),
         ylab=sprintf("%d bp allele %s",locus.details$lengths[y.row], axis.label.info), pch=spot.pch, 
         xlim=c(0, max.val), ylim=c(0, max.val), main=main.title, col=as.character(pt.color))
    abline(a=0, b=1, col="black", lwd=2, lty=2)
    legend(x="topright", legend=legend.info$genotype, col=as.character(legend.info$color), pch=spot.pch,
           title="Genotypes", bg="white")
  
    if(!missing(highlight.family)) {
      fam.colors <- data.frame(member=c("mother","father","proband","sibling"),
                               color=c(rgb(1, 0, 0, 0.7), rgb(0, 0, 1, 0.7), rgb(0, 1, 0, 0.7), rgb(0, 0, 0, 0.7)))
      point.lwd <- ifelse(file.type == "png", 6, 2)
      points(data[x.row, highlight.family], data[y.row, highlight.family], pch=21,
             col=fam.colors$color, bg=NA, lwd=point.lwd)
      legend(x="top", legend=fam.colors$member, col=fam.colors$color, lwd=point.lwd, pch=21, xjust=0.5, lty=NA, horiz=TRUE)
    }
  } else {
    require(scatterplot3d)
    if(nrow(data) > 2) {
      zdata <- {
        if(nrow(data) > 3){
          colSums(data[-c(x.row, y.row), ])
        } else {
          data[-c(x.row, y.row), ]
        }
      }
            
      s3d <- scatterplot3d(x=data[x.row, ], y=data[y.row, ], z=zdata,
                  xlab=sprintf("%d bp allele %s", locus.details$lengths[x.row], axis.label.info),
                  ylab=sprintf("%d bp allele %s", locus.details$lengths[y.row], axis.label.info), 
                  zlab=sprintf("other allele %s", axis.label.info), pch=spot.pch, xlim=c(0, max.val), 
                  ylim=c(0, max.val), zlim=c(0,max.val), main=main.title, color=as.character(pt.color),
                  cex.symbols=1, cex.axis=4, cex.lab=4)
      s3d$points(x=data[x.row, ], y=data[y.row, ], z=zdata, type="h", pch= " ",
                 col=rgb(0, 0, 0, 0.3))
      legend(s3d$xyz(max.val * 0.8, max.val * 0.8, max.val * 0.8), legend=legend.info$genotype, col=as.character(legend.info$color),
             pch=spot.pch, title="Genotypes", bg="white", xjust=0)
    
      if(!missing(highlight.family)) {
        fam.colors <- data.frame(member=c("mother", "father", "proband", "sibling"),
                                color=c(rgb(1, 0, 0, 0.7), rgb(0, 0, 1, 0.7), rgb(0, 1, 0, 0.7), rgb(0, 0, 0, 0.7)))
      
        s3d$points3d(x=data[x.row, highlight.family], y=data[y.row, highlight.family], 
                     z=zdata[highlight.family], pch=21, lwd=6,
                     col=fam.colors$color, bg=NA)
        legend(s3d$xyz(max.val, max.val, max.val), legend=fam.colors$member, col=fam.colors$color, pch=21, lwd=6, lty=NA,
               horiz=TRUE, xjust=1, yjust=0)
      }
    } else {
      plot(data[x.row, ],data[y.row, ], xlab=sprintf("%d bp allele %s",locus.details$lengths[x.row], axis.label.info),
           ylab=sprintf("%d bp allele %s",locus.details$lengths[y.row], axis.label.info), pch=spot.pch, 
           xlim=c(0, max.val), ylim=c(0, max.val), main=main.title, col=as.character(pt.color))
      abline(a=0, b=1, col="black", lwd=2, lty=2)
      legend(x="topright", legend=legend.info$genotype, col=as.character(legend.info$color), pch=spot.pch,
             title="Genotypes", bg="white")
      
      if(!missing(highlight.family)) {
        fam.colors <- data.frame(member=c("mother","father","proband","sibling"),
                                 color=c(rgb(1, 0, 0, 0.7), rgb(0, 0, 1, 0.7), rgb(0, 1, 0, 0.7), rgb(0, 0, 0, 0.7)))
        points(data[x.row, highlight.family], data[y.row, highlight.family], pch=18, cex=1.5,
               col=fam.colors$color)
        legend(x="top", legend=fam.colors$member, col=fam.colors$color, pch=18, xjust=0.5, horiz=TRUE)
      }
    }
  }
}

get.fam.member.indices <- function(individual.ind=NULL, family.ind=NULL, individual.id=NULL,
                                 family.id=NULL, person.info) {
  fam.ind <- {
    if(!is.null(individual.ind) && is.integer(individual.ind)) {
      if(individual.ind <= nrow(person.info) & individual.ind > 0) {
        fam.ind <- which(person.info$family.id == person.info$family.id[individual.ind])
      } else {
        stop(sprintf("Invalid person index, index must be between 1 and %d",nrow(person.info)))
      }
    } else if(!is.null(individual.id) && is.character(individual.id)) {
      if(individual.id %in% person.info$individual.id) {
        fam.ind <- which(person.info$family.id == person.info$family.id[person.info$individual.id == individual.id])
      } else {
        stop(sprintf("%s could not be found in provided person.info matrix",individual.id))
      }
    } else if(!is.null(family.ind) && all(family.ind == floor(family.ind))) {
      if(length(family.ind) == 4 && max(family.ind) <= nrow(person.info) && min(family.ind) >= 1) {
        fam.ind <- family.ind
      } else {
        stop("Invalid family indices")
      }
    } else if(!is.null(family.id) && is.character(family.id)) {
      if(family.id %in% person.info$family.id) {
        fam.ind <- which(person.info$family.id == family.id)
      }else {
        print(sprintf("%s could not be found in provided person.info matrix",family.id))
      }
    } else {
      stop("Unrecognized person or family specifier!  Please provide a an individual index, a set of 4 family indices, an individual's ID, or a family ID")
    }
  }
  return(fam.ind)
}

get.genotype.colors <- function(genotypes, locus.ind, alpha=0.6, palette="Set1") {
  require(RColorBrewer)
  genotype.calls      <- apply(genotypes[locus.ind, , 1:2], 1, function(x) paste(x, collapse="|"))
  locus.genotypes     <- sort(unique(genotype.calls))
  num.colors          <- ifelse(length(locus.genotypes) >  brewer.pal.info[palette, 'maxcolors'], 
                                brewer.pal.info[palette, 'maxcolors'], length(locus.genotypes)) 
  genotype.color.func <- colorRampPalette(brewer.pal(num.colors, palette))
  genotype.colors     <- data.frame(genotype=locus.genotypes, 
                                    color=do.call("rbind", lapply(genotype.color.func(length(locus.genotypes)), 
                                                                  addAlphaToHex, alpha=0.6)))
  person.genotype.colors <- do.call("rbind", lapply(genotype.calls, function(x) as.character(genotype.colors[which(genotype.colors[, 1] == x), 2])))
  
  return(list(legend=genotype.colors, pt.colors=person.genotype.colors))
}
  
get.locus.title<-function(locus.info, locus.ind) {
  title<-sprintf("%s:%s %s repeat, reference length %s", locus.info$chr[locus.ind], prettyNum(locus.info$pos[locus.ind], big.mark=","),
                 locus.info$unit[locus.ind], locus.info$ref.length[locus.ind])
  return(title)
}

append.dn.info <- function(main.title, who.dn, fam.genotypes, genotypes, fam.id) {
  pro.score <- get.dn.score(1, fam.genotypes)
  sib.score <- get.dn.score(2, fam.genotypes)
  
  idxs          <- NULL
  dn.string     <- NULL
  if(who.dn == 1) {
    idxs <- 1:3
    dn.string <- "proband"
  } else if(who.dn == 2) {
    idxs <- c(1, 2, 4)
    dn.string <- "sibling"
  } else if(who.dn == 3) {
    idxs <- 1:4
    dn.string <- "both"
  } else {
    stop("Unrecognized who.dn code:", who.dn)    
  }
  
  main.title <- sprintf("%s\n%s %s denovo, obedience scores: %.1f (pro); %.1f (sib)",
                      main.title, fam.id, dn.string, pro.score, sib.score)
  
  return(main.title)
}

add.trio.metrics <- function(who.dn, fam.genotypes, genotypes, fam.id, palette="Blues",file.type=c("png","pdf")) {
  require(RColorBrewer)
  tab.col <- colorRampPalette(brewer.pal(brewer.pal.info[palette,'maxcolors'], palette))

  file.type <- match.arg(file.type)
    
  pro.score <- get.dn.score(1, fam.genotypes)
  sib.score <- get.dn.score(2, fam.genotypes)
  
  idxs          <- NULL
  swapped.score <- NULL
  if(who.dn == 1) {
    idxs <- 1:3
    swapped.score <- fam.genotypes[13]
  } else if(who.dn == 2) {
    idxs <- c(1, 2, 4)
    swapped.score <- fam.genotypes[14]
  } else if(who.dn == 3) {
    idxs <- 1:4
    swapped.score <- min(fam.genotypes[13],fam.genotypes[13])
  } else {
    stop("Unrecognized who.dn code:", who.dn)    
  }
  
  tab.matrix      <- matrix(0,nc=2,nr=5)
  tab.matrix[, 1] <- 1
  par(mar=c(10, 0, 13, 2) + 0.1)
  if(file.type == "png"){
    par(cex=3)
  }
  image(x=1:2, y=1:5, t(tab.matrix), col=tab.col(2), xlab="", ylab="", xaxt="n",
        yaxt="n", bty="n")
  title(main="Trio metrics", line=0.1)
  abline(h=seq(0.5, 4.5, 1), col="white", lwd=4)
  abline(v=1.5, col="white", lwd=4)
  text(y=1:5, x=1.4, labels=c("Null\nprobability", "Confidence", "Allele fit", "Noise fit", "Kinship score"),
       col="gold", font=2, adj=c(1, 0.5))

  text(y=1:5, x=1.9, labels=sapply(c(1 - prod(1 - genotypes[idxs, 11]), prod(genotypes[idxs, 3]),
                              prod(genotypes[idxs, 4]), prod(genotypes[idxs, 5]), swapped.score),
                              function(x) sprintf("%.2E",x)))
       
}

append.annotation <- function(main.title, anno) {
  anno.string <- annotation.string(anno)
  if(nchar(anno.string) > 0) {
    main.title <- sprintf("%s\n%s", anno.string, main.title)
  }
  return(main.title)
}

get.seq.axis <- function(unit, max.length) {
  split.seq <- unlist(strsplit(unit, ""))
  unit.length <- nchar(unit)
  
  complete <- max.length %/% unit.length
  incomplete <- max.length %% unit.length
  
  seq <- rep(split.seq, complete)
  seq <- c(seq, seq[seq_len(incomplete)])
  
  return(seq)
} 

get.dn.score <- function(who.dn, fam.genotypes) {
  if(who.dn == 1) {
    score <- log10(fam.genotypes[11]) * -10
  } else if(who.dn == 2) {
    score <- log10(fam.genotypes[12]) * -10
  } else if(who.dn == 3) {
    score<-min(log10(fam.genotypes[11:12])) * -10
  }
  return(score)
}

get.plot.name.preamble <- function(print.dir, locus.info, locus.ind, dn.score, zero.padding=3) {
  fname <- sprintf("%s/", print.dir)
  if(!missing(dn.score)) {
    fname <- sprintf("%s%0*d_", fname, zero.padding, ceiling(dn.score))
  }
  fname <- sprintf("%s%s_%s", fname, locus.info$chr[locus.ind], locus.info$pos[locus.ind])
  return(fname)
}

create.dn.figure.dirs <- function(denovo.dir, additional.child.cat=NULL) {
  mutation.type   <- c("omissions", "commissions")
  violation.child <- c("proband", "sibling", "both", additional.child.cat)
  context         <- c("exon", "intron", "UTR", "miRNA", "intergenic")
  
  conditional.make.dir(denovo.dir)
  
  directories <- list()
  
  for(i in mutation.type) {
    directories[[i]] <- list()
    dirname <- sprintf("%s/%s/", denovo.dir, i)
    conditional.make.dir(dirname)
    for(j in violation.child) {
      directories[[i]][[j]] <- list()
      dirname <- sprintf("%s/%s/%s/", denovo.dir, i, j)
      conditional.make.dir(dirname)
      for(k in context) {
        directories[[i]][[j]][[k]] <- sprintf("%s/%s/%s/%s/", denovo.dir, i, j, k)
        conditional.make.dir(directories[[i]][[j]][[k]])
      }
    }
  }
  return(directories)
}

get.pop.subgroup <- function(people, subgroup) {
  subgroup.info <- list(ind=(1:nrow(people)), anno=NULL, size=nrow(people))
  if(!missing(subgroup)) {
    if(grepl("parents", subgroup, ignore.case=TRUE)){
      subgroup.info$ind  <- which(people$rel.code %in% 1:2)
      subgroup.info$anno <- "parents"
      subgroup.info$size <- sum(people$rel.code %in% 1:2)
    } else if(grepl("(kids|children)", subgroup, ignore.case=TRUE)){
      subgroup.info$ind  <- which(people$rel.code %in% 3:4)
      subgroup.info$anno <- "children"
      subgroup.info$size <- sum(people$rel.code %in% 3:4)
    } else if(grepl("mother", subgroup, ignore.case=TRUE)){
      subgroup.info$ind  <- which(people$rel.code == 1)
      subgroup.info$anno <- "mothers"
      subgroup.info$size <- sum(people$rel.code == 1)
    } else if(grepl("father", subgroup, ignore.case=TRUE)) {
      subgroup.info$ind  <- which(people$rel.code == 2)
      subgroup.info$anno <- "fathers"
      subgroup.info$size <- sum(people$rel.code == 2)
    } else if(grepl("(proband|aut)", subgroup, ignore.case=TRUE)) {
      subgroup.info$ind  <- which(people$rel.code == 3)
      subgroup.info$anno <- "probands" 
      subgroup.info$size <- sum(people$rel.code == 3)
    } else if(grepl("sib", subgroup, ignore.case=TRUE)) {
      subgroup.info$ind  <- which(people$rel.code == 4)
      subgroup.info$anno <- "siblings"
      subgroup.info$size <- sum(people$rel.code == 4)
    } 
  }
  return(subgroup.info)
}

append.subgroup <- function(main.title, group.info) {
  if(length(group.info$anno) > 0) {
    main.title<-sprintf("%s %s only",main.title,group.info$anno)
  }
  return(main.title)
}

pop.total.coverage.hist <- function(cov.matrix, dirname, fig.height=2400, fig.width=2400) {
  figure.dir <- get.figure.dir(dirname)
  png(sprintf("%s/totalPerPersonCoverage.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  total.cov <- colSums(cov.matrix)
  max.cov   <- (log10(max(total.cov)) %/% 1) + 1
  h <- hist(log10(total.cov), breaks=seq(0, max.cov, 0.1))
  h$counts <- h$counts / sum(h$counts)
  plot(h, col="grey", border=NA, main=sprintf("Total per person coverage\n%d total people", ncol(cov.matrix)),
     xlab="Coverage", xaxt="n")
  axis(side=1, at=seq(0, max.cov), labels=10^(seq(0, max.cov)))
  abline(v=log10(mean(total.cov)), col="red", lwd=2)
  mtext(sprintf("%.2f", mean(total.cov)), side=3, at=log10(mean(total.cov)),
        cex=1.5, font=2)
  abline(v=log10(median(total.cov)), col="green", lwd=2)
  mtext(sprintf("%.2f", median(total.cov)), side=3, at=log10(mean(total.cov)),
        cex=1.5, font=2, line=-0.5)  
  dev.off()
  return(total.cov)
}

pop.mean.coverage.hist <- function(cov.matrix, dirname, fig.height=2400, fig.width=2400) {
  figure.dir <- get.figure.dir(dirname)
  png(sprintf("%s/meanPerPersonLocusCoverage.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  mean.cov <- colMeans(cov.matrix)
  max.cov  <- ((max(mean.cov) %/% 10) * 10) + 10
  h <- hist(mean.cov, breaks=seq(0, max.cov, 1))
  h$counts <- h$counts / sum(h$counts)
  plot(h, col="grey", border=NA, main=sprintf("Mean locus coverage per person\n%d total people", ncol(cov.matrix)),
       xlab="Coverage", xaxt="n")
  axis(side=1, at=seq(0, max.cov, 20))
  abline(v=mean(mean.cov), col="red", lwd=2)
  mtext(sprintf("%.2f", mean(mean.cov)), side=3, at=mean(mean.cov),
        cex=1.5, font=2)
  abline(v=median(mean.cov), col="green", lwd=2)
  mtext(sprintf("%.2f", median(mean.cov)), side=3, at=median(mean.cov),
        line=-0.5, cex=1.5, font=2)
  dev.off()
  return(mean.cov)
}

graphBiasInfo <- function(dirname, locus.info, fig.height=2400, fig.width=2400) {
  figure.dir <- get.figure.dir(dirname)
  
  png(sprintf("%salleleBiases.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  biases <- locus.info$bias[locus.info$allele.no > 0 & !is.na(locus.info$bias) &
                            locus.info$num.called > 0]
  max.bias <- ((max(biases) %/% 0.1) * 0.1) + 0.1
  max.plotted.bias <- if(max.bias > 3) 3 else ((max.bias %/% 0.25) * 0.25) + 0.25
  
  h <- hist(biases, breaks=seq(0, max.bias, 0.025), plot=FALSE)
  h$counts <- h$counts / sum(h$counts)
  max.y    <- ((max(h$counts) %/% 0.1) * 0.1) + 0.1
  
  plot(h, col="grey", border=NA, main="All allele biases", xlab="Allele biases",
       ylab="Frequency", xlim=c(0, max.plotted.bias), ylim=c(0, max.y), xaxt="n")
  axis(side=1, at=seq(0, max.plotted.bias, 0.25))
  dev.off()
  
  png(sprintf("%srefAlleleBias.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  
  ref.biases <- locus.info$bias[locus.info$allele.no == locus.info$ref.length & !is.na(locus.info$bias) &
                                locus.info$num.called > 0]
  max.bias   <- ((max(ref.biases) %/% 0.1) * 0.1) + 0.1
  max.plotted.bias <- ((max.bias %/% 0.25) * 0.25) + 0.25

  h <- hist(ref.biases, breaks=seq(0, max.bias, 0.025), plot=FALSE)
  h$counts <- h$counts / sum(h$counts)
  max.y    <- ((max(h$counts) %/% 0.1) * 0.1) + 0.1
  
  plot(h, col="grey", xlim=c(0, max.plotted.bias), ylim=c(0, max.y),
       main="Allele biases for reference alleles", xlab="Allele bias", 
       ylab="Frequency", border=NA, xaxt="n")
  axis(side=1, at=seq(0, max.plotted.bias, 0.25))
  dev.off()
  
  png(sprintf("%snonRefAlleleBias.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  
  non.ref.biases <- locus.info$bias[locus.info$allele.no != locus.info$ref.length &
                                    locus.info$allele.no > 0 & locus.info$num.called > 0 & 
                                    !is.na(locus.info$bias)]
  max.bias <- ((max(non.ref.biases) %/% 0.1) * 0.1) + 0.1
  max.plotted.bias <- ((max.bias %/% 0.25) * 0.25) + 0.25
  
  h <- hist(non.ref.biases, breaks=seq(0, max.bias, 0.025), plot=FALSE)
  h$counts <- h$counts / sum(h$counts)
  max.y    <- ((max(h$counts) %/% 0.1) * 0.1) + 0.1
  
  plot(h, col="grey", xlim=c(0, max.plotted.bias), ylim=c(0, max.y),
       main="Allele biases for non-reference alleles",xlab="Allele bias",
       ylab="Frequency", border=NA, xaxt="n")
  axis(side=1, at=seq(0, max.plotted.bias, 0.25))
  dev.off()
}

graphGenotypeOverviewStats<-function(dirname, f.geno) {
  graphNullFrequency(dirname, f.geno)
  graphConfidenceAndFits(dirname, f.geno)
  graphHetInfo(dirname, f.geno)
  graphEMInfo(dirname, f.geno)
}

graphNullFrequency <- function(dirname, f.geno, fig.height=2400, fig.width=2400) {
  figure.dir <- get.figure.dir(dirname)
  png(sprintf("%ssingleNullFreq.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  h <- hist(apply(f.geno$genotypes, 1, function(x) sum(x[, 2] == -1 & x[, 1] != -1)) / pop,
            breaks=seq(0, 1, 0.01), plot=FALSE)
  h$counts <- h$counts / sum(h$counts)
  max.y    <- ((max(h$counts) %/% 0.1) * 0.1) + 0.1
  plot(h, col="grey", border=NA, main="Histogram of single null frequency per locus",
       xlab="Single null frequency", xlim=c(0, 1), ylim=c(0, max.y), ylab="Frequency")
  dev.off()
  
  png(sprintf("%sdoubleNullFreq.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  h <- hist(apply(f.geno$genotypes, 1, function(x) sum(x[, 1] == -1)) / pop, 
            breaks=seq(0, 1, 0.01), plot=FALSE)
  h$counts <- h$counts / sum(h$counts)
  max.y    <- ((max(h$counts) %/% 0.1) * 0.1) + 0.1
  plot(h, col="grey", border=NA, main="Histogram of double null frequency per locus",
       xlab="Double null frequency", xlim=c(0, 1),ylim=c(0, max.y), ylab="Frequency")
  dev.off()
}

graphConfidenceAndFits <- function(dirname, f.geno, fig.height=1000, fig.width=1000) {
  figure.dir <- get.figure.dir(dirname)
  png(sprintf("%sgenotypeConfidence.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  h <- hist(f.geno$genotypes[, , 3], breaks=seq(0, 1 , 0.01), plot=FALSE)
  h$counts <- h$counts / sum(h$counts)
  max.y    <-((max(h$counts) %/% 0.1) * 0.1) + 0.1
  plot(h, col="grey", border=NA, main="Histogram of confidence scores",
       xlab="Confidence scores", ylab="Frequency", xlim=c(0, 1), ylim=c(0, max.y))
  dev.off()
  
  png(sprintf("%salleleGoodnessOfFit.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  h <- hist(f.geno$genotypes[, , 4], breaks=seq(0, 1, 0.01), plot=FALSE)
  h$counts <- h$counts / sum(h$counts)
  max.y <- ((max(h$counts) %/% 0.1) * 0.1) + 0.1
  plot(h, col="grey", border=NA, main="Histogram of coverage goodness-of-fit",
       xlab="Coverage goodness-of-fit", ylab="Frequency", xlim=c(0, 1), ylim=c(0, max.y))
  dev.off()
  
  png(sprintf("%snoiseGoodnessOfFit.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  h <- hist(f.geno$genotypes[, , 5],breaks=seq(0, 1, 0.01), plot=FALSE)
  h$counts <- h$counts / sum(h$counts)
  max.y    <- ((max(h$counts) %/% 0.1) * 0.1) + 0.1
  plot(h, col="grey", border=NA, main="Histogram of noise goodness-of-fit",
       xlab="Noise goodness-of-fit", ylab="Frequency", xlim=c(0, 1), ylim=c(0, max.y))
  dev.off()
}

graphHetInfo <- function(dirname, f.geno, fig.height=2400, fig.width=2400) {
  figure.dir <- get.figure.dir(dirname)
  png(sprintf("%snumAllelesPerLocus.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  num.alleles.per.locus <- apply(f.geno$genotypes, 1, function(x) sum(unique(unique(x[, 1]), unique(x[, 2])) > 0))
  freqs <- tabulate(num.alleles.per.locus) / length(num.alleles.per.locus)
  max.y <- ((max(freqs) %/% 0.1) * 0.1) + 0.1
  bp    <- barplot(freqs, border=NA, ylim=c(0, max.y), main="Histogram of number of alleles per locus",
                   xlab="Number of alleles", ylab="Frequency", xaxt="n", col=rgb(228, 108, 10, 200, max=255))
  axis(side=1, at=bp, labels=1:max(num.alleles.per.locus), lwd=0)
  dev.off()
  
  png(sprintf("%spctHetPerLocus.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  pct.het <- apply(f.geno$genotypes[, , 1:2], 1, function(x) sum(x[, 1] != x[, 2] & x[, 2] > 0)/pop)
  h <- hist(pct.het[pct.het > 0], breaks=seq(0, 1, 0.01), plot=FALSE)
  h$counts <- h$counts / sum(h$counts)
  max.y    <- ((max(h$counts) %/% 0.1) * 0.1) + 0.1
  plot(h, col="grey", border=NA, main="Histogram of % people heterozygous per locus\nloci with at least one heterozygous call",
       xlab="Percent heterozygous per locus", ylab="Frequency", xlim=c(0, 1), ylim=c(0, max.y))
  dev.off()
}

graphEMInfo <- function(dirname, f.geno, fig.height=2400, fig.width=2400) {
  figure.dir <- get.figure.dir(dirname)

  png(sprintf("%snumItersEM.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  freqs <- tabulate(f.geno$n.em.iter) / length(f.geno$n.em.iter)
  max.y <- ((max(freqs) %/% 0.1) * 0.1) + 0.1
  bp    <-  barplot(freqs, col="grey", border=NA, main="Histogram of EM iterations per locus",
                    xlab="Number of EM iterations", ylab="Frequency", ylim=c(0, max.y), xlim=c(0, 30))
  axis(side=1, at=bp, labels=1:max(f.geno$n.em.iter))
  dev.off()

  png(sprintf("%serrorRates.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  h <- hist(f.geno$locus.error.rates, breaks=seq(0, 1, 0.01), plot=FALSE)
  h$counts <- h$counts / sum(h$counts)
  max.y <- ((max(h$counts) %/% 0.1) * 0.1) + 0.1
  plot(h, col="grey", border=NA, main="Histogram of EM-derived error rate",
       xlab="Error rate", ylab="Frequency", xlim=c(0, 1), ylim=c(0, max.y))
  dev.off()
  
  png(sprintf("%slog10ErrorRates.png", figure.dir), height=fig.height, width=fig.width)
  par(cex=4)
  min.log.err <- ((min(log10(f.geno$locus.error.rates)) %/% 10) * 10) - 10
  h <- hist(log10(f.geno$locus.error.rates), breaks=seq(min.log.err, 0, 0.1), plot=FALSE)
  h$counts <- h$counts / sum(h$counts)
  min.x    <- if(min.log.err < -10) -10 else min.log.err
  max.y    <- ((max(h$counts) %/% 0.1) * 0.1) + 0.1
  plot(h, col="grey", border=NA, main="Histogram of log10 EM-derived error rate",
       xlab="log10 error rate", ylab="Frequency", xlim=c(min.x, 0), ylim=c(0, max.y), xaxt="n")
  axis(side=1, at=seq(-10, 0, 2), labels=10^seq(-10, 0, 2))
  dev.off()  
}

get.figure.dir <- function(dirname) {
  figure.dir <- sprintf("%s/figures/", dirname)
  conditional.make.dir(figure.dir)
  return(figure.dir)
}

#################################################
# OBSOLETE PLOTTING FUNCTIONS                   #
#################################################

# plot.pop.scatter<-function(genotypes,pop.size,locus.ind,locus.info,allele.info,print.dir=NULL,new.display=TRUE,HWE.p.val=NULL,fam.genotypes=NULL,who.dn=NULL,anno=NULL)
# {
#   require(ggplot2)
#   require(fields)
#   require(graphics)
#   require(colorRamps)
#   
#   la.ind<-which(locus.info$locus.row.ind == locus.ind & locus.info$allele.no == 0)
#   
#   locus<-data.frame(al.1=genotypes[locus.ind,,6],rest=rowSums(genotypes[locus.ind,,6:8]) - genotypes[locus.ind,,6],
#                     Genotype=sapply(1:pop.size,function(x) sprintf("%d|%d",genotypes[locus.ind,x,1],genotypes[locus.ind,x,2])))
#   
#   max.count<- 1.025*max(locus$al.1,locus$rest)
#   min.count<- -0.025*max.count
#   dn.score<-NULL
#   
#   if(length(who.dn) > 0)
#   {
#     dn.score<-get.dn.score(who.dn,fam.genotypes)
#   }
#   
#   p<-ggplot(locus,aes(x=rest,y=al.1))
#   p1<-p + geom_jitter(aes(colour=Genotype))
#   p1<-p1 + coord_cartesian(ylim=c(min.count,max.count),xlim=c(min.count,max.count))
#   p1<-p1 + ylab("Allele 1 coverage") + xlab("Remaining coverage")
# 
#   main.title<-get.locus.title(locus.info,la.ind)
#   if(length(anno) > 0){main.title<-append.annotation(main.title,anno)}
#   if(length(who.dn) > 0){main.title<-append.dn.info(main.title,who.dn,fam.genotypes)}
#   
#   if(length(HWE.p.val) > 0)
#   {
#     main.title<-sprintf("%s\nHWE p-value: %.3E",main.title,HWE.p.val)
#   }
#   
#   p1<-p1 + opts(title=main.title)
#   
#   if(is.character(print.dir))
#   {
#     fname<-sprintf("%s_genoFreq.png",get.plot.name.preamble(print.dir,locus.info,la.ind,dn.score))
#     png(fname,width=1400,height=1400)
#     p1<-p1 + opts(theme_text=20)
#   } else if(new.display==TRUE)
#   {
#     quartz()
#   }
#   print(p1)
#   if(is.character(print.dir))
#   {
#     dev.off()
#   }
# }
