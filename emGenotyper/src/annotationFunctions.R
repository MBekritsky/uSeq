get.annotations <- function(exon.file, intron.file=NULL, mirna.file=NULL, utr.file=NULL,
                            locus.info, verbose=1, cran.repos){
  # return microsatellite loci annotations for exon information, optionally also intron,
  # miRNA and UTR information files must be in BED format with the feature ID in column 4
  # NB: takes a matrix in the locus.info format as input
  if(verbose >= 1) cat("Annotating loci...\n")
  ptm <- proc.time()
  zero.indices <- which(locus.info$allele.no == 0)
  chr.pos.locus.info <-do.call("paste", list(as.character(locus.info[zero.indices, 1]),
                                             locus.info[zero.indices, 2], sep="_"))
  
  # find microsatellite loci that overlap with exons
  annotations <- data.frame(exon=feature.annotate(exon.file, chr.pos.locus.info))
  
  if(length(intron.file) > 0) annotations$intron <- feature.annotate(intron.file,
                                                                     chr.pos.locus.info)
  if(length(mirna.file) > 0)  annotations$mirna  <- feature.annotate(mirna.file,
                                                                     chr.pos.locus.info)
  if(length(utr.file) > 0)    annotations$utr    <- feature.annotate(utr.file,
                                                                     chr.pos.locus.info)
  
  
  if(verbose >= 1) cat("Finished annotating loci in", (proc.time() - ptm)[3], "seconds\n")
  return(annotations)
}

feature.annotate <- function(bed.filename, chr.pos.locus.info, verbose=1) {
  # annotate microsatellites loci with a particular feature
  # files must be in BED format with the feature ID in column 4
  ptm <- proc.time()
  anno.feature <- read.table(bed.filename, stringsAsFactors = FALSE)
  names(anno.feature) <- c("chr", "start", "stop", "info")
  num.loci            <- length(chr.pos.locus.info)
  ms.and.feature      <- rep(NA, num.loci)

  # create strings of (chr)_(position) to match between BED annotation file
  # (provided by bed.filename) and locus.info
  chr.pos.anno.feature <- do.call("paste", list(as.character(anno.feature[, 1]),
                                                anno.feature[, 2], sep="_"))

  # get the indices of loci that overlap annotation
  has.feature <- which(chr.pos.locus.info %in% chr.pos.anno.feature)
  #get the annotation that each locus has each locus should only have one
  # entry in the BED file)
  feature.contained <- which(chr.pos.anno.feature %in% chr.pos.locus.info)
  
  # in case a locus has more than one feature, merge features by microsatellite
  # chromosome and position, so that each microsatellite locus has only one
  # feature matching it
  final.features <- anno.feature[feature.contained, ]
  
  # assign features to microsatellite loci
  ms.and.feature[has.feature] <- final.features[, 3]

  if(verbose >= 2) cat("Annotated with", bed.filename, "in", (proc.time() - ptm)[3], "seconds\n")
  rm(anno.feature)
  return(ms.and.feature)
}

load.gene.ids <- function(gene.id.file) {
  # load gene ID table to translate from CCDS and refGene IDs to descriptive names
  # the rows of gene.id.file should have the format <refGene ID>, <geneID>, <CCDS ID>
  # and should be tab-delimited.  The file should not have a header
  gene.ids <- read.table(gene.id.file)
  names(gene.ids) <- c("refGene", "geneID", "CCDS")
  return(gene.ids)
}

get.locus.annotations <-function(anno.idx, annotations, gene.ids) {
  # get all translated annotations for a locus
  l.anno <- list(exon=translate.exon(annotations[anno.idx, ], gene.ids))
  if(length(annotations$intron) > 0) l.anno$intron <- translate.intron(annotations[anno.idx, ],
                                                                       gene.ids)
  # miRNAs do not have entries in gene.ids
  if(length(annotations$mirna) > 0)  l.anno$mirna  <- translate.mirna(annotations[anno.idx, ])
  if(length(annotations$utr) > 0)    l.anno$utr    <- translate.utr(annotations[anno.idx, ],
                                                                    gene.ids)
  return(l.anno)
}

translate.exon   <- function(locus.annotation, gene.ids) {
  # translate exon refGene and CCDS IDs to descriptive IDs
  gene.info <- get.anno.info(as.character(locus.annotation$exon), gene.ids)
  return(unique(gene.info))
}
  
translate.intron <- function(locus.annotation, gene.ids) {  
  # translate intron refGene and CCDS IDs to descriptive IDs
  gene.info <- get.anno.info(as.character(locus.annotation$intron), gene.ids)
  return(gene.info)
}

translate.mirna  <- function(locus.annotation) {
  # translate miRNA ID to descriptive ID
  if(!is.na(locus.annotation$mirna)) {
    mirna.names <- sapply(strsplit(locus.annotation$mirna, ";|,"),
                          function(x) grep("Name", x, value=TRUE, fixed=TRUE))
    mirna.names <- unlist(sapply(mirna.names, function(x) gsub("Name=hsa-", "", x)))
    gene.info   <- data.frame(geneID=mirna.names, feature.no=rep(1, length(mirna.names)),
                              row.names=NULL)
    return(gene.info)
  } else {
    return(NA)
  }
}

translate.utr    <- function(locus.annotation, gene.ids) {
  # translate UTR refGene IDs to descriptive IDs
  # NB teh version of CCDS used for these annotations
  # does NOT contain any UTR annotations, therefore, all
  # annotations are refGene annotations
  
  if(!is.na(locus.annotation$utr)) {
    genes     <- unlist(strsplit(locus.annotation$utr, ";|,"))
    gene.info <- data.frame(geneID=rep(NA, length(genes)), 
                            feature.no=rep(NA, length(genes)),
                            which=rep(NA, length(genes)))
    
    for(i in seq_along(genes)) {
      gene.vec            <- unlist(strsplit(genes[i], "_"))
      rg.gene.id          <- sprintf("%s_%s", gene.vec[1], gene.vec[2])
      gene.name.idx       <- which(gene.ids$refGene == rg.gene.id)[1]
      gene.info$geneID[i] <- as.character(gene.ids$geneID[gene.name.idx])
      
      # refGene IDs with an NR prefix are non-coding RNAs
      if(gene.vec[1] == "NR") gene.info$geneID[i] <- sprintf("%s (ncRNA)", 
                                                             gene.info$geneID[i])
      
      gene.info$feature.no[i] <- as.integer(gene.vec[4])
      # which reports which UTR is overlapped (5' or 3')
      gene.info$which[i]      <- substr(gene.vec[3], 4, 4)
    }
    return(unique(gene.info)) 
  } else {
    return(NA)
  }
}

get.anno.info <- function(feature.annotation, gene.ids) {
  # function that handles annotation when both CCDS and refGene IDs are
  # possible (i.e. exons and introns)
  if(!(is.na(feature.annotation))) {
    genes     <- unlist(strsplit(feature.annotation, ";|,"))
    gene.info <- data.frame(geneID=rep(NA, length(genes)),
                            feature.no=rep(NA, length(genes)))
    
    for(i in 1:length(genes)) {
      if(!(is.na(pmatch("CCDS", genes[i])))) {
        # if any of the annotations begin with CCDS
        gene.vec                <- unlist(strsplit(genes[i], "_"))

        # force uniqueness if there are more than 2 gene names matching the CCDS ID
        gene.name.idx           <- which(gene.ids$CCDS == gene.vec[1])[1] 
        gene.info$geneID[i]     <- as.character(gene.ids[gene.name.idx,"geneID"])
        gene.info$feature.no[i] <- as.integer(gene.vec[3]) + 1 
        # CCDS feature numbers are zero-indexed
      }
      if(!(is.na(pmatch("NM", genes[i]))) | 
         !(is.na(pmatch("NR", genes[i])))) {
        # if any of the annotations begin with an NM or an NR
        # NB: NR annotations are only found in introns or UTRs
        # not exons
        gene.vec                <- unlist(strsplit(genes[i], "_"))
        rg.gene.id              <- sprintf("%s_%s", gene.vec[1], gene.vec[2])
        gene.name.idx           <- which(gene.ids$refGene == rg.gene.id)[1]
        gene.info$geneID[i]     <- as.character(gene.ids$geneID[gene.name.idx])
        if(gene.vec[1] == "NR") {
          gene.info$geneID[i]   <- sprintf("%s (ncRNA)", gene.info$geneID[i])
        }
        gene.info$feature.no[i] <- as.integer(gene.vec[4])
      }
    }
    return(unique(gene.info))
  } else {
    return(NA)
  }
}

annotation.string <- function(locus.annotation, max.anno=3) {
  # return a string of all annotations for a locus up to a maximum
  # of max.anno
  total.anno  <- 0
  anno.string <- ""
  if(sum(is.na(locus.annotation$exon)) == 0) {
    names <- unique(locus.annotation$exon$geneID)
    for(i in names) {
      if(total.anno <= max.anno) {
        if(total.anno > 0) {
          anno.string <- sprintf("%s; ", anno.string)
        }
        anno.string <- sprintf("%s%s exon", anno.string, i)
        exons       <- sort(locus.annotation$exon$feature.no[which(locus.annotation$exon == i)])
        ifelse(length(exons) == 1, anno.string <- sprintf("%s %d", anno.string, exons),
                                   anno.string <- sprintf("%ss %s", anno.string, paste(exons, collapse="/"))
               )
        total.anno <- total.anno + 1
      }
    }
  }
  
  if(sum(is.na(locus.annotation$intron)) == 0) {
    names <- unique(locus.annotation$intron$geneID)
    for(i in names) {
      if(total.anno <= max.anno) {
        if(total.anno > 0) {
          anno.string <- sprintf("%s; ", anno.string)
        }
        anno.string <- sprintf("%s%s intron", anno.string, i)
        introns     <- sort(locus.annotation$intron$feature.no[which(locus.annotation$intron == i)])
        ifelse(length(introns) == 1, anno.string <- sprintf("%s %d", anno.string, introns),
                                     anno.string <- sprintf("%ss %s", anno.string, paste(introns, collapse="/"))
        )
        total.anno<-total.anno + 1
      }
    }
  }
  
  if(sum(is.na(locus.annotation$utr)) == 0) {
    names <- unique(locus.annotation$utr$geneID)
    for(i in names) {
      if(total.anno <= max.anno) {
        if(total.anno > 0) {
          anno.string <- sprintf("%s; ", anno.string)
        }
        end <- locus.annotation$utr$which[which(locus.annotation$utr == i)][1]
        anno.string <- sprintf("%s%s %s\' UTR", anno.string, i, end)
        total.anno  <- total.anno + 1
      }
    }
  }
  
  if(sum(is.na(locus.annotation$mirna)) == 0) {
    names <- unique(locus.annotation$mirna$geneID)
    for(i in names) {
      if(total.anno <= max.anno) {
        if(total.anno > 0) {
          anno.string <- sprintf("%s; ", anno.string)
        }
        anno.string <- sprintf("%s%s ", anno.string, i)
        total.anno  <- total.anno + 1
      }
    }
  }
  
  if(total.anno == 0) anno.string <- "Intergenic"
  return(anno.string)
}

annotation.precedence<-function(locus.annotation) {
  #get annotation precedence
  # exon > intron > UTR > miRNA > intergenic
  ifelse(sum(is.na(locus.annotation$exon))   == 0, return("exon"),
  ifelse(sum(is.na(locus.annotation$intron)) == 0, return("intron"),
  ifelse(sum(is.na(locus.annotation$utr))    == 0, return("UTR"),
  ifelse(sum(is.na(locus.annotation$mirna))  == 0, return("miRNA"),
                                                   return("intergenic")
  ))))
}
