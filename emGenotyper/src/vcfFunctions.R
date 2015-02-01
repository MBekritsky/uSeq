info.fields <- function(){
  vcf.info <- data.frame(ID=c("END","RL","RU","POP","SUM","TOP","NR"), stringsAsFactors=F)
  vcf.info$Number <- 1
  vcf.info$Type   <- c("Integer", "Integer", "String", "Integer", "Integer", "Integer", "Float")
  vcf.info$Description <- c("End of microsatellite tract", "Reference length of microsatellite tract", 
                            "Reference motif of microsatellite tract", "Number of people with locus in population",
                            "Total locus coverage across population", "Maximum locus coverage across population",
                            "log10 of locus-wide error rate")
  return(vcf.info)
}

format.fields <- function(){
  format.info <- data.frame(ID=c("GT", "GB", "GC", "EC", "AF", "NF", "AL1", "AL2", "ALN", "PNULL", "DP", "SEEN"),
                            stringsAsFactors=F)
  format.info$Number <- 1
  format.info$Type <- c("String", "String", "Float", "Float", "Float", "Float", "Integer", "Integer", "Integer",
                     "Float", "Integer", "String")
  format.info$Number[format.info$ID == "SEEN"] <- 2
  format.info$Description <- c("Genotype","Genotype in bp (-1 means null allele)", "Genotype confidence",
                            "Expected per-allele coverage", "Exact Poisson allele fit",
                            "Exact binomial noise fit", "Coverage for allele 1",
                            "Coverage for allele 2","Noise coverage", 
                            "-10log10 of marginal likelihood for null allele", "Read depth", "Genotyped alleles")
  return(format.info)
}

print.info.or.fields <- function(meta.table, vcf.file="", desc=c("INFO","FORMAT")){
  if(missing(meta.table))
    stop("No table was passed to print to the VCF file")

  desc <- match.arg(desc)
  
  apply(meta.table, 1, function(x) cat(sprintf("##%s=<ID=%s,Number=%d,Type=%s,Description=\"%s\">\n", desc,
                                                  x['ID'], as.integer(x['Number']), x['Type'], x['Description']),
                                       file=vcf.file, append=T))
}

print.vcf.meta.information <- function(vcf.file,person.id,family.id) {
  cat("##fileformat=VCFv4.1\n",file=vcf.file)
  cat("##fileDate=",strftime(Sys.time()),"\n",file=vcf.file,sep="",append=T)
  cat("##source=uSeq_v1\n",file=vcf.file,append=T)
  cat("##personInfo=<IndividualID=",person.id,",FamilyID=",family.id,file=vcf.file,">\n",sep="",
      append=T)
  print.info.or.fields(info.fields(), vcf.file, desc="INFO")
  print.info.or.fields(format.fields(), vcf.file, desc="FORMAT")
  # print header line
  cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tuSeq\n",file=vcf.file, append=T)
}

get.format.string <- function(){
  format.fields <- format.fields()
  return(paste(format.fields$ID,collapse=":"))
}

get.ref.string <- function(locus.info){
  ref.strings <- apply(locus.info[,c("unit","ref.length")], 1, function(x) sprintf("%s%s",
                       paste(rep(x["unit"], as.integer(x["ref.length"]) %/% nchar(x["unit"])), collapse=""), 
                       substr(x["unit"], 1, as.integer(x["ref.length"]) %% nchar(x["unit"]))))
  return(ref.strings)
}

info.table.2.string <- function(info.table){
  info.table$END <- sub("^","END=",as.character(info.table$END))
  info.table$RL  <- sub("^","RL=",as.character(info.table$RL))
  info.table$RU  <- sub("^","RU=",info.table$RU)
  info.table$POP <- sub("^","POP=",as.character(info.table$POP))
  info.table$SUM <- sub("^","SUM=",as.character(info.table$SUM))
  info.table$TOP <- sub("^","TOP=",as.character(info.table$TOP))
  info.table$NR  <- sub("^","NR=",sapply(info.table$NR, function(x) sprintf("%0.2f",log10(x))))
  
  strings <- apply(info.table, 1, function(x) paste(x, collapse=";"))
  
  return(strings)
}

useq.table.2.string <- function(useq.table){
  useq.table$GT  <- as.character(useq.table$GT)
  useq.table$GB  <- as.character(useq.table$GB)
  useq.table$GC  <- sapply(useq.table$GC, function(x) sprintf("%.3f",x))
  useq.table$EC  <- sapply(useq.table$EC, function(x) sprintf("%.3f",x))
  useq.table$AF  <- sapply(useq.table$AF, function(x) sprintf("%.3f",x))
  useq.table$NF  <- sapply(useq.table$NF, function(x) sprintf("%.3f",x))
  useq.table$AL1 <- as.character(useq.table$AL1)
  useq.table$AL2 <- as.character(useq.table$AL2)
  useq.table$ALN <- as.character(useq.table$ALN)
  useq.table$PNULL <- sapply(useq.table$PNULL, function(x) sprintf("%.3f",x))
  useq.table$DP <- as.character(useq.table$DP)
  useq.table$SEEN <- as.character(useq.table$SEEN)
  
  strings <- apply(useq.table[,c("GT","GB","GC","EC","AF","NF","AL1","AL2","ALN","PNULL","DP","SEEN")], 1, function(x) paste(x, collapse=":"))
  
  return(strings)
}