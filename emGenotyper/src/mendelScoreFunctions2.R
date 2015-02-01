# Requires genotyperFunctions2.R as a source file
get.mendel.trio.indices <- function(n.alleles) {
  # get indices for all possible Mendelian inheritance patterns up to n.alleles
  
  # enumerate all possible genotypes
  genos <- expand.grid(1:n.alleles, 1:n.alleles)
  
  # order genotypes s.t. most common allele (allele with higher index number) is
  # always first
  genos <- t(apply(genos, 1, sort))
  # switch genotype structure from list to matrix
  genos <- matrix(unlist(unique(genos)), ncol=2)
  
  # enumerate all unique parental genotype combinations
  p.geno.pairs <- expand.grid(1:nrow(genos), 1:nrow(genos))
  p.geno.pairs <- t(apply(p.geno.pairs, 1, sort))
  p.geno.pairs <- unique(p.geno.pairs)

  p.geno.alleles <- cbind(genos[p.geno.pairs[, 1], 1], genos[p.geno.pairs[, 1], 2],
                          genos[p.geno.pairs[, 2], 1], genos[p.geno.pairs[, 2], 2])
  p.geno.alleles <- rbind(p.geno.alleles, cbind(p.geno.alleles[, 3:4], p.geno.alleles[, 1:2]))
  p.geno.alleles <- p.geno.alleles[order(p.geno.alleles[, 1], p.geno.alleles[, 2],
                                         p.geno.alleles[, 3], p.geno.alleles[, 4]), ]
  p.geno.alleles <- unique(p.geno.alleles)
  
  # enumerate all possible unique Mendelian offspring for all parental genotype combinations
  mendel.geno.ind <- numeric(0)
  for(i in seq_len(nrow(p.geno.alleles))) {
    num.par.1.alleles  <- length(unique(p.geno.alleles[i, 1:2]))
    num.par.2.alleles  <- length(unique(p.geno.alleles[i, 3:4]))
    mendel.pattern.ind <- matrix(unlist(expand.grid(1:num.par.1.alleles, 1:num.par.2.alleles)),ncol=2)
    mendel.patterns    <- cbind(p.geno.alleles[i, mendel.pattern.ind[, 1]], 
                                p.geno.alleles[i, mendel.pattern.ind[, 2] + 2])
    mendel.patterns    <- t(apply(mendel.patterns, 1, sort))
      
    mendel.geno.ind    <- rbind(mendel.geno.ind, cbind(matrix(rep(p.geno.alleles[i, ],nrow(mendel.patterns)),
                                                              ncol=4, byrow=T), mendel.patterns))
  }
  mendel.geno.ind <- mendel.geno.ind[order(mendel.geno.ind[, 1], mendel.geno.ind[, 2],
                                           mendel.geno.ind[, 3], mendel.geno.ind[, 4],
                                           mendel.geno.ind[, 5], mendel.geno.ind[, 6]), ]
  mendel.geno.ind <- unique(mendel.geno.ind)
  
  return(mendel.geno.ind)
}

get.possible.genos<-function(n.alleles, max.alleles) {
  possible.alleles <- ifelse(n.alleles > max.alleles, max.alleles,
                             n.alleles)
  possible.genos      <- matrix(0, nrow=sum(seq_len(possible.alleles)), ncol=2)
  possible.genos[, 1] <- rep.int(1:possible.alleles, possible.alleles:1)
  possible.genos[, 2] <- unlist(lapply(seq_len(possible.alleles), function(i) i:possible.alleles))
  return(possible.genos)
}

get.fam.geno.prob <- function(i, p.geno, by.fam, likely, n.alleles) {
    return(c(prod(p.geno[cbind(by.fam[i, seq_len(3)], matrix(likely[i, seq_len(6)], ncol=2, byrow=T))]),
             prod(p.geno[cbind(by.fam[i, c(1, 2, 4)], matrix(likely[i, c(seq_len(4), 7, 8)], ncol=2, byrow=T))])))
}

get.fam.sorted.alleles<-function(l.geno, by.fam) {
  n.fams    <- nrow(by.fam)
  n.alleles <- dim(l.geno)[2]
  # get each person's probability for having every allele we observe
  allele.freqs <- (apply(l.geno, 1:2, sum) + apply(l.geno, c(1,3), sum)) / 2
  
  pro.sorted <- apply(by.fam, 1, function(x) order(colMeans(allele.freqs[x[1:3],      ]),
                                                   decreasing=TRUE))
  sib.sorted <- apply(by.fam, 1, function(x) order(colMeans(allele.freqs[x[c(1,2,4)], ]),
                                                   decreasing=TRUE))

  return(list(pro=pro.sorted, sib=sib.sorted))
}

get.max.ind<-function(x) {
  if(any(!is.nan(x))) {
    return(which(x == max(x), arr.ind=T)[1, ])
  } else {
    return(c(-2, -2))
  }
}

get.likely.geno.codes <- function(p.geno, n.fams, max.children, by.fam) {
  likely.geno.codes        <- array(0, dim=c(n.fams, max.children * 2))
  likely.geno.codes[, 1:2] <- t(apply(p.geno[by.fam[, 1], , ], 1, get.max.ind))
  likely.geno.codes[, 3:4] <- t(apply(p.geno[by.fam[, 2], , ], 1, get.max.ind))
  likely.geno.codes[, 5:6] <- t(apply(p.geno[by.fam[, 3], , ], 1, get.max.ind))
  likely.geno.codes[, 7:8] <- t(apply(p.geno[by.fam[, 4], , ], 1, get.max.ind))
  return(likely.geno.codes)
}

get.mendel.scores <- function(allele.info, locus.info, geno.list, n.fams, by.fam,
                              pop.size, coverage.model, mendel.ind, max.locus=NULL,
                              start.locus=1, debug=FALSE) {
  n.ind <- nrow(mendel.ind)
  ptm   <- proc.time()
  if(missing(max.locus)) {
    max.locus     <- allele.info$num.loci
    fam.genotypes <- array(0, c(allele.info$num.loci, n.fams, 14))
  } else {
    fam.genotypes <- array(0, c(max.locus - start.locus + 1, n.fams, 14))
  }
  
  field.names<-c("mom.allele.1", "mom.allele.2", "dad.allele.1", "dad.allele.2",
                 "pro.allele.1", "pro.allele.2", "sib.allele.1", "sib.allele.2",
                 "pro.trio.likelihood", "sib.trio.likelihood", "pro.mendel.obedience.score",
                 "sib.mendel.obedience.score", "pro.kinship.1", "sib.kinship.1",
                 "pro.kinship.2", "sib.kinship.2")
  dimnames(fam.genotypes) <- list(NULL, NULL, field.names)
  
  for(i in start.locus:max.locus) {
    if(debug) print(i)

    # get indices in allele and locus matrices for locus and allele coverage, as well as number of alleles
    locus.details <- get.locus.details(i, locus.info)
    # get coverage information for locus and for each allele
    coverage.info <- get.coverage(allele.info, locus.details, pop.size)
    
    # get matrices for genotyper
    mats <- get.obs.exp.mat(pop=pop.size, n.alleles=locus.details$n.alleles, coverage.info=coverage.info, 
                            locus.coverage.model=coverage.model$exp.cov[i, ], allele.biases=geno.list$allele.biases[[i]][, 2])
    geno.likelihood <- get.geno.probs(mats, geno.list$locus.error.rates[i], geno.list$weights[[i]])
    obs.alleles     <- c(locus.info$allele.no[locus.details$ordered.lengths], -1)
    
    # Get the likely genotype codes for each family
    max.children <- ncol(by.fam)
    likely.geno.codes <- get.likely.geno.codes(geno.likelihood, n.fams, max.children, by.fam)
    fam.genotypes[(i - start.locus + 1), , 1:8] <- matrix(obs.alleles[likely.geno.codes], ncol=8)

    # Get joint probability for all genotypes in family (currently, probability for the child is not conditional on the parents' genotypes)
    fam.genotypes[(i - start.locus + 1), , 9:10] <- matrix(unlist(lapply(seq_len(n.fams), get.fam.geno.prob, p.geno=geno.likelihood, 
                                                                         by.fam=by.fam, likely=likely.geno.codes, 
                                                                         n.alleles=locus.details$n.alleles)),
                                                           nrow=n.fams,byrow=TRUE)
    
    possible.genos <- get.possible.genos(locus.details$n.alleles, max.alleles)
    n.genotypes    <- nrow(possible.genos)
    #Get alleles sorted by frequency within family
    fam.sorted.alleles <- get.fam.sorted.alleles(geno.likelihood, by.fam)
    pro.sorted.alleles <- fam.sorted.alleles$pro
    sib.sorted.alleles <- fam.sorted.alleles$sib
    
    # get possible genotypes for trio including proband or sibling
    # based on the most common alleles within the trio
    pro.trio.genos <- array(0, dim=c(n.fams, n.genotypes, 2))
    sib.trio.genos <- array(0, dim=c(n.fams, n.genotypes, 2))
    for(j in seq_len(n.fams)) {
      pro.trio.genos[j, , ] <- t(apply(possible.genos, 1, function(x) c(pro.sorted.alleles[x[1], j], pro.sorted.alleles[x[2] ,j])))
      pro.trio.genos[j, , ] <- t(apply(pro.trio.genos[j, , ], 1, sort))
      
      sib.trio.genos[j, , ] <- t(apply(possible.genos, 1, function(x) c(sib.sorted.alleles[x[1], j], sib.sorted.alleles[x[2], j])))
      sib.trio.genos[j, , ] <- t(apply(sib.trio.genos[j, , ], 1, sort))
    }
   
    
    # calculate the probands' mendel scores
    mendel.pro <- array(0, dim=c(n.fams, 3, n.ind))
    fam.array  <- array(0, dim=c(n.fams, 3, max.alleles, max.alleles))
    for(k in seq_len(3)) {
      #for every family member, set likelihoods in fam.array in order of allele's prevalence in that family
      #for each family member, index into genotype array for each possible Mendelian genotype they can have
      source.ind <- cbind(rep.int(by.fam[, k], rep.int(n.genotypes, n.fams)), 
                          t(matrix(t(pro.trio.genos[, , 1]), nrow=1)),
                          t(matrix(t(pro.trio.genos[, , 2]), nrow=1)))
      target.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.genotypes, n.fams)),
                          rep.int(k, n.fams * n.genotypes), 
                          rep.int(possible.genos[, 1], n.fams),
                          rep.int(possible.genos[, 2], n.fams))
      fam.array[target.ind] <- geno.likelihood[source.ind]
    }  
    
    for(k in seq_len(3)) {
      mendel.col    <- c((2 * k) - 1, (2 * k))
      mendel.t.ind  <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)),
                             rep.int(k, n.fams * n.ind), 
                             rep.int(seq_len(n.ind), n.fams))
      fam.array.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)),
                             rep.int(k, n.fams * n.ind), 
                             rep.int(mendel.ind[, mendel.col[1]], n.fams),
                             rep.int(mendel.ind[, mendel.col[2]], n.fams))
      mendel.pro[mendel.t.ind] <- fam.array[fam.array.ind]
    }

    pro.mendel.prob <- mendel.pro[, 1, ] * mendel.pro[, 2, ] * mendel.pro[, 3, ]
    fam.genotypes[(i - start.locus + 1), , 11] <- rowSums(pro.mendel.prob) 
    
    rm(pro.mendel.prob, mendel.pro)
    
    # switch mom and proband, get mendel score
    fam.array  <- array(0, dim=c(n.fams, 3, max.alleles, max.alleles))
    for(k in seq_len(3)) {
      #for every family member, set likelihoods in fam.array in order of allele's prevalence in that family
      #for each family member, index into genotype array for each possible Mendelian genotype they can have
      target.adj <- if(k == 1) 3 else if (k == 3) 1 else k
      source.ind <- cbind(rep.int(by.fam[, k], rep.int(n.genotypes, n.fams)), 
                          t(matrix(t(pro.trio.genos[, , 1]), nrow=1)),
                          t(matrix(t(pro.trio.genos[, , 2]), nrow=1)))
      target.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.genotypes, n.fams)),
                          rep.int(target.adj, n.fams * n.genotypes), 
                          rep.int(possible.genos[, 1], n.fams),
                          rep.int(possible.genos[, 2], n.fams))
      fam.array[target.ind] <- geno.likelihood[source.ind]
    }
    
    switch.mom.for.pro <- array(0, dim=c(n.fams, 3, n.ind))
    for(k in seq_len(3)) {
      mendel.col    <- c((2 * k) - 1, (2 * k))
      mendel.t.ind  <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)),
                             rep.int(k, n.fams * n.ind), 
                             rep.int(seq_len(n.ind), n.fams))
      fam.array.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)),
                             rep.int(k, n.fams * n.ind), 
                             rep.int(mendel.ind[, mendel.col[1]], n.fams),
                             rep.int(mendel.ind[, mendel.col[2]], n.fams))
      switch.mom.for.pro[mendel.t.ind] <- fam.array[fam.array.ind]
    
    }
    mom.switched.prob<-switch.mom.for.pro[, 1, ] * switch.mom.for.pro[, 2, ] * switch.mom.for.pro[, 3, ]
    
    # switch dad and proband, get mendel score
    fam.array  <- array(0, dim=c(n.fams, 3, max.alleles, max.alleles))
    for(k in seq_len(3)) {
      #for every family member, set likelihoods in fam.array in order of allele's prevalence in that family
      #for each family member, index into genotype array for each possible Mendelian genotype they can have
      target.adj <- if(k == 2) 3 else if (k == 3) 2 else k
      source.ind <- cbind(rep.int(by.fam[, k], rep.int(n.genotypes, n.fams)), 
                          t(matrix(t(pro.trio.genos[, , 1]), nrow=1)),
                          t(matrix(t(pro.trio.genos[, , 2]), nrow=1)))
      target.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.genotypes, n.fams)),
                          rep.int(target.adj, n.fams * n.genotypes), 
                          rep.int(possible.genos[, 1], n.fams),
                          rep.int(possible.genos[, 2], n.fams))
      fam.array[target.ind] <- geno.likelihood[source.ind]
    }
    
    switch.dad.for.pro <- array(0, dim=c(n.fams, 3, n.ind))
    for(k in seq_len(3)) {
      mendel.col    <- c((2 * k) - 1, (2 * k))
      mendel.t.ind  <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)),
                             rep.int(k, n.fams * n.ind), 
                             rep.int(seq_len(n.ind), n.fams))
      fam.array.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)),
                             rep.int(k, n.fams * n.ind), 
                             rep.int(mendel.ind[, mendel.col[1]], n.fams),
                             rep.int(mendel.ind[, mendel.col[2]], n.fams))
      switch.dad.for.pro[mendel.t.ind] <- fam.array[fam.array.ind]
      
    }
    dad.switched.prob<-switch.dad.for.pro[, 1, ] * switch.dad.for.pro[, 2, ] * switch.dad.for.pro[, 3, ]
    
    fam.genotypes[(i - start.locus + 1), , 13] <- apply(cbind(rowSums(dad.switched.prob),rowSums(mom.switched.prob)),1,max)
    rm(dad.switched.prob, mom.switched.prob, switch.mom.for.pro, switch.dad.for.pro)

    #calculate the siblings' mendel scores
    mendel.sib <- array(0, dim=c(n.fams, 3, n.ind))
    fam.array  <- array(0, dim=c(n.fams, 3, max.alleles, max.alleles))
    for(k in c(1,2,4)) {
      target.adj <- if(k == 4) 3 else k
      #for every family member, set likelihoods in fam.array in order of allele's prevalence in that family
      source.ind <- cbind(rep.int(by.fam[, k], rep.int(n.genotypes, n.fams)),
                          t(matrix(t(sib.trio.genos[, , 1]), nrow=1)),
                          t(matrix(t(sib.trio.genos[, , 2]), nrow=1)))

      target.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.genotypes, n.fams)),
                          rep.int(target.adj, n.fams * n.genotypes),
                          rep.int(possible.genos[, 1], n.fams), 
                          rep.int(possible.genos[, 2], n.fams))
      fam.array[target.ind]<-geno.likelihood[source.ind]
    }  
    
    for(k in seq_len(3)) {
      mendel.col <- c((2 * k) - 1, (2 * k))
      mendel.t.ind  <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)), 
                             rep.int(k, n.fams * n.ind), rep.int(seq_len(nrow(mendel.ind)), n.fams))
      fam.array.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)),
                             rep.int(k, n.fams * n.ind),
                             rep.int(mendel.ind[, mendel.col[1]], n.fams),
                             rep.int(mendel.ind[, mendel.col[2]], n.fams))
      mendel.sib[mendel.t.ind] <- fam.array[fam.array.ind]
    }
        
    sib.mendel.prob<-mendel.sib[, 1, ] * mendel.sib[, 2,] * mendel.sib[, 3, ]
    fam.genotypes[(i - start.locus + 1), , 12] <- rowSums(sib.mendel.prob) 
    
    rm(sib.mendel.prob, mendel.sib)
    
    fam.array  <- array(0, dim=c(n.fams, 3, max.alleles, max.alleles))
    for(k in c(1, 2, 4)) {
      #for every family member, set likelihoods in fam.array in order of allele's prevalence in that family
      #for each family member, index into genotype array for each possible Mendelian genotype they can have
      target.adj <- if(k == 1) 3 else if (k == 4) 1 else k
      source.ind <- cbind(rep.int(by.fam[, k], rep.int(n.genotypes, n.fams)), 
                          t(matrix(t(sib.trio.genos[, , 1]), nrow=1)),
                          t(matrix(t(sib.trio.genos[, , 2]), nrow=1)))
      target.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.genotypes, n.fams)),
                          rep.int(target.adj, n.fams * n.genotypes), 
                          rep.int(possible.genos[, 1], n.fams),
                          rep.int(possible.genos[, 2], n.fams))
      fam.array[target.ind] <- geno.likelihood[source.ind]
    }
    
    # switch mom and sibling, get mendel score
    switch.mom.for.sib <- array(0, dim=c(n.fams, 3, n.ind))
    for(k in seq_len(3)) {
      mendel.col    <- c((2 * k) - 1, (2 * k))
      mendel.t.ind  <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)),
                             rep.int(k, n.fams * n.ind), 
                             rep.int(seq_len(n.ind), n.fams))
      fam.array.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)),
                             rep.int(k, n.fams * n.ind), 
                             rep.int(mendel.ind[, mendel.col[1]], n.fams),
                             rep.int(mendel.ind[, mendel.col[2]], n.fams))
      switch.mom.for.sib[mendel.t.ind] <- fam.array[fam.array.ind]
      
    }
    mom.switched.prob<-switch.mom.for.sib[, 1, ] * switch.mom.for.sib[, 2, ] * switch.mom.for.sib[, 3, ]
    
    # switch dad and sibling, get mendel score
    for(k in c(1, 2, 4)) {
      #for every family member, set likelihoods in fam.array in order of allele's prevalence in that family
      #for each family member, index into genotype array for each possible Mendelian genotype they can have
      target.adj <- if(k == 2) 3 else if (k == 4) 2 else k
      source.ind <- cbind(rep.int(by.fam[, k], rep.int(n.genotypes, n.fams)), 
                          t(matrix(t(sib.trio.genos[, , 1]), nrow=1)),
                          t(matrix(t(sib.trio.genos[, , 2]), nrow=1)))
      target.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.genotypes, n.fams)),
                          rep.int(target.adj, n.fams * n.genotypes), 
                          rep.int(possible.genos[, 1], n.fams),
                          rep.int(possible.genos[, 2], n.fams))
      fam.array[target.ind] <- geno.likelihood[source.ind]
    }
    
    switch.dad.for.sib <- array(0, dim=c(n.fams, 3, n.ind))
    for(k in seq_len(3)) {
      mendel.col    <- c((2 * k) - 1, (2 * k))
      mendel.t.ind  <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)),
                             rep.int(k, n.fams * n.ind), 
                             rep.int(seq_len(n.ind), n.fams))
      fam.array.ind <- cbind(rep.int(seq_len(n.fams), rep.int(n.ind, n.fams)),
                             rep.int(k, n.fams * n.ind), 
                             rep.int(mendel.ind[, mendel.col[1]], n.fams),
                             rep.int(mendel.ind[, mendel.col[2]], n.fams))
      switch.dad.for.sib[mendel.t.ind] <- fam.array[fam.array.ind]
      
    }
    dad.switched.prob<-switch.dad.for.sib[, 1, ] * switch.dad.for.sib[, 2, ] * switch.dad.for.sib[, 3, ]
    
    fam.genotypes[(i - start.locus + 1), , 14] <- apply(cbind(rowSums(dad.switched.prob),rowSums(mom.switched.prob)),1,max)
  }
  cat("Completed calculating Mendel scores in",(proc.time() - ptm)[3],"\n")  
  return(fam.genotypes)
}
