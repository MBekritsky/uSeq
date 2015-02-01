#Returns a n.fams x 4 matrix of person indices grouped by family ID
get.fam.mat <- function(person.info) {
  fam.ids      <- unique(person.info$family.id)
  max.fam.size <- max(tabulate(person.info$family.id))
  pop          <- nrow(person.info)
  by.fam       <- array(0,dim=c(length(fam.ids),max.fam.size))
  
  if(!("rel.code" %in% colnames(person.info))) {
    # if there is no column called "rel.code", code families
    # such that mother is first, father second, proband third
    # and sibling fourth
    person.info$rel.code[person.info[,'relation'] == "mother"]  <- 1
    person.info$rel.code[person.info[,'relation'] == "father"]  <- 2
    person.info$rel.code[person.info[,'relation'] == "proband"] <- 3
    person.info$rel.code[person.info[,'relation'] == "sibling"] <- 4
  }
  person.info$rel.code[person.info[,'relation'] == "s2"]      <- 4
  # "s2" is a designation for some siblings in the Eichler and State
  # SSC datasets
  
  for(f in fam.ids) {
    fam.idx     <- which(fam.ids == f)
    fam.idxs    <- which(person.info$family.id == f)
    fam.columns <- person.info$rel.code[fam.idxs]
#     fam.idxs    <- fam.idxs[fam.columns]
    by.fam[fam.idx, fam.columns] <- fam.idxs
  }
  colnames(by.fam) <- c("mother","father","proband","sibling")
  rownames(by.fam) <- fam.ids
  
  return(by.fam)
}

default.fields <- function() {
  return(c("allele.1", "allele.2", "confidence", "allele.fit",
           "noise.fit", "allele.1.cov", "allele.2.cov", "rest.cov",
           "exp.allele.cov", "exp.err", "p.null"))
}

# initialize genotype array
init.genotype.array <- function(n.loci, pop) {
  field.names <- default.fields()
  n.fields  <- length(field.names)
  genotypes <- array(0, dim=c(n.loci, pop, n.fields))
  dimnames(genotypes)<-list(NULL, NULL, field.names)
  return(genotypes)
}

init.genotype.array.slice <- function(pop) {
  field.names <- default.fields()
  n.fields    <- length(field.names)
  slice       <- array(0, dim=c(pop, n.fields))
  dimnames(slice) <- list(NULL, field.names)
  return(slice)
}

# get information about alleles at a locus, their row numbers
# in the allele matrix, and orders them by total coverage per
# allele throughout the population
get.locus.details <- function(locus.index, locus.info) {
  locus.ind  <- which(locus.info$locus.row.ind == locus.index & locus.info$allele.no == 0)
  allele.ind <- locus.ind + seq_len(locus.info$num.lengths[locus.ind])

  ordered.lengths <- allele.ind[order(locus.info$sum.count[allele.ind],decreasing=TRUE)]
  lengths         <- locus.info$allele.no[ordered.lengths]
  n.alleles       <- locus.info$num.lengths[locus.ind] + 1 #total observed alleles and null allele
  
  return(list(locus.ind=locus.ind,allele.ind=allele.ind,
              ordered.lengths=ordered.lengths,lengths=lengths,
              n.alleles=n.alleles))
}

#get allele coverage and total locus coverage for every person
# at a locus
get.coverage <- function(allele.info, locus.details, pop) {
  allele.cov <- matrix(allele.info$alleles[locus.details$ordered.lengths, ], ncol=pop)
  locus.cov <- colSums(allele.cov)
  return(list(allele.cov=allele.cov, locus.cov=locus.cov))
}

#returns a list of matrices used for genotyping (pop size x n x n)
get.obs.exp.mat <- function(pop, n.alleles, coverage.info, locus.coverage.model, allele.biases) {
  n.obs.alleles <- n.alleles - 1 #n.alleles is total number of alleles, including the null allele
  
  allele.cov <- coverage.info$allele.cov
  locus.cov  <- coverage.info$locus.cov
  
  obs.allele.1  <- array(0,  dim=c(pop, n.alleles, n.alleles))
  obs.allele.2  <- array(0,  dim=c(pop, n.alleles, n.alleles))
  obs.total     <- array(0,  dim=c(pop, n.alleles, n.alleles))
  cn.allele.1   <- array(0,  dim=c(pop, n.alleles, n.alleles))
  cn.allele.2   <- array(0,  dim=c(pop, n.alleles, n.alleles))
  expected.cov  <- array(NA, dim=c(pop, n.alleles, n.alleles))
  bias.allele.1 <- array(0,  dim=c(pop, n.alleles, n.alleles))
  bias.allele.2 <- array(0,  dim=c(pop, n.alleles, n.alleles))
  
  for(i in seq_len(n.obs.alleles)) {
    obs.allele.1[, i, i:n.alleles]    <- rep(allele.cov[i, ], n.alleles - i + 1)
    obs.total[, i, i:n.alleles]       <- rep(locus.cov, n.alleles - i + 1)
    bias.allele.1[, i, i:n.alleles]   <- allele.biases[i]
    expected.cov[, i, i:n.alleles]    <- rep(locus.coverage.model, n.alleles - i + 1)
    
    # copy number of homozygous alleles is 2
    # copy number of heterozygous alleles is 1
    cn.allele.1[, i, (i+1):n.alleles] <- 1
    cn.allele.1[, i, i]               <- 2
  }
  
  
  # every locus will have at least one allele, so seq_len will
  # never be called with a negative number
  for(i in seq_len(n.obs.alleles - 1)) {
    obs.allele.2[, 1:i, (i + 1)]  <- rep(allele.cov[(i + 1), ], i)
    bias.allele.2[, 1:i, (i + 1)] <- allele.biases[(i + 1)]
    cn.allele.2[, 1:i, (i + 1)]   <- 1
  }
  
  # set total locus cov for double null case
  obs.total[, n.alleles, n.alleles] <- locus.cov
  # set observed error
  obs.error <- obs.total - obs.allele.1 - obs.allele.2
  
  # when people have expected coverage below 1, set their expected coverage to 1
  expected.cov[expected.cov < 1]    <- 1
  expected.cov[is.na(expected.cov)] <- 0
  
  return(list(obs.allele.1=obs.allele.1,obs.allele.2=obs.allele.2,obs.total=obs.total,obs.error=obs.error,
              cn.allele.1=cn.allele.1,cn.allele.2=cn.allele.2,bias.allele.1=bias.allele.1,bias.allele.2=bias.allele.2,
              expected.cov=expected.cov))
}

# estimate locus noise paramater and bias parameters for each allele
# at the locus, then assign genotypes for a particular locus, including
# a genotype confidence, expected coverage and noise, exact p-values for
# alelle fit (exact two-sided Poisson p-value) and noise fit (exact
# binomial p-value for p(observed) > p(expected)), and probability that
# this person has any null genotype (p.null)
assign.genotypes <- function(ag.mats, pop, locus.info, locus.details, coverage.model, locus.num,
                             coverage.info, initial.noise.estimate=0.001, initial.bias.estimates=1,
                             em.threshold=10, max.em.iter=25) {
  locus.cov       <- coverage.info$locus.cov
  allele.cov      <- coverage.info$allele.cov
  ordered.lengths <- locus.details$ordered.lengths
  n.alleles       <- locus.details$n.alleles
  genotypes       <- init.genotype.array.slice(pop)
    
  # Assign bias terms for each observed allele (i.e. not the null allele)
  lb <- length(initial.bias.estimates)
  if( lb == 1L) {
    init.bias.estimates <- rep(initial.bias.estimates, n.alleles - 1)
  } else if(lb == (locus.details$n.alleles - 1))
  {
    init.bias.estimates<-initial.bias.estimates
  } else {
    stop("Initial bias estimate is neither a scalar nor a vector with one bias parameter per allele")
  }
  
  weights <- initialize.weights(n.alleles)
  p.geno  <- get.geno.probs(ag.mats, initial.noise.estimate, weights)
  
  #get the initial model likelihood and set the initial maximum likelihood estimators
  #for bias and error rate
  err.rate.mle       <- initial.noise.estimate
  bias.estimates.mle <- init.bias.estimates
  ll.geno            <- get.model.likelihood(ag.mats, initial.noise.estimate, weights)
  
  old.ll.geno <- ll.geno
  new.ll.geno <- ll.geno + (2 * em.threshold)
  n.iter      <- 0
  
  while(abs(new.ll.geno - old.ll.geno) > em.threshold & n.iter < max.em.iter & err.rate.mle < 0.4)
  {
    n.iter <- n.iter + 1
    #expectation step--get the likelihood for each genotype for each member in the population
    p.geno <- get.geno.probs(ag.mats, err.rate.mle, weights)
    
    #maximization step--get the maximum likelihood estimators for each genotype weight and bias, and the locus error rate
    weights            <- get.genotype.weights(p.geno)    
    err.rate.mle       <- mle.error.prob(p.geno, ag.mats, weights)
    bias.estimates.mle <- mle.allele.biases(p.geno, ag.mats, weights)
    ag.mats            <- update.bias.mats(ag.mats, bias.estimates.mle)

    #calculate the new model's likelihood
    old.ll.geno <- new.ll.geno
    new.ll.geno <- get.model.likelihood(ag.mats, err.rate.mle, weights)
  }
  
  # get the most likely genotype in each person at this locus
  allele.ind  <- t(apply(p.geno, 1, function(x) which(x == max(x),arr.ind=TRUE)[1,]))
  allele.ind  <- cbind(allele.ind,seq_len(pop)) 
  # create a list of all possible alleles
  all.alleles <- c(locus.details$lengths, -1)
  #assign each person's alleles
  alleles     <- t(apply(allele.ind, 1, function(x) c(all.alleles[x[1]], all.alleles[x[2]])))
  
  good.for.nf <- rowSums(is.finite(allele.ind[,1:2])) == 2 & locus.cov > 0
  #1-sided exact binomial test giving probability that observed error rate is greater
  # than expected error rate
  likeliest.error<-ag.mats$obs.error[cbind(which(good.for.nf), allele.ind[good.for.nf, 1], allele.ind[good.for.nf, 2])]
  genotypes[good.for.nf, 5] <- vectorized.binom.p.val(x=likeliest.error, n=locus.cov[good.for.nf], p=err.rate.mle, 
                                                      alternative="greater")
  
  # most likely bi-allelic genotype
  genotypes[, 1:2] <- alleles
  
  # confidence for likeliest genotype among all binomial genotypes
  genotypes[, 3]   <- apply(p.geno, 1, max)
  
  # expected per-allele coverage
  genotypes[, 9]   <- coverage.model$exp.cov[locus.num, ]
  
  # expected number of error reads
  genotypes[, 10]  <- err.rate.mle * locus.cov
  
  #total probability that an allele is missing (marginal probability of null allele)
  genotypes[, 11]  <- apply(p.geno, 1, function(x) sum(x[, n.alleles])) 
  
  hom         <- alleles[, 1] == alleles[, 2] & alleles[, 1] > 0
  het         <- alleles[, 1] != alleles[, 2] & rowSums(alleles > 0) == 2
  null        <- alleles[, 1] > 0 & alleles[, 2] == -1
  double.null <- alleles[, 1] == -1 & alleles[, 2] == -1

  #Assign coverage for each allele and noise
  if(any(hom) | any(het) | any(null)) {
    has.one.allele    <- which(hom | het | null)
    allele.one.ind    <- allele.ind[has.one.allele, 1]
    genotypes[has.one.allele, 6] <- allele.cov[cbind(allele.one.ind, has.one.allele)]
  }
  if(any(het)) {
    genotypes[het, 7] <- allele.cov[cbind(allele.ind[which(het),2],which(het))]
  }
  genotypes[, 8] <- locus.cov - rowSums(genotypes[, 6:7])
  
  # determine the allele fit by a two-tailed exact Poisson test
  exp.allele.1 <- ag.mats$expected.cov[allele.ind[, c(3, 1, 2)]] * ag.mats$cn.allele.1[allele.ind[, c(3, 1, 2)]] * 
                  ag.mats$bias.allele.1[allele.ind[, c(3, 1, 2)]]
  exp.allele.2 <- ag.mats$expected.cov[allele.ind[, c(3, 1, 2)]] * ag.mats$cn.allele.2[allele.ind[, c(3, 1, 2)]] * 
                  ag.mats$bias.allele.2[allele.ind[, c(3, 1, 2)]]

  genotypes[, 4] <- vectorized.poisson.p.val(x=ag.mats$obs.allele.1[allele.ind[, c(3, 1, 2)]],r=exp.allele.1, alternative="two.sided") * 
                   vectorized.poisson.p.val(x=ag.mats$obs.allele.2[allele.ind[, c(3, 1, 2)]],r=exp.allele.2, alternative="two.sided")
  
  return(list(genotypes=genotypes, p.geno=p.geno, biases=data.frame(length=locus.details$lengths, bias=bias.estimates.mle),
              err.rate=err.rate.mle, n.iter=n.iter, weights=weights))
}

get.parental.genotype.frequency<-function(l.geno,mom.and.pop)
{
  return(apply(l.geno[mom.and.pop,,],2:3,sum))
}

call.genotypes<-function(allele.info,locus.info,coverage.model,pop.size,max.locus=NULL,start.locus=1,debug=FALSE,verbose=1,max.em.iter=25)
{
  if(is.null(max.locus))
  {
    max.locus<-allele.info$num.loci
  }

  n.loci<-max.locus - start.locus + 1
  
  called.genotypes<-init.genotype.array(n.loci,pop.size)
  allele.biases<-vector("list",n.loci)
  locus.error.rates<-rep(NA,n.loci)
  n.em.iter<-rep(NA,n.loci)
  weights<-vector("list",n.loci)
  
  ptm<-proc.time()
  for(i in start.locus:max.locus)
  {
    if(i %% 1000 == 0 & verbose >= 1)
    {
      cat("Processed locus ",i,". Total elapsed time: ",(proc.time() - ptm)[3]," seconds\n",sep="")
    }
    if(debug)
    {
      print(i)
    }
    #get indices in allele and locus matrices for locus and allele coverage, as well as number of alleles
    locus.details<-get.locus.details(i,locus.info)
    
    #get coverage information for locus and for each allele
    coverage.info<-get.coverage(allele.info,locus.details,pop.size)

    #get matrices for genotyper
    mats<-get.obs.exp.mat(pop.size,locus.details$n.alleles,coverage.info,coverage.model$exp.cov[i,],rep(1,locus.details$n.alleles - 1))

    #get genotypes, genotype confidences, and coverage info
    genotype.output<-assign.genotypes(mats,pop.size,locus.info,locus.details,coverage.model,i,coverage.info,max.em.iter=max.em.iter)
    called.genotypes[(i - start.locus + 1),,]<-genotype.output$genotypes
    n.em.iter[(i - start.locus + 1)]<-genotype.output$n.iter
    locus.error.rates[(i - start.locus + 1)]<-genotype.output$err.rate
    allele.biases[[(i - start.locus + 1)]]<-genotype.output$biases
    weights[[i - start.locus + 1]]<-genotype.output$weights
  }
  cat("Called genotypes in",(proc.time() - ptm)[3],"seconds\n")
  return(list(genotypes=called.genotypes,
              n.em.iter=n.em.iter,locus.error.rates=locus.error.rates,allele.biases=allele.biases,weights=weights))
}

initialize.weights <- function(num.alleles) {
  weights <- matrix(1, ncol=num.alleles, nrow=num.alleles)
  weights[lower.tri(weights, diag=FALSE)] <- 0
  return(weights)
}

make.weight.array <- function(weights, num.people) {
  # makes a weight array of pop x n.alleles x n.alleles
  # input weights is a square matrix
  num.alleles <- nrow(weights)
  weight.arr <- array(rep(weights, num.people), dim=c(num.alleles, num.alleles, num.people))
  weight.arr <- aperm(weight.arr, c(3, 1, 2))
  return(weight.arr)
}

get.geno.probs <- function(f.mats, error.estimate, weights) {
  # given a locus error estimate, genotype weights, and a list of matrices
  # with observed and expected counts, allele copy number, expected per allele
  # coverage, and per allele bias, call genotypes
  
  num.alleles  <- dim(f.mats$obs.allele.1)[2]
  num.people   <- dim(f.mats$obs.allele.1)[1]
  weight.array <- make.weight.array(weights, num.people)
  probs        <- weight.array * 
                  dpois(x=f.mats$obs.allele.1, lambda=(f.mats$cn.allele.1 * f.mats$expected.cov * f.mats$bias.allele.1)) * # allele 1
                  dpois(x=f.mats$obs.allele.2, lambda=(f.mats$cn.allele.2 * f.mats$expected.cov * f.mats$bias.allele.2)) * # allele 2
                  dbinom(x=f.mats$obs.error, size=f.mats$obs.total, prob=error.estimate)                                   # noise
  
  for(i in seq_len(num.alleles - 1))
  {
    probs[,(i + 1),seq_len(i)] <- 0
  }
  
  nan.sums<-apply(probs,1,function(x) sum(is.na(x)))
  zero.sums<-apply(probs,1,function(x) sum(x == 0))

  # find people who either have no likely genotypes or have all NaN genotypes
  # and set them to be double null
  no.probs<-which(nan.sums == num.alleles^2 | zero.sums == num.alleles^2)
  probs[no.probs,num.alleles,num.alleles]<-1

  # any remaining probabilities that are NaN are set to 0
  still.nan<-which(is.nan(probs))
  probs[still.nan]<-0
  
  probs<-probs/apply(probs,1,sum)
  
  return(probs)
}

get.model.likelihood <- function(f.mats, error.estimate, weights) {
  num.people   <- dim(f.mats$obs.allele.1)[1]
  weight.array <- make.weight.array(weights, num.people)
  likelihood   <- sum(weight.array * 
                      (dpois(x=f.mats$obs.allele.1, lambda=(f.mats$cn.allele.1 * f.mats$expected.cov * f.mats$bias.allele.1), log=TRUE) +
                       dpois(x=f.mats$obs.allele.2, lambda=(f.mats$cn.allele.2 * f.mats$expected.cov * f.mats$bias.allele.2), log=TRUE) +
                       dbinom(x=f.mats$obs.error, size=f.mats$obs.total, prob=error.estimate, log=TRUE)))
  return(likelihood)
}

get.genotype.weights <- function(membership.probs) {
  n.people <- dim(membership.probs)[1]
  weights  <- apply(membership.probs, 2:3, sum) / n.people
  return(weights)
}

mle.error.prob <- function(membership.probs, f.mats, weights) {
  prob <- sum(weights * apply(membership.probs * f.mats$obs.error, 2:3, sum)) / 
          sum(weights * apply(membership.probs * f.mats$obs.total, 2:3, sum))
      
  if(prob == 0) {
    # if the error probability is zero, the error probability become 1/(total coverage at locus)
    prob <- 1 / sum(f.mats$obs.total[, 1, 1])
  }
  return(prob)
}

mle.allele.biases <- function(membership.probs, f.mats, weights) {
  #allele bias terms are only estimated for observed alleles, not the null allele
  num.alleles <- dim(membership.probs)[2] - 1
  
  #total observed coverage of each allele when it is allele 1
  numer.a1    <- apply(membership.probs * f.mats$obs.allele.1, 2:3, sum)
  #total expected coverage of each allele when it is allele 1
  denom.a1    <- apply(membership.probs * f.mats$cn.allele.1 * f.mats$expected.cov, 2:3, sum) 
  #total observed coverage of each allele when it is allele 2
  numer.a2    <- apply(membership.probs * f.mats$obs.allele.2, 2:3, sum)
  #total expected coverage of each allele when it is allele 2
  denom.a2    <- apply(membership.probs * f.mats$cn.allele.2 * f.mats$expected.cov, 2:3, sum)
  
  new.bias.estimates    <- rep.int(NA, num.alleles)
  new.bias.estimates[1] <- log(sum(weights[1, ] * numer.a1[1, ])) - log(sum(weights[1, ] * denom.a1[1, ]))
  
  for(i in seq_len(num.alleles - 1)) {
    # positions where this allele occurs as allele 1
    a1.inds <- cbind((i + 1),((i+1):(num.alleles + 1)))
    # positions where this allele occurs as allele 2
    a2.inds <- cbind(seq_len(i),(i + 1))
    
    a1.allele.numer <- weights[a1.inds] * numer.a1[a1.inds]
    a2.allele.numer <- weights[a2.inds] * numer.a2[a2.inds]
    a1.allele.denom <- weights[a1.inds] * denom.a1[a1.inds]
    a2.allele.denom <- weights[a2.inds] * denom.a2[a2.inds]
    
    new.bias.estimates[i + 1] <- log(sum(a1.allele.numer, a2.allele.numer)) - log(sum(a1.allele.denom, a2.allele.denom))
  }
  
  new.bias.estimates[(new.bias.estimates <= log(1e-07) |
                      is.na(new.bias.estimates))] <- log(1e-07) #set the bias for an allele with bias 0 to be exceedingly small
  return(exp(new.bias.estimates))
}

update.bias.mats <- function(f.mats, new.allele.biases) {
  num.people      <- dim(f.mats$obs.1)[1]
  num.obs.alleles <- length(new.allele.biases)
  num.alleles     <- num.obs.alleles + 1 #all alleles, includes null allele
  
  f.mats$bias.allele.1[,1 , 1:num.alleles] <- new.allele.biases[1]
  if(num.alleles > 2) {
    for(i in seq_len(num.obs.alleles - 1)) {
      f.mats$bias.allele.1[, (i+1), seq_len(num.alleles - i) + i] <- new.allele.biases[i]
      f.mats$bias.allele.2[, seq_len(i), (i + 1)] <- new.allele.biases[i + 1]
    }
  }
  return(f.mats)
}