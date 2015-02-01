#Get a diagonal matrix with the top k eigenvalues from SVD
svd.diag<-function(i.svd,k=3)
{
  d<-i.svd$d
  if(k < length(i.svd$d))
  {
    d[(k+1):length(d)]<-0
  }
  return(diag(d))
}

#Graph eigenvalues of SVD vs random data and identify how many significant eigenvalues there are,
#optionally set a threshold from random that an eigenvalue must exceed to be considered significant
graph.eigenvalue.vs.random.pick.sig<-function(svd,original.data,n.sim=5,print.dir=NULL,file.name=NULL,svd.type=NULL,new.display=TRUE,min.eigen=1,max.eigen=30)
{
  require(irlba)
  if(is.character(print.dir))
  {
    fname<-NULL
    if(is.character(file.name))
    {
      fname<-sprintf("%s/%s.png",print.dir,file.name)
    } else
    {
      fname<-sprintf("%s/ObsVsShuffledEigenvalues.png",print.dir,file.name)
    }
    png(fname,width=800,height=800)
    par(cex=2)
  } else if(new.display==TRUE)
  {
    quartz(width=14.5,height=8)
  }
  
  min.dim<-min(ncol(original.data),nrow(original.data))
  eigenvalues.to.consider<-1:max.eigen
  
  plot(svd$d,type="p",col="blue",xlab="Eigenvalue",ylab="Magnitude",
       main=sprintf("Eigenvalues of observed vs randomized %s\n%d simulations",svd.type,n.sim),xlim=c(0,max.eigen),cex=0.5,pch=16)
  legend(x="topright",c("Actual","Randomized"),col=c("blue","red"),pch=c(16,16),pt.cex=0.5)
  eigenvalues<-array(NA,dim=c(n.sim,length(svd$d)))
  for(i in 1:n.sim)
  {
    shuf<-matrix(sample(original.data,length(original.data),replace=TRUE),ncol=ncol(original.data))
    shuf.svd<-irlba(shuf,nu=max.eigen,nv=max.eigen)
    eigenvalues[i,]<-shuf.svd$d
    points(shuf.svd$d,col="red",cex=0.5,pch=16)
  }
  if(is.character(print.dir))
  {
    dev.off()
  }

  num.non.random<-max(which(sapply(1:max.eigen, function(j) sum(svd$d[j] > eigenvalues[,j]) == n.sim)))

  num.non.random<-max(min.eigen,num.non.random)
  return(list(k=num.non.random,eigenvalues=eigenvalues))
}

#simulate Poisson coverage given a matrix of lambdas
sim.poisson.cov<-function(lambdas)
{
  #simulate bi-allelic  coverage
  sim.cov<-array(rpois(length(lambdas),2*lambdas),dim=c(nrow(lambdas),ncol(lambdas)))
  return(sim.cov)
}

get.sim.correlation.gain<-function(sim.cov,sim.svd,max.eigen)
{
  cor.gain<-sapply(1:max.eigen, function(i)
  {
        cor(matrix(get.low.rank.from.svd(sim.svd,i),ncol=1),matrix(sim.cov,ncol=1))
  })
  return(cor.gain)
}

get.low.rank.from.svd<-function(svd,rank)
{
  d<-svd.diag(svd,rank)
  return(svd$u %*% d %*% t(svd$v))
}

pick.eigen.poisson.bootstrap.cm<-function(input.svd,num.components=50,max.eigen=50,cor.threshold=1e-3,print.dir=NULL,file.name=NULL,svd.type=NULL,new.display=TRUE)
{
  require(irlba)
  d<-svd.diag(input.svd,num.components)
  cov.rate<-input.svd$u %*% d %*% t(input.svd$v)
  cov.rate[cov.rate < 0]<-min(cov.rate[cov.rate > 0])
  
  #simulate copy number 2 coverage and get low-rank estimation of expected coverage for every eigenvalue up to max.eigen
  ptm<-proc.time()
  sim.cov<-sim.poisson.cov(cov.rate)
  print("Simulated coverage")
  print(proc.time() - ptm)
  print(max.eigen)
  ptm<-proc.time()
  sim.cov.svd<-irlba(sim.cov,nu=max.eigen,nv=max.eigen)
  print("Got SVD of simulated coverage")
  print(proc.time() - ptm)
  cor.to.obs<-get.sim.correlation.gain(sim.cov,sim.cov.svd,max.eigen)
  
  #measure "gain" -- increase in correlation from adding every extra eigenvalue to low-rank approximation
  gain<-sapply(2:max.eigen,function(i) {
      cor.to.obs[i]-cor.to.obs[i-1]
    })
  
  if(is.character(print.dir))
  {
    fname<-sprintf("%s/EigenvalueGain.png",print.dir)
    png(fname,width=800,height=800)
    par(cex=2)
  } else if(new.display==TRUE)
  {
    quartz()
  }
  plot(gain,main="Average correlation gain by eigenvalue",ylab="Correlation gain",xlab="Eigenvalue",xaxt="n",pch=16,cex=0.75)
  axis(side=1,at=c(1,seq(9,max.eigen,10)),labels=c(2,seq(10,max.eigen,10)))
  if(is.character(print.dir))
  {
    dev.off()
  }
  print("Calculated eigenvalue gain")
  
  #first estimate of how many eigenvalues to use is all eigenvalues with gain above cor.threshold (1e-3 by default)
  new.eigen.guess<-max(which(gain > cor.threshold)) + 1
  print(paste("First guess of significant eigenvalues: ",new.eigen.guess,sep=""))
  old.eigen.guess<-new.eigen.guess + 1
  num.iter<-0
  
  #create low-rank approximation, simulate poisson random variables from it, and ask how many eigenvalues are significant (exceed randomness)
  #until the number of significant eigenvalues converges
  while(new.eigen.guess != old.eigen.guess)
  {
    d<-svd.diag(input.svd,new.eigen.guess)
    cov.rate<-input.svd$u %*% d %*% t(input.svd$v)
    cov.rate[cov.rate < 0]<-min(cov.rate[cov.rate > 0])

    sim.cov<-sim.poisson.cov(cov.rate)
    sim.cov.svd<-irlba(sim.cov/2,nu=50,nv=50)
    
    new.randomization<-graph.eigenvalue.vs.random.pick.sig(sim.cov.svd,sim.cov,max.eigen=50)
    
    old.eigen.guess<-new.eigen.guess
    new.eigen.guess<-new.randomization$k
    num.iter<-num.iter + 1
  }
  print(sprintf("# of significant eigenvalues: %d",new.eigen.guess))
  print(sprintf("Converged in %d iterations",num.iter))
  return(new.eigen.guess)
}

linear.coverage.model<-function(locus.cov,n.sim=5,print.dir=NULL,suffix=NULL,num.components=NULL,verbose=1)
{
  ptm<-proc.time()
  require(irlba)
  if(is.character(print.dir))
  {
    to.file<-TRUE
  } else
  {
    to.file<-FALSE
  }
  fig.dir<-NULL
  if(is.character(print.dir))
  {
    fig.dir<-sprintf("%s/LinearCoverageModel",print.dir)
    if(!is.null(suffix))
    {
      fig.dir<-sprintf("%s%s",fig.dir,suffix)
    }
    conditional.make.dir(fig.dir)
  }
  
  allele.cov<-locus.cov/2
	t.nu <- ifelse(any(dim(allele.cov) < 30), min(dim(allele.cov)), 30)
	t.nv <- t.nu
  allele.cov.svd<-irlba(allele.cov,nu=t.nu,nv=t.nv) #because we'll only be using the top eigenvalues of the matrix anyway, there's no reason to take a full SVD of the matrix
  
  if(!is.integer(num.components))
  {
    sig.eigen<-graph.eigenvalue.vs.random.pick.sig(allele.cov.svd,allele.cov,print.dir=fig.dir,file.name="coverageSVD",svd.type="coverage",n.sim=n.sim)
    num.components<-sig.eigen$k
#     num.components<-pick.eigen.poisson.bootstrap.cm(allele.cov.svd,print.dir=fig.dir,file.name="coverageSVD",svd.type="coverage")
  }
  
  if(verbose >= 1) cat("Using",num.components,"SVD components to model coverage\n")
  k<-num.components
  d<-svd.diag(allele.cov.svd,k)
  exp.allele.cov<-allele.cov.svd$u %*% d %*% t(allele.cov.svd$v)
  exp.allele.cov[exp.allele.cov < 0]<-min(exp.allele.cov[exp.allele.cov > 0])
  rm(allele.cov.svd)
 
  cov.cor<-cor(matrix(exp.allele.cov,ncol=1),matrix(allele.cov,ncol=1))
  rmse <- sqrt(mean((exp.allele.cov - allele.cov)^2))
  if(to.file)
  {
    fname<-sprintf("%s/observedVsEstimatedCoverage.png",fig.dir)
    png(fname,width=1800,height=1800,res=300)
  } else
  {
    quartz()
  }
  random.sample<-sample(length(exp.allele.cov),5e5)
  plot(exp.allele.cov[random.sample],locus.cov[random.sample],pch=16,xlab="Expected coverage per allele",ylab="Observed coverage",
       main=sprintf("Expected vs observed allele coverage\ncor: %.2f; RMSE: %.2f",cov.cor,rmse),col=rgb(0,0,0,0.1),xlim=c(0,500),ylim=c(0,500))
  abline(a=0,b=2,col="red",lwd=2)
  if(to.file)
  {
    dev.off()
  }
  
  #in case there are people with no coverage, only calculate correlation for people who have coverage
  has.cov<-which(colSums(locus.cov) > 0)
  cor<-rep(NA,ncol(exp.allele.cov))
  cor<-sapply(has.cov,function(i) cor(exp.allele.cov[,i],allele.cov[,i]))
  if(to.file)
  {
    fname<-sprintf("%s/individualCorrelations.png",fig.dir)
    png(fname,width=1200,height=1200)
    par(cex=2)
  } else
  {
    quartz()
  } 
  hist(cor,breaks=seq(0,1,0.01),xlab="expected and observed per-allele coverage correlation",border=NA,
       main=sprintf("Expected vs observed allele coverage correlation histogram\n overall cor: %.3f",cov.cor),col=rgb(0,0,0,0.4))
  if(to.file)
  {
    dev.off()
  }
  return(list(exp.cov=exp.allele.cov,cor=cor,num.components=k))
  
}
