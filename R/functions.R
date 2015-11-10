



#######################################################
### Functions

.rateClonemates1G=function(counts1){
  return(sum(counts1/sum(counts1)*pmax(0,counts1-1)/(sum(counts1)-1)))
}

.rateClonemates2G=function(counts1,counts2){
  return(sum(counts1/sum(counts1)*counts2/sum(counts2)))
}
	
.count2profil=function(R){
  out=NULL
  for(i in 1:length(R)){ out=c(out,rep(i,R[i])) }
  return(out)
}
	
.profil2count=function(Z){
  out=NULL
  for(i in 1:length(unique(Z))){ out=c(out,sum(Z==unique(Z)[i])) }
  return(out)
}

.cycle=function(R0,s){
  Z0=.count2profil(R0)
  N=length(Z0)
  S=rbinom(1,N,s)
  C=N-S
  Z1C=sample(Z0,C,replace=T)
  R1=NULL
  for(i in 1:length(R0)){ R1=c(R1,sum(Z1C==i)) }
  return(cbind(c(R0,rep(0,S)),c(R1,rep(1,S))))
}
	
.cycleMult=function(R0,s,n1,n2){
  Ri=R0
  if(n1>0){
    for(i in 1:n1){
      Ri=.cycle(Ri,0)[,2]
    }
  }	
  Rn1plus1=.cycle(Ri,s)[,2]
  R0=c(R0,rep(0,length(Rn1plus1)-length(R0)))
  Ri=Rn1plus1
  if(n2>0){
    for(i in 1:n2){
      Ri=.cycle(Ri,0)[,2]
    }
  }
  return(cbind(R0,Ri))
}

.logit=function(p){ return(log(p/(1-p))) }
.logitinv=function(x){ return(exp(x)/(1+exp(x))) }

## Contrast
.critMult=function(f0raw,f1raw,p0,p1,g01,n1,n2,logN,logits){
  N=round(exp(logN))
  s=.logitinv(logits)
  f0=N*f0raw/(N-1)-p0/(N-1)
  Fobs=N*f1raw/(N-1)-p1/(N-1)
  Gobs=g01
  a=1/N
  b=1-1/N
  G=(1-s)*(a+b*f0)
  fprime=(1-s)^2*(a*sum(b^(0:n1))+b^(n1+1)*f0)
  if(n2>0){
    fprime=a*sum(b^(0:(n2-1)))+b^n2*fprime
  }
  return(abs(fprime-Fobs)+abs(G-Gobs))
}


## Optimization function
.optim.Ne.s=function(f0raw,f1raw,p0,p1,g01,n1,n2,initial,maxit){
  h=function(par){ 
    .critMult(f0raw,f1raw,p0,p1,g01,n1=n1,n2=n2,par[1],par[2])
  }
  oo=NULL
  oo$value=Inf
  for(i in 1:nrow(initial)){
    ooNew=optim(c(log(initial[i,1]),.logit(initial[i,2])),h,
      control=list(maxit=maxit))
    if(ooNew$value<oo$value){ oo=ooNew }
  }
  oo$partransf=c(round(exp(oo$par[1])),.logitinv(oo$par[2]))
  return(oo)
}

## Simulation
simul=function(start.counts,param,nb1ClonalCycles,nb2ClonalCycles,nb1Sample,nb2Sample,
	nbPrelCycles){
	Ne=param[1]
	s=param[2]
  ## initial conditions
  if(sum(start.counts)<Ne){
    start.counts=start.counts[start.counts>0]
    start.countsstar=c(sum(start.counts==1),start.counts[start.counts>1])	
    R0=as.numeric(rmultinom(1,Ne,start.countsstar))
    R0=c(R0[-1],rep(1,R0[1]))
  }
  if(sum(start.counts)==Ne){
    R0=start.counts
  }
  if(sum(start.counts)>Ne){
    R0=as.numeric(rmultinom(1,Ne,start.counts))
  }
  ## R1 simulation
  temp=.cycleMult(R0,s,nb1ClonalCycles,nb2ClonalCycles)
  if(nbPrelCycles>0){
    for(i in 1:nbPrelCycles){
      temp=.cycleMult(temp[temp[,2]>0,2],s,nb1ClonalCycles,nb2ClonalCycles)
    }
  }
  R0=temp[,1]
  R1=temp[,2]
  ## sampling 
  z0=sample(.count2profil(R0),nb1Sample,replace=TRUE)
  z1=sample(.count2profil(R1),nb2Sample,replace=TRUE)
  r0=NULL
  r1=NULL
  for(i in sort(unique(c(z0,z1)))){
    r0=c(r0,sum(z0==i))
    r1=c(r1,sum(z1==i))
  }
  ## statistics
  f0raw=.rateClonemates1G(r0)
  f1raw=.rateClonemates1G(r1)
  p0=sum(r0/nb1Sample)
  p1=sum(r1/nb2Sample)
  g01=.rateClonemates2G(r0,r1)
  
  #f0raw=(.rateClonemates1G(R0)+1/(Ne-1))*(Ne-1)/Ne
  #f1raw=(.rateClonemates1G(R1)+1/(Ne-1))*(Ne-1)/Ne
  #p0=sum(r0/nb1Sample)
  #p1=sum(r1/nb2Sample)
  #g01=.rateClonemates2G(R0,R1)
  population=temp[apply(temp,1,function(x) sum(x)>0),]
  colnames(population)=c("t0","t1")
  samples=cbind(r0,r1)
  colnames(samples)=c("t0","t1")
  return(list(population=population,samples=samples,
  		statistics=c(f0raw,f1raw,p0,p1,g01)))
}
	
## parametric bootstrap	
.parboot.parfix=function(r00,nech1,nech2,param,n1,n2,nbB,nbPrelCycles,param0,maxit,trace){
  temp0=NULL
  for(k in 1:nbB){
    if(k/100==round(k/100) & trace){ print(paste("Iteration in parametric bootstrap:",k),quote=F) }
    simk=simul(r00,param,n1,n2,nech1,nech2,nbPrelCycles)
    f0raw=simk$statistics[1]
    f1raw=simk$statistics[2]
    p0=simk$statistics[3]
    p1=simk$statistics[4]
    g01=simk$statistics[5]
    ## parameter estimation    
    oo=.optim.Ne.s(f0raw,f1raw,p0,p1,g01,n1,n2,param0,maxit)
    temp0=rbind(temp0,c(param,oo$partransf))
  }
  return(temp0)
}

.estim.CI=function(param,paramB,confidenceLevel){
  ## without bias correction
  endpoints=c((1-confidenceLevel)/2,(1+confidenceLevel)/2)
  est1=cbind(param,t(apply(paramB[,3:4],2,quantile,endpoints)))
  ## with bias correction
  bias=c(mean(param[1]-paramB[,3]),mean(param[2]-paramB[,4]))
  sdhat=apply(paramB[,3:4],2,sd)
  param2=param+bias
  est2=cbind(param2,rbind(qnorm(endpoints,param[1]+bias[1],sdhat[1]),
    qnorm(endpoints,param[2]+bias[2],sdhat[2])))
  est=cbind(est1,est2)
  est[1,]=round(est[1,])
  est[2,]=round(est[2,],digits=4)
  rownames(est)=c("Ne","s")
  colnames(est)=c("Estim","IC1","IC2","EstimBiasCorr","IC1BiasCorr",
            "IC2BiasCorr")
  return(est)
}


.estim.CItransf=function(param,paramB,confidenceLevel){
  ## without bias correction
  endpoints=c((1-confidenceLevel)/2,(1+confidenceLevel)/2)
  est1=cbind(param,t(apply(paramB[,3:4],2,quantile,endpoints)))
  ## with bias correction
  ##bias=c(mean(param[1]-paramB[,3]),mean(param[2]-paramB[,4]))
  ##sdhat=apply(paramB[,3:4],2,sd)
  ##param2=param+bias
  ##est2=cbind(param2,rbind(qnorm(endpoints,param[1]+bias[1],sdhat[1]),
  ##  qnorm(endpoints,param[2]+bias[2],sdhat[2])))
  ## with bias correction
  ## after log transf for Ne and logit transformation for s
  lparam=c(log(param[1]),log(param[2]/(1-param[2])))
  lparamB=cbind(log(paramB[,3]),log(paramB[,4]/(1-paramB[,4])))
  bias=c(mean(lparam[1]-lparamB[,1]),mean(lparam[2]-lparamB[,2]))
  sdhat=apply(lparamB,2,sd)
  lparam3=lparam+bias
  browser()
  est3=cbind(lparam3,rbind(qnorm(endpoints,lparam[1]+bias[1],sdhat[1]),
    qnorm(endpoints,lparam[2]+bias[2],sdhat[2])))
  est3=rbind(exp(est3[1,]),exp(est3[2,])/(1+exp(est3[2,])))
  est=cbind(est1,est3)
  est[1,]=round(est[1,])
  est[2,]=round(est[2,],digits=4)
  rownames(est)=c("Ne","s")
  colnames(est)=c("Estim","IC1","IC2","EstimBiasCorr",
  	"IC1BiasCorr","IC2BiasCorr")
  est
}


.coverage=function(param,nbCoverage,confidenceLevel,nbBootstrap,
  r00,nech0,nech1,n1,n2,param0,maxit,trace){
  cover=cbind(matrix(0,2,2),NA)
  confint=NULL
  colnames(cover)=c("Cover","CoverBiasCorr","nbIteration")
  rownames(cover)=c("Ne","s")
  for(k in 1:nbCoverage){
    if(k/10==round(k/10) & trace){ 
      print(paste("Iteration in coverage assessment:",k),quote=F)
      cover[1,3]=k^2
    }
    oostar=.parboot.parfix(r00,nech0,nech1,param,
      n1=n1,n2=n2,1,0,param0,maxit,trace=FALSE)
    PBstar=.parboot.parfix(r00,nech0,nech1,oostar[3:4],
      n1=n1,n2=n2,nbBootstrap,0,param0,maxit,trace=FALSE)
    est=.estim.CI(oostar[3:4],PBstar,confidenceLevel)
    cover[1,1]=cover[1,1]+((param[1]>=est[1,2]) & (param[1]<=est[1,3]))
    cover[2,1]=cover[2,1]+((param[2]>=est[2,2]) & (param[2]<=est[2,3]))
    cover[1,2]=cover[1,2]+((param[1]>=est[1,5]) & (param[1]<=est[1,6]))
    cover[2,2]=cover[2,2]+((param[2]>=est[2,5]) & (param[2]<=est[2,6]))
    confint=rbind(confint,c(est[1,2:3],est[2,2:3],est[1,5:6],est[2,5:6]))
  }
  cover[,1:2]=cover[,1:2]/nbCoverage
  cover[1,3]=nbCoverage
  return(cover)
}

###################################################

ncase=function(counts,nb1ClonalCycles,nb2ClonalCycles,start.param,maxit=500,
	nbBootstrap=0,conf.level=0.95,nbCoverage=0,trace=TRUE){
	counts1=counts[,1]
	counts2=counts[,2]
	## Computation of clonemate rates	
	f=.rateClonemates1G(counts1)	
	fprime=.rateClonemates1G(counts2)
	p=sum(counts1/sum(counts1))
	pprime=sum(counts2/sum(counts2))
	g=.rateClonemates2G(counts1,counts2)
	## Estimation
	estim1=.optim.Ne.s(f,fprime,p,pprime,g,nb1ClonalCycles,nb2ClonalCycles,start.param,
	  maxit)
	param1=estim1$partransf
	est=param1
	names(est)=c("Ne","s")
	if(nbBootstrap>0){
  		print("PARAMETRIC BOOTSTRAP: BIAS CORRECTION AND CONFIDENCE INTERVALS (may be time consuming)",
        	quote=F)
		estimB=.parboot.parfix(counts1,sum(counts1),sum(counts2),param1,
    		nb1ClonalCycles,nb2ClonalCycles,nbBootstrap,nbPrelCycles=0,start.param,maxit,
    		trace)
  		est=.estim.CI(param1,estimB,conf.level)
    	param2=est[,4]
  		if(nbCoverage>0){
    			print("ASSESSMENT OF THE COVERAGE OF CONFIDENCE INTERVALS (may be time consuming)",quote=F)
    			cover1=.coverage(param1,nbCoverage,conf.level,nbBootstrap,counts1,
      			sum(counts1),sum(counts2),nb1ClonalCycles,nb2ClonalCycles,start.param,
      			maxit,trace)
      		return(list(bootstrap.estim=estimB[,3:4],cover=cover1,
      			optim.output=estim1,estimates=est))
      	} else {
      		return(list(bootstrap.estim=estimB[,3:4],optim.output=estim1,estimates=est))
      	}
	} else {
		return(est)
	}
}


ncase4simul=function(param,nb1ClonalCycles,nb2ClonalCycles,nb1Sample,nb2Sample,
	start.counts,nbPrelCycles,nbSimul,start.param,maxit=500,trace=TRUE){
  	print("SIMULATIONS UNDER FIXED PARAMETERS AND RAW ESTIMATIONS (may be time consuming)",quote=F)
    estimB=.parboot.parfix(start.counts,nb1Sample,nb2Sample,param,
      		nb1ClonalCycles,nb2ClonalCycles,nbB=nbSimul,nbPrelCycles,start.param,maxit,trace)  
    colnames(estimB)=c("Ne","s","Ne.hat","s.hat")
    return(estimB)
}



