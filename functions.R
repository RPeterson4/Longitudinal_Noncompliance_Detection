# File of objects and functions needed to implement method

## Section 1: Compound Symmetry

### Load required packages gtools, mvtnorm, and pROC
lapply(c("gtools","mvtnorm","matrixStats","pROC"),require,character.only=TRUE) 

### Set the true parameter values
alpha0_t = 4; alpha1_t = -2.1; tau_t = 0.7; sigma_t = 0.3
beta0_t = 0; gamma_t = 1

### Define objects needed for data generation, parameter estimation, and AUC
z = rep(0,t); zk1 = rep(0,t-1)
ps = permutations(2,t,v=0:1,repeats.allowed=TRUE)
psk1 = permutations(2,t-1,v=0:1,repeats.allowed=TRUE)
cam = ps[rep(1:2^t,n),]
psr = combinations(2,t,v=0:1,repeats.allowed=TRUE)
psr1 = combinations(2,t-1,v=0:1,repeats.allowed=TRUE)
lu_mat = cbind(apply(ps,2,function(x){ifelse(x==1,-Inf,0)}),
	apply(ps,2,function(x){ifelse(x==1,0,Inf)}))
lu_matr = cbind(apply(psr,2,function(x){ifelse(x==1,-Inf,0)}),
	apply(psr,2,function(x){ifelse(x==1,0,Inf)}))
lu_matr1 = cbind(apply(psr1,2,function(x){ifelse(x==1,-Inf,0)}),
	apply(psr1,2,function(x){ifelse(x==1,0,Inf)}))
rsums = rowSums(ps) + 1
rsums1 = rowSums(psk1) + 1

### Ancillary function needed to estimate marginal compliance probabilities
pfun = function(vec,beta0,s_mat,time){
	return(max(.Machine$double.xmin,pmvnorm(vec[1:time],vec[(time+1):(2*time)],
	-beta0,sigma=s_mat,algorithm=Miwa()))[1])
}

### Compliance log likelihood function, defaulted to Compound Symmetry
pnlog_lik = function(x,w,rsms=rsums){
	return(-sum(w*log(apply(lu_matr,1,pfun,beta0=x[1],
	s_mat=matrix(c(1+exp(x[2]),rep(c(rep(exp(x[2]),t),1+exp(x[2])),t-1)),
	nrow=t),time=t)[rsms])))
}

### Biomarker log likelihood function, defaulted to Compound Symmetry
bnlog_lik = function(x,bset,w){
	return(-sum(w*dmvnorm(bset-x[1]-x[2]*cam,z, 
	matrix(c(exp(x[3])+exp(x[4]),rep(c(rep(exp(x[3]),t),exp(x[3])+exp(x[4])),t-1)),
	nrow=t),TRUE)))
}

### Function that arranges all the integrals needed to integrate out the future biomarker
omf = function(subm,a0,a1,ta,sa){
	return(Vectorize(function(bk){
	m_f = subm[,t]*dmvnorm(cbind(subm[,1:(t-1)],bk)-a0-a1*ps,z,
			matrix(c(ta+sa,rep(c(rep(ta,t),ta+sa),t-1)),nrow=t))
	ball = sum(m_f)
	ltp = sum((m_f/ball)[seq(2,2^t,2)])
	ltp[is.nan(ltp)] = 0
	return(ltp*(ball/subm[1,t+1]))
}))
}

### Integration function
int = function(x){integrate(x,-Inf,Inf)[1]$value}

## Section 2: AR(1)

### Set the AR(1) covariance structure for compliance
rho = 0.8; H <- abs(outer(1:t, 1:t, "-"))
cl_sq = 1 + gamma_t;
cv <- cl_sq * rho^H

### Set the AR(1) covariance structure for the biomarker
bl_sq = tau_t + sigma_t;
cb <- bl_sq * rho^H

### Calculate the marginal compliance probabilities under AR(1)
probs_AR1 = apply(lu_mat,1,pfun,beta0=beta0_t,s_mat=cv,time=t)

### Biomarker data generation function under AR(1)
bgen = function(row){rmvnorm(1,mean=alpha0_t+alpha1_t*row,sigma=cb)}

## Section 3: Compound Symmetry for Left-Censored Biomarker Data

### Set the limit of quantification (LOQ)
loq = 1.2

### Helper functions for left-censored data
f = function(x,minimum=loq){return(which(x==minimum))}
f2 = function(x){return(length(f(x)))}

### Helper function for calculating multivariate normal CDF 
### for each row of a matrix for left-censored biomarker
pfun2 = function(vec,s_mat,time){
	set.seed(1)
	return(pmvnorm(vec[1:time],vec[(time+1):(2*time)],
	vec[(2*time+1):(3*time)],sigma=s_mat)[1])
}

### Helper function for calculating multivariate normal PDF for biomarker
bfun = function(vec,s_mat){
	return(dmvnorm(vec,sigma=s_mat))
}

### Biomarker log likelihood function with Compound Symmetry and left-censoring
bno = function(x,bset,lg=TRUE,ft=log,tot=sum,neg=-1,wv=1,time,kset){
	a0 = x[1]; a1 = x[2]; tau = exp(x[3]); sig = exp(x[4])
	sm = matrix(c(tau+sig,rep(c(rep(tau,time),tau+sig),time-1)),nrow=time)
	subao = bset[which(rowSums(inds_mat)==time),]
	cam = ps[rep(1:nrow(ps),nrow(subao)/2^time),]
	lik = dmvnorm(subao-a0-a1*cam,sigma=sm,log=lg)
	for(k in kset){
		subm = bset[which(rowSums(inds_mat)==time-k),]
		bs = cbind(subm[,(k+1):time])
		nc = nrow(subm)
		psk = permutations(2,time-k,v=0:1,repeats.allowed=TRUE)
		camk = cbind(psk[rep(1:nrow(psk),2^k*nrow(subm)/2^time),])
		vals = bs-a0-a1*camk
		tt = matrix(c(tau+sig,rep(c(rep(tau,time-k),tau+sig),time-k-1)),nrow=time-k)
		fs = apply(vals,1,bfun,s_mat=rbind(tt))
	
		oo = matrix(c(tau+sig,rep(c(rep(tau,k),tau+sig),k-1)),nrow=k)
		ot = matrix(tau,nrow=k,ncol=time-k)
		tti = solve(tt)
		dm = ot%*%tti[1:nrow(tti),1:ncol(tti)]%*%t(vals[1:nrow(vals),1:ncol(vals)])
		sc = oo-ot%*%tti%*%t(ot)
		psk2 = cbind(ps[,1:k])
		camk2 = cbind(psk2[rep(1:nrow(psk2),nrow(subm)/2^time),])
		mu = a0+a1*camk2+t(dm)
		bounds = rep(c(-Inf,loq),each=k)
		a_m = cbind(matrix(rep(bounds,nc),nrow=nc,byrow=TRUE),mu)
		fc = apply(a_m,1,pfun2,s_mat=sc,time=k)
		lik = c(lik,ft(fs*fc))
	}
	ai = rbind(rep(c(-Inf,loq),each=time))
	a_a = cbind(ai[rep(1,2^time),],a0+a1*ps)
	ps = apply(a_a,1,pfun2,s_mat=sm,time=time)
	lik = c(lik,rep(ft(ps),nrow(bset[which(rowSums(inds_mat)==0),])/2^time))
	return(neg*tot(wv*lik))
}