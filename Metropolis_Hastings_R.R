
# This is from the discussion of Hierarchical Bayesian Modelswhere 
# y_i~N(mu,1)
# mu~t(0,1,1)
# P(mu|y_1,y_2,...y_n) propto exp[n(ybar*mu-mu^2/2)]/(1+mu^2)
# we call this g(mu)
# log(g(mu))=n(ybar*mu-mu^2/2)-log(1+mu^2)
# 

#Likelihood function
lg<-function(ybar,n,mu){
  mu2=mu^2
  n*(ybar*mu-mu2/2)-log(1+mu^2)
}

#Metropolis Hastings function
mh<- function(n, ybar, n_iter, mu_init, cand_sd){
  mu_out=numeric(n_iter) #output mu
  accpt=0 #to know our acceptance rate
  mu_now=mu_init
  lg_now=lg(ybar, mu=mu_now,n=n)
  #Step2
  for (i in 1:n_iter){
    mu_cand=rnorm(1,mean=mu_now,sd=cand_sd)#a
    lg_cand=lg(ybar,n=n,mu=mu_cand)#b
    lalpha=lg_cand-lg_now
    alpha=exp(lalpha)
    
    u=runif(1)
    if (u<alpha){
      mu_now=mu_cand
      accpt=accpt+1
      lg_now=lg_cand
    }
    
    mu_out[i]=mu_now
  }
  list(mu=mu_out,accpt=accpt/n_iter)
}


### Set up the problem
y=c(1.2,1.4,1.5,1.6,0.8,0.9,1.1,-0.5,0,2,0.3,0.5,0.4,0.45,0.3,0.2,0.8,0.75,0.6)
hist(y,freq=F,xlim = c(-1,3))
points(y, rep(0.0,n))
ybar=mean(y)
n=length(y)
curve(dt(x,df=1),lty=2,add=TRUE)


#posterior sampling
set.seed(43)
post=mh(n=n, ybar = ybar, n_iter=1e3,mu_init = 0,cand_sd=0.7)
str(post)


#Diagnostics
library("coda")
traceplot(as.mcmc(post$mu))


#postanalysis
post$mu_keep=post$mu[-c(1:100)]
plot(density(post$mu_keep),xlim=c(-1,3))
curve(dt(x,df=1),lty=2,add=TRUE)
