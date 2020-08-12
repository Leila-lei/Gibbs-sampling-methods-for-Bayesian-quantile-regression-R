# install.packages("matrixNormal")
# install.packages("GeneralizedHyperbolic")
# install.packages("mcmcse")
# install.packages("invgamma")
library(MASS)
library(matrixNormal)
library(GeneralizedHyperbolic)
library(mcmcse)
library(invgamma)
rm(list=ls())
set.seed(0)
setwd('C:\\Users\\juyil\\Dropbox\\Study\\GitHub\\Bayesian-gibbs-sampler\\R')
source("data_sim1.R")

## Setup
quantiles=c(0.1,0.5,0.9)
dist_eps_list=c("stdN","heteroN","student")  ## heteroN requires d>=3
# dist_eps_list=c("stdN","student")
itr_total = 12000
n=100 
d=3 #dimension of x
beta_p=matrix(c(1,1,1))
n0=3
s0=0.1

for (dist_eps in dist_eps_list){
  for (p in quantiles){
    ## Generate Data
    data = generate_data(n=n, d=d, dist_eps=dist_eps, beta_p=beta_p)
    x = data[[1]]; y=data[[2]]
    
    sigma=rinvgamma(1,n0/2,s0/2)
    u=rnorm(n,0,1)
    z=rexp(n)
    v=sigma*z
    
    theta = (1-2*p)/(p*(1-p))
    tau2 = 2/(p*(1-p))
    gamma_i_hat = 2/sigma + theta^2/(tau2*sigma)
    
    ## Gibbs sampler with scale
    ptm = proc.time()
    beta_p0=matrix(rep(0,d),d,1)
    B_p0=diag(d)*100
    beta_p = matrix(mvrnorm(1, beta_p0, B_p0))
    beta_p_record = rep(list(), itr_total)
    v_record = rep(list(), itr_total)
    for (itr in 1:itr_total){
      XtX = crossprod(t(x))
      XtY = crossprod(t(x),y)
      Q_theta = solve(B_p0)
      for (i in 1:n){
        Q_theta = Q_theta + (x[,i]%*%t(x[,i]) / (tau2*z[i]))
      }
      
      ell_theta = solve(B_p0)%*%beta_p0
      for (i in 1:n){
        ell_theta = ell_theta + matrix(x[,i],d,1)*(y[i]-theta*z[i]) / (tau2*z[i])
      }
      ch_Q = chol(Q_theta)
      beta_p = backsolve(ch_Q, forwardsolve(t(ch_Q), ell_theta) + rnorm(d))
      beta_p_record[[itr]] = beta_p
      
      ## v
      for (i in 1:n){
        delta_i_hat = (y[i]-(t(x[,1])%*%beta_p))^2/(tau2*sigma)
        v[i] = rgig(1, lambda=0.5, chi=delta_i_hat, psi=gamma_i_hat)
      }
      v_record[[itr]] = v
      
      if (itr%%1200==0){print(paste0("Dist:",dist_eps," p:",p," itr: ",itr,"/",itr_total))}
      
    }
    t = as.numeric(proc.time()-ptm)[3]
    
    save(beta_p_record, file=paste0("cholesky\\Gibbs with Scale\\","GWS_",p,"_",dist_eps,"_betap.RData"))
    save(v_record, file=paste0("cholesky\\Gibbs with Scale\\","GWS_",p,"_",dist_eps,"_v.RData"))
    save(t, file=paste0("cholesky\\Gibbs with Scale\\","GWS_",p,"_",dist_eps,"_t.RData"))
  }
}
