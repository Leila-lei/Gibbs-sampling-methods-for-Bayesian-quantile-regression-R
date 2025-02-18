# install.packages("GeneralizedHyperbolic")
# install.packages("mcmcse")
# install.packages("purrr")
# install.packages("coda")
library(MASS)
library(GeneralizedHyperbolic)
library(mcmcse)
library(purrr)
library(coda)
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


for (dist_eps in dist_eps_list){
  ## Generate Data
  data = generate_data(n=n, d=d, dist_eps=dist_eps, beta_p=beta_p)
  x = data[[1]]; y=data[[2]]
  y[which(y < 0)] = 0
  
  z=rexp(n)
  
  ## Tobit without scale
  for (p in quantiles){
    ptm = proc.time()
    ## Initialization
    beta_p0=matrix(rep(0,d),d,1)
    B_p0=diag(d)*100
    beta_p = matrix(mvrnorm(1, beta_p0, B_p0))
    theta = (1-2*p)/(p*(1-p))
    tau2 = 2/(p*(1-p))
    
    beta_p_record = rep(list(), itr_total)
    z_record = rep(list(), itr_total)
    
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
      
      ## z
      for (i in 1:n){
        delta_i_hat = sqrt((y[i]-(t(x[,i])%*%beta_p))^2/tau2)
        gamma_i_hat = sqrt(theta^2/tau2 +2)
        z[i] = rgig(1, param=c(0.5,delta_i_hat,gamma_i_hat))
      }
      z_record[[itr]] = z
      
      if (itr%%1200==0){print(paste0("Dist:",dist_eps," p:",p," itr: ",itr,"/",itr_total, " beta_p: ",beta_p))}
    }
    t = as.numeric(proc.time()-ptm)[3]
    
    save(beta_p_record, file=paste0("cholesky\\Tobit without Scale\\","TWOS_",p,"_",dist_eps,"_betap.RData"))
    save(z_record, file=paste0("cholesky\\Tobit without Scale\\","TWOS_",p,"_",dist_eps,"_z.RData"))
    save(t, file=paste0("cholesky\\Tobit without Scale\\","TWOS_",p,"_",dist_eps,"_t.RData"))
  }
}


