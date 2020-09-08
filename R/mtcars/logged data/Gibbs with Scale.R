# install.packages("matrixNormal")
# install.packages("GeneralizedHyperbolic")
# install.packages("mcmcse")
# install.packages("invgamma")
library(MASS)
library(matrixNormal)
library(GeneralizedHyperbolic)
library(mcmcse)
library(invgamma)
library(quantreg)
rm(list=ls())
setwd("C:\\Users\\juyil\\Dropbox\\Study\\GitHub\\Bayesian-gibbs-sampler\\R\\mtcars\\logged data")
set.seed(0)
data(mtcars)
mtcars = mtcars[c(1,2,5,8)]
## Log transform
mtcars[2:3] = log(mtcars[2:3])

## Setup
# quantiles=c(0.5)
quantiles=c(0.1,0.5,0.9)
dist_eps_list=c("stdN")  ## heteroN requires d>=3
# dist_eps_list=c("stdN","student")
itr_total = 12000
n=nrow(mtcars) 
d=ncol(mtcars) #dimension of x
n0=3
s0=0.1

for (dist_eps in dist_eps_list){
  for (p in quantiles){
    ## Import Data
    x = t(cbind(matrix(rep(1,n),n,1),matrix(mtcars$cyl,n,1),matrix(mtcars$drat,n,1),matrix(mtcars$vs,n,1)))
    y = matrix(mtcars$mpg)
    
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
      ## B_p_hat
      B_p_hat = solve(B_p0)
      for (i in 1:n){
        B_p_hat = B_p_hat + (matrix(x[,i])%*%t(x[,i]) / (tau2*sigma*v[i]))
      }
      B_p_hat = solve(B_p_hat)
      
      ## beta_p_hat
      beta_p_hat = solve(B_p0)%*%beta_p0
      for (i in 1:n){
        beta_p_hat = beta_p_hat + matrix(x[,i],d,1)*(y[i]-theta*v[i]) / (tau2*sigma*v[i])
      }
      beta_p_hat = B_p_hat %*% beta_p_hat
      
      beta_p = matrix(mvrnorm(1,beta_p_hat,B_p_hat),d,1)
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
    
    save(beta_p_record, file=paste0("Gibbs with Scale\\","GWS_",p,"_",dist_eps,"_betap.RData"))
    save(v_record, file=paste0("Gibbs with Scale\\","GWS_",p,"_",dist_eps,"_v.RData"))
    save(t, file=paste0("Gibbs with Scale\\","GWS_",p,"_",dist_eps,"_t.RData"))
  }
}

