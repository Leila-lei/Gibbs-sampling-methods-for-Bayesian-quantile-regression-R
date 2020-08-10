# install.packages("matrixNormal")
# install.packages("GeneralizedHyperbolic")
# install.packages("mcmcse")
# install.packages("nimble")
# install.packages("purrr")
# install.packages("coda")
library(MASS)
library(matrixNormal)
library(GeneralizedHyperbolic)
library(mcmcse)
library(nimble)
library(coda)
library(quantreg)
rm(list=ls())
set.seed(0)
data(mtcars)
mtcars = mtcars[c(1,2,5,8)]

## Setup
# quantiles=c(0.5)
quantiles=c(0.1,0.5,0.9)
dist_eps_list=c("stdN","heteroN","student")  ## heteroN requires d>=3
# dist_eps_list=c("stdN","student")
itr_total = 12000
n=nrow(mtcars) 
d=ncol(mtcars) #dimension of x
lambda0 = 0.14

## Generate Data
for (dist_eps in dist_eps_list){
    for (p in quantiles){
        ## Import Data
        x = t(cbind(matrix(rep(1,n),n,1),matrix(mtcars$cyl,n,1),matrix(mtcars$drat,n,1),matrix(mtcars$vs,n,1)))
        y = matrix(mtcars$mpg)
        
        ## Gibbs sampler with Double Exponential Prior
        ptm = proc.time()
        beta_p0=matrix(rep(0,d),d,1)
        omega=rexp(d, 2/(lambda0^2))
        Omega=diag(d)*omega
        beta_p = matrix(mvrnorm(1,beta_p0,Omega))
        beta_p_record = rep(list(), itr_total)
        omega_record = rep(list(), itr_total)
        for (itr in 1:itr_total){
          ## omega
          for (i in 1:d){
            omega[i] = rgig(1, 0.5, abs(beta_p-beta_p0)[i,], lambda0)
          }
          omega_record[[itr]] = omega
          Omega=diag(d)*omega
          
          ## beta_p
          beta_p = matrix(mvrnorm(1,beta_p0,Omega),d,1)
          beta_p_record[[itr]] = beta_p
          
          
          if (itr%%1200==0){print(paste0("Dist:",dist_eps," p:",p," itr: ",itr,"/",itr_total))}
        }
        t = as.numeric(proc.time()-ptm)[3]
        save(beta_p_record, file=paste0("C:\\Users\\juyil\\Documents\\la\\R\\mtcars\\Double Exponential Prior\\","DEP_",p,"_",dist_eps,"_betap.RData"))
        save(omega_record, file=paste0("C:\\Users\\juyil\\Documents\\la\\R\\mtcars\\Double Exponential Prior\\","DEP_",p,"_",dist_eps,"_v.RData"))
        save(t, file=paste0("C:\\Users\\juyil\\Documents\\la\\R\\mtcars\\Double Exponential Prior\\","DEP_",p,"_",dist_eps,"_t.RData"))
    }
}






