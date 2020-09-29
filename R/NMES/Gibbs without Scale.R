# install.packages("GeneralizedHyperbolic")
# install.packages("mcmcse")
# install.packages("coda")
rm(list=ls())
set.seed(0)

library(MASS)
library(GeneralizedHyperbolic)
library(mcmcse)
library(coda)
library(AER)
data("NMES1988")

dat = NMES1988[c(1,8,11,15,16)] # Remove the responses & select numeric variables

# Center and scale the continuous predictors:
# dat[c(2,4)]= scale(dat[c(2,4)])



## Setup
# quantiles=c(0.5)
quantiles=c(0.1,0.5,0.9)
# dist_eps_list=c("stdN","heteroN","student")  ## heteroN requires d>=3
dist_eps_list=c("stdN")
itr_total = 1200
n=nrow(dat) 
d=ncol(dat) #dimension of x


for (dist_eps in dist_eps_list){
    ## Import Data
    x = t(cbind(matrix(rep(1,n),n,1),
                matrix(dat$chronic,n,1),
                matrix(dat$age,n,1),
                matrix(dat$school,n,1),
                matrix(dat$income,n,1)))
    y = matrix(dat$visits)
    z=rexp(n)
    
    ## Generic Gibbs sampler
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
            ## B_p_hat
            B_p_hat = solve(B_p0)
            for (i in 1:n){
              B_p_hat = B_p_hat + (x[,i]%*%t(x[,i]) / (tau2*z[i]))
            }
            B_p_hat = solve(B_p_hat)
            
            ## beta_p_hat
            beta_p_hat = solve(B_p0)%*%beta_p0
            for (i in 1:n){
              beta_p_hat = beta_p_hat + matrix(x[,i],d,1)*(y[i]-theta*z[i]) / (tau2*z[i])
            }
            beta_p_hat = B_p_hat %*% beta_p_hat
            
            beta_p = matrix(mvrnorm(1,beta_p_hat,B_p_hat),d,1)
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
        
        save(beta_p_record, file=paste0("NMES\\Gibbs without Scale\\","GWOS_",p,"_",dist_eps,"_betap.RData"))
        save(z_record, file=paste0("NMES\\Gibbs without Scale\\","GWOS_",p,"_",dist_eps,"_z.RData"))
        save(t, file=paste0("NMES\\Gibbs without Scale\\","GWOS_",p,"_",dist_eps,"_t.RData"))
    }
}

# plot(unlist(lapply(beta_p_record, function(x){return(x[1,])})))
