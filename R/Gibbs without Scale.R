rm(list=ls())
set.seed(0)
source("~/Study/Github/Bayesian-gibbs-sampler/R/data_sim1.R")
source("~/Study/Github/Bayesian-gibbs-sampler/R/Gibbs without Scale function.R")

## Setup
quantiles=c(0.1,0.5,0.9)
# dist_eps_list=c("stdN","heteroN","student")  ## heteroN requires d>=3
dist_eps_list=c("stdN")
itr_total = 12
n=100 
d=3 #dimension of x
beta_p=matrix(c(1,1,1))


for (dist_eps in dist_eps_list){
    ## Generate Data
    data = generate_data(n=n, d=d, dist_eps=dist_eps, beta_p=beta_p)
    x = data[[1]]; y=data[[2]]
    z=rexp(n)
    
    for (p in quantiles){
        ## Gibbs sampler without scale 
        result = GWOS(quantile=p, dist_eps=dist_eps, itr_total=itr_total, 
                      n=n, d=d, beta_p=beta_p)
        
        # save(result[[1]], file=paste0("~/Study/Github/Bayesian-gibbs-sampler/R/Gibbs without Scale/","GWOS_",p,"_",dist_eps,"_betap.RData"))
        # save(result[[2]], file=paste0("~/Study/Github/Bayesian-gibbs-sampler/R/Gibbs without Scale/","GWOS_",p,"_",dist_eps,"_z.RData"))
        # save(result[[3]], file=paste0("~/Study/Github/Bayesian-gibbs-sampler/R/Gibbs without Scale/","GWOS_",p,"_",dist_eps,"_t.RData"))
    }
}

