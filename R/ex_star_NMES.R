rm(list = ls())
#----------------------------------------------------------------------------
# Application for STAR-MCMC: NMES data

# Source the necessary files:
library(rSTAR)

# Relevant packages:
library(coda)
library(rstanarm)
# library(devtools); devtools::install_github("vdorie/dbarts")
library(dbarts)

# For computing multivariate ESS:
library(mcmcse)

# Include stan?
use_stan = TRUE
nsave = 1000; nburn = 1000; nskip = 2
#----------------------------------------------------------------------------
# Source: US National Medical Expenditure Survey [NMES] data for 1987/88
  # Available as NMES1988 in package AER (Kleiber and Zeileis 2008).
  # Originally taken from Deb and Trivedi (J. Applied Econometrics 1997)
library(AER)
data("NMES1988")

# IDEA: use the first 4 as examples of Deb and Trivedi (w/ increasing zero inflation)
#response = 'visits'   # Number of physician office visits.
response = 'nvisits'  # Number of non-physician office visits.
#response = 'ovisits'  # Number of physician hospital outpatient visits.
#response = 'novisits' # Number of non-physician hospital outpatient visits.

# Plot the variables together:
# dev.new(); par(mfrow = c(2,2))
stickplot = function(y, ...){js = 0:max(y); plot(js, sapply(js, function(js) mean(js == y)), type='h', lwd=2, ...)}
# stickplot(NMES1988[,'visits'], main = 'Physician office visits (visits)', xlab = 'Visits', ylab = 'Probability mass', cex.axis = 1.5, cex.lab = 1.5)
# stickplot(NMES1988[,'nvisits'], main = 'Non-physician office visits (nvisits)', xlab = 'Visits', ylab = 'Probability mass', cex.axis = 1.5, cex.lab = 1.5)
# stickplot(NMES1988[,'ovisits'], main = 'Physician hospital outpatient visits (ovisits)', xlab = 'Visits', ylab = 'Probability mass', cex.axis = 1.5, cex.lab = 1.5)
# stickplot(NMES1988[,'novisits'], main = 'Non-physician hospital outpatient visits (novisits)', xlab = 'Visits', ylab = 'Probability mass', cex.axis = 1.5, cex.lab = 1.5)

# Remove the values with negative incomes and log-transform:
#NMES1988 = NMES1988[NMES1988$income >=0 ,]; NMES1988$age = log(NMES1988$age); NMES1988$income = log(NMES1988$income + 1);

# Center and scale the continuous predictors:
dat = NMES1988[,-(1:6)] # Remove the responses
dat[,c('age', 'school', 'income')] = scale(dat[,c('age', 'school', 'income')])
X = model.matrix( ~ ., data = dat)

# Response variable:
y = NMES1988[, response]

# Emergency room visits:
# y = NMES1988$emergency;
# X = cbind(X, as.matrix(NMES1988[,c('visits', 'nvisits', 'ovisits', 'novisits', 'hospital')]))

# Number of hospital stays:
#y = NMES1988$hospital;
#X = cbind(X, as.matrix(NMES1988[,c('visits', 'nvisits', 'ovisits', 'novisits', 'emergency')]))

# sapply(1:6, function(j){plot(NMES1988[,j], main = colnames(NMES1988)[j]); hist(NMES1988[,j], main = colnames(NMES1988)[j]); paste(colnames(NMES1988)[j], 'prop zero:', round(mean(NMES1988[,j]==0), 2))})

# Sample size:
n = nrow(X); p = ncol(X)

# Summary:
stickplot(y)
#----------------------------------------------------------------------------
# Linear Models
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Gaussian LM for y:
fit_gauss = Gauss_MCMC(y = y,
                       sample_params = function(y, params) sample_params_lm_hs(y, X, params),
                       init_params = function(y) init_params_lm_hs(y, X),
                       nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#multiESS(fit_gauss$post.coefficients[,1:p])
#----------------------------------------------------------------------------
# Gaussian LM for sqrt(y):
fit_gauss_sqrt = Gauss_MCMC(y = sqrt(y),
                           sample_params = function(y, params) sample_params_lm_hs(y, X, params),
                           init_params = function(y) init_params_lm_hs(y, X),
                           nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#multiESS(fit_gauss_sqrt$post.coefficients[,1:p])

# WAIC for the sqrt-transformed data:
fit_gauss_sqrt$post.log.like.point = fit_gauss_sqrt$post.log.like.point +
  log(0.5*(rep(y, each = length(fit_gauss_sqrt$post.sigma)))^(-1/2))
lppd = sum(log(colMeans(exp(fit_gauss_sqrt$post.log.like.point))))
fit_gauss_sqrt$p_waic = sum(apply(fit_gauss_sqrt$post.log.like.point, 2, function(x) sd(x)^2))
fit_gauss_sqrt$WAIC = -2*(lppd - fit_gauss_sqrt$p_waic)
#----------------------------------------------------------------------------
# Gaussian LM for log(y+1):
fit_gauss_log = Gauss_MCMC(y = log(y+1),
                           sample_params = function(y, params) sample_params_lm_hs(y, X, params),
                           init_params = function(y) init_params_lm_hs(y, X),
                           nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#multiESS(fit_gauss_log$post.coefficients[,1:p])

# WAIC for the log-transformed data:
fit_gauss_log$post.log.like.point = fit_gauss_log$post.log.like.point +
  log(abs(1/(rep(y, each = length(fit_gauss_log$post.sigma)) + 1)))
lppd = sum(log(colMeans(exp(fit_gauss_log$post.log.like.point))))
fit_gauss_log$p_waic = sum(apply(fit_gauss_log$post.log.like.point, 2, function(x) sd(x)^2))
fit_gauss_log$WAIC = -2*(lppd - fit_gauss_log$p_waic)
#----------------------------------------------------------------------------
# Poisson regression
if(use_stan){
  stan_pois = stan_glm(y ~ X-1, family = poisson(link="log"),
                       prior = normal(0, 2.5), #prior = hs(),
                       prior_intercept = normal(0, 5),
                       data = data.frame(y = y, X = X),
                       chains = 1, iter = nburn+(nskip+1)*(nsave), warmup = nburn, thin = nskip + 1)
  #multiESS(as.matrix(stan_pois)[,1:p])
}
#----------------------------------------------------------------------------
# Negative binomial regression
if(use_stan){
  stan_nb = stan_glm.nb(y ~ X-1, link = "log", prior_aux = exponential(1),
                        prior = normal(0, 2.5), #prior = hs(),
                        prior_intercept = normal(0, 5),
                        data = data.frame(y = y, X = X),
                        chains = 1, iter = nburn+(nskip+1)*(nsave), warmup = nburn, thin = nskip + 1)
  #multiESS(as.matrix(stan_nb)[,1:p])
}
#----------------------------------------------------------------------------
# STAR: log-transformation:
fit_star_log = star_MCMC(y = y,
                         sample_params = function(y, params) sample_params_lm_hs(y, X, params),
                         init_params = function(y) init_params_lm_hs(y, X),
                         transformation = 'log',
                         nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#multiESS(fit_star_log$post.coefficients[,1:p])
#----------------------------------------------------------------------------
# STAR: sqrt-transformation:
fit_star_sqrt = star_MCMC(y = y,
                          sample_params = function(y, params) sample_params_lm_hs(y, X, params),
                          init_params = function(y) init_params_lm_hs(y, X),
                          transformation = 'sqrt',
                          nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#multiESS(fit_star_sqrt$post.coefficients[,1:p])
#----------------------------------------------------------------------------
# STAR: identity-transformation:
fit_star_id = star_MCMC(y = y,
                        sample_params = function(y, params) sample_params_lm_hs(y, X, params),
                        init_params = function(y) init_params_lm_hs(y, X),
                        transformation = 'identity',
                        nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#multiESS(fit_star_id$post.coefficients[,1:p])
#----------------------------------------------------------------------------
# STAR: unknown Box-Cox transformation
temp0 = proc.time()[3]
fit_star_bc = star_MCMC(y = y,
                        sample_params = function(y, params) sample_params_lm_hs(y, X, params),
                        init_params = function(y) init_params_lm_hs(y, X),
                        transformation = 'box-cox',
                        nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#multiESS(fit_star_bc$post.coefficients[,1:p])
time0 = proc.time()[3] - temp0
multiESS(fit_star_bc$post.coefficients[,1:p])/time0
#----------------------------------------------------------------------------
# STAR: unknown I-spline transformation
fit_star_np = star_np_MCMC(y = y,
                         sample_params = function(y, params) sample_params_lm_hs(y, X, params),
                         init_params = function(y) init_params_lm_hs(y, X),
                         nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#multiESS(fit_star_np$post.coefficients[,1:p])
#----------------------------------------------------------------------------
# Plot the LM results:
#----------------------------------------------------------------------------
fit_mcmc = fit_star_log # Pick a model
#----------------------------------------------------------------------------
#dev.new();
if(use_stan){jit = 0.05} else jit = 0 # Jitter for Stan
ci = t(apply(fit_mcmc$post.coefficients[,1:p], 2, quantile,
             c(0.05/2, 1-0.05/2)))
ylim = range(ci[-1,])
if(use_stan) ylim = range(ylim,  posterior_interval(stan_nb, prob = 0.95)[2:p,])
plot(2:p, 2:p, ylim = ylim, type='n', ylab='', xaxt='n', xlab = '',
     main = paste('Regression Coefficients: NMES Data (', response, ')', sep=""))
axis(1, at = 2:p, labels=FALSE)
text(2:p, par("usr")[3] - 0.08, labels = colnames(X)[-1], srt = 45, pos = 1, xpd = TRUE)
#axis(1, at = 2:p, colnames(X)[-1])
if(use_stan){
  lines(1:p - jit, stan_nb$coefficients[1:p], type='p', pch = 1, lwd=6, cex = 2, col='darkgray')
  arrows(1:p - jit, posterior_interval(stan_nb, prob = 0.95)[,1], 1:p - jit, posterior_interval(stan_nb, prob = 0.95)[,2],
         length=0.08, angle=90, code=3, lwd=8, col='darkgray')
}
lines(1:p + jit, colMeans(fit_mcmc$post.coefficients[,1:p]),
      type='p', pch=4, lwd = 6, cex = 2, col='black')
arrows(1:p + jit, ci[,1], 1:p + jit, ci[,2],
       length=0.08, angle=90, code=3, lwd=8, col='black')
abline(h = 0, lwd=3, col='gray', lty=2)
legend('bottomright', c('NegBin (Stan)', 'LM-STAR-log'), lwd=8, col=c('darkgray', 'black'))
#----------------------------------------------------------------------------
# BART Models
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# BART for y:
fit_bart_raw = bart(x.train = X, y.train = y,
                    ntree = 200, ndpost = nsave, nskip = nburn, verbose = TRUE)

# WAIC for the raw data:
fit_bart_raw$post.log.like.point = array(NA, c(length(fit_bart_raw$sigma), n))
for(nsi in 1:length(fit_bart_raw$sigma))
  fit_bart_raw$post.log.like.point[nsi,] = dnorm(y,
                                                 mean = fit_bart_raw$yhat.train[nsi,],
                                                 sd = fit_bart_raw$sigma[nsi], log=TRUE)
lppd_raw = sum(log(colMeans(exp(fit_bart_raw$post.log.like.point))))
fit_bart_raw$p_waic = sum(apply(fit_bart_raw$post.log.like.point, 2, function(x) sd(x)^2))
fit_bart_raw$WAIC = -2*(lppd_raw - fit_bart_raw$p_waic)
#----------------------------------------------------------------------------
# BART for log(y+1)
fit_bart_log = bart(x.train = X, y.train = log(y+1),
                    ntree = 200, ndpost = nsave, nskip = nburn, verbose = TRUE)

# WAIC for the log-transformed data:
fit_bart_log$post.log.like.point = array(NA, c(length(fit_bart_log$sigma), n))
for(nsi in 1:length(fit_bart_log$sigma))
  fit_bart_log$post.log.like.point[nsi,] = log(abs(1/(y + 1))) + dnorm(log(y+1),
                                                                       mean = fit_bart_log$yhat.train[nsi,],
                                                                       sd = fit_bart_log$sigma[nsi], log=TRUE)
lppd_log = sum(log(colMeans(exp(fit_bart_log$post.log.like.point))))
fit_bart_log$p_waic = sum(apply(fit_bart_log$post.log.like.point, 2, function(x) sd(x)^2))
fit_bart_log$WAIC = -2*(lppd_log - fit_bart_log$p_waic)
#----------------------------------------------------------------------------
# BART: log-transformation:
fit_bart_star_log = bart_star_MCMC(y = y, X = X,
                                   transformation = 'log',
                                   nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#----------------------------------------------------------------------------
# BART: sqrt-transformation:
fit_bart_star_sqrt = bart_star_MCMC(y = y, X = X,
                                    transformation = 'sqrt',
                                    nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#----------------------------------------------------------------------------
# BART: identity-transformation:
fit_bart_star_id = bart_star_MCMC(y = y, X = X,
                                  transformation = 'identity',
                                  nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#----------------------------------------------------------------------------
# BART: unknown Box-Cox transformation
fit_bart_star_bc = bart_star_MCMC(y = y, X = X,
                                  transformation = 'box-cox',
                                  sigest = median(fit_star_bc$post.sigma),
                                  nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#----------------------------------------------------------------------------
# BART: unknown I-spline transformation:
fit_bart_star_np = bart_star_np_MCMC(y = y, X = X,
                                   verbose = TRUE)
#----------------------------------------------------------------------------
# Additive Models
#----------------------------------------------------------------------------
# Define the linear/nonlinear predictors:
id.nonlin = which(apply(X, 2, function(x) length(unique(x)) > 10))
X_nonlin = as.matrix(X[,id.nonlin]); X_lin = as.matrix(X[,-id.nonlin])

# Basis matrices for all nonlinear predictors:
library(spikeSlabGAM)
#B_all = lapply(1:ncol(X_nonlin), function(j) splineBasis(X_nonlin[,j], sumToZero = TRUE))
B_all = lapply(1:ncol(X_nonlin), function(j) {B0 = sm(X_nonlin[,j]); B0/sqrt(sum(diag(crossprod(B0))))})

# And a recurring term (one-time cost):
# BtB_all = lapply(B_all, crossprod)
diagBtB_all = lapply(1:ncol(X_nonlin), function(j) colSums(B_all[[j]]^2))
#----------------------------------------------------------------------------
# Gaussian AM for y:
fit_am_gauss = Gauss_MCMC(y = y,
                          sample_params = function(y,params)
                            sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, diagBtB_all = diagBtB_all, A = 10^6),
                            #sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, BtB_all = BtB_all, A = 10^6),
                          init_params = function(y,params)
                            init_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, B_all = B_all),
                          nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#----------------------------------------------------------------------------
# Gaussian AM for log(y+1):
fit_am_gauss_log = Gauss_MCMC(y = log(y+1),
                              sample_params = function(y,params)
                                sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, diagBtB_all = diagBtB_all, A = 10^6),
                              #sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, BtB_all = BtB_all, A = 10^6),
                              init_params = function(y,params)
                                init_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, B_all = B_all),
                              nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)

# WAIC for the log-transformed data:
fit_am_gauss_log$post.log.like.point = fit_am_gauss_log$post.log.like.point +
  log(abs(1/(rep(y, each = length(fit_am_gauss_log$post.sigma)) + 1)))
lppd = sum(log(colMeans(exp(fit_am_gauss_log$post.log.like.point))))
fit_am_gauss_log$p_waic = sum(apply(fit_am_gauss_log$post.log.like.point, 2, function(x) sd(x)^2))
fit_am_gauss_log$WAIC = -2*(lppd - fit_am_gauss_log$p_waic)
#----------------------------------------------------------------------------
# STAR AM: log-transformation
fit_am_star_log = star_MCMC(y = y,
                            sample_params = function(y,params)
                              sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, diagBtB_all = diagBtB_all, A = 10^6),
                            #sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, BtB_all = BtB_all, A = 10^6),
                            init_params = function(y,params)
                              init_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, B_all = B_all),
                            transformation = 'log',
                            nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#----------------------------------------------------------------------------
# STAR AM: sqrt-transformation
# fit_am_star_sqrt = star_MCMC(y = y,
#                              sample_params = function(y,params)
#                                sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, diagBtB_all = diagBtB_all, A = 10^6),
#                              #sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, BtB_all = BtB_all, A = 10^6),
#                              init_params = function(y,params)
#                                init_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, B_all = B_all),
#                              transformation = 'sqrt',
#                              nsave = nsave, nburn = nburn, verbose = TRUE)
#----------------------------------------------------------------------------
# STAR AM: identity-transformation:
# fit_am_star_id = star_MCMC(y = y,
#                            sample_params = function(y,params)
#                              sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, diagBtB_all = diagBtB_all, A = 10^6),
#                            #sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, BtB_all = BtB_all, A = 10^6),
#                            init_params = function(y,params)
#                              init_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, B_all = B_all),
#                            transformation = 'identity',
#                            nsave = nsave, nburn = nburn, verbose = TRUE)
#----------------------------------------------------------------------------
# STAR AM: unknown Box-Cox transformation
fit_am_star_bc = star_MCMC(y = y,
                           sample_params = function(y,params)
                             sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, diagBtB_all = diagBtB_all, A = 10^6),
                           #sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, BtB_all = BtB_all, A = 10^6),
                           init_params = function(y,params)
                             init_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, B_all = B_all),
                           transformation = 'box-cox',
                           nsave = nsave, nburn = nburn, nskip = nskip, verbose = TRUE)
#----------------------------------------------------------------------------
# STAR AM: unknown I-spline transformation
fit_am_star_np = star_np_MCMC(y = y,
                            sample_params = function(y,params)
                              sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, diagBtB_all = diagBtB_all, A = 10^6),
                            #sample_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, params = params, B_all = B_all, BtB_all = BtB_all, A = 10^6),
                            init_params = function(y,params)
                              init_params_additive(y = y, X_lin = X_lin, X_nonlin = X_nonlin, B_all = B_all),
                            verbose = TRUE)

# Clean up:
rm(B_all)#, BtB_all)
#----------------------------------------------------------------------------
# Plot the AM results:
#----------------------------------------------------------------------------
fit_mcmc_am = fit_am_star_log # Pick a model
#----------------------------------------------------------------------------
# Variable names:
vnames = colnames(fit_mcmc_am$post.coefficients)

# Nonlinear curves:
post_fj = array(fit_mcmc_am$post.coefficients[,which('f_j' == substr(vnames, 0, 3))],
                c(length(fit_mcmc_am$post.sigma), n, ncol(X_nonlin)))
# Plot the resulting additive functions:
for(j in 1:ncol(X_nonlin)){
  x = X_nonlin[,j];
  j.ind = order(x); x = x[j.ind];
  ci = t(apply(post_fj[,j.ind,j], 2, quantile, c(.025, .975)))
  #ci = HPDinterval(as.mcmc(post_fj[,j.ind,j]))
  plot(x, ci[,1], ylim = range(ci),  main = colnames(X_nonlin)[j], type='n')
  polygon(c(x, rev(x)), c(ci[,2], rev(ci[,1])), col='grey', border=NA)
  lines(x, colMeans(post_fj[,j.ind,j]), lwd=6, col='blue')
  lines(x, rep(min(ci), length(x)), type= 'p', pch = '|')
  abline(h = 0, lty=6, lwd=2)
}
getEffSize(post_fj)

# Clean up:
rm(post_fj)
#----------------------------------------------------------------------------
# Posterior predictive diagnostics:
#----------------------------------------------------------------------------
if(FALSE){
  # Select a model:
  fit_mcmc = fit_star_log

  # Posterior predictive diagnostics:
  dev.new();
  par(mfrow = c(2,3))
  post_pred_summary = as.matrix(apply(fit_mcmc$post.pred, 1, summary))
  for(j in 2:5) {hist(post_pred_summary[j,], main = rownames(post_pred_summary)[j], xlab=response); abline(v = summary(as.numeric(y))[j], lwd=4, col ='blue')}
  hist(apply(fit_mcmc$post.pred, 1, sd), main = 'SD', xlab=response); abline(v = sd(as.numeric(y)), lwd=4, col ='blue')
  hist(apply(fit_mcmc$post.pred, 1, function(x) mean(x==0)), main = 'Proportion of Zeros', xlab=response); abline(v = mean(y==0), lwd=4, col ='blue')

  # PMF comparison:
  plot_pmf(y, fit_mcmc$post.pred)
}
#----------------------------------------------------------------------------
# Model comparisons: WAIC and LOO
#----------------------------------------------------------------------------
# Linear models: WAIC and LOO
WAIC_lm = c(
  # LM-Gauss:
  fit_gauss$WAIC,
  fit_gauss_log$WAIC,

  # LM-STAR:
  fit_star_log$WAIC,
  fit_star_sqrt$WAIC,
  fit_star_id$WAIC,
  fit_star_bc$WAIC,
  fit_star_np$WAIC
); names(WAIC_lm) = c('LM-Gauss-id', 'LM-Gauss-log', 'LM-STAR-log', 'LM-STAR-sqrt', 'LM-STAR-id', 'LM-STAR-bc', 'LM-STAR-np')
if(use_stan) {WAIC_lm = c(WAIC_lm, waic(stan_pois)$estimates['waic',1], waic(stan_nb)$estimates['waic',1]); names(WAIC_lm)[(length(WAIC_lm)-1):length(WAIC_lm)] =  c('LM-Pois', 'LM-NegBin')}

# BART: WAIC and LOO
WAIC_bart = c(
  # BART-Gaussian:
  fit_bart_raw$WAIC,
  fit_bart_log$WAIC,

  # BART-STAR
  fit_bart_star_log$WAIC,
  fit_bart_star_sqrt$WAIC,
  fit_bart_star_id$WAIC,
  fit_bart_star_bc$WAIC,
  fit_bart_star_np$WAIC
); names(WAIC_bart) = c('BART-id', 'BART-log', 'BART-STAR-log', 'BART-STAR-sqrt', 'BART-STAR-id', 'BART-STAR-bc', 'BART-STAR-np')

# Additive models: WAIC and LOO
WAIC_am = c(
  # AM-Gauss:
  fit_am_gauss$WAIC,
  fit_am_gauss_log$WAIC,

  # AM-STAR:
  fit_am_star_log$WAIC,
  NA, #fit_am_star_sqrt$WAIC,
  NA, #fit_am_star_id$WAIC,
  fit_am_star_bc$WAIC,
  fit_am_star_np$WAIC
); names(WAIC_am) = c('AM-Gauss-id', 'AM-Gauss-log', 'AM-STAR-log', 'AM-STAR-sqrt', 'AM-STAR-id', 'AM-STAR-bc', 'AM-STAR-np')

library(xtable)
WAIC_all = rbind(WAIC_lm,
                 c(WAIC_am, rep(NA,2)),
                 c(WAIC_bart, rep(NA, 2)));
colnames(WAIC_all) = c("Gauss-id", "Gauss-log",
                       "STAR-log", "STAR-sqrt", "STAR-id", "STAR-bc", "STAR-np",
                       'Pois', 'NB');
rownames(WAIC_all) = c('Linear', 'Additive', 'BART')
xtable(WAIC_all, digits = 0)

# LOO:
library(loo)
LOOs = compare(
  # LM-Gauss:
  loo(fit_gauss$post.log.like.point),
  loo(fit_gauss_log$post.log.like.point),

  # LM-STAR:
  loo(fit_star_log$post.log.like.point),
  loo(fit_star_sqrt$post.log.like.point),
  loo(fit_star_id$post.log.like.point),
  loo(fit_star_bc$post.log.like.point),
  loo(fit_star_np$post.log.like.point),

  # Stan:
  loo(stan_pois),
  loo(stan_nb),

  # BART-Gaussian:
  loo(fit_bart_raw$post.log.like.point),
  loo(fit_bart_log$post.log.like.point),

  # BART-STAR
  loo(fit_bart_star_log$post.log.like.point),
  loo(fit_bart_star_sqrt$post.log.like.point),
  loo(fit_bart_star_id$post.log.like.point),
  loo(fit_bart_star_bc$post.log.like.point),
  loo(fit_bart_star_np$post.log.like.point),

  # AM-Gauss:
  loo(fit_am_gauss$post.log.like.point),
  loo(fit_am_gauss_log$post.log.like.point),

  # AM-STAR:
  loo(fit_am_star_log$post.log.like.point),
  # loo(fit_am_star_sqrt$post.log.like.point),
  # loo(fit_am_star_id$post.log.like.point),
  loo(fit_am_star_bc$post.log.like.point),
  loo(fit_am_star_np$post.log.like.point)
); print(LOOs)

# All WAICs (Ordered)
print(WAIC_lm[order(WAIC_lm)])
print(WAIC_am[order(WAIC_am)])
print(WAIC_bart[order(WAIC_bart)])

# Save the workspace:
#save.image(paste("~/Dropbox/Projects/Count Data Analysis/R code/Results/Applications/nmes_", response, ".RData", sep = ""))
