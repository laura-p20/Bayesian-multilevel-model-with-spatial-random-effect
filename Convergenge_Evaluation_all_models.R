#============================================================================---
#               Model Evaluation 
# This script corresponds to the evaluation using the samples generated with 
# the model Gibbs Sampler of the model without any shrinkage applied, for each of the
# 3 simulation scenarios
#============================================================================---


#===========Libraries and settings====

setwd("/Users/macbookpro/Desktop/Tesis/Results")

#===============Libraries and preparation----

library(dplyr)
library(tidyverse)
library(ggplot2)
library(sf)
library(spdep)
library(dplyr)
library(stringi)
library(stringr)
suppressMessages(suppressWarnings(library(coda)))

#Set the work directory to get a previus result of the MCMC using the real data

setwd("/Users/macbookpro/Desktop/Tesis/Results")
load("/Users/macbookpro/Desktop/Tesis/Procesed_Data/data_ready.RData")

#Upload the covariables matrices and some numbers
X_ijk<-data_ready$X_ijk
Z_jk<-data_ready$Z_jk
Z_jk<-scale(Z_jk)
W_k<-data_ready$W_k
W_k<-scale(W_k)
n<-data_ready$n
n_jk_1<-data_ready$n_jk
m<-data_ready$m
d<-data_ready$d

#Upload the response variable (global score)

y_ijk<-data_ready$punt_global

#============= Functions -----

CV<-function(x){
  cv<-sd(x)/abs(mean(x))*100
  return(cv)
}

bayesian_summary<-function(x,alpha=0.05){
  posterior_mean<-mean(x)
  credibility_interval<-quantile(x,c(alpha/2,1-(alpha/2)))
  CV_samples<-sd(x)/abs(mean(x))*100
  info<-c(posterior_mean,credibility_interval,CV_samples)
  return(info)
  cat("Posterior mean: ",posterior_mean,"\n","Credibility interval ",alpha," level: ",credibility_interval,"\n","CV samples (%): ",CV_samples)
}

bayesian_summary2<-function(x,alpha=0.05){
  posterior_mean<-mean(x)
  credibility_interval<-quantile(x,c(alpha/2,1-(alpha/2)))
  CV_samples<-sd(x)/abs(mean(x))*100
  cat("Posterior mean: ",posterior_mean,"\n","Credibility interval ",alpha," level: ",credibility_interval,"\n","CV samples (%): ",CV_samples,"\n","\n")
}


stat_mode<-function(x){
  freq<-hist(x,plot=FALSE)
  mode<-freq$mids[which(freq$counts==max(freq$counts))]
  return(mode)
}

##===========Extracting the samples of each model----
#Basis model 
#load("~/Desktop/Tesis/Results/MCMC_model2_real_data_standarized_.RData")
load("~/Desktop/Tesis/Results/MCMC_model2_real_data_standarized.RData")

beta<-THETA[[1]]
beta_E<-THETA[[2]]
beta_M<-THETA[[3]]
beta_D<-THETA[[4]]

sigma2_beta<-THETA[[5]]
sigma2_E<-THETA[[6]]
sigma2_M<-THETA[[7]]
sigma2_D<-THETA[[8]]

phi<-THETA[[9]]
tau_phi<-THETA[[10]]

kappa2_jk<-THETA[[11]]
kappa2_k<-THETA[[12]]

alpha_kappa<-THETA[[13]]
beta_kappa<-THETA[[14]]


#Ridge samples
#load("~/Desktop/Tesis/Results/MCMC_ridge_real_data_non_standarized.RData")
load("~/Desktop/Tesis/Results/MCMC_ridge_real_data.RData")

betar<-THETA[[1]]
beta_Er<-THETA[[2]]
beta_Mr<-THETA[[3]]
beta_Dr<-THETA[[4]]

sigma2_betar<-THETA[[5]]
lambda2_Er<-THETA[[6]]
lambda2_Mr<-THETA[[7]]
lambda2_Dr<-THETA[[8]]

phir<-THETA[[9]]
tau_phir<-THETA[[10]]

kappa2_jkr<-THETA[[11]]
kappa2_kr<-THETA[[12]]

alpha_kappar<-THETA[[13]]
beta_kappar<-THETA[[14]]

acr_r<-THETA[[15]]




# Lasso samples
load("~/Desktop/Tesis/Results/MCMC_lasso_real_data_standarized.RData")

betal<-THETA[[1]]
beta_El<-THETA[[2]]
beta_Ml<-THETA[[3]]
beta_Dl<-THETA[[4]]

sigma2_betal<-THETA[[5]]
tau2_El<-THETA[[6]]
tau2_Ml<-THETA[[7]]
tau2_Dl<-THETA[[8]]

lambda2_El<-THETA[[9]]
lambda2_Ml<-THETA[[10]]
lambda2_Dl<-THETA[[11]]

phil<-THETA[[12]]
tau_phil<-THETA[[13]]


kappa2_jkl<-THETA[[14]]
kappa2_kl<-THETA[[15]]

alpha_kappal<-THETA[[16]]
beta_kappal<-THETA[[17]]

acr_l<-THETA[[18]]

nsams<-length(betal)


#============= Convergence Analysis  ========

##========= Log-Likelyhood=========


ll_basis<-rep(NA,nsams)
ll_ridge<-rep(NA,nsams)
ll_lasso<-rep(NA,nsams)

#Calculus of the log-likelyhoods based on the samples

#Baasis model
for(i in c(1:nsams)){
  
  mean<-rep(beta[i],n)+X_ijk%*%beta_E[i,]+Z_jk%*%beta_M[i,]+W_k%*%beta_D[i,]+rep(phi[i,],n_jk)
  v<-rep(sqrt(kappa2_jk[i,]),n_jk)
  ll_basis[i]<-sum(dnorm(x=y_ijk,mean=mean,sd=v,log = TRUE))
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
}

#Ridge model
for(i in c(1:nsams)){
  
  mean<-rep(betar[i],n)+X_ijk%*%beta_Er[i,]+Z_jk%*%beta_Mr[i,]+W_k%*%beta_Dr[i,]+rep(phir[i,],n_jk)
  v<-rep(sqrt(kappa2_jkr[i,]),n_jk)
  ll_ridge[i]<-sum(dnorm(x=y_ijk,mean=mean,sd=v,log = TRUE))
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
}

#Lasso model
for(i in c(1:nsams)){
  
  mean<-rep(betal[i],n)+X_ijk%*%beta_El[i,]+Z_jk%*%beta_Ml[i,]+W_k%*%beta_Dl[i,]+rep(phil[i,],n_jk)
  v<-rep(sqrt(kappa2_jkl[i,]),n_jk)
  ll_lasso[i]<-sum(dnorm(x=y_ijk,mean=mean,sd=v,log = TRUE))
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
}

ll_mean<-c(mean(ll_basis),mean(ll_ridge),mean(ll_lasso))

#Plot

par(mfrow=c(1,3))
plot(ll_basis,pch = 20, cex = 0.45,xlab = "",ylab = "Log likelihood",main="Basic",col="#70284a",col.main="#70284a")
abline(h=mean(ll_mean[1]),col="#70284a",lwd="2")
plot(ll_ridge,pch = 20, cex = 0.45,xlab = "Iteration",ylab = "",main="Ridge",col="#dc7176",col.main="#70284a")
abline(h=mean(ll_mean[2]),col="#dc7176",lwd="2")
plot(ll_lasso,pch = 20, cex = 0.45,xlab = "",ylab = "",main="Lasso",col="#f2a28a",col.main="#70284a")
abline(h=mean(ll_mean[3]),col="#f2a28a",lwd="2")

ll<-data.frame(
  Loglikelihood_basic=ll_basis,
  Loglikelihood_ridge=ll_ridge,
  Loglikelihood_lasso=ll_lasso
)

##=========DIC==========
###-------Basis---------

mean_est_bayes_basis<-mean(beta)+X_ijk%*%apply(beta_E,2,mean)+Z_jk%*%apply(beta_M,2,mean)+W_k%*%apply(beta_D,2,mean)+rep(apply(phi,2,mean),n_jk)
var_est_bayes_basis<-apply(kappa2_jk,2,mean)

lp_basis<-sum(dnorm(x=y_ijk,mean=mean_est_bayes_basis,sd=sqrt(var_est_bayes_basis)),log=TRUE)
p_DIC_basis<-2*(lp_basis-mean(ll_basis))

DIC_basis<--2*lp_basis+2*p_DIC_basis

###-------Ridge---------

mean_est_bayes_ridge<-mean(betar)+X_ijk%*%apply(beta_Er,2,stat_mode)+Z_jk%*%apply(beta_Mr,2,stat_mode)+W_k%*%apply(beta_Dr,2,stat_mode)+rep(apply(phir,2,mean),n_jk)
var_est_bayes_ridge<-apply(kappa2_jkr,2,mean)

lp_ridge<-sum(dnorm(x=y_ijk,mean=mean_est_bayes_ridge,sd=sqrt(var_est_bayes_ridge)),log=TRUE)
p_DIC_ridge<-2*(lp_ridge-mean(ll_ridge))

DIC_ridge<--2*lp_ridge+2*p_DIC_ridge

###------- Lasso ---------

mean_est_bayes_lasso<-mean(betal)+X_ijk%*%apply(beta_El,2,stat_mode)+Z_jk%*%apply(beta_Ml,2,stat_mode)+W_k%*%apply(beta_Dl,2,stat_mode)+rep(apply(phil,2,mean),n_jk)
var_est_bayes_lasso<-apply(kappa2_jkl,2,mean)

lp_lasso<-sum(dnorm(x=y_ijk,mean=mean_est_bayes_lasso,sd=sqrt(var_est_bayes_lasso)),log=TRUE)
p_DIC_lasso<-2*(lp_lasso-mean(ll_lasso))

DIC_lasso<--2*lp_lasso+2*p_DIC_lasso

##=====WAIC=========
index<-rep(seq_len(m),n_jk)

###-------Basis---------
tmp1<-0
tmp2<-0
lppd_basis<-0
pWAIC_basis<-0

for (j in c(1:n)) {
  mean<-beta+X_ijk[j,]%*%t(beta_E)+Z_jk[j,]%*%t(beta_M)+W_k[j,]%*%t(beta_D)+phi[,index[j]]
  v<-kappa2_jk[,index[j]]
  
  tmp1<-dnorm(x=y_ijk[j],mean=mean,sd=sqrt(v))
  lppd_basis<- lppd_basis + log(mean(tmp1))
  
  tmp2<-dnorm(x=y_ijk[j],mean=mean,sd=sqrt(v),log=TRUE)
  pWAIC_basis<-pWAIC_basis+ 2 * (log(mean(tmp1)) - mean(tmp2))
  
  # Advance
  if (j %% ceiling(n / 10) == 0) {
    cat(paste0("Iteración ", j, " de ", n, " (", round(100 * j / n), "%)\n"))
  }
}

WAIC_basis <- -2 * lppd_basis + 2 * pWAIC_basis


###-------Ridge---------
tmp1<-0
tmp2<-0
lppd_ridge<-0
pWAIC_ridge<-0


for (j in c(1:n)) {
  mean<-betar+X_ijk[j,]%*%t(beta_Er)+Z_jk[j,]%*%t(beta_Mr)+W_k[j,]%*%t(beta_Dr)+phir[,index[j]]
  v<-kappa2_jkr[,index[j]]
  
  tmp1<-dnorm(x=y_ijk[j],mean=mean,sd=sqrt(v))
  lppd_ridge<- lppd_ridge + log(mean(tmp1))
  
  tmp2<-dnorm(x=y_ijk[j],mean=mean,sd=sqrt(v),log=TRUE)
  pWAIC_ridge<-pWAIC_ridge+ 2 * (log(mean(tmp1)) - mean(tmp2))
  
  # Advance
  if (j %% ceiling(n / 10) == 0) {
    cat(paste0("Iteración ", j, " de ", n, " (", round(100 * j / n), "%)\n"))
  }
}

WAIC_ridge <- -2 * lppd_ridge + 2 * pWAIC_ridge

###-------Lasso---------
tmp1<-0
tmp2<-0
lppd_lasso<-0
pWAIC_lasso<-0


for (j in c(1:n)) {
  mean<-betal+X_ijk[j,]%*%t(beta_El)+Z_jk[j,]%*%t(beta_Ml)+W_k[j,]%*%t(beta_Dl)+phil[,index[j]]
  v<-kappa2_jkl[,index[j]]
  
  tmp1<-dnorm(x=y_ijk[j],mean=mean,sd=sqrt(v))
  lppd_lasso<- lppd_lasso + log(mean(tmp1))
  
  tmp2<-dnorm(x=y_ijk[j],mean=mean,sd=sqrt(v),log=TRUE)
  pWAIC_lasso<-pWAIC_lasso+ 2 * (log(mean(tmp1)) - mean(tmp2))
  
  # Advance
  if (j %% ceiling(n / 10) == 0) {
    cat(paste0("Iteración ", j, " de ", n, " (", round(100 * j / n), "%)\n"))
  }
}

WAIC_lasso <- -2 * lppd_lasso + 2 * pWAIC_lasso



DIC_WAIC<-data.frame(
  Model=c("Basic","Ridge","Lasso"),
  lp=c(lp_basis,lp_ridge,lp_lasso),
  pDIC=c(p_DIC_basis,p_DIC_ridge,p_DIC_lasso),
  DIC=c(DIC_basis,DIC_ridge,DIC_lasso),
  lppd=c(lppd_basis,lppd_ridge,lppd_lasso),
  pWAIC=c(pWAIC_basis,pWAIC_ridge,pWAIC_lasso),
  WAIC=c(WAIC_basis,WAIC_ridge,WAIC_lasso)
)

## =======Basis =======
### ===== Effective sample size summary for each parameter family====

neff_beta <- round(coda::effectiveSize(beta), 0)
neff_beta_E<-apply(beta_E, 2,effectiveSize)
neff_beta_M<-apply(beta_M, 2,effectiveSize)
neff_beta_D<-apply(beta_D, 2,effectiveSize)

neff_beta_summary<-summary(c(neff_beta,neff_beta_E,neff_beta_M,neff_beta_D))

neff_sigma2_beta<-round(coda::effectiveSize(sigma2_beta), 0)
neff_sigma2_E<-round(coda::effectiveSize(sigma2_E), 0)
neff_sigma2_M<-round(coda::effectiveSize(sigma2_M), 0)
neff_sigma2_D<-round(coda::effectiveSize(sigma2_D), 0)

neff_phi<-apply(phi, 2,effectiveSize)
neff_phi_summary<-summary(neff_phi)

neff_tau_phi<-round(coda::effectiveSize(tau_phi), 0)

neff_kappa2_jk<-apply(kappa2_jk, 2,effectiveSize)
neff_kappa2_jk_summary<-summary(neff_kappa2_jk)
neff_kappa2_k<-apply(kappa2_k, 2,effectiveSize)
neff_kappa2_k_summary<-summary(neff_kappa2_jk)

neff_alpha_kappa<-round(coda::effectiveSize(alpha_kappa), 0)
neff_beta_kappa<-round(coda::effectiveSize(beta_kappa), 0)

neff_summary_param_families_basis<-data.frame(
  parameter_family=c("Beta","Phi","Kappa2_jk","Kappa2_k"),
  min=c(neff_beta_summary[1],neff_phi_summary[1],neff_kappa2_jk_summary[1],neff_kappa2_k_summary[1]),
  Q_1=c(neff_beta_summary[2],neff_phi_summary[2],neff_kappa2_jk_summary[2],neff_kappa2_k_summary[2]),
  Median=c(neff_beta_summary[3],neff_phi_summary[3],neff_kappa2_jk_summary[3],neff_kappa2_k_summary[3]),
  Mean= c(neff_beta_summary[4],neff_phi_summary[4],neff_kappa2_jk_summary[4],neff_kappa2_k_summary[4]),
  Q_3=c(neff_beta_summary[5],neff_phi_summary[5],neff_kappa2_jk_summary[5],neff_kappa2_k_summary[5]),
  max=c(neff_beta_summary[6],neff_phi_summary[6],neff_kappa2_jk_summary[6],neff_kappa2_k_summary[6])
)

neff_summary_others<-data.frame(
  parameter=c("sigma2 beta","sigma2_E","sigma2_M","sigma2_D","lambda2_E","lambda2_M","lambda2_D","tau_phi","alpha_kappa","beta_kappa"),
  basis_model=c(neff_sigma2_beta,neff_sigma2_E,neff_sigma2_M,neff_sigma2_D,rep("-",3),neff_tau_phi,neff_alpha_kappa,neff_beta_kappa)
)

###======Monte Carlo Standard Error======
EMC_beta <- sd(beta)/sqrt(neff_beta)
EMC_beta_E<-apply(beta_E, 2, sd)/sqrt(neff_beta_E)
EMC_beta_M<-apply(beta_M, 2, sd)/sqrt(neff_beta_M)
EMC_beta_D<-apply(beta_D, 2, sd)/sqrt(neff_beta_D)

EMC_beta_summary<-summary(c(EMC_beta,EMC_beta_E,EMC_beta_M,EMC_beta_D))

EMC_sigma2_beta<-sd(sigma2_beta)/sqrt(neff_sigma2_beta)
EMC_sigma2_E<-sd(sigma2_E)/sqrt(neff_sigma2_E)
EMC_sigma2_M<-sd(sigma2_M)/sqrt(neff_sigma2_M)
EMC_sigma2_D<-sd(sigma2_D)/sqrt(neff_sigma2_D)

EMC_phi<-apply(phi, 2, sd)/sqrt(neff_phi)
EMC_phi_summary<-summary(EMC_phi)

EMC_tau_phi<-sd(tau_phi)/sqrt(neff_tau_phi)

EMC_kappa2_jk<-apply(kappa2_jk, 2,sd)/sqrt(neff_kappa2_jk)
EMC_kappa2_jk_summary<-summary(EMC_kappa2_jk)
EMC_kappa2_k<-apply(kappa2_k, 2,sd)/sqrt(neff_kappa2_k)
EMC_kappa2_k_summary<-summary(EMC_kappa2_k)

EMC_alpha_kappa<-sd(alpha_kappa)/sqrt(neff_alpha_kappa)
EMC_beta_kappa<-sd(beta_kappa)/sqrt(neff_beta_kappa)

EMC_summary_param_families_basis<-data.frame(
  parameter_family=c("Beta","Phi","Kappa2_jk","Kappa2_k"),
  min=c(EMC_beta_summary[1],EMC_phi_summary[1],EMC_kappa2_jk_summary[1],EMC_kappa2_k_summary[1]),
  Q_1=c(EMC_beta_summary[2],EMC_phi_summary[2],EMC_kappa2_jk_summary[2],EMC_kappa2_k_summary[2]),
  Median=c(EMC_beta_summary[3],EMC_phi_summary[3],EMC_kappa2_jk_summary[3],EMC_kappa2_k_summary[3]),
  Mean= c(EMC_beta_summary[4],EMC_phi_summary[4],EMC_kappa2_jk_summary[4],EMC_kappa2_k_summary[4]),
  Q_3=c(EMC_beta_summary[5],EMC_phi_summary[5],EMC_kappa2_jk_summary[5],EMC_kappa2_k_summary[5]),
  max=c(EMC_beta_summary[6],EMC_phi_summary[6],EMC_kappa2_jk_summary[6],EMC_kappa2_k_summary[6])
)

EMC_summary_others<-data.frame(
  parameter=c("sigma2 beta","sigma2_E","sigma2_M","sigma2_D","lambda2_E","lambda2_M","lambda2_D","tau_phi","alpha_kappa","beta_kappa"),
  basis_model=c(EMC_sigma2_beta,EMC_sigma2_E,EMC_sigma2_M,EMC_sigma2_D,rep("-",3),EMC_tau_phi,EMC_alpha_kappa,EMC_beta_kappa)
)

###========= Monte Carlo Variation Coefficient % ====

CVMC_beta <- 100*EMC_beta/abs(mean(beta))

CVMC_beta_E <- 100*EMC_beta_E/abs(apply(beta_E,2,mean))
CVMC_beta_M <- 100*EMC_beta_M/abs(apply(beta_M,2,mean))
CVMC_beta_D <- 100*EMC_beta_D/abs(apply(beta_D,2,mean))

CVMC_beta_summary<-summary(c(CVMC_beta,CVMC_beta_E,CVMC_beta_M,CVMC_beta_D))

CVMC_sigma2_beta<-100*EMC_sigma2_beta/abs(mean(sigma2_beta))
CVMC_sigma2_E<-100*EMC_sigma2_E/abs(mean(sigma2_E))
CVMC_sigma2_M<-100*EMC_sigma2_M/abs(mean(sigma2_M))
CVMC_sigma2_D<-100*EMC_sigma2_D/abs(mean(sigma2_D))



CVMC_phi <- 100*EMC_phi/abs(apply(phi,2,mean))
CVMC_phi_summary<-summary(CVMC_phi)
CVMC_tau_phi <- 100*EMC_tau_phi/abs(mean(tau_phi))

CVMC_kappa2_jk<- 100*EMC_kappa2_jk/abs(apply(kappa2_jk,2,mean))
CVMC_kappa2_jk_summary<-summary(CVMC_kappa2_jk)
CVMC_kappa2_k<- 100*EMC_kappa2_k/abs(apply(kappa2_k,2,mean))
CVMC_kappa2_k_summary<-summary(CVMC_kappa2_k)

CVMC_alpha_kappa <- 100*EMC_alpha_kappa/abs(mean(alpha_kappa))
CVMC_beta_kappa <- 100*EMC_beta_kappa/abs(mean(beta_kappa))


CVMC_summary_param_families_basis<-data.frame(
  parameter_family=c("Beta","Phi","Kappa2_jk","Kappa2_k"),
  min=c(CVMC_beta_summary[1],CVMC_phi_summary[1],CVMC_kappa2_jk_summary[1],CVMC_kappa2_k_summary[1]),
  Q_1=c(CVMC_beta_summary[2],CVMC_phi_summary[2],CVMC_kappa2_jk_summary[2],CVMC_kappa2_k_summary[2]),
  Median=c(CVMC_beta_summary[3],CVMC_phi_summary[3],CVMC_kappa2_jk_summary[3],CVMC_kappa2_k_summary[3]),
  Mean= c(CVMC_beta_summary[4],CVMC_phi_summary[4],CVMC_kappa2_jk_summary[4],CVMC_kappa2_k_summary[4]),
  Q_3=c(CVMC_beta_summary[5],CVMC_phi_summary[5],CVMC_kappa2_jk_summary[5],CVMC_kappa2_k_summary[5]),
  max=c(CVMC_beta_summary[6],CVMC_phi_summary[6],CVMC_kappa2_jk_summary[6],CVMC_kappa2_k_summary[6])
)

CVMC_summary_others<-data.frame(
  parameter=c("sigma2 beta","sigma2_E","sigma2_M","sigma2_D","lambda2_E","lambda2_M","lambda2_D","tau_phi","alpha_kappa","beta_kappa"),
  basis_model=c(CVMC_sigma2_beta,CVMC_sigma2_E,CVMC_sigma2_M,CVMC_sigma2_D,rep("-",3),CVMC_tau_phi,CVMC_alpha_kappa,CVMC_beta_kappa)
)



##===== Ridge====


### ===== Effective sample size summary for each parameter family====
neff_beta <- round(coda::effectiveSize(betar), 0)
neff_beta_E<-apply(beta_Er, 2,effectiveSize)
neff_beta_M<-apply(beta_Mr, 2,effectiveSize)
neff_beta_D<-apply(beta_Dr, 2,effectiveSize)

neff_beta_summary<-summary(c(neff_beta,neff_beta_E,neff_beta_M,neff_beta_D))

neff_sigma2_beta<-round(coda::effectiveSize(sigma2_betar), 0)
neff_lambda2_E<-round(coda::effectiveSize(lambda2_Er), 0)
neff_lambda2_M<-round(coda::effectiveSize(lambda2_Mr), 0)
neff_lambda2_D<-round(coda::effectiveSize(lambda2_Dr), 0)

neff_phi<-apply(phir, 2,effectiveSize)
neff_phi_summary<-summary(neff_phi)

neff_tau_phi<-round(coda::effectiveSize(tau_phir), 0)

neff_kappa2_jk<-apply(kappa2_jkr, 2,effectiveSize)
neff_kappa2_jk_summary<-summary(neff_kappa2_jk)
neff_kappa2_k<-apply(kappa2_kr, 2,effectiveSize)
neff_kappa2_k_summary<-summary(neff_kappa2_jk)

neff_alpha_kappa<-round(coda::effectiveSize(alpha_kappar), 0)
neff_beta_kappa<-round(coda::effectiveSize(beta_kappar), 0)

neff_summary_param_families_ridge<-data.frame(
  parameter_family=c("Beta","Phi","Kappa2_jk","Kappa2_k"),
  min=c(neff_beta_summary[1],neff_phi_summary[1],neff_kappa2_jk_summary[1],neff_kappa2_k_summary[1]),
  Q_1=c(neff_beta_summary[2],neff_phi_summary[2],neff_kappa2_jk_summary[2],neff_kappa2_k_summary[2]),
  Median=c(neff_beta_summary[3],neff_phi_summary[3],neff_kappa2_jk_summary[3],neff_kappa2_k_summary[3]),
  Mean= c(neff_beta_summary[4],neff_phi_summary[4],neff_kappa2_jk_summary[4],neff_kappa2_k_summary[4]),
  Q_3=c(neff_beta_summary[5],neff_phi_summary[5],neff_kappa2_jk_summary[5],neff_kappa2_k_summary[5]),
  max=c(neff_beta_summary[6],neff_phi_summary[6],neff_kappa2_jk_summary[6],neff_kappa2_k_summary[6])
)

neff_summary_others$ridge_model<-c(neff_sigma2_beta,rep("-",3),neff_lambda2_E,neff_lambda2_M,neff_lambda2_D,neff_tau_phi,neff_alpha_kappa,neff_beta_kappa)


###======Monte Carlo Standard Error======
EMC_beta <- sd(betar)/sqrt(neff_beta)
EMC_beta_E<-apply(beta_Er, 2, sd)/sqrt(neff_beta_E)
EMC_beta_M<-apply(beta_Mr, 2, sd)/sqrt(neff_beta_M)
EMC_beta_D<-apply(beta_Dr, 2, sd)/sqrt(neff_beta_D)

EMC_beta_summary<-summary(c(EMC_beta,EMC_beta_E,EMC_beta_M,EMC_beta_D))

EMC_sigma2_beta <- sd(sigma2_betar)/sqrt(neff_sigma2_beta)
EMC_lambda2_E <- sd(lambda2_Er)/sqrt(neff_lambda2_E)
EMC_lambda2_M <- sd(lambda2_Mr)/sqrt(neff_lambda2_M)
EMC_lambda2_D <- sd(lambda2_Dr)/sqrt(neff_lambda2_D)

EMC_phi<-apply(phir, 2, sd)/sqrt(neff_phi)
EMC_phi_summary<-summary(EMC_phi)

EMC_tau_phi<-sd(tau_phir)/sqrt(neff_tau_phi)

EMC_kappa2_jk<-apply(kappa2_jkr, 2,sd)/sqrt(neff_kappa2_jk)
EMC_kappa2_jk_summary<-summary(EMC_kappa2_jk)
EMC_kappa2_k<-apply(kappa2_kr, 2,sd)/sqrt(neff_kappa2_k)
EMC_kappa2_k_summary<-summary(EMC_kappa2_k)

EMC_alpha_kappa<-sd(alpha_kappar)/sqrt(neff_alpha_kappa)
EMC_beta_kappa<-sd(beta_kappar)/sqrt(neff_beta_kappa)

EMC_summary_param_families_ridge<-data.frame(
  parameter_family=c("Beta","Phi","Kappa2_jk","Kappa2_k"),
  min=c(EMC_beta_summary[1],EMC_phi_summary[1],EMC_kappa2_jk_summary[1],EMC_kappa2_k_summary[1]),
  Q_1=c(EMC_beta_summary[2],EMC_phi_summary[2],EMC_kappa2_jk_summary[2],EMC_kappa2_k_summary[2]),
  Median=c(EMC_beta_summary[3],EMC_phi_summary[3],EMC_kappa2_jk_summary[3],EMC_kappa2_k_summary[3]),
  Mean= c(EMC_beta_summary[4],EMC_phi_summary[4],EMC_kappa2_jk_summary[4],EMC_kappa2_k_summary[4]),
  Q_3=c(EMC_beta_summary[5],EMC_phi_summary[5],EMC_kappa2_jk_summary[5],EMC_kappa2_k_summary[5]),
  max=c(EMC_beta_summary[6],EMC_phi_summary[6],EMC_kappa2_jk_summary[6],EMC_kappa2_k_summary[6])
)

EMC_summary_others$ridge_model<-c(EMC_sigma2_beta,rep("-",3),EMC_lambda2_E,EMC_lambda2_M,EMC_lambda2_D,EMC_tau_phi,EMC_alpha_kappa,EMC_beta_kappa)


###========= Monte Carlo Variation Coefficient % ====

CVMC_beta <- 100*EMC_beta/abs(mean(betar))

CVMC_beta_E <- 100*EMC_beta_E/abs(apply(beta_Er,2,mean))
CVMC_beta_M <- 100*EMC_beta_M/abs(apply(beta_Mr,2,mean))
CVMC_beta_D <- 100*EMC_beta_D/abs(apply(beta_Dr,2,mean))

CVMC_beta_summary<-summary(c(CVMC_beta,CVMC_beta_E,CVMC_beta_M,CVMC_beta_D))

CVMC_sigma2_beta <- 100*EMC_sigma2_beta/abs(mean(sigma2_betar))
CVMC_lambda2_E <- 100*EMC_lambda2_E/abs(mean(lambda2_Er))
CVMC_lambda2_M <- 100*EMC_lambda2_M/abs(mean(lambda2_Mr))
CVMC_lambda2_D <- 100*EMC_lambda2_D/abs(mean(lambda2_Dr))


CVMC_phi <- 100*EMC_phi/abs(apply(phir,2,mean))
CVMC_phi_summary<-summary(CVMC_phi)
CVMC_tau_phi <- 100*EMC_tau_phi/abs(mean(tau_phir))

CVMC_kappa2_jk<- 100*EMC_kappa2_jk/abs(apply(kappa2_jkr,2,mean))
CVMC_kappa2_jk_summary<-summary(CVMC_kappa2_jk)
CVMC_kappa2_k<- 100*EMC_kappa2_k/abs(apply(kappa2_kr,2,mean))
CVMC_kappa2_k_summary<-summary(CVMC_kappa2_k)

CVMC_alpha_kappa <- 100*EMC_alpha_kappa/abs(mean(alpha_kappar))
CVMC_beta_kappa <- 100*EMC_beta_kappa/abs(mean(beta_kappar))


CVMC_summary_param_families_ridge<-data.frame(
  parameter_family=c("Beta","Phi","Kappa2_jk","Kappa2_k"),
  min=c(CVMC_beta_summary[1],CVMC_phi_summary[1],CVMC_kappa2_jk_summary[1],CVMC_kappa2_k_summary[1]),
  Q_1=c(CVMC_beta_summary[2],CVMC_phi_summary[2],CVMC_kappa2_jk_summary[2],CVMC_kappa2_k_summary[2]),
  Median=c(CVMC_beta_summary[3],CVMC_phi_summary[3],CVMC_kappa2_jk_summary[3],CVMC_kappa2_k_summary[3]),
  Mean= c(CVMC_beta_summary[4],CVMC_phi_summary[4],CVMC_kappa2_jk_summary[4],CVMC_kappa2_k_summary[4]),
  Q_3=c(CVMC_beta_summary[5],CVMC_phi_summary[5],CVMC_kappa2_jk_summary[5],CVMC_kappa2_k_summary[5]),
  max=c(CVMC_beta_summary[6],CVMC_phi_summary[6],CVMC_kappa2_jk_summary[6],CVMC_kappa2_k_summary[6])
)

CVMC_summary_others$ridge_model<-c(CVMC_sigma2_beta,rep("-",3),CVMC_lambda2_E,CVMC_lambda2_M,CVMC_lambda2_D,CVMC_tau_phi,CVMC_alpha_kappa,CVMC_beta_kappa)

##====Lasso======


### ===== Effective sample size summary for each parameter family====
neff_beta <- round(coda::effectiveSize(betal), 0)
neff_beta_E<-apply(beta_El, 2,effectiveSize)
neff_beta_M<-apply(beta_Ml, 2,effectiveSize)
neff_beta_D<-apply(beta_Dl, 2,effectiveSize)

neff_beta_summary<-summary(c(neff_beta,neff_beta_E,neff_beta_M,neff_beta_D))

neff_sigma2_beta<-round(coda::effectiveSize(sigma2_betal), 0)

neff_tau2_E<-apply(tau2_El, 2,effectiveSize)
neff_tau2_M<-apply(tau2_Ml, 2,effectiveSize)
neff_tau2_D<-apply(tau2_Dl, 2,effectiveSize)

neff_tau_summary<-summary(c(neff_tau2_E,neff_tau2_M,neff_tau2_D))

neff_lambda2_E<-round(coda::effectiveSize(lambda2_El), 0)
neff_lambda2_M<-round(coda::effectiveSize(lambda2_Ml), 0)
neff_lambda2_D<-round(coda::effectiveSize(lambda2_Dl), 0)

neff_phi<-apply(phil, 2,effectiveSize)
neff_phi_summary<-summary(neff_phi)

neff_tau_phi<-round(coda::effectiveSize(tau_phil), 0)

neff_kappa2_jk<-apply(kappa2_jkl, 2,effectiveSize)
neff_kappa2_jk_summary<-summary(neff_kappa2_jk)
neff_kappa2_k<-apply(kappa2_kl, 2,effectiveSize)
neff_kappa2_k_summary<-summary(neff_kappa2_jk)

neff_alpha_kappa<-round(coda::effectiveSize(alpha_kappal), 0)
neff_beta_kappa<-round(coda::effectiveSize(beta_kappal), 0)

neff_summary_param_families_lasso<-data.frame(
  parameter_family=c("Beta","Tau","Phi","Kappa2_jk","Kappa2_k"),
  min=c(neff_beta_summary[1],neff_tau_summary[1],neff_phi_summary[1],neff_kappa2_jk_summary[1],neff_kappa2_k_summary[1]),
  Q_1=c(neff_beta_summary[2],neff_tau_summary[2],neff_phi_summary[2],neff_kappa2_jk_summary[2],neff_kappa2_k_summary[2]),
  Median=c(neff_beta_summary[3],neff_tau_summary[3],neff_phi_summary[3],neff_kappa2_jk_summary[3],neff_kappa2_k_summary[3]),
  Mean= c(neff_beta_summary[4],neff_tau_summary[4],neff_phi_summary[4],neff_kappa2_jk_summary[4],neff_kappa2_k_summary[4]),
  Q_3=c(neff_beta_summary[5],neff_tau_summary[5],neff_phi_summary[5],neff_kappa2_jk_summary[5],neff_kappa2_k_summary[5]),
  max=c(neff_beta_summary[6],neff_tau_summary[6],neff_phi_summary[6],neff_kappa2_jk_summary[6],neff_kappa2_k_summary[6])
)


neff_summary_others$lasso_model<-c(neff_sigma2_beta,rep("-",3),neff_lambda2_E,neff_lambda2_M,neff_lambda2_D,neff_tau_phi,neff_alpha_kappa,neff_beta_kappa)


###======Monte Carlo Standard Error======
EMC_beta <- sd(betal)/sqrt(neff_beta)
EMC_beta_E<-apply(beta_El, 2, sd)/sqrt(neff_beta_E)
EMC_beta_M<-apply(beta_Ml, 2, sd)/sqrt(neff_beta_M)
EMC_beta_D<-apply(beta_Dl, 2, sd)/sqrt(neff_beta_D)

EMC_beta_summary<-summary(c(EMC_beta,EMC_beta_E,EMC_beta_M,EMC_beta_D))

EMC_tau2_E<-apply(tau2_El, 2, sd)/sqrt(neff_tau2_E)
EMC_tau2_M<-apply(tau2_Ml, 2, sd)/sqrt(neff_tau2_M)
EMC_tau2_D<-apply(tau2_Dl, 2, sd)/sqrt(neff_tau2_D)

EMC_tau_summary<-summary(c(EMC_tau2_E,EMC_tau2_M,EMC_tau2_D))

EMC_sigma2_beta<-sd(sigma2_betal)/sqrt(neff_sigma2_beta)
EMC_lambda2_E<-sd(lambda2_El)/sqrt(neff_lambda2_E)
EMC_lambda2_M<-sd(lambda2_Ml)/sqrt(neff_lambda2_M)
EMC_lambda2_D<-sd(lambda2_Dl)/sqrt(neff_lambda2_D)

EMC_phi<-apply(phil, 2, sd)/sqrt(neff_phi)
EMC_phi_summary<-summary(EMC_phi)

EMC_tau_phi<-sd(tau_phil)/sqrt(neff_tau_phi)

EMC_kappa2_jk<-apply(kappa2_jkl, 2,sd)/sqrt(neff_kappa2_jk)
EMC_kappa2_jk_summary<-summary(EMC_kappa2_jk)
EMC_kappa2_k<-apply(kappa2_kl, 2,sd)/sqrt(neff_kappa2_k)
EMC_kappa2_k_summary<-summary(EMC_kappa2_k)

EMC_alpha_kappa<-sd(alpha_kappal)/sqrt(neff_alpha_kappa)
EMC_beta_kappa<-sd(beta_kappal)/sqrt(neff_beta_kappa)

EMC_summary_param_families_lasso<-data.frame(
  parameter_family=c("Beta","Tau","Phi","Kappa2_jk","Kappa2_k"),
  min=c(EMC_beta_summary[1],EMC_tau_summary[1],EMC_phi_summary[1],EMC_kappa2_jk_summary[1],EMC_kappa2_k_summary[1]),
  Q_1=c(EMC_beta_summary[2],EMC_tau_summary[2],EMC_phi_summary[2],EMC_kappa2_jk_summary[2],EMC_kappa2_k_summary[2]),
  Median=c(EMC_beta_summary[3],EMC_tau_summary[3],EMC_phi_summary[3],EMC_kappa2_jk_summary[3],EMC_kappa2_k_summary[3]),
  Mean= c(EMC_beta_summary[4],EMC_tau_summary[4],EMC_phi_summary[4],EMC_kappa2_jk_summary[4],EMC_kappa2_k_summary[4]),
  Q_3=c(EMC_beta_summary[5],EMC_tau_summary[5],EMC_phi_summary[5],EMC_kappa2_jk_summary[5],EMC_kappa2_k_summary[5]),
  max=c(EMC_beta_summary[6],EMC_tau_summary[6],EMC_phi_summary[6],EMC_kappa2_jk_summary[6],EMC_kappa2_k_summary[6])
)

EMC_summary_others$lasso_model<-c(EMC_sigma2_beta,rep("-",3),EMC_lambda2_E,EMC_lambda2_M,EMC_lambda2_D,EMC_tau_phi,EMC_alpha_kappa,EMC_beta_kappa)


###========= Monte Carlo Variation Coefficient % ====

CVMC_beta <- 100*EMC_beta/abs(mean(betal))

CVMC_beta_E <- 100*EMC_beta_E/abs(apply(beta_El,2,mean))
CVMC_beta_M <- 100*EMC_beta_M/abs(apply(beta_Ml,2,mean))
CVMC_beta_D <- 100*EMC_beta_D/abs(apply(beta_Dl,2,mean))

CVMC_beta_summary<-summary(c(CVMC_beta,CVMC_beta_E,CVMC_beta_M,CVMC_beta_D))

CVMC_tau2_E <- 100*EMC_tau2_E/abs(apply(tau2_El,2,mean))
CVMC_tau2_M <- 100*EMC_tau2_M/abs(apply(tau2_Ml,2,mean))
CVMC_tau2_D <- 100*EMC_tau2_D/abs(apply(tau2_Dl,2,mean))

CVMC_tau2_summary<-summary(c(CVMC_tau2_E,CVMC_tau2_M,CVMC_tau2_D))

CVMC_sigma2_beta<-100*EMC_sigma2_beta/abs(mean(sigma2_betal))

CVMC_lambda2_E<-100*EMC_lambda2_E/abs(mean(lambda2_El))
CVMC_lambda2_M<-100*EMC_lambda2_M/abs(mean(lambda2_Ml))
CVMC_lambda2_D<-100*EMC_lambda2_D/abs(mean(lambda2_Dl))

CVMC_phi <- 100*EMC_phi/abs(apply(phil,2,mean))
CVMC_phi_summary<-summary(CVMC_phi)
CVMC_tau_phi <- 100*EMC_tau_phi/abs(mean(tau_phil))

CVMC_kappa2_jk<- 100*EMC_kappa2_jk/abs(apply(kappa2_jkl,2,mean))
CVMC_kappa2_jk_summary<-summary(CVMC_kappa2_jk)
CVMC_kappa2_k<- 100*EMC_kappa2_k/abs(apply(kappa2_kl,2,mean))
CVMC_kappa2_k_summary<-summary(CVMC_kappa2_k)

CVMC_alpha_kappa <- 100*EMC_alpha_kappa/abs(mean(alpha_kappal))
CVMC_beta_kappa <- 100*EMC_beta_kappa/abs(mean(beta_kappal))


CVMC_summary_param_families_lasso<-data.frame(
  parameter_family=c("Beta","Tau","Phi","Kappa2_jk","Kappa2_k"),
  min=c(CVMC_beta_summary[1],CVMC_tau2_summary[1],CVMC_phi_summary[1],CVMC_kappa2_jk_summary[1],CVMC_kappa2_k_summary[1]),
  Q_1=c(CVMC_beta_summary[2],CVMC_tau2_summary[2],CVMC_phi_summary[2],CVMC_kappa2_jk_summary[2],CVMC_kappa2_k_summary[2]),
  Median=c(CVMC_beta_summary[3],CVMC_tau2_summary[3],CVMC_phi_summary[3],CVMC_kappa2_jk_summary[3],CVMC_kappa2_k_summary[3]),
  Mean= c(CVMC_beta_summary[4],CVMC_tau2_summary[4],CVMC_phi_summary[4],CVMC_kappa2_jk_summary[4],CVMC_kappa2_k_summary[4]),
  Q_3=c(CVMC_beta_summary[5],CVMC_tau2_summary[5],CVMC_phi_summary[5],CVMC_kappa2_jk_summary[5],CVMC_kappa2_k_summary[5]),
  max=c(CVMC_beta_summary[6],CVMC_tau2_summary[6],CVMC_phi_summary[6],CVMC_kappa2_jk_summary[6],CVMC_kappa2_k_summary[6])
)

CVMC_summary_others$lasso_model<-c(CVMC_sigma2_beta,rep("-",3),CVMC_lambda2_E,CVMC_lambda2_M,CVMC_lambda2_D,CVMC_tau_phi,CVMC_alpha_kappa,CVMC_beta_kappa)


information_criteria<-list(Log_likelihoods=ll,NEFF_basic=neff_summary_param_families_basis,NEFF_ridge=neff_summary_param_families_ridge,
                           NEFF_Lasso=neff_summary_param_families_lasso,EMC_basic=EMC_summary_param_families_basis,EMC_ridge=EMC_summary_param_families_ridge,
                           EMC_Lasso=EMC_summary_param_families_lasso,CVMC_basic=CVMC_summary_param_families_basis,CVMC_ridge=CVMC_summary_param_families_ridge,
                           CVMC_Lasso=CVMC_summary_param_families_lasso,neff_others=neff_summary_others,EMC_others=EMC_summary_others,CVMC_others=CVMC_summary_others,
                           DIC_WAIC_table=DIC_WAIC)

save(information_criteria,file="Informacion_criteria_real_data.RData")


# ============ Covariates effects posterior mean =============


df_basic<- data.frame(
  variable = c("beta",colnames(beta_E),colnames(beta_M),colnames(beta_D)),
  media = c(mean(beta),apply(beta_E,2,stat_mode),apply(beta_M,2,stat_mode),apply(beta_D,2,stat_mode)),
  lim_inf=c(bayesian_summary(beta)[2],apply(beta_E, 2, bayesian_summary)[2,],apply(beta_M, 2, bayesian_summary)[2,],apply(beta_D,2, bayesian_summary)[2,]),
  lim_sup=c(bayesian_summary(beta)[3],apply(beta_E, 2, bayesian_summary)[3,],apply(beta_M, 2, bayesian_summary)[3,],apply(beta_D, 2, bayesian_summary)[3,]),
  Nivel=c("Intercepto", rep("Student",ncol(beta_E)),rep("Municipal",ncol(beta_M)),rep("Departamental",ncol(beta_D))),
  Model=rep("Basic",34)
)
df_basic$Signif<-ifelse(0>=df_basic$lim_inf & 0<=df_basic$lim_sup,0,1)

df_basic$variable <- recode(df_basic$variable,
                            "beta"                        = "Intercepto",
                            "fami_educacionmadre_modif"   = "Mother's Education",
                            "computador"                  = "Computer access",
                            "internet"                    = "Internet access",
                            "etnia"                       = "Ethinicity",
                            "libros_11_25"                = "Books (11–25)",
                            "libros_26_100"               = "Books (26–100)",
                            "libros_mas100"               = "Books (>100)",
                            "estrato_1"                   = "Socio-Econ level 1",
                            "estrato_2"                   = "Socio-Econ level 2",
                            "estrato_3"                   = "Socio-Econ level 3",
                            "estrato_4"                   = "Socio-Econ level 4",
                            "estrato_5"                   = "Socio-Econ level 5",
                            "estrato_6"                   = "Socio-Econ level 6",
                            "Genero_mujer"                = "Gender (women)",
                            "Calendario_A"                = "Academic calendar A",
                            "Calendario_B"                = "Academic calendar B",
                            "cole_privado"                = "Private school",
                            "trabaja_menos_de_10_horas"   = "Work (<10h)",
                            "trabaja__11_a_20_horas"      = "Work (11–20h)",
                            "trabaja__21_a_30_horas"      = "Work (21–30h)",
                            "trabaja_mas_de_30_horas"     = "Work (>30h)",
                            "docenttotal_alumtotal"       = "Teachers-Student ratio",
                            "RISK_VICTIM_2022"            = " Victimization Risk",
                            "Homi_x_100k_habitantes"      = "Homicides (x100k hab.)",
                            "porc_alumn_col_publico"      = "% Students in public school",
                            "terrorismot"                 = "Terrorism index",
                            "hurto"                       = "Theft rate",
                            "secuestros"                  = "Kidnapping",
                            "discapital"                  = "Distance to capital",
                            "PIB_percapita_DPTO"          = "GDP per capita",
                            "proporcio_pob_rural"         = "% Rural population",
                            "X..municipios.con.riesgo"    = "% Municipalities at risk",
                            "Homicidios_ponderado_x_100k" = "Weighted homicides (x100k)"
)


df_ridge<- data.frame(
  variable = c("beta",colnames(beta_Er),colnames(beta_Mr),colnames(beta_Dr)),
  media = c(mean(betar),apply(beta_Er,2,stat_mode),apply(beta_Mr,2,stat_mode),apply(beta_Dr,2,stat_mode)),
  lim_inf=c(bayesian_summary(betar)[2],apply(beta_Er, 2, bayesian_summary)[2,],apply(beta_Mr, 2, bayesian_summary)[2,],apply(beta_Dr,2, bayesian_summary)[2,]),
  lim_sup=c(bayesian_summary(betar)[3],apply(beta_Er, 2, bayesian_summary)[3,],apply(beta_Mr, 2, bayesian_summary)[3,],apply(beta_Dr, 2, bayesian_summary)[3,]),
  Nivel=c("Intercepto", rep("Student",ncol(beta_Er)),rep("Municipal",ncol(beta_Mr)),rep("Departamental",ncol(beta_Dr))),
  Model=rep("Ridge",34)
)
df_ridge$Signif<-ifelse(0>=df_ridge$lim_inf & 0<=df_ridge$lim_sup,0,1)

df_ridge$variable <- recode(df_ridge$variable,
                            "beta"                        = "Intercepto",
                            "fami_educacionmadre_modif"   = "Mother's Education",
                            "computador"                  = "Computer access",
                            "internet"                    = "Internet access",
                            "etnia"                       = "Ethinicity",
                            "libros_11_25"                = "Books (11–25)",
                            "libros_26_100"               = "Books (26–100)",
                            "libros_mas100"               = "Books (>100)",
                            "estrato_1"                   = "Socio-Econ level 1",
                            "estrato_2"                   = "Socio-Econ level 2",
                            "estrato_3"                   = "Socio-Econ level 3",
                            "estrato_4"                   = "Socio-Econ level 4",
                            "estrato_5"                   = "Socio-Econ level 5",
                            "estrato_6"                   = "Socio-Econ level 6",
                            "Genero_mujer"                = "Gender (women)",
                            "Calendario_A"                = "Academic calendar A",
                            "Calendario_B"                = "Academic calendar B",
                            "cole_privado"                = "Private school",
                            "trabaja_menos_de_10_horas"   = "Work (<10h)",
                            "trabaja__11_a_20_horas"      = "Work (11–20h)",
                            "trabaja__21_a_30_horas"      = "Work (21–30h)",
                            "trabaja_mas_de_30_horas"     = "Work (>30h)",
                            "docenttotal_alumtotal"       = "Teachers-Student ratio",
                            "RISK_VICTIM_2022"            = " Victimization Risk",
                            "Homi_x_100k_habitantes"      = "Homicides (x100k hab.)",
                            "porc_alumn_col_publico"      = "% Students in public school",
                            "terrorismot"                 = "Terrorism index",
                            "hurto"                       = "Theft rate",
                            "secuestros"                  = "Kidnapping",
                            "discapital"                  = "Distance to capital",
                            "PIB_percapita_DPTO"          = "GDP per capita",
                            "proporcio_pob_rural"         = "% Rural population",
                            "X..municipios.con.riesgo"    = "% Municipalities at risk",
                            "Homicidios_ponderado_x_100k" = "Weighted homicides (x100k)"
)



df_lasso<- data.frame(
  variable = c("beta",colnames(beta_El),colnames(beta_Ml),colnames(beta_Dl)),
  media = c(mean(betal),apply(beta_El,2,stat_mode),apply(beta_Ml,2,stat_mode),apply(beta_Dl,2,stat_mode)),
  lim_inf=c(bayesian_summary(betal)[2],apply(beta_El, 2, bayesian_summary)[2,],apply(beta_Ml, 2, bayesian_summary)[2,],apply(beta_Dl,2, bayesian_summary)[2,]),
  lim_sup=c(bayesian_summary(betal)[3],apply(beta_El, 2, bayesian_summary)[3,],apply(beta_Ml, 2, bayesian_summary)[3,],apply(beta_Dl, 2, bayesian_summary)[3,]),
  Nivel=c("Intercepto", rep("Student",ncol(beta_El)),rep("Municipal",ncol(beta_Ml)),rep("Departamental",ncol(beta_Dl))),
  Model=rep("Lasso",34)
)

df_lasso$Signif<-ifelse(0>=df_lasso$lim_inf & 0<=df_lasso$lim_sup,0,1)

df_lasso$variable <- recode(df_lasso$variable,
                       "beta"                        = "Intercepto",
                       "fami_educacionmadre_modif"   = "Mother's Education",
                       "computador"                  = "Computer access",
                       "internet"                    = "Internet access",
                       "etnia"                       = "Ethinicity",
                       "libros_11_25"                = "Books (11–25)",
                       "libros_26_100"               = "Books (26–100)",
                       "libros_mas100"               = "Books (>100)",
                       "estrato_1"                   = "Socio-Econ level 1",
                       "estrato_2"                   = "Socio-Econ level 2",
                       "estrato_3"                   = "Socio-Econ level 3",
                       "estrato_4"                   = "Socio-Econ level 4",
                       "estrato_5"                   = "Socio-Econ level 5",
                       "estrato_6"                   = "Socio-Econ level 6",
                       "Genero_mujer"                = "Gender (women)",
                       "Calendario_A"                = "Academic calendar A",
                       "Calendario_B"                = "Academic calendar B",
                       "cole_privado"                = "Private school",
                       "trabaja_menos_de_10_horas"   = "Work (<10h)",
                       "trabaja__11_a_20_horas"      = "Work (11–20h)",
                       "trabaja__21_a_30_horas"      = "Work (21–30h)",
                       "trabaja_mas_de_30_horas"     = "Work (>30h)",
                       "docenttotal_alumtotal"       = "Teachers-Student ratio",
                       "RISK_VICTIM_2022"            = " Victimization Risk",
                       "Homi_x_100k_habitantes"      = "Homicides (x100k hab.)",
                       "porc_alumn_col_publico"      = "% Students in public school",
                       "terrorismot"                 = "Terrorism index",
                       "hurto"                       = "Theft rate",
                       "secuestros"                  = "Kidnapping",
                       "discapital"                  = "Distance to capital",
                       "PIB_percapita_DPTO"          = "GDP per capita",
                       "proporcio_pob_rural"         = "% Rural population",
                       "X..municipios.con.riesgo"    = "% Municipalities at risk",
                       "Homicidios_ponderado_x_100k" = "Weighted homicides (x100k)"
)


df_basic<-df_basic%>%
  arrange(desc(media))
df_ridge<-df_ridge%>%
  arrange(desc(media))
df_lasso<-df_lasso%>%
  arrange(desc(media))

df<-rbind(df_basic,df_ridge,df_lasso)
df<-df%>%
  arrange(desc(media))

#A plot in which is possible to visualice posterior means and the crebilibity intervals

ggplot(df[-(1:3),], aes(y = reorder(variable, media), x = media, color = Nivel)) +
  geom_point(size = 1.5) +
  geom_errorbarh(aes(xmin = lim_inf, xmax = lim_sup), height = 0.3)  +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "Estimation", y = "Covariate", color = "Level", title ="Regression Coefficient Estimation Analysis" ,subtitle = "Puntual estimation and Credibility Intervals") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right") +
  scale_color_manual(values =  c("#dc7176", "#70284a", "#f2a28a"))+
  theme_minimal(base_size = 13,base_family = "Helvetica") +
  facet_wrap(~Model,ncol=3,nrow=1,strip.position="top",dir="h",scales="fixed") +
  theme(
    plot.title = element_text(color = "#7C3B5E",face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40",size = 12, hjust = 0.5),
    axis.title = element_text(color = "black",face = "plain"),
    legend.position = "bottom"
  ) 



#============ PPP ==========


index<-data_ready$id_dep_mun[,2]
index_dep<-data_ready$id_dep_mun[,1]
# Statistic of the data

global_statistic<-data.frame(
  mean=mean(y_ijk),
  CV=CV(y_ijk),
  median=quantile(y_ijk,c(0.5))
)



municipal_statistic<-data.frame(
  mean=tapply(y_ijk,index,mean),
  var=tapply(y_ijk,index,var)
)

departamental_statistic<-data.frame(
  mean=tapply(y_ijk,index_dep,mean),
  var=tapply(y_ijk,index_dep,var)
)

y_ppp<-0

general_statistics_basic<-matrix(NA,nrow = nsams,ncol=3)
municipal_mean_basic<-matrix(NA,nrow = nsams,ncol=length(n_jk))
variance_mean_basic<-matrix(NA,nrow = nsams,ncol=length(n_jk))
departamental_mean_basic<-matrix(NA,nrow = nsams,ncol=d)
departamental_variance_mean_basic<-matrix(NA,nrow = nsams,ncol=d)


general_statistics_ridge<-matrix(NA,nrow = nsams,ncol=3)
municipal_mean_ridge<-matrix(NA,nrow = nsams,ncol=length(n_jk))
variance_mean_ridge<-matrix(NA,nrow = nsams,ncol=length(n_jk))
departamental_mean_ridge<-matrix(NA,nrow = nsams,ncol=d)
departamental_variance_mean_ridge<-matrix(NA,nrow = nsams,ncol=d)

general_statistics_lasso<-matrix(NA,nrow = nsams,ncol=3)
municipal_mean_lasso<-matrix(NA,nrow = nsams,ncol=length(n_jk))
variance_mean_lasso<-matrix(NA,nrow = nsams,ncol=length(n_jk))
departamental_mean_lasso<-matrix(NA,nrow = nsams,ncol=d)
departamental_variance_mean_lasso<-matrix(NA,nrow = nsams,ncol=d)

set.seed(2008)

for(i in c(1:nsams)){
  
  
  #Generation of the data with basic
  mean<-rep(beta[i],n)+X_ijk%*%beta_E[i,]+Z_jk%*%beta_M[i,]+W_k%*%beta_D[i,]+rep(phi[i,],n_jk)
  v<-rep(sqrt(kappa2_jk[i,]),n_jk)
  
  y_ppp<-rnorm(n,mean=mean,sd=v)
  
  #Calculus of the general statistic
  general_statistics_basic[i,]<-c(mean(y_ppp),sd(y_ppp)/abs(mean(y_ppp)),quantile(y_ppp,c(0.5)))
  
  municipal_mean_basic[i,]<-tapply(y_ppp,index,mean)
  variance_mean_basic[i,]<-tapply(y_ppp,index,var)
  departamental_mean_basic[i,]<-tapply(y_ppp,index_dep,mean)
  departamental_variance_mean_basic[i,]<-tapply(y_ppp,index_dep,var)
  
  #------------------------------------------------------------ #
  
  #Generation of the data lasso
  
  mean<-rep(betal[i],n)+X_ijk%*%beta_El[i,]+Z_jk%*%beta_Ml[i,]+W_k%*%beta_Dl[i,]+rep(phil[i,],n_jk)
  v<-rep(sqrt(kappa2_jkl[i,]),n_jk)
  
  y_ppp<-rnorm(n,mean=mean,sd=v)
  
  #Calculus of the general statistic
  general_statistics_lasso[i,]<-c(mean(y_ppp),sd(y_ppp)/abs(mean(y_ppp)),quantile(y_ppp,c(0.5)))
  
  municipal_mean_lasso[i,]<-tapply(y_ppp,index,mean)
  variance_mean_lasso[i,]<-tapply(y_ppp,index,var)
  departamental_mean_lasso[i,]<-tapply(y_ppp,index_dep,mean)
  departamental_variance_mean_lasso[i,]<-tapply(y_ppp,index_dep,var)
  
  #------------------------------------------------------------ #
  
  #Generation of the data ridge
  mean<-rep(betar[i],n)+X_ijk%*%beta_Er[i,]+Z_jk%*%beta_Mr[i,]+W_k%*%beta_Dr[i,]+rep(phir[i,],n_jk)
  v<-rep(sqrt(kappa2_jkr[i,]),n_jk)
  
  y_ppp<-rnorm(n,mean=mean,sd=v)
  
  #Calculus of the general statistic
  general_statistics_ridge[i,]<-c(mean(y_ppp),sd(y_ppp)/abs(mean(y_ppp)),quantile(y_ppp,c(0.5)))
  
  municipal_mean_ridge[i,]<-tapply(y_ppp,index,mean)
  variance_mean_ridge[i,]<-tapply(y_ppp,index,var)
  departamental_mean_ridge[i,]<-tapply(y_ppp,index_dep,mean)
  departamental_variance_mean_ridge[i,]<-tapply(y_ppp,index_dep,var)
  
  
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
}


# PPP for municipal mean, and variance

municipal_PPP<-data.frame(
  PPP_mean_basic=rep(NA,m),
  PPP_var_basic=rep(NA,m),
  PPP_mean_ridge=rep(NA,m),
  PPP_var_ridge=rep(NA,m),
  PPP_mean_lasso=rep(NA,m),
  PPP_var_lasso=rep(NA,m)
)

departamental_PPP<-data.frame(
  PPP_mean_basic=rep(NA,d),
  PPP_var_basic=rep(NA,d),
  PPP_mean_ridge=rep(NA,d),
  PPP_var_ridge=rep(NA,d),
  PPP_mean_lasso=rep(NA,d),
  PPP_var_lasso=rep(NA,d)
)

for (i in c(1:m)) {
  municipal_PPP$PPP_mean_basic[i]<-mean(municipal_mean_basic[,i]<= municipal_statistic$mean[1])
  municipal_PPP$PPP_var_basic[i]<-mean(variance_mean_basic[,i]<= municipal_statistic$var[1])
}

for (i in c(1:m)) {
  municipal_PPP$PPP_mean_ridge[i]<-mean(municipal_mean_ridge[,i]<=municipal_statistic$mean[1])
  municipal_PPP$PPP_var_ridge[i]<-mean(variance_mean_ridge[,i]<=municipal_statistic$var[1])
}

for (i in c(1:m)) {
  municipal_PPP$PPP_mean_lasso[i]<-mean(municipal_mean_lasso[,i]<= municipal_statistic$mean[1])
  municipal_PPP$PPP_var_lasso[i]<-mean(variance_mean_lasso[,i]<= municipal_statistic$var[1])
}


for (i in c(1:d)) {
  departamental_PPP$PPP_mean_basic[i]<-mean(departamental_mean_basic[,i]<= departamental_statistic$mean[1])
  departamental_PPP$PPP_var_basic[i]<-mean(departamental_variance_mean_basic[,i]<= departamental_statistic$var[1])
}

for (i in c(1:d)) {
  departamental_PPP$PPP_mean_ridge[i]<-mean(departamental_mean_ridge[,i]<=departamental_statistic$mean[1])
  departamental_PPP$PPP_var_ridge[i]<-mean(departamental_variance_mean_ridge[,i]<=departamental_statistic$var[1])
}

for (i in c(1:d)) {
  departamental_PPP$PPP_mean_lasso[i]<-mean(departamental_mean_lasso[,i]<= departamental_statistic$mean[1])
  departamental_PPP$PPP_var_lasso[i]<-mean(departamental_variance_mean_lasso[,i]<= departamental_statistic$var[1])
}


#Global PPP

Global_PPP<-data.frame(
  Statistic=c("Mean","CV","Median"),
  basic=c(mean(general_statistics_basic[,1]<global_statistic$mean),mean(general_statistics_basic[,2]<global_statistic$CV),mean(general_statistics_basic[,3]<global_statistic$median)),
  rigde=c(mean(general_statistics_ridge[,1]<global_statistic$mean),mean(general_statistics_ridge[,2]<global_statistic$CV),mean(general_statistics_ridge[,3]<global_statistic$median)),
  lasso=c(mean(general_statistics_lasso[,1]<global_statistic$mean),mean(general_statistics_lasso[,2]<global_statistic$CV),mean(general_statistics_lasso[,3]<global_statistic$median))
)


municipal_PPP_boxplot_df<-data.frame(
  scenario=c(rep("Basic",m),rep("Ridge",m),rep("Lasso",m)),
  PPP_mean=c(municipal_PPP$PPP_mean_basic,municipal_PPP$PPP_mean_ridge,municipal_PPP$PPP_mean_lasso),
  PPP_var=c(municipal_PPP$PPP_var_basic,municipal_PPP$PPP_var_ridge,municipal_PPP$PPP_var_lasso)
)

departamental_PPP_boxplot_df<-data.frame(
  scenario=c(rep("Basic",d),rep("Ridge",d),rep("Lasso",d)),
  PPP_mean=c(departamental_PPP$PPP_mean_basic,departamental_PPP$PPP_mean_ridge,departamental_PPP$PPP_mean_lasso),
  PPP_var=c(departamental_PPP$PPP_var_basic,departamental_PPP$PPP_var_ridge,departamental_PPP$PPP_var_lasso)
)

All_model_PPP<-list(municipal_PPP_boxplot_df=municipal_PPP_boxplot_df,departamental_PPP_boxplot_df=departamental_PPP_boxplot_df)
save(All_model_PPP,file="All_model_real_data_PPP.RData")


ggplot(municipal_PPP_boxplot_df, aes(x = scenario, y = PPP_mean, fill = scenario, colour = scenario)) +
  geom_boxplot(
    alpha = 0.7,          # Transparencia media
    outlier.shape = NA   # Ocultamos outliers  # Borde gris oscuro uniforme
  ) +
  geom_jitter(
    aes(color = scenario), 
    width = 0.1, 
    alpha = 0.25
  ) +
  scale_color_discrete_sequential(
    palette = "BurgYl",
    rev = FALSE,
    name = "Model"
  )+
  scale_fill_discrete_sequential(
    palette="BurgYl",
    rev=FALSE,
    name="Model"
  )+
  theme_minimal()+
  labs(title = " Posterior Predictive P-value",subtitle = "Municipal posterior mean", x = NULL, y = NULL,fill="Probability")+theme(
    plot.title = element_text(
      color = "#7C3B5E",   
      size  = 18,          
      face  = "bold",     
      hjust = 0.5          
    ),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
  )


ggplot(municipal_PPP_boxplot_df, aes(x = scenario, y = PPP_var, fill = scenario, colour = scenario)) +
  geom_boxplot(
    alpha = 0.7,          # Transparencia media
    outlier.shape = NA   # Ocultamos outliers  # Borde gris oscuro uniforme
  ) +
  geom_jitter(
    aes(color = scenario), 
    width = 0.1, 
    alpha = 0.25
  ) +
  scale_color_discrete_sequential(
    palette = "BurgYl",
    rev = FALSE,
    name = "Model"
  )+
  scale_fill_discrete_sequential(
    palette="BurgYl",
    rev=FALSE,
    name="Model"
  )+
  theme_minimal()+
  labs(title = " Posterior Predictive P-value",subtitle = "Municipal posterior variance", x = NULL, y = NULL,fill="Probability")+theme(
    plot.title = element_text(
      color = "#7C3B5E",   
      size  = 18,          
      face  = "bold",     
      hjust = 0.5          
    ),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
  )



ggplot(municipal_PPP_boxplot_df, aes(x = scenario, y = PPP_var, fill = scenario, colour = scenario)) +
  geom_boxplot(
    alpha = 0.7,          # Transparencia media
    outlier.shape = NA   # Ocultamos outliers  # Borde gris oscuro uniforme
  ) +
  geom_jitter(
    aes(color = scenario), 
    width = 0.1, 
    alpha = 0.25
  ) +
  scale_color_discrete_sequential(
    palette = "BurgYl",
    rev = FALSE,
    name = "Model"
  )+
  scale_fill_discrete_sequential(
    palette="BurgYl",
    rev=FALSE,
    name="Model"
  )+
  theme_minimal()+
  labs(title = " Posterior Predictive P-value",subtitle = "Municipal posterior variance", x = NULL, y = NULL,fill="Probability")+theme(
    plot.title = element_text(
      color = "#7C3B5E",   
      size  = 18,          
      face  = "bold",     
      hjust = 0.5          
    ),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
  )


tapply(departamental_PPP_boxplot_df, departamental_PPP_boxplot_df$scenario, summary)

