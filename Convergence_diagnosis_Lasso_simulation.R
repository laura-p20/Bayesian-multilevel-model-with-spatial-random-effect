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
library(RColorBrewer)
library(colorspace)

#Set the work directory to get a previus result of the MCMC using the real data

load("data_ready.RData")

#Upload the covariables matrices and some numbers
X_ijk<-data_ready$X_ijk
Z_jk<-data_ready$Z_jk
Z_jk<-scale(Z_jk)
W_k<-data_ready$W_k
W_k<-scale(W_k)
n<-data_ready$n
n_jk<-data_ready$n_jk
m<-data_ready$m

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

#Extracting the global score simulated in each scenario

load("simulation_response_scene_1.RData")
load("simulation_response_scene_2.RData")
load("simulation_response_scene_3.RData")

y_ijk_scenes<-cbind(y_ijk_1,y_ijk_2,y_ijk_3)
#Extracting the scenarios values 

load("simulation_values.RData")

##===========Extracting the parameter of the samples----

# Scenario 1
load("MCMC_lasso_simulated_e1.RData")

beta1<-THETA[[1]]
beta_E1<-THETA[[2]]
beta_M1<-THETA[[3]]
beta_D1<-THETA[[4]]

sigma2_beta1<-THETA[[5]]
tau2_E1<-THETA[[6]]
tau2_M1<-THETA[[7]]
tau2_D1<-THETA[[8]]

lamda2_E1<-THETA[[9]]
lamda2_M1<-THETA[[10]]
lamda2_D1<-THETA[[11]]

phi1<-THETA[[12]]
tau_phi1<-THETA[[13]]


kappa2_jk1<-THETA[[14]]
kappa2_k1<-THETA[[15]]

alpha_kappa1<-THETA[[16]]
beta_kappa1<-THETA[[17]]

acr1<-THETA[[18]]

#Scenario 2
load("MCMC_lasso_simulated_e2.RData")

beta2<-THETA[[1]]
beta_E2<-THETA[[2]]
beta_M2<-THETA[[3]]
beta_D2<-THETA[[4]]

sigma2_beta2<-THETA[[5]]
tau2_E2<-THETA[[6]]
tau2_M2<-THETA[[7]]
tau2_D2<-THETA[[8]]

lamda2_E2<-THETA[[9]]
lamda2_M2<-THETA[[10]]
lamda2_D2<-THETA[[11]]

phi2<-THETA[[12]]
tau_phi2<-THETA[[13]]


kappa2_jk2<-THETA[[14]]
kappa2_k2<-THETA[[15]]

alpha_kappa2<-THETA[[16]]
beta_kappa2<-THETA[[17]]

acr2<-THETA[[18]]


load("MCMC_lasso_simulated_e3.RData")

beta3<-THETA[[1]]
beta_E3<-THETA[[2]]
beta_M3<-THETA[[3]]
beta_D3<-THETA[[4]]

sigma2_beta3<-THETA[[5]]
tau2_E3<-THETA[[6]]
tau2_M3<-THETA[[7]]
tau2_D3<-THETA[[8]]

lamda2_E3<-THETA[[9]]
lamda2_M3<-THETA[[10]]
lamda2_D3<-THETA[[11]]

phi3<-THETA[[12]]
tau_phi3<-THETA[[13]]


kappa2_jk3<-THETA[[14]]
kappa2_k3<-THETA[[15]]

alpha_kappa3<-THETA[[16]]
beta_kappa3<-THETA[[17]]

acr3<-THETA[[18]]


nsams<-length(beta1)

##============= Convergence Analysis  ========

###========= Log-Likelihood=========
ll_1<-rep(NA,nsams)
ll_2<-rep(NA,nsams)
ll_3<-rep(NA,nsams)

ll_true<-data.frame(
  ll_true_1=0,
  ll_true_2=0,
  ll_true_3=0
)


#Calculus of the real log-likelihoods

for(i in c(2:4)){
  mean_true<-scenarios_parameter$student[[i]][1]+ X_ijk%*%scenarios_parameter$student[[i]][-1]+Z_jk%*%scenarios_parameter$Municipal[[i]]+W_k%*%scenarios_parameter$Departamental[[i]]+rep(scenarios_parameter$phi_post,n_jk)
  ll_true[1,i-1]<-sum(dnorm(x=y_ijk_scenes[,i-1],mean=mean_true,sd=sqrt(rep(scenarios_parameter$kappa2_jk_post,n_jk)),log = TRUE))
}

print(ll_true)

#Calculus of the log-likelihoods based on the samples

for(i in c(1:nsams)){
  
  mean<-rep(beta1[i],n)+X_ijk%*%beta_E1[i,]+Z_jk%*%beta_M1[i,]+W_k%*%beta_D1[i,]+rep(phi1[i,],n_jk)
  v<-rep(sqrt(kappa2_jk1[i,]),n_jk)
  ll_1[i]<-sum(dnorm(x=y_ijk_1,mean=mean,sd=v,log = TRUE))
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
}


for(i in c(1:nsams)){
  
  mean<-rep(beta2[i],n)+X_ijk%*%beta_E2[i,]+Z_jk%*%beta_M2[i,]+W_k%*%beta_D2[i,]+rep(phi2[i,],n_jk)
  v<-rep(sqrt(kappa2_jk2[i,]),n_jk)
  ll_2[i]<-sum(dnorm(x=y_ijk_2,mean=mean,sd=v,log = TRUE))
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
}

for(i in c(1:nsams)){
  
  mean<-rep(beta3[i],n)+X_ijk%*%beta_E3[i,]+Z_jk%*%beta_M3[i,]+W_k%*%beta_D3[i,]+rep(phi3[i,],n_jk)
  v<-rep(sqrt(kappa2_jk3[i,]),n_jk)
  ll_3[i]<-sum(dnorm(x=y_ijk_3,mean=mean,sd=v,log = TRUE))
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
}



#Plot


par(mfrow=c(1,3))
plot(ll_1,pch = 20, cex = 0.45,xlab = "",ylab = "Log likelihood",main="Scenario 1",col="#70284a",col.main="#70284a",ylim = c(-1279319 ,-1232250))
abline(h=mean(ll_true[1,1]),col="#A0153E",lwd="2")
plot(ll_2,pch = 20, cex = 0.45,xlab = "Iteration",ylab = "",main="Scenario 2",col="#dc7176",col.main="#70284a",ylim = c(-1279319 ,-1232250))
abline(h=mean(ll_true[1,2]),col="#A0153E",lwd="2")
plot(ll_3,pch = 20, cex = 0.45,xlab = "",ylab = "",main="Scenario 3",col="#f2a28a",col.main="#70284a",ylim = c(-1279319 ,-1232250))
abline(h=mean(ll_true[1,3]),col="#A0153E",lwd="2")



Lasso_ll<-data.frame(
  ll_1=ll_1,
  ll_2=ll_2,
  ll_3=ll_3
)

save(Lasso_ll,file = "Lasso_log_like.RData")
##======== Information criteria====

### =======Scenario 1=======
#### ===== Effective sample size summary for each parameter family====
neff_beta <- round(coda::effectiveSize(beta1), 0)
neff_beta_E<-apply(beta_E1, 2,effectiveSize)
neff_beta_M<-apply(beta_M1, 2,effectiveSize)
neff_beta_D<-apply(beta_D1, 2,effectiveSize)

neff_beta_summary<-summary(c(neff_beta,neff_beta_E,neff_beta_M,neff_beta_D))

neff_sigma2_beta<-round(coda::effectiveSize(sigma2_beta1), 0)

neff_tau2_E<-apply(tau2_E1, 2,effectiveSize)
neff_tau2_M<-apply(tau2_M1, 2,effectiveSize)
neff_tau2_D<-apply(tau2_D1, 2,effectiveSize)

neff_tau2_summary<-summary(c(neff_tau2_E,neff_tau2_M,neff_tau2_D))


neff_phi<-apply(phi1, 2,effectiveSize)
neff_phi_summary<-summary(neff_phi)

neff_tau_phi<-round(coda::effectiveSize(tau_phi1), 0)

neff_kappa2_jk<-apply(kappa2_jk1, 2,effectiveSize)
neff_kappa2_jk_summary<-summary(neff_kappa2_jk)
neff_kappa2_k<-apply(kappa2_k1, 2,effectiveSize)
neff_kappa2_k_summary<-summary(neff_kappa2_jk)

neff_alpha_kappa<-round(coda::effectiveSize(alpha_kappa1), 0)
neff_beta_kappa<-round(coda::effectiveSize(beta_kappa1), 0)

neff_summary1<-data.frame(
  parameter_family=c("Beta","Tau2","Phi","Kappa2_jk","Kappa2_k"),
  min=c(neff_beta_summary[1],neff_tau2_summary[1],neff_phi_summary[1],neff_kappa2_jk_summary[1],neff_kappa2_k_summary[1]),
  Q_1=c(neff_beta_summary[2],neff_tau2_summary[2],neff_phi_summary[2],neff_kappa2_jk_summary[2],neff_kappa2_k_summary[2]),
  Median=c(neff_beta_summary[3],neff_tau2_summary[3],neff_phi_summary[3],neff_kappa2_jk_summary[3],neff_kappa2_k_summary[3]),
  Mean= c(neff_beta_summary[4],neff_tau2_summary[4],neff_phi_summary[4],neff_kappa2_jk_summary[4],neff_kappa2_k_summary[4]),
  Q_3=c(neff_beta_summary[5],neff_tau2_summary[5],neff_phi_summary[5],neff_kappa2_jk_summary[5],neff_kappa2_k_summary[5]),
  max=c(neff_beta_summary[6],neff_tau2_summary[6],neff_phi_summary[6],neff_kappa2_jk_summary[6],neff_kappa2_k_summary[6])
)

####======Monte Carlo Standard Error======
EMC_beta <- sd(beta1)/sqrt(neff_beta)
EMC_beta_E<-apply(beta_E1, 2, sd)/sqrt(neff_beta_E)
EMC_beta_M<-apply(beta_M1, 2, sd)/sqrt(neff_beta_M)
EMC_beta_D<-apply(beta_D1, 2, sd)/sqrt(neff_beta_D)

EMC_beta_summary<-summary(c(EMC_beta,EMC_beta_E,EMC_beta_M,EMC_beta_D))

EMC_tau2_E<-apply(tau2_E1, 2, sd)/sqrt(neff_tau2_E)
EMC_tau2_M<-apply(tau2_M1, 2, sd)/sqrt(neff_tau2_M)
EMC_tau2_D<-apply(tau2_D1, 2, sd)/sqrt(neff_tau2_D)

EMC_tau2_summary<-summary(c(EMC_tau2_E,EMC_tau2_M,EMC_tau2_D))


EMC_phi<-apply(phi1, 2, sd)/sqrt(neff_phi)
EMC_phi_summary<-summary(EMC_phi)

EMC_tau_phi<-sd(tau_phi1)/sqrt(neff_tau_phi)

EMC_kappa2_jk<-apply(kappa2_jk1, 2,sd)/sqrt(neff_kappa2_jk)
EMC_kappa2_jk_summary<-summary(EMC_kappa2_jk)
EMC_kappa2_k<-apply(kappa2_k1, 2,sd)/sqrt(neff_kappa2_k)
EMC_kappa2_k_summary<-summary(EMC_kappa2_k)

EMC_alpha_kappa<-sd(alpha_kappa1)/sqrt(neff_alpha_kappa)
EMC_beta_kappa<-sd(beta_kappa1)/sqrt(neff_beta_kappa)

EMC_summary1<-data.frame(
  parameter_family=c("Beta","Tau2","Phi","Kappa2_jk","Kappa2_k"),
  min=c(EMC_beta_summary[1],EMC_tau2_summary[1],EMC_phi_summary[1],EMC_kappa2_jk_summary[1],EMC_kappa2_k_summary[1]),
  Q_1=c(EMC_beta_summary[2],EMC_tau2_summary[2],EMC_phi_summary[2],EMC_kappa2_jk_summary[2],EMC_kappa2_k_summary[2]),
  Median=c(EMC_beta_summary[3],EMC_tau2_summary[3],EMC_phi_summary[3],EMC_kappa2_jk_summary[3],EMC_kappa2_k_summary[3]),
  Mean= c(EMC_beta_summary[4],EMC_tau2_summary[4],EMC_phi_summary[4],EMC_kappa2_jk_summary[4],EMC_kappa2_k_summary[4]),
  Q_3=c(EMC_beta_summary[5],EMC_tau2_summary[5],EMC_phi_summary[5],EMC_kappa2_jk_summary[5],EMC_kappa2_k_summary[5]),
  max=c(EMC_beta_summary[6],EMC_tau2_summary[6],EMC_phi_summary[6],EMC_kappa2_jk_summary[6],EMC_kappa2_k_summary[6])
)

####========= Monte Carlo Variation Coefficient % ====

CVMC_beta <- 100*EMC_beta/abs(mean(beta1))

CVMC_beta_E <- 100*EMC_beta_E/abs(apply(beta_E1,2,mean))
CVMC_beta_M <- 100*EMC_beta_M/abs(apply(beta_M1,2,mean))
CVMC_beta_D <- 100*EMC_beta_D/abs(apply(beta_D1,2,mean))

CVMC_beta_summary<-summary(c(CVMC_beta,CVMC_beta_E,CVMC_beta_M,CVMC_beta_D))


CVMC_tau2_E <- 100*EMC_tau2_E/abs(apply(tau2_E1,2,mean))
CVMC_tau2_M <- 100*EMC_tau2_M/abs(apply(tau2_M1,2,mean))
CVMC_tau2_D <- 100*EMC_tau2_D/abs(apply(tau2_D1,2,mean))

CVMC_tau2_summary<-summary(c(CVMC_tau2_E,CVMC_tau2_M,CVMC_tau2_D))


CVMC_phi <- 100*EMC_phi/abs(apply(phi1,2,mean))
CVMC_phi_summary<-summary(CVMC_phi)
CVMC_tau_phi <- 100*EMC_tau_phi/abs(mean(tau_phi1))

CVMC_kappa2_jk<- 100*EMC_kappa2_jk/abs(apply(kappa2_jk1,2,mean))
CVMC_kappa2_jk_summary<-summary(CVMC_kappa2_jk)
CVMC_kappa2_k<- 100*EMC_kappa2_k/abs(apply(kappa2_k1,2,mean))
CVMC_kappa2_k_summary<-summary(CVMC_kappa2_k)

CVMC_alpha_kappa <- 100*EMC_alpha_kappa/abs(mean(alpha_kappa1))
CVMC_beta_kappa <- 100*EMC_beta_kappa/abs(mean(beta_kappa1))


CVMC_summary1<-data.frame(
  parameter_family=c("Beta","Tau2","Phi","Kappa2_jk","Kappa2_k"),
  min=c(CVMC_beta_summary[1],CVMC_tau2_summary[1],CVMC_phi_summary[1],CVMC_kappa2_jk_summary[1],CVMC_kappa2_k_summary[1]),
  Q_1=c(CVMC_beta_summary[2],CVMC_tau2_summary[2],CVMC_phi_summary[2],CVMC_kappa2_jk_summary[2],CVMC_kappa2_k_summary[2]),
  Median=c(CVMC_beta_summary[3],CVMC_tau2_summary[3],CVMC_phi_summary[3],CVMC_kappa2_jk_summary[3],CVMC_kappa2_k_summary[3]),
  Mean= c(CVMC_beta_summary[4],CVMC_tau2_summary[4],CVMC_phi_summary[4],CVMC_kappa2_jk_summary[4],CVMC_kappa2_k_summary[4]),
  Q_3=c(CVMC_beta_summary[5],CVMC_tau2_summary[5],CVMC_phi_summary[5],CVMC_kappa2_jk_summary[5],CVMC_kappa2_k_summary[5]),
  max=c(CVMC_beta_summary[6],CVMC_tau2_summary[6],CVMC_phi_summary[6],CVMC_kappa2_jk_summary[6],CVMC_kappa2_k_summary[6])
)



###======Scenario 2=====

#### ===== Effective sample size summary for each parameter family====


neff_beta <- round(coda::effectiveSize(beta2), 0)
neff_beta_E<-apply(beta_E2, 2,effectiveSize)
neff_beta_M<-apply(beta_M2, 2,effectiveSize)
neff_beta_D<-apply(beta_D2, 2,effectiveSize)

neff_beta_summary<-summary(c(neff_beta,neff_beta_E,neff_beta_M,neff_beta_D))

neff_sigma2_beta<-round(coda::effectiveSize(sigma2_beta2), 0)
neff_tau2_E<-apply(tau2_E2, 2,effectiveSize)
neff_tau2_M<-apply(tau2_M2, 2,effectiveSize)
neff_tau2_D<-apply(tau2_D2, 2,effectiveSize)

neff_tau2_summary<-summary(c(neff_tau2_E,neff_tau2_M,neff_tau2_D))


neff_phi<-apply(phi2, 2,effectiveSize)
neff_phi_summary<-summary(neff_phi)

neff_tau_phi<-round(coda::effectiveSize(tau_phi2), 0)

neff_kappa2_jk<-apply(kappa2_jk2, 2,effectiveSize)
neff_kappa2_jk_summary<-summary(neff_kappa2_jk)
neff_kappa2_k<-apply(kappa2_k2, 2,effectiveSize)
neff_kappa2_k_summary<-summary(neff_kappa2_jk)

neff_alpha_kappa<-round(coda::effectiveSize(alpha_kappa2), 0)
neff_beta_kappa<-round(coda::effectiveSize(beta_kappa2), 0)

neff_summary2<-data.frame(
  parameter_family=c("Beta","Tau2","Phi","Kappa2_jk","Kappa2_k"),
  min=c(neff_beta_summary[1],neff_tau2_summary[1],neff_phi_summary[1],neff_kappa2_jk_summary[1],neff_kappa2_k_summary[1]),
  Q_1=c(neff_beta_summary[2],neff_tau2_summary[2],neff_phi_summary[2],neff_kappa2_jk_summary[2],neff_kappa2_k_summary[2]),
  Median=c(neff_beta_summary[3],neff_tau2_summary[3],neff_phi_summary[3],neff_kappa2_jk_summary[3],neff_kappa2_k_summary[3]),
  Mean= c(neff_beta_summary[4],neff_tau2_summary[4],neff_phi_summary[4],neff_kappa2_jk_summary[4],neff_kappa2_k_summary[4]),
  Q_3=c(neff_beta_summary[5],neff_tau2_summary[5],neff_phi_summary[5],neff_kappa2_jk_summary[5],neff_kappa2_k_summary[5]),
  max=c(neff_beta_summary[6],neff_tau2_summary[6],neff_phi_summary[6],neff_kappa2_jk_summary[6],neff_kappa2_k_summary[6])
)


####======Monte Carlo Standard Error======


EMC_beta <- sd(beta2)/sqrt(neff_beta)
EMC_beta_E<-apply(beta_E2, 2, sd)/sqrt(neff_beta_E)
EMC_beta_M<-apply(beta_M2, 2, sd)/sqrt(neff_beta_M)
EMC_beta_D<-apply(beta_D2, 2, sd)/sqrt(neff_beta_D)


EMC_tau2_E<-apply(tau2_E2, 2, sd)/sqrt(neff_tau2_E)
EMC_tau2_M<-apply(tau2_M2, 2, sd)/sqrt(neff_tau2_M)
EMC_tau2_D<-apply(tau2_D2, 2, sd)/sqrt(neff_tau2_D)

EMC_tau2_summary<-summary(c(EMC_tau2_E,EMC_tau2_M,EMC_tau2_D))


EMC_beta_summary<-summary(c(EMC_beta,EMC_beta_E,EMC_beta_M,EMC_beta_D))

EMC_phi<-apply(phi2, 2, sd)/sqrt(neff_phi)
EMC_phi_summary<-summary(EMC_phi)

EMC_tau_phi<-sd(tau_phi2)/sqrt(neff_tau_phi)

EMC_kappa2_jk<-apply(kappa2_jk2, 2,sd)/sqrt(neff_kappa2_jk)
EMC_kappa2_jk_summary<-summary(EMC_kappa2_jk)
EMC_kappa2_k<-apply(kappa2_k2, 2,sd)/sqrt(neff_kappa2_k)
EMC_kappa2_k_summary<-summary(EMC_kappa2_k)

EMC_alpha_kappa<-sd(alpha_kappa2)/sqrt(neff_alpha_kappa)
EMC_beta_kappa<-sd(beta_kappa2)/sqrt(neff_beta_kappa)

EMC_summary2<-data.frame(
  parameter_family=c("Beta","Tau2","Phi","Kappa2_jk","Kappa2_k"),
  min=c(EMC_beta_summary[1],EMC_tau2_summary[1],EMC_phi_summary[1],EMC_kappa2_jk_summary[1],EMC_kappa2_k_summary[1]),
  Q_1=c(EMC_beta_summary[2],EMC_tau2_summary[2],EMC_phi_summary[2],EMC_kappa2_jk_summary[2],EMC_kappa2_k_summary[2]),
  Median=c(EMC_beta_summary[3],EMC_tau2_summary[3],EMC_phi_summary[3],EMC_kappa2_jk_summary[3],EMC_kappa2_k_summary[3]),
  Mean= c(EMC_beta_summary[4],EMC_tau2_summary[4],EMC_phi_summary[4],EMC_kappa2_jk_summary[4],EMC_kappa2_k_summary[4]),
  Q_3=c(EMC_beta_summary[5],EMC_tau2_summary[5],EMC_phi_summary[5],EMC_kappa2_jk_summary[5],EMC_kappa2_k_summary[5]),
  max=c(EMC_beta_summary[6],EMC_tau2_summary[6],EMC_phi_summary[6],EMC_kappa2_jk_summary[6],EMC_kappa2_k_summary[6])
)

####========= Monte Carlo Variation Coefficient % ====
CVMC_beta <- 100*EMC_beta/abs(mean(beta2))

CVMC_beta_E <- 100*EMC_beta_E/abs(apply(beta_E2,2,mean))
CVMC_beta_M <- 100*EMC_beta_M/abs(apply(beta_M2,2,mean))
CVMC_beta_D <- 100*EMC_beta_D/abs(apply(beta_D2,2,mean))

CVMC_beta_summary<-summary(c(CVMC_beta,CVMC_beta_E,CVMC_beta_M,CVMC_beta_D))


CVMC_tau2_E <- 100*EMC_tau2_E/abs(apply(tau2_E2,2,mean))
CVMC_tau2_M <- 100*EMC_tau2_M/abs(apply(tau2_M2,2,mean))
CVMC_tau2_D <- 100*EMC_tau2_D/abs(apply(tau2_D2,2,mean))

CVMC_tau2_summary<-summary(c(CVMC_tau2_E,CVMC_tau2_M,CVMC_tau2_D))


CVMC_phi <- 100*EMC_phi/abs(apply(phi2,2,mean))
CVMC_phi_summary<-summary(CVMC_phi)
CVMC_tau_phi <- 100*EMC_tau_phi/abs(mean(tau_phi2))

CVMC_kappa2_jk<- 100*EMC_kappa2_jk/abs(apply(kappa2_jk2,2,mean))
CVMC_kappa2_jk_summary<-summary(CVMC_kappa2_jk)
CVMC_kappa2_k<- 100*EMC_kappa2_k/abs(apply(kappa2_k2,2,mean))
CVMC_kappa2_k_summary<-summary(CVMC_kappa2_k)

CVMC_alpha_kappa <- 100*EMC_alpha_kappa/abs(mean(alpha_kappa2))
CVMC_beta_kappa <- 100*EMC_beta_kappa/abs(mean(beta_kappa2))


CVMC_summary2<-data.frame(
  parameter_family=c("Beta","Tau2","Phi","Kappa2_jk","Kappa2_k"),
  min=c(CVMC_beta_summary[1],CVMC_tau2_summary[1],CVMC_phi_summary[1],CVMC_kappa2_jk_summary[1],CVMC_kappa2_k_summary[1]),
  Q_1=c(CVMC_beta_summary[2],CVMC_tau2_summary[2],CVMC_phi_summary[2],CVMC_kappa2_jk_summary[2],CVMC_kappa2_k_summary[2]),
  Median=c(CVMC_beta_summary[3],CVMC_tau2_summary[3],CVMC_phi_summary[3],CVMC_kappa2_jk_summary[3],CVMC_kappa2_k_summary[3]),
  Mean= c(CVMC_beta_summary[4],CVMC_tau2_summary[4],CVMC_phi_summary[4],CVMC_kappa2_jk_summary[4],CVMC_kappa2_k_summary[4]),
  Q_3=c(CVMC_beta_summary[5],CVMC_tau2_summary[5],CVMC_phi_summary[5],CVMC_kappa2_jk_summary[5],CVMC_kappa2_k_summary[5]),
  max=c(CVMC_beta_summary[6],CVMC_tau2_summary[6],CVMC_phi_summary[6],CVMC_kappa2_jk_summary[6],CVMC_kappa2_k_summary[6])
)

###====Scenario 3=======

#### ===== Effective sample size summary for each parameter family====

neff_beta <- round(coda::effectiveSize(beta3), 0)
neff_beta_E<-apply(beta_E3, 2,effectiveSize)
neff_beta_M<-apply(beta_M3, 2,effectiveSize)
neff_beta_D<-apply(beta_D3, 2,effectiveSize)

neff_beta_summary<-summary(c(neff_beta,neff_beta_E,neff_beta_M,neff_beta_D))

neff_tau2_E<-apply(tau2_E3, 2,effectiveSize)
neff_tau2_M<-apply(tau2_M3, 2,effectiveSize)
neff_tau2_D<-apply(tau2_D3, 2,effectiveSize)

neff_tau2_summary<-summary(c(neff_tau2_E,neff_tau2_M,neff_tau2_D))


neff_sigma2_beta<-round(coda::effectiveSize(sigma2_beta3), 0)
neff_tau2_E<-round(coda::effectiveSize(tau2_E3), 0)
neff_tau2_M<-round(coda::effectiveSize(tau2_M3), 0)
neff_tau2_D<-round(coda::effectiveSize(tau2_D3), 0)

neff_phi<-apply(phi3, 2,effectiveSize)
neff_phi_summary<-summary(neff_phi)

neff_tau_phi<-round(coda::effectiveSize(tau_phi3), 0)

neff_kappa2_jk<-apply(kappa2_jk3, 2,effectiveSize)
neff_kappa2_jk_summary<-summary(neff_kappa2_jk)
neff_kappa2_k<-apply(kappa2_k3, 2,effectiveSize)
neff_kappa2_k_summary<-summary(neff_kappa2_jk)

neff_alpha_kappa<-round(coda::effectiveSize(alpha_kappa3), 0)
neff_beta_kappa<-round(coda::effectiveSize(beta_kappa3), 0)

neff_summary3<-data.frame(
  parameter_family=c("Beta","Tau2","Phi","Kappa2_jk","Kappa2_k"),
  min=c(neff_beta_summary[1],neff_tau2_summary[1],neff_phi_summary[1],neff_kappa2_jk_summary[1],neff_kappa2_k_summary[1]),
  Q_1=c(neff_beta_summary[2],neff_tau2_summary[2],neff_phi_summary[2],neff_kappa2_jk_summary[2],neff_kappa2_k_summary[2]),
  Median=c(neff_beta_summary[3],neff_tau2_summary[3],neff_phi_summary[3],neff_kappa2_jk_summary[3],neff_kappa2_k_summary[3]),
  Mean= c(neff_beta_summary[4],neff_tau2_summary[4],neff_phi_summary[4],neff_kappa2_jk_summary[4],neff_kappa2_k_summary[4]),
  Q_3=c(neff_beta_summary[5],neff_tau2_summary[5],neff_phi_summary[5],neff_kappa2_jk_summary[5],neff_kappa2_k_summary[5]),
  max=c(neff_beta_summary[6],neff_tau2_summary[6],neff_phi_summary[6],neff_kappa2_jk_summary[6],neff_kappa2_k_summary[6])
)

####======Monte Carlo Standard Error======



EMC_beta <- sd(beta3)/sqrt(neff_beta)
EMC_beta_E<-apply(beta_E3, 2, sd)/sqrt(neff_beta_E)
EMC_beta_M<-apply(beta_M3, 2, sd)/sqrt(neff_beta_M)
EMC_beta_D<-apply(beta_D3, 2, sd)/sqrt(neff_beta_D)

EMC_beta_summary<-summary(c(EMC_beta,EMC_beta_E,EMC_beta_M,EMC_beta_D))

EMC_tau2_E<-apply(tau2_E3, 2, sd)/sqrt(neff_tau2_E)
EMC_tau2_M<-apply(tau2_M3, 2, sd)/sqrt(neff_tau2_M)
EMC_tau2_D<-apply(tau2_D3, 2, sd)/sqrt(neff_tau2_D)

EMC_tau2_summary<-summary(c(EMC_tau2_E,EMC_tau2_M,EMC_tau2_D))


EMC_phi<-apply(phi3, 2, sd)/sqrt(neff_phi)
EMC_phi_summary<-summary(EMC_phi)

EMC_tau_phi<-sd(tau_phi3)/sqrt(neff_tau_phi)

EMC_kappa2_jk<-apply(kappa2_jk3, 2,sd)/sqrt(neff_kappa2_jk)
EMC_kappa2_jk_summary<-summary(EMC_kappa2_jk)
EMC_kappa2_k<-apply(kappa2_k3, 2,sd)/sqrt(neff_kappa2_k)
EMC_kappa2_k_summary<-summary(EMC_kappa2_k)

EMC_alpha_kappa<-sd(alpha_kappa3)/sqrt(neff_alpha_kappa)
EMC_beta_kappa<-sd(beta_kappa3)/sqrt(neff_beta_kappa)

EMC_summary3<-data.frame(
  parameter_family=c("Beta","Tau2","Phi","Kappa2_jk","Kappa2_k"),
  min=c(EMC_beta_summary[1],EMC_tau2_summary[1],EMC_phi_summary[1],EMC_kappa2_jk_summary[1],EMC_kappa2_k_summary[1]),
  Q_1=c(EMC_beta_summary[2],EMC_tau2_summary[2],EMC_phi_summary[2],EMC_kappa2_jk_summary[2],EMC_kappa2_k_summary[2]),
  Median=c(EMC_beta_summary[3],EMC_tau2_summary[3],EMC_phi_summary[3],EMC_kappa2_jk_summary[3],EMC_kappa2_k_summary[3]),
  Mean= c(EMC_beta_summary[4],EMC_tau2_summary[4],EMC_phi_summary[4],EMC_kappa2_jk_summary[4],EMC_kappa2_k_summary[4]),
  Q_3=c(EMC_beta_summary[5],EMC_tau2_summary[5],EMC_phi_summary[5],EMC_kappa2_jk_summary[5],EMC_kappa2_k_summary[5]),
  max=c(EMC_beta_summary[6],EMC_tau2_summary[6],EMC_phi_summary[6],EMC_kappa2_jk_summary[6],EMC_kappa2_k_summary[6])
)
####========= Monte Carlo Variation Coefficient % ====

CVMC_beta <- 100*EMC_beta/abs(mean(beta3))

CVMC_beta_E <- 100*EMC_beta_E/abs(apply(beta_E3,2,mean))
CVMC_beta_M <- 100*EMC_beta_M/abs(apply(beta_M3,2,mean))
CVMC_beta_D <- 100*EMC_beta_D/abs(apply(beta_D3,2,mean))

CVMC_beta_summary<-summary(c(CVMC_beta,CVMC_beta_E,CVMC_beta_M,CVMC_beta_D))

CVMC_tau2_E <- 100*EMC_tau2_E/abs(apply(tau2_E3,2,mean))
CVMC_tau2_M <- 100*EMC_tau2_M/abs(apply(tau2_M3,2,mean))
CVMC_tau2_D <- 100*EMC_tau2_D/abs(apply(tau2_D3,2,mean))

CVMC_tau2_summary<-summary(c(CVMC_tau2_E,CVMC_tau2_M,CVMC_tau2_D))


CVMC_phi <- 100*EMC_phi/abs(apply(phi3,2,mean))
CVMC_phi_summary<-summary(CVMC_phi)
CVMC_tau_phi <- 100*EMC_tau_phi/abs(mean(tau_phi3))

CVMC_kappa2_jk<- 100*EMC_kappa2_jk/abs(apply(kappa2_jk3,2,mean))
CVMC_kappa2_jk_summary<-summary(CVMC_kappa2_jk)
CVMC_kappa2_k<- 100*EMC_kappa2_k/abs(apply(kappa2_k3,2,mean))
CVMC_kappa2_k_summary<-summary(CVMC_kappa2_k)

CVMC_alpha_kappa <- 100*EMC_alpha_kappa/abs(mean(alpha_kappa3))
CVMC_beta_kappa <- 100*EMC_beta_kappa/abs(mean(beta_kappa3))


CVMC_summary3<-data.frame(
  parameter_family=c("Beta","Tau2","Phi","Kappa2_jk","Kappa2_k"),
  min=c(CVMC_beta_summary[1],CVMC_tau2_summary[1],CVMC_phi_summary[1],CVMC_kappa2_jk_summary[1],CVMC_kappa2_k_summary[1]),
  Q_1=c(CVMC_beta_summary[2],CVMC_tau2_summary[2],CVMC_phi_summary[2],CVMC_kappa2_jk_summary[2],CVMC_kappa2_k_summary[2]),
  Median=c(CVMC_beta_summary[3],CVMC_tau2_summary[3],CVMC_phi_summary[3],CVMC_kappa2_jk_summary[3],CVMC_kappa2_k_summary[3]),
  Mean= c(CVMC_beta_summary[4],CVMC_tau2_summary[4],CVMC_phi_summary[4],CVMC_kappa2_jk_summary[4],CVMC_kappa2_k_summary[4]),
  Q_3=c(CVMC_beta_summary[5],CVMC_tau2_summary[5],CVMC_phi_summary[5],CVMC_kappa2_jk_summary[5],CVMC_kappa2_k_summary[5]),
  max=c(CVMC_beta_summary[6],CVMC_tau2_summary[6],CVMC_phi_summary[6],CVMC_kappa2_jk_summary[6],CVMC_kappa2_k_summary[6])
)


##============ Goodness of fit ========



df1<- data.frame(
  variable = c("beta",colnames(beta_E1),colnames(beta_M1),colnames(beta_D1)),
  media = c(mean(beta1),apply(beta_E1,2,mean),apply(beta_M1,2,mean),apply(beta_D1,2,mean)),
  lim_inf=c(bayesian_summary(beta1)[2],apply(beta_E1, 2, bayesian_summary)[2,],apply(beta_M1, 2, bayesian_summary)[2,],apply(beta_D1,2, bayesian_summary)[2,]),
  lim_sup=c(bayesian_summary(beta1)[3],apply(beta_E1, 2, bayesian_summary)[3,],apply(beta_M1, 2, bayesian_summary)[3,],apply(beta_D1, 2, bayesian_summary)[3,]),
  Nivel=c("Intercepto", rep("Student",ncol(beta_E1)),rep("Municipal",ncol(beta_M1)),rep("Departamental",ncol(beta_D1))),
  true_value=c(scenarios_parameter$student[,2],scenarios_parameter$Municipal[,2],scenarios_parameter$Departamental[,2]),
  scenario=rep("scenario 1",34)
)
df1$true_value_in<-ifelse(df1$true_value>=df1$lim_inf & df1$true_value<=df1$lim_sup,1,0)

df1$variable <- recode(df1$variable,
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


df2<- data.frame(
  variable = c("beta",colnames(beta_E2),colnames(beta_M2),colnames(beta_D2)),
  media = c(mean(beta2),apply(beta_E2,2,mean),apply(beta_M2,2,mean),apply(beta_D2,2,mean)),
  lim_inf=c(bayesian_summary(beta2)[2],apply(beta_E2, 2, bayesian_summary)[2,],apply(beta_M2, 2, bayesian_summary)[2,],apply(beta_D2,2, bayesian_summary)[2,]),
  lim_sup=c(bayesian_summary(beta2)[3],apply(beta_E2, 2, bayesian_summary)[3,],apply(beta_M2, 2, bayesian_summary)[3,],apply(beta_D2, 2, bayesian_summary)[3,]),
  Nivel=c("Intercepto", rep("Student",ncol(beta_E2)),rep("Municipal",ncol(beta_M2)),rep("Departamental",ncol(beta_D2))),
  true_value=c(scenarios_parameter$student[,3],scenarios_parameter$Municipal[,3],scenarios_parameter$Departamental[,3]),
  scenario=rep("scenario 2",34)
)
df2$true_value_in<-ifelse(df2$true_value>=df2$lim_inf & df2$true_value<=df2$lim_sup,1,0)

df2$variable <- recode(df2$variable,
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



df3<- data.frame(
  variable = c("beta",colnames(beta_E3),colnames(beta_M3),colnames(beta_D3)),
  media = c(mean(beta3),apply(beta_E3,2,mean),apply(beta_M3,2,mean),apply(beta_D3,2,mean)),
  lim_inf=c(bayesian_summary(beta3)[2],apply(beta_E3, 2, bayesian_summary)[2,],apply(beta_M3, 2, bayesian_summary)[2,],apply(beta_D3,2, bayesian_summary)[2,]),
  lim_sup=c(bayesian_summary(beta3)[3],apply(beta_E3, 2, bayesian_summary)[3,],apply(beta_M3, 2, bayesian_summary)[3,],apply(beta_D3, 2, bayesian_summary)[3,]),
  Nivel=c("Intercepto", rep("Student",ncol(beta_E3)),rep("Municipal",ncol(beta_M3)),rep("Departamental",ncol(beta_D3))),
  true_value=c(scenarios_parameter$student[,4],scenarios_parameter$Municipal[,4],scenarios_parameter$Departamental[,4]),
  scenario=rep("scenario 3",34)
)
df3$true_value_in<-ifelse(df3$true_value>=df3$lim_inf & df3$true_value<=df3$lim_sup,1,0)

df3$variable <- recode(df3$variable,
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


df1<-df1%>%
  arrange(desc(media))
df2<-df2%>%
  arrange(desc(media))
df3<-df3%>%
  arrange(desc(media))

df<-rbind(df1,df2,df3)
df<-df%>%
  arrange(desc(mean))

#A plot in which is possible to visualice posterior means and the crebilibity intervals

ggplot(df[-c(1,35,69),], aes(y = reorder(variable, media), x = media, color = Nivel)) +
  geom_point(size = 1.5) +
  geom_errorbarh(aes(xmin = lim_inf, xmax = lim_sup), height = 0.3) +
  geom_point(aes(x = true_value, y = variable),
             shape = 18,        # Forma de diamante para diferenciar
             size = 1,
             color = "black",   # Color negro para que resalte
             inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "Posterior Mean", y = "Parameter", color = "Level", title ="Covariates posterior mean ranking" ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right") +
  scale_color_manual(values = c("#dc7176", "#70284a", "#f2a28a"))+
  theme_minimal(base_family = "Helvetica") +
  facet_wrap(~scenario,ncol=3,nrow=1,strip.position="top",dir="h",scales="fixed")  +
  theme(
    plot.title = element_text(color = "#7C3B5E",face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40",size = 12, hjust = 0.5),
    axis.title = element_text(color = "black",face = "plain"),
    legend.position = "bottom"
  ) 

# A set of data is going to be generated which each sample of the chain
# then, some test statistics, like the mean, median and variation coefficient
# of the global score.

percentage_in_true_val<-c(mean(df1$true_value_in),mean(df2$true_value_in),mean(df3$true_value_in))


###=======ppp=====

index<-data_ready$id_dep_mun[,2]
index_dep<-data_ready$id_dep_mun[,1]
# Statistic of the data

global_statistic<-data.frame(
  mean_s1=mean(y_ijk_1),
  mean_s2=mean(y_ijk_2),
  mean_s3=mean(y_ijk_3),
  CV_1=CV(y_ijk_1),
  CV_2=CV(y_ijk_2),
  CV_3=CV(y_ijk_3),
  median1=quantile(y_ijk_1,c(0.5)),
  median2=quantile(y_ijk_2,c(0.5)),
  median3=quantile(y_ijk_3,c(0.5))
)



municipal_statistic<-data.frame(
  mean_s1=tapply(y_ijk_1,index,mean),
  var_s1=tapply(y_ijk_1,index,var),
  mean_s2=tapply(y_ijk_2,index,mean),
  var_s2=tapply(y_ijk_2,index,var),
  mean_s3=tapply(y_ijk_3,index,mean),
  var_s3=tapply(y_ijk_3,index,var)
)

departamental_statistic<-data.frame(
  mean_s1=tapply(y_ijk_1,index_dep,mean),
  var_s1=tapply(y_ijk_1,index_dep,var),
  mean_s2=tapply(y_ijk_2,index_dep,mean),
  var_s2=tapply(y_ijk_2,index_dep,var),
  mean_s3=tapply(y_ijk_3,index_dep,mean),
  var_s3=tapply(y_ijk_3,index_dep,var)
)

y_ppp<-0

general_statistics1<-matrix(NA,nrow = nsams,ncol=3)
municipal_mean1<-matrix(NA,nrow = nsams,ncol=length(n_jk))
variance_mean1<-matrix(NA,nrow = nsams,ncol=length(n_jk))
departamental_mean1<-matrix(NA,nrow = nsams,ncol=d)
departamental_variance_mean1<-matrix(NA,nrow = nsams,ncol=d)


general_statistics2<-matrix(NA,nrow = nsams,ncol=3)
municipal_mean2<-matrix(NA,nrow = nsams,ncol=length(n_jk))
variance_mean2<-matrix(NA,nrow = nsams,ncol=length(n_jk))
departamental_mean2<-matrix(NA,nrow = nsams,ncol=d)
departamental_variance_mean2<-matrix(NA,nrow = nsams,ncol=d)

general_statistics3<-matrix(NA,nrow = nsams,ncol=3)
municipal_mean3<-matrix(NA,nrow = nsams,ncol=length(n_jk))
variance_mean3<-matrix(NA,nrow = nsams,ncol=length(n_jk))
departamental_mean3<-matrix(NA,nrow = nsams,ncol=d)
departamental_variance_mean3<-matrix(NA,nrow = nsams,ncol=d)

set.seed(2008)

for(i in c(1:nsams)){
  
  
  #Generation of the data scenario 1
  mean<-rep(beta1[i],n)+X_ijk%*%beta_E1[i,]+Z_jk%*%beta_M1[i,]+W_k%*%beta_D1[i,]+rep(phi1[i,],n_jk)
  v<-rep(sqrt(kappa2_jk1[i,]),n_jk)
  
  y_ppp<-rnorm(n,mean=mean,sd=v)
  
  #Calculus of the general statistic
  general_statistics1[i,]<-c(mean(y_ppp),sd(y_ppp)/abs(mean(y_ppp)),quantile(y_ppp,c(0.5)))
  
  municipal_mean1[i,]<-tapply(y_ppp,index,mean)
  variance_mean1[i,]<-tapply(y_ppp,index,var)
  departamental_mean1[i,]<-tapply(y_ppp,index_dep,mean)
  departamental_variance_mean1[i,]<-tapply(y_ppp,index_dep,var)
  
  #------------------------------------------------------------ #
  
  #Generation of the data scenario 2
  
  mean<-rep(beta2[i],n)+X_ijk%*%beta_E2[i,]+Z_jk%*%beta_M2[i,]+W_k%*%beta_D2[i,]+rep(phi2[i,],n_jk)
  v<-rep(sqrt(kappa2_jk2[i,]),n_jk)
  
  y_ppp<-rnorm(n,mean=mean,sd=v)
  
  #Calculus of the general statistic
  general_statistics2[i,]<-c(mean(y_ppp),sd(y_ppp)/abs(mean(y_ppp)),quantile(y_ppp,c(0.5)))
  
  municipal_mean2[i,]<-tapply(y_ppp,index,mean)
  variance_mean2[i,]<-tapply(y_ppp,index,var)
  departamental_mean2[i,]<-tapply(y_ppp,index_dep,mean)
  departamental_variance_mean2[i,]<-tapply(y_ppp,index_dep,var)
  
  #------------------------------------------------------------ #
  
  #Generation of the data scenario 3
  mean<-rep(beta3[i],n)+X_ijk%*%beta_E3[i,]+Z_jk%*%beta_M3[i,]+W_k%*%beta_D3[i,]+rep(phi3[i,],n_jk)
  v<-rep(sqrt(kappa2_jk3[i,]),n_jk)
  
  y_ppp<-rnorm(n,mean=mean,sd=v)
  
  #Calculus of the general statistic
  general_statistics3[i,]<-c(mean(y_ppp),sd(y_ppp)/abs(mean(y_ppp)),quantile(y_ppp,c(0.5)))
  
  municipal_mean3[i,]<-tapply(y_ppp,index,mean)
  variance_mean3[i,]<-tapply(y_ppp,index,var)
  departamental_mean3[i,]<-tapply(y_ppp,index_dep,mean)
  departamental_variance_mean3[i,]<-tapply(y_ppp,index_dep,var)
  
  
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
}


# PPP for municipal mean, and variance

municipal_PPP<-data.frame(
  PPP_mean_1=rep(NA,m),
  PPP_var_1=rep(NA,m),
  PPP_mean_2=rep(NA,m),
  PPP_var_2=rep(NA,m),
  PPP_mean_3=rep(NA,m),
  PPP_var_3=rep(NA,m)
)

departamental_PPP<-data.frame(
  PPP_mean_1=rep(NA,d),
  PPP_var_1=rep(NA,d),
  PPP_mean_2=rep(NA,d),
  PPP_var_2=rep(NA,d),
  PPP_mean_3=rep(NA,d),
  PPP_var_3=rep(NA,d)
)

for (i in c(1:m)) {
  municipal_PPP$PPP_mean_1[i]<-mean(municipal_mean1[,i]<= municipal_statistic$mean_s1[1])
  municipal_PPP$PPP_var_1[i]<-mean(variance_mean1[,i]<= municipal_statistic$var_s1[1])
}

for (i in c(1:m)) {
  municipal_PPP$PPP_mean_2[i]<-mean(municipal_mean2[,i]<=municipal_statistic$mean_s2[1])
  municipal_PPP$PPP_var_2[i]<-mean(variance_mean2[,i]<=municipal_statistic$var_s2[1])
}

for (i in c(1:m)) {
  municipal_PPP$PPP_mean_3[i]<-mean(municipal_mean3[,i]<= municipal_statistic$mean_s3[1])
  municipal_PPP$PPP_var_3[i]<-mean(variance_mean3[,i]<= municipal_statistic$var_s3[1])
}


for (i in c(1:d)) {
  departamental_PPP$PPP_mean_1[i]<-mean(departamental_mean1[,i]<= departamental_statistic$mean_s1[1])
  departamental_PPP$PPP_var_1[i]<-mean(departamental_variance_mean1[,i]<= departamental_statistic$var_s1[1])
}

for (i in c(1:d)) {
  departamental_PPP$PPP_mean_2[i]<-mean(departamental_mean2[,i]<=departamental_statistic$mean_s2[1])
  departamental_PPP$PPP_var_2[i]<-mean(departamental_variance_mean2[,i]<=departamental_statistic$var_s2[1])
}

for (i in c(1:d)) {
  departamental_PPP$PPP_mean_3[i]<-mean(departamental_mean3[,i]<= departamental_statistic$mean_s3[1])
  departamental_PPP$PPP_var_3[i]<-mean(departamental_variance_mean3[,i]<= departamental_statistic$var_s3[1])
}


#Global PPP

Global_PPP<-data.frame(
  Statistic=c("Mean","CV","Median"),
  scenario1=c(mean(general_statistics1[,1]<global_statistic$mean_s1),mean(general_statistics1[,2]<global_statistic$CV_1),mean(general_statistics1[,3]<global_statistic$median1)),
  scenario2=c(mean(general_statistics2[,1]<global_statistic$mean_s2),mean(general_statistics2[,2]<global_statistic$CV_2),mean(general_statistics2[,3]<global_statistic$median2)),
  scenario3=c(mean(general_statistics3[,1]<global_statistic$mean_s3),mean(general_statistics3[,2]<global_statistic$CV_3),mean(general_statistics3[,3]<global_statistic$median3))
)


municipal_PPP_boxplot_df<-data.frame(
  scenario=c(rep("1",m),rep("2",m),rep("3",m)),
  PPP_mean=c(municipal_PPP$PPP_mean_1,municipal_PPP$PPP_mean_2,municipal_PPP$PPP_mean_3),
  PPP_var=c(municipal_PPP$PPP_var_1,municipal_PPP$PPP_var_2,municipal_PPP$PPP_var_3)
)

departamental_PPP_boxplot_df<-data.frame(
  scenario=c(rep("1",d),rep("2",d),rep("3",d)),
  PPP_mean=c(departamental_PPP$PPP_mean_1,departamental_PPP$PPP_mean_2,departamental_PPP$PPP_mean_3),
  PPP_var=c(departamental_PPP$PPP_var_1,departamental_PPP$PPP_var_2,departamental_PPP$PPP_var_3)
)

Lasso_PPP<-list(municipal_PPP_boxplot_df=municipal_PPP_boxplot_df,departamental_PPP_boxplot_df=departamental_PPP_boxplot_df)
save(Lasso_PPP,file="Lasso_PPP.RData")

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
    name = "scenario"
  )+
  scale_fill_discrete_sequential(
    palette="BurgYl",
    rev=FALSE,
    name="scenario"
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


ggplot(departament_PPP_boxplot_df, aes(x = scenario, y = PPP_var, fill = scenario, colour = scenario)) +
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
    name = "scenario"
  )+
  scale_fill_discrete_sequential(
    palette="BurgYl",
    rev=FALSE,
    name="scenario"
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


###=======RMSE- MAE =======


RMSE_MAE_summary<-list()

RMSE_MAE_summary$betas<-data.frame(
  parameter=c("Student_level","Municipal_level","Departamental_level"),
  MAE_1=rep(NA,3),
  RMSE_1=rep(NA,3),
  MAE_2=rep(NA,3),
  RMSE_2=rep(NA,3),
  MAE_3=rep(NA,3),
  RMSE_3=rep(NA,3)
)

prueba<-matrix(NA,nrow = 3,ncol=1)

for(i in c(2:4)){
  
  MAE<-c(mean(abs(c(mean(beta1),apply(beta_E1,2,mean))-scenarios_parameter$student[,i])),mean(abs(apply(beta_M1,2,mean)-scenarios_parameter$Municipal[,i])),
         mean(abs(apply(beta_D1,2,mean)-scenarios_parameter$Departamental[,i])))
  RMSE<-sqrt(c(mean((c(mean(beta1),apply(beta_E1,2,mean))-scenarios_parameter$student[,i])^2),mean((apply(beta_M1,2,mean)-scenarios_parameter$Municipal[,i])^2),
               mean((apply(beta_D1,2,mean)-scenarios_parameter$Departamental[,i])^2)))
  
  prueba<-cbind(prueba,MAE,RMSE)
}

RMSE_MAE_summary$betas[,c(2:7)]<-prueba[,-1]
rm(prueba)

RMSE_MAE_summary$others<-data.frame(
  metric=c("MAE","RMSE"),
  municipal_variance=c(mean(abs(apply(kappa2_jk1,2,mean)-scenarios_parameter$kappa2_jk_post)),
                       sqrt(mean((apply(kappa2_jk1,2,mean)-scenarios_parameter$kappa2_jk_post)^2))),
  spatial_effect=c(mean(abs(apply(phi1,2,mean)-scenarios_parameter$phi_post)),
                   sqrt(mean((apply(phi1,2,mean)-scenarios_parameter$phi_post)^2)))
)

##=========DIC==========
###-------scenario 1---------

mean_est_bayes_1<-mean(beta1)+X_ijk%*%apply(beta_E1,2,stat_mode)+Z_jk%*%apply(beta_M1,2,stat_mode)+W_k%*%apply(beta_D1,2,stat_mode)+rep(apply(phi1,2,mean),n_jk)
var_est_bayes_1<-apply(kappa2_jk1,2,mean)

lp_1<-sum(dnorm(x=y_ijk_1,mean=mean_est_bayes_1,sd=sqrt(var_est_bayes_1)),log=TRUE)
p_DIC_1<-2*(lp_1-mean(ll_1))

DIC_1<--2*lp_1+2*p_DIC_1

###-------Scenario 2--------

mean_est_bayes_2<-mean(beta2)+X_ijk%*%apply(beta_E2,2,stat_mode)+Z_jk%*%apply(beta_M2,2,stat_mode)+W_k%*%apply(beta_D2,2,stat_mode)+rep(apply(phi2,2,mean),n_jk)
var_est_bayes_2<-apply(kappa2_jk2,2,mean)

lp_2<-sum(dnorm(x=y_ijk_2,mean=mean_est_bayes_2,sd=sqrt(var_est_bayes_2)),log=TRUE)
p_DIC_2<-2*(lp_2-mean(ll_2))

DIC_2<--2*lp_2+2*p_DIC_2

###------- Scenario 3 ---------

mean_est_bayes_3<-mean(beta3)+X_ijk%*%apply(beta_E3,2,stat_mode)+Z_jk%*%apply(beta_M3,2,stat_mode)+W_k%*%apply(beta_D3,2,stat_mode)+rep(apply(phi3,2,mean),n_jk)
var_est_bayes_3<-apply(kappa2_jk3,2,mean)

lp_3<-sum(dnorm(x=y_ijk_3,mean=mean_est_bayes_3,sd=sqrt(var_est_bayes_3)),log=TRUE)
p_DIC_3<-2*(lp_3-mean(ll_3))

DIC_3<--2*lp_3+2*p_DIC_3

##=====WAIC=========
index<-rep(seq_len(m),n_jk)

###-------Scenario 1---------
tmp1<-0
tmp2<-0
lppd_1<-0
pWAIC_1<-0

for (j in c(1:n)) {
  mean<-beta1+X_ijk[j,]%*%t(beta_E1)+Z_jk[j,]%*%t(beta_M1)+W_k[j,]%*%t(beta_D1)+phi1[,index[j]]
  v<-kappa2_jk1[,index[j]]
  
  tmp1<-dnorm(x=y_ijk_1[j],mean=mean,sd=sqrt(v))
  lppd_1<- lppd_1 + log(mean(tmp1))
  
  tmp2<-dnorm(x=y_ijk_1[j],mean=mean,sd=sqrt(v),log=TRUE)
  pWAIC_1<-pWAIC_1+ 2 * (log(mean(tmp1)) - mean(tmp2))
  
  # Advance
  if (j %% ceiling(n / 10) == 0) {
    cat(paste0("Iteración ", j, " de ", n, " (", round(100 * j / n), "%)\n"))
  }
}

WAIC_1 <- -2 * lppd_1 + 2 * pWAIC_1


###-------Scenario 2---------
tmp1<-0
tmp2<-0
lppd_2<-0
pWAIC_2<-0


for (j in c(1:n)) {
  mean<-beta2+X_ijk[j,]%*%t(beta_E2)+Z_jk[j,]%*%t(beta_M2)+W_k[j,]%*%t(beta_D2)+phi2[,index[j]]
  v<-kappa2_jk2[,index[j]]
  
  tmp1<-dnorm(x=y_ijk_2[j],mean=mean,sd=sqrt(v))
  lppd_2<- lppd_2 + log(mean(tmp1))
  
  tmp2<-dnorm(x=y_ijk_2[j],mean=mean,sd=sqrt(v),log=TRUE)
  pWAIC_2<-pWAIC_2+ 2 * (log(mean(tmp1)) - mean(tmp2))
  
  # Advance
  if (j %% ceiling(n / 10) == 0) {
    cat(paste0("Iteración ", j, " de ", n, " (", round(100 * j / n), "%)\n"))
  }
}

WAIC_2<- -2 * lppd_2 + 2 * pWAIC_2

###-------Scenario 3---------
tmp1<-0
tmp2<-0
lppd_3<-0
pWAIC_3<-0


for (j in c(1:n)) {
  mean<-beta3+X_ijk[j,]%*%t(beta_E3)+Z_jk[j,]%*%t(beta_M3)+W_k[j,]%*%t(beta_D3)+phi3[,index[j]]
  v<-kappa2_jk3[,index[j]]
  
  tmp1<-dnorm(x=y_ijk_3[j],mean=mean,sd=sqrt(v))
  lppd_3<- lppd_3 + log(mean(tmp1))
  
  tmp2<-dnorm(x=y_ijk_3[j],mean=mean,sd=sqrt(v),log=TRUE)
  pWAIC_3<-pWAIC_3+ 2 * (log(mean(tmp1)) - mean(tmp2))
  
  # Advance
  if (j %% ceiling(n / 10) == 0) {
    cat(paste0("Iteración ", j, " de ", n, " (", round(100 * j / n), "%)\n"))
  }
}

WAIC_3 <- -2 * lppd_3 + 2 * pWAIC_3

#Corrige DIC_WAIC

DIC_WAIC<-data.frame(
  Scenario=c("1","2","3"),
  lp=c(lp_1,lp_2,lp_3),
  pDIC=c(p_DIC_1,p_DIC_2,p_DIC_3),
  DIC=c(DIC_1,DIC_2,DIC_3),
  lppd=c(lppd_1,lppd_2,lppd_3),
  pWAIC=c(pWAIC_1,pWAIC_2,pWAIC_3),
  WAIC=c(WAIC_1,WAIC_2,WAIC_3)
)

save(DIC_WAIC,file="DIC_WAIC_lasso.RData")





set.seed(2008)

for(i in c(1:nsams)){
  
  
  #Generation of the data scenario 1
  mean<-rep(beta1[i],n)+X_ijk%*%beta_E1[i,]+Z_jk%*%beta_M1[i,]+W_k%*%beta_D1[i,]+rep(phi1[i,],n_jk)
  v<-rep(sqrt(kappa2_jk1[i,]),n_jk)
  
  y_ppp<-rnorm(n,mean=mean,sd=v)
  
  #Calculus of the general statistic
  #general_statistics1[i,]<-c(mean(y_ppp),sd(y_ppp)/abs(mean(y_ppp)),quantile(y_ppp,c(0.5)))
  
  #municipal_mean1[i,]<-tapply(y_ppp,index,mean)
  #variance_mean1[i,]<-tapply(y_ppp,index,var)
  departamental_mean1[i,]<-tapply(y_ppp,index_dep,mean)
  departamental_variance_mean1[i,]<-tapply(y_ppp,index_dep,var)
  
  #------------------------------------------------------------ #
  
  #Generation of the data scenario 2
  
  mean<-rep(beta2[i],n)+X_ijk%*%beta_E2[i,]+Z_jk%*%beta_M2[i,]+W_k%*%beta_D2[i,]+rep(phi2[i,],n_jk)
  v<-rep(sqrt(kappa2_jk2[i,]),n_jk)
  
  y_ppp<-rnorm(n,mean=mean,sd=v)
  
  #Calculus of the general statistic
  #general_statistics2[i,]<-c(mean(y_ppp),sd(y_ppp)/abs(mean(y_ppp)),quantile(y_ppp,c(0.5)))
  
  #municipal_mean2[i,]<-tapply(y_ppp,index,mean)
  #variance_mean2[i,]<-tapply(y_ppp,index,var)
  departamental_mean2[i,]<-tapply(y_ppp,index_dep,mean)
  departamental_variance_mean2[i,]<-tapply(y_ppp,index_dep,var)
  
  #------------------------------------------------------------ #
  
  #Generation of the data scenario 3
  mean<-rep(beta3[i],n)+X_ijk%*%beta_E3[i,]+Z_jk%*%beta_M3[i,]+W_k%*%beta_D3[i,]+rep(phi3[i,],n_jk)
  v<-rep(sqrt(kappa2_jk3[i,]),n_jk)
  
  y_ppp<-rnorm(n,mean=mean,sd=v)
  
  #Calculus of the general statistic
  #general_statistics3[i,]<-c(mean(y_ppp),sd(y_ppp)/abs(mean(y_ppp)),quantile(y_ppp,c(0.5)))
  
  #municipal_mean3[i,]<-tapply(y_ppp,index,mean)
  #variance_mean3[i,]<-tapply(y_ppp,index,var)
  departamental_mean3[i,]<-tapply(y_ppp,index_dep,mean)
  departamental_variance_mean3[i,]<-tapply(y_ppp,index_dep,var)
  
  
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
}

# CUANDO TERMINE, ARREGLA WAIC Y PPP