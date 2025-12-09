#=======================================================================================
#. Description:  
#. These script corresponds to the 1st version of Gibbs Sampler from the Ridge model
#=======================================================================================


#===============Libraries and preparation----

library(dplyr)
library(stringi)
library(stringr)
library(mvtnorm)
library(tictoc)


#=== Settings and preparation====

### EN LA SIGUIENTE LÍNEA, PON LA RUTA DE DONDE VAS A DEJAR LOS ARCHIVOS 
#.  DESCARGADOS 

setwd("/Users/macbookpro/Desktop/Tesis/Procesed_Data")
load("data_ready.RData")
load("Spatial_Component.RData")

#====Variable and covariables====
datos<-data_ready$data_complete
y_ijk<-data_ready$punt_global
y_ijk_b<-cbind(y_ijk,data_ready$id_dep_mun)
names(y_ijk_b)<-c("punt_global","estu_cod_reside_depto" ,"estu_cod_reside_mcpio")
data_prior<-data_ready$desing_matrix
X_ijk<-data_ready$X_ijk
Z_jk<-data_ready$Z_jk
W_k<-data_ready$W_k

W<-neighbors_municipios$neighbors
d_jk<-neighbors_municipios$d_jk

#==== Numbers

pe<-data_ready$pe
pm<-data_ready$pm
pd<-data_ready$pd

n<-data_ready$n
n_jk<-data_ready$n_jk
n_k<-data_ready$n_k
m_k<-data_ready$m_k
m<-data_ready$m
d<-data_ready$d


#========== Hiperparámeters ==========----

nu_phi<-2
gamma_phi<-1

nu_beta<- 2
nu_E<-2
nu_M<-2
nu_D<-2
gamma_beta <-1
gamma_E <- 1
gamma_M <-1
gamma_D <- 1

nu_kappa <-2
a_alpha_kappa <- 2
b_alpha_kappa <- 5
a_beta_kappa <- 2
b_beta_kappa <- 5




##========= Unitary prior ===============----

u_prior <- lm(y_ijk ~ ., data = data_prior)

#Extract the beta_ols to set them as the mu hiperparameters 

mu_beta<-u_prior$coefficients[[1]]
mu_E <-u_prior$coefficients[seq(2,pe+1)]
mu_M <-u_prior$coefficients[seq(pe+2,pe+pm+1)]
mu_D <-u_prior$coefficients[seq(pe+pm+2,34)]








#======================== Actualization functions ==============================




phi_sampler<-function(phi,kappa2_jk,d_jk,tau_phi,n_jk,W,y_ijk,beta,n,data_prior,beta_E,beta_M,beta_D,m){
  r_jk<-y_ijk-rep(beta,times=n)-as.matrix(data_prior)%*%c(beta_E,beta_M,beta_D)
  R_jk<-rowsum(r_jk,rep(seq(1:m),n_jk))/kappa2_jk  
  
  b2<-d_jk/tau_phi+n_jk/kappa2_jk 
  a <- (sapply(W, function(idx) sum(phi[idx]))/tau_phi)+R_jk
  
  phi<-rnorm(a/b2,sqrt(1/b2))
  phi<-phi-mean(phi)
  
  return(phi)
}


beta_sampler<-function(kappa2_jk,mu_beta,sigma2_beta,n_jk,y_ijk,data_prior,beta_E,beta_M,beta_D,m,phi){
  
  r_jk<-y_ijk-as.matrix(data_prior)%*%c(beta_E,beta_M,beta_D)-rep(phi,n_jk)
  R_jk<-rowsum(r_jk,rep(seq(1:m),n_jk))/kappa2_jk 
  R<-sum(R_jk)
  #R<-ifelse(R>500000,749620,r_jk_groupped)
  
  a<-R+(mu_beta/sigma2_beta)
  
  b2<-1/(sum(n_jk/kappa2_jk)+(1/sigma2_beta))
  if (!is.finite(b2) || b2 <= 0) {
    b2 <- 1e-9
  }
  
  beta<-rnorm(1,a*b2,sqrt(b2))
  return(beta)
}

betaE_sampler<-function(mu_E,lambda2_E,kappa2_jk,n_jk,pe,y_ijk,beta,n,X_ijk,Z_jk,W_k,beta_M,beta_D,phi){
  
  r_jk_E<-y_ijk-rep(beta,times=n)-as.matrix(cbind(Z_jk,W_k))%*%c(beta_M,beta_D)-rep(phi,n_jk)
  
  inv_kappa<-rep(1/kappa2_jk,n_jk)
  inv_kappa<-ifelse(is.na(inv_kappa) | !is.finite(inv_kappa) | inv_kappa <= 0, 1e-9, inv_kappa)
  
  a<-t(X_ijk)%*%(r_jk_E*inv_kappa)
  b<-(diag(rep(lambda2_E,pe))+t(X_ijk)%*%(X_ijk*inv_kappa))
  
  mu<-as.vector(solve(b)%*%a)
  Sigma<-solve(b)
  
  betaE<-c(mvtnorm::rmvnorm(1,mu,Sigma))
  return(betaE)
}


betaM_sampler<-function(mu_M,lambda2_M,kappa2_jk,n_jk,pm,y_ijk,beta,n,X_ijk,Z_jk,W_k,beta_E,beta_D,phi){
  
  r_jk_M<-y_ijk-rep(beta,times=n)-as.matrix(cbind(X_ijk,W_k))%*%c(beta_E,beta_D)-rep(phi,n_jk)
  
  inv_kappa<-rep(1/kappa2_jk,n_jk)
  inv_kappa<-ifelse(is.na(inv_kappa) | !is.finite(inv_kappa) | inv_kappa <= 0, 1e-9, inv_kappa)
  
  a<-t(Z_jk)%*%(r_jk_M*inv_kappa)
  b<-(diag(rep(lambda2_M,pm))+t(Z_jk)%*%(Z_jk*inv_kappa))
  
  mu<-as.vector(solve(b)%*%a)
  Sigma<-solve(b)
  
  betaM<-c(mvtnorm::rmvnorm(1,mu,Sigma)) 
  
  return(betaM)
}



betaD_sampler<-function(mu_D,lambda2_D,kappa2_jk,n_jk,pd,y_ijk,beta,n,X_ijk,Z_jk,W_k,beta_E,beta_M,phi){
  
  r_jk_D<-y_ijk-rep(beta,times=n)-as.matrix(cbind(X_ijk,Z_jk))%*%c(beta_E,beta_M)-rep(phi,n_jk)
  
  inv_kappa<-rep(1/kappa2_jk,n_jk)
  inv_kappa<-ifelse(is.na(inv_kappa) | !is.finite(inv_kappa) | inv_kappa <= 0, 1e-9, inv_kappa)
  
  a<-t(W_k)%*%(r_jk_D*(inv_kappa))
  b<-(diag(rep(lambda2_D,pd))+t(W_k)%*%(W_k*inv_kappa))
  
  mu<-as.vector(solve(b)%*%a)
  Sigma<-solve(b)
  
  betaD<-c(mvtnorm::rmvnorm(1,mu,Sigma))
  return(betaD)
}


tauPhi_sampler<-function(phi,d_jk,nu_phi,gamma_phi,m,W){
  
  s<-sum((phi-(sapply(W, function(idx) sum(phi[idx]))/d_jk))^2)
  
  tauPhi<-1/rgamma(1,(nu_phi+m)/2,(nu_phi*gamma_phi+s)/2)
  
  return(tauPhi)
}


kappa2jk_sampler<-function(nu_kappa,n_jk,r_jk,kappa2_k,m_k,m){
  kappa2jk<-rep(NA,m)
  r_jk_groupped<-rowsum((r_jk*r_jk),rep(seq(1:m),n_jk))
  #r_jk_groupped<-ifelse(r_jk_groupped>50000000,429423,r_jk_groupped)
  a_kappa2_jk<-(nu_kappa+n_jk)/2
  b_kappa2_jk<-(nu_kappa*rep(kappa2_k,m_k)+r_jk_groupped)/2
  
  for (idx in seq_len(m)) {
    
    kappa2jk[idx]<-pmax(1e-9, 1 / rgamma(1,a_kappa2_jk[idx], b_kappa2_jk[idx]))
    
  }
  
  return(kappa2jk)
}

kappa2k_sampler<-function(alpha_kappa,nu_kappa,beta_kappa,kappa2_jk,d,m_k){
  kappa2k<-rep(NA,d)
  kappa2_jk_in_grouped<-rowsum((1/kappa2_jk),rep(seq_len(d),m_k))
  a_kappa2_k<-(alpha_kappa+nu_kappa*m_k)/2
  b_kappa2_k<-(beta_kappa+nu_kappa*kappa2_jk_in_grouped)/2
  for (idx in seq_len(d)) {
    
    kappa2k[idx]<-pmax(1e-9,rgamma(1,a_kappa2_k[idx], b_kappa2_k[idx]))
    
  }
  
  
  return(kappa2k)
}



sigma2Beta_sampler<-function(nu_beta,gamma_beta,beta,mu_beta){
  
  sigma2Beta<-1/rgamma(1,(nu_beta+1)/2,(nu_beta*gamma_beta+(beta-mu_beta)^2)/2)
  return(sigma2Beta)
}


lambda2E_sampler<-function(nu_E,gamma_E,beta_E,mu_E,pe){
  
  lambda2E<-rgamma(1,(nu_E+pe)/2,(nu_E*gamma_E+sum(beta_E*beta_E))/2)
  return(lambda2E)
}


lambda2M_sampler<-function(nu_M,gamma_M,beta_M,mu_M,pm){
  
  lambda2M<-rgamma(1,(nu_M+pm)/2,(nu_M*gamma_M+sum(beta_M*beta_M))/2)
  return(lambda2M)
}

lambda2D_sampler<-function(nu_D,gamma_D,beta_D,mu_D,pd){
  
  lambda2D<-rgamma(1,(nu_D+pd)/2,(nu_D*gamma_D+sum(beta_D*beta_D))/2)
  return(lambda2D)
}







MCMC<-function(B,nburn,jumps,y_ijk,y_ijk_b,data_prior,X_ijk,Z_jk,W_k,n,n_jk,m_k,m,pe,pm,pd,
               W,d_jk,nu_phi,gamma_phi,mu_beta,nu_beta,gamma_beta,mu_E,nu_E,gamma_E,
               mu_M,nu_M,gamma_M,mu_D,nu_D,gamma_D,nu_kappa,a_alpha_kappa,b_alpha_kappa,
               a_beta_kappa,b_beta_kappa,delta2,d){
  
  #==============Chain settings
  # B: Iterations desired
  # nburn: Iterations percentage that will be burn i.e. (nburn=0.1)
  # jumps: The number of iterations that are going to be omitted 
  #-------------------------------
  burn_in<-B*nburn
  tot_it<-burn_in+B*jumps
  
  #===================
  
  #======== Parameter Initialization ----
  
  set.seed(123)
  
  ##Sigmas----
  
  sigma2_beta<-2
  lambda2_E<-1/2
  lambda2_M<-1/2
  lambda2_D<-1/2
  
  ##Betas----
  
  beta<-mu_beta
  beta_E<-mu_E
  beta_M<-mu_M
  beta_D<-mu_D
  
  ## alpha and beta sub kappa----
  
  acr<-0 #Includede to verify the chain acceptation rate
  alpha_kappa<-3
  beta_kappa<-3
  
  ## Kappas----
  #Initialization of kappa's was made up from sufficient statistics
  
  kappa2_k<-y_ijk_b%>%
    group_by(estu_cod_reside_depto)%>%
    summarise(s2=var(punt_global))
  
  kappa2_k<-kappa2_k$s2
  
  kappa2_jk<-y_ijk_b%>%
    group_by(estu_cod_reside_mcpio)%>%
    summarise(s2=var(punt_global,na.rm = TRUE))%>%
    mutate(s2 = case_when(
      is.na(s2) ~ 0.1,
      TRUE ~ s2
    ))
  
  kappa2_jk<-kappa2_jk$s2
  
  
  ##phis----
  #The initialization of phi's accomplish the restriction (sum (phi's)=0)
  
  phi<-rep(0,times=m)
  
  ##Tau_phi----
  
  tau_phi<-1/rgamma(1,nu_phi/2,nu_phi*gamma_phi/2)
  
  
  
  #==================  Creating the space for the chain results============================
  
  BETA<-vector("numeric",length = B)
  
  BETA_E<-matrix(NA,nrow = B,ncol = pe)
  colnames(BETA_E)<-colnames(X_ijk)
  
  BETA_M<-matrix(NA,nrow = B,ncol = pm)
  colnames(BETA_M)<-colnames(Z_jk)
  
  BETA_D<-matrix(NA,nrow = B,ncol = pd)
  colnames(BETA_D)<-colnames(W_k)
  
  SIGMA2_BETA<-vector("numeric",length = B)
  LAMBDA2_E<-vector("numeric",length = B)
  LAMBDA2_M<-vector("numeric",length = B)
  LAMBDA2_D<-vector("numeric",length = B)
  
  PHI<-matrix(NA,nrow = B, ncol = m)
  TAU<-vector("numeric",length = B)
  
  KAPPA2_JK<-matrix(NA,nrow = B, ncol = m)
  
  KAPPA2_K<-matrix(NA,nrow = B, ncol = d)
  
  ALPHA_KAPPA<-vector("numeric",length = B)
  BETA_KAPPA<-vector("numeric",length = B)
  
  ACR<-NA
  
  #====================================================================================
  
  
  #=============== Some counts and functions before actualization================================
  
  r_jk<-y_ijk-(rep(beta,times=n)+as.matrix(data_prior)%*%c(beta_E,beta_M,beta_D)+rep(phi,n_jk))
  
  p_cond_complete<-function(x){
    p<-d*(x/2*log(beta_kappa)-log(gamma(x/2)))+(x/2)*sum(log(kappa2_k))+(a_alpha_kappa-1)*log(x)-b_alpha_kappa*x
    return(p)
  }
  
  alphaKappa_sampler<-function(delta2,alpha_kappa,beta_kappa,a_alpha_kappa,b_alpha_kappa,kappa2_k,acr){
    
    # Propossal
    alpha_prop<-rnorm(1,alpha_kappa,sqrt(delta2))
    alpha_prop<-ifelse(alpha_prop<=0,abs(alpha_prop),alpha_prop)
    
    #Log Acceptation rate
    log_r<-p_cond_complete(alpha_prop)+dnorm(alpha_kappa,alpha_prop,sqrt(delta2),log = TRUE)-p_cond_complete(alpha_kappa)-dnorm(alpha_prop,alpha_kappa,sqrt(delta2),log = TRUE)
    
    #Actualización
    
    if(runif(1)< exp(log_r)){
      
      alpha_kappa<-alpha_prop
      acr<-acr+1/B #To calculate acceptation rate
    }
    
    alpha<-list(alpha_kappa,acr)
    return(alpha)
  }
  
  
  betaKappa_sampler<-function(alpha_kappa,beta_kappa,a_beta_kappa,b_beta_kappa,kappa2_k,d){
    
    betaKappa<-rgamma(1,d*alpha_kappa/2+a_beta_kappa,b_beta_kappa+(sum(kappa2_k)/2))
    
    return(betaKappa)
  }
  
  
  #================ Actualization ====================================================----
  set.seed(123)
  
  for (b in c(1:tot_it)) {
    
    phi<-phi_sampler(phi,kappa2_jk,d_jk,tau_phi,n_jk,W,y_ijk,beta,n,data_prior,beta_E,beta_M,beta_D,m)
    tau_phi<-tauPhi_sampler(phi,d_jk,nu_phi,gamma_phi,m,W)
    
    beta<-beta_sampler(kappa2_jk,mu_beta,sigma2_beta,n_jk,y_ijk,data_prior,beta_E,beta_M,beta_D,m,phi)
    beta_E<-betaE_sampler(mu_E,lambda2_E,kappa2_jk,n_jk,pe,y_ijk,beta,n,X_ijk,Z_jk,W_k,beta_M,beta_D,phi)
    beta_M<-betaM_sampler(mu_M,lambda2_M,kappa2_jk,n_jk,pm,y_ijk,beta,n,X_ijk,Z_jk,W_k,beta_E,beta_D,phi)
    beta_D<-betaD_sampler(mu_D,lambda2_D,kappa2_jk,n_jk,pd,y_ijk,beta,n,X_ijk,Z_jk,W_k,beta_E,beta_M,phi)
    
    kappa2_jk<-kappa2jk_sampler(nu_kappa,n_jk,r_jk,kappa2_k,m_k,m)
    kappa2_jk<-ifelse(is.na(kappa2_jk) | !is.finite(kappa2_jk) | kappa2_jk <= 0, 1e-9, kappa2_jk)
    kappa2_k<-kappa2k_sampler(alpha_kappa,nu_kappa,beta_kappa,kappa2_jk,d,m_k)
    
    sigma2_beta<-sigma2Beta_sampler(nu_beta,gamma_beta,beta,mu_beta)
    lambda2_E<-lambda2E_sampler(nu_E,gamma_E,beta_E,mu_E,pe)
    lambda2_M<-lambda2M_sampler(nu_M,gamma_M,beta_M,mu_M,pm)
    lambda2_D<-lambda2D_sampler(nu_D,gamma_D,beta_D,mu_D,pd)
    
    
    alpha_kappa_res<-alphaKappa_sampler(delta2,alpha_kappa,beta_kappa,a_alpha_kappa,b_alpha_kappa,kappa2_k,acr)
    alpha_kappa<-alpha_kappa_res[[1]]
    
    acr<-alpha_kappa_res[[2]]
    
    beta_kappa<-betaKappa_sampler(alpha_kappa,beta_kappa,a_beta_kappa,b_beta_kappa,kappa2_k,d)
    
    #Actualization of 
    r_jk<-y_ijk-(rep(beta,times=n)+as.matrix(data_prior)%*%c(beta_E,beta_M,beta_D)+rep(phi,n_jk))
    
    
    #=================== To store ===============================================
    
    
    if(b>burn_in && b %% jumps ==0 ){
      
      #Iterations to be stored
      it_stored<-(b-burn_in)/jumps
      
      #Betas
      BETA[it_stored]<-beta
      BETA_E[it_stored,]<-beta_E
      BETA_M[it_stored,]<-beta_M
      BETA_D[it_stored,]<-beta_D
      
      #Sigmas
      
      SIGMA2_BETA[it_stored]<-sigma2_beta
      LAMBDA2_E[it_stored]<-lambda2_E
      LAMBDA2_M[it_stored]<-lambda2_M
      LAMBDA2_D[it_stored]<-lambda2_D
      
      
      #Phi
      PHI[it_stored,]<-phi
      
      #Tau_phi
      TAU[it_stored]<-tau_phi
      
      #Kappas
      KAPPA2_JK[it_stored,]<-kappa2_jk
      KAPPA2_K[it_stored,]<-kappa2_k
      
      #Alpha_Kappa Beta_Kappa
      ALPHA_KAPPA[it_stored]<-alpha_kappa
      ACR<-acr
      BETA_KAPPA[it_stored]<-beta_kappa
      
      
    }
    
    # Advance
    if (b %% ceiling(tot_it / 10) == 0) {
      cat(paste0("Iteración ", b, " de ", tot_it, " (", round(100 * b / tot_it), "%)\n"))
    }
    
    
  }
  
  THETA<-list(BETA,BETA_E,BETA_M,BETA_D,SIGMA2_BETA,LAMBDA2_E,LAMBDA2_M,LAMBDA2_D,PHI,TAU,KAPPA2_JK,KAPPA2_K,ALPHA_KAPPA,BETA_KAPPA,ACR)
  
  return(THETA)
}        


#================ Running the Gibbs Sampler =================================

tictoc::tic()

THETA<-MCMC(B=25000,nburn=0.1,jumps=5,y_ijk,y_ijk_b,data_prior,X_ijk,Z_jk,W_k,n,n_jk,m_k,m,pe,pm,pd,
            W,d_jk,nu_phi,gamma_phi,mu_beta,nu_beta,gamma_beta,mu_E,nu_E,gamma_E,
            mu_M,nu_M,gamma_M,mu_D,nu_D,gamma_D,nu_kappa,a_alpha_kappa,b_alpha_kappa,
            a_beta_kappa,b_beta_kappa,delta2=1,d)

tiempo = tictoc::toc()

####### 
# EL ARCHIVO VA A QUEDAR GUARDADO EN LA RUTA QUE ESPECIFICASTE AL PRINCIPIO

save(THETA,file = "MCMC_ridge_real_data_non_standarized.RData")

