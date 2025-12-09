
#===============Libraries and preparation----

library(dplyr)
library(tidyverse)
library(stringi)
library(stringr)
library(mvtnorm)
library(tictoc)
library(Rcpp)
library(RcppArmadillo)

#Set the work directory to get de clean data
setwd("/Users/macbookpro/Desktop/Tesis/Procesed_data")

#These files correspond to the depurated data which is going to be used
datos <- read.csv('Datos_ProcesadosExamen_Saber_11_2022_2_Entrenamiento.csv')
datos <- datos[!(datos$estu_cod_reside_depto == 99999), ] %>%
  drop_na()

datos <- datos[!(datos$estu_cod_reside_depto == 68 & datos$estu_cod_reside_mcpio == 68264), ]
datos <- datos[!(datos$estu_cod_reside_depto == 88 ), ]
datos <- datos[!(datos$estu_cod_reside_mcpio==27493 ), ]
datos <- datos[, !names(datos) %in% c("homicidios", "codprovincia"
                                      , "ano")]%>%
  arrange(estu_cod_reside_depto,estu_cod_reside_mcpio)%>%
  na.omit(datos) 

table(datos$estu_depto_reside)

x_ijk_names <- c("fami_educacionmadre_modif", "computador", "internet", "etnia",
                 "libros_11_25", "libros_26_100", "libros_mas100", 
                 "estrato_1", "estrato_2", "estrato_3", "estrato_4", "estrato_5", "estrato_6",
                 "Genero_mujer", "Calendario_A", "Calendario_B", "cole_privado", "trabaja_menos_de_10_horas",  
                 "trabaja__11_a_20_horas", "trabaja__21_a_30_horas", "trabaja_mas_de_30_horas")

z_jk_names <-  c("docenttotal_alumtotal", "RISK_VICTIM_2022", "Homi_x_100k_habitantes",
                 "porc_alumn_col_publico", "terrorismot", "hurto", "secuestros", "discapital")

w_k_names <-  c("PIB_percapita_DPTO", "proporcio_pob_rural", "X..municipios.con.riesgo", 
                "Homicidios_ponderado_x_100k")

X_ijk <- as.matrix(subset(datos, select = x_ijk_names))
Z_jk  <- as.matrix(subset(datos, select = z_jk_names))
W_k   <- as.matrix(subset(datos, select = w_k_names))

Z_jk<-scale(Z_jk)
W_k<-scale(W_k)

load("Spatial_Component.RData")

#Extract the trainning info

#The response vector
y_ijk<-datos$punt_global
y_ijk_b<-datos[,c(3,5,16)]
y_ijk<-as.vector(y_ijk)
id_dep_mun_punt<-y_ijk_b[,c(1,2)]

#Full matrix

data_prior<-as.data.frame(cbind(X_ijk,Z_jk,W_k))


#Spatial information

# The number of each municipio's neighbors 
d_jk<-neighbors_municipios$d_jk

# A list with the index of every municipio's neighbor

W<-neighbors_municipios$neighbors


#==============Numbers=======================----
#Total of students
n<-length(y_ijk)
#Number of students from j municpio of the k departament 
n_jk<-datos%>%
  group_by(estu_cod_reside_mcpio)%>%
  summarise(n=n())%>%
  arrange(estu_cod_reside_mcpio)

n_jk<-n_jk$n

#Number of students from k departament 

n_k<-datos%>%
  group_by(estu_cod_reside_depto)%>%
  summarise(n=n())%>%
  arrange(estu_cod_reside_depto)

n_k<-n_k$n

#Number of municipios at every departament

m_k<-datos%>%
  group_by(estu_cod_reside_mcpio,estu_cod_reside_depto)%>%
  summarise(n=n())

m_k<-m_k%>%
  group_by(estu_cod_reside_depto)%>%
  summarise(m_k=n())
m_k<-m_k$m_k

#Number of municipios
m<-length(n_jk)

#Number of departaments

d<-length(m_k)



#Number of regression parameter
pe<-ncol(X_ijk)
pm<-ncol(Z_jk)
pd<-ncol(W_k)

data_ready<-list(data_complete=datos,punt_global=y_ijk,desing_matrix=data_prior,id_dep_mun=id_dep_mun_punt,X_ijk=X_ijk,Z_jk=Z_jk,W_k=W_k,pe=pe,pm=pm,pd=pd,n=n,n_jk=n_jk,n_k=n_k,m_k=m_k,m=m,d=d)

save(data_ready,file="data_ready.RData")

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



# --- To create necessary index vectors for C++ ----

# 1. Municipio Index Vector (mun_idx_vec): Maps each student to a municipio index (1 to m)
# Required for efficient rowsum in C++
mun_map_df <- datos %>%
  group_by(estu_cod_reside_mcpio) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(idx = seq_along(estu_cod_reside_mcpio)) %>%
  arrange(estu_cod_reside_mcpio) # Ensure correct ordering
mun_idx_vec <- rep(mun_map_df$idx, times = n_jk)

# 2. Kappa2_k Mapping Vector (kappa2_k_mun_map): Maps the department variance (kappa2_k) 
# to each of the 'm' municipalities. (This replaces rep(kappa2_k, m_k))
# NOTE: This requires kappa2_k to be initialized or ordered by department index (1:d).

# Create a mapping of department indices for each of the 'm' municipalities
dep_mun_map <- datos %>%
  group_by(estu_cod_reside_mcpio, estu_cod_reside_depto) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(estu_cod_reside_depto, estu_cod_reside_mcpio)

dep_map_k <- dep_mun_map %>%
  group_by(estu_cod_reside_depto) %>%
  mutate(dep_idx = cur_group_id()) %>% # Assigns index 1 to d to each department
  ungroup() %>%
  dplyr::select(dep_idx)


#Setting the work directory to 
# Load the C++ functions

setwd("/Users/macbookpro/Desktop/Tesis/Codes")

Rcpp::sourceCpp("Samplers.cpp") 

#
setwd("/Users/macbookpro/Desktop/Tesis/Procesed_Data")

#======================== Remaining R Sampler Functions ==============================

tauPhi_sampler<-function(phi,d_jk,nu_phi,gamma_phi,m,W){
  s<-sum((phi-(sapply(W, function(idx) sum(phi[idx]))/d_jk))^2)
  tauPhi<-1/rgamma(1,(nu_phi+m)/2,(nu_phi*gamma_phi+s)/2)
  return(tauPhi)
}

beta_sampler<-function(kappa2_jk,mu_beta,sigma2_beta,n_jk,y_ijk,data_prior,beta_E,beta_M,beta_D,m,phi){
  r_jk<-y_ijk-as.matrix(data_prior)%*%c(beta_E,beta_M,beta_D)-rep(phi,n_jk)
  R_jk<-rowsum(r_jk,rep(seq(1:m),n_jk))/kappa2_jk
  R<-sum(R_jk)
  
  a<-R+(mu_beta/sigma2_beta)
  b2<-1/(sum(n_jk/kappa2_jk)+(1/sigma2_beta))
  if (!is.finite(b2) || b2 <= 0) {
    b2 <- 1e-9
  }
  
  beta<-rnorm(1,a*b2,sqrt(b2))
  return(beta)
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

sigma2E_sampler<-function(nu_E,gamma_E,beta_E,mu_E,pe){
  sigma2E<-1/rgamma(1,(nu_E+pe)/2,(nu_E*gamma_E+sum((beta_E-mu_E)*(beta_E-mu_E)))/2)
  return(sigma2E)
}

sigma2M_sampler<-function(nu_M,gamma_M,beta_M,mu_M,pm){
  sigma2M<-1/rgamma(1,(nu_M+pm)/2,(nu_M*gamma_M+sum((beta_M-mu_M)*(beta_M-mu_M)))/2)
  return(sigma2M)
}

sigma2D_sampler<-function(nu_D,gamma_D,beta_D,mu_D,pd){
  sigma2D<-1/rgamma(1,(nu_D+pd)/2,(nu_D*gamma_D+sum((beta_D-mu_D)*(beta_D-mu_D)))/2)
  return(sigma2D)
}

# --- End of Preparation ---


MCMC<-function(B,nburn,jumps,y_ijk,y_ijk_b,data_prior,X_ijk,Z_jk,W_k,n,n_jk,m_k,m,pe,pm,pd,
               W,d_jk,nu_phi,gamma_phi,mu_beta,nu_beta,gamma_beta,mu_E,nu_E,gamma_E,
               mu_M,nu_M,gamma_M,mu_D,nu_D,gamma_D,nu_kappa,a_alpha_kappa,b_alpha_kappa,
               a_beta_kappa,b_beta_kappa,delta2,d, mun_idx_vec, dep_map_k){ # NEW ARGS ADDED
  
  #==============Chain settings
  # B: Iterations desired
  # nburn: Iterations percentage that will be burn i.e. (nburn=0.1)
  # jumps: The number of iterations that are going to be omitted 
  #-------------------------------
  burn_in<-B*nburn
  tot_it<-burn_in+B*jumps
  
  
  
  # ... (Initialization section) ...
  
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
  
  # -------------------------------------------------------------
  # Calculate kappa2_k_mun_map once after initialization
  # -------------------------------------------------------------
  kappa2_k_mun_map <- kappa2_k[dep_map_k$dep_idx]
  
  set.seed(123)
  
  ##Sigmas----
  
  sigma2_beta<-2
  sigma2_E<-2
  sigma2_M<-2
  sigma2_D<-2
  
  ##Betas----
  
  beta<-mu_beta
  beta_E<-mu_E
  beta_M<-mu_M
  beta_D<-mu_D
  
  ## alpha and beta sub kappa----
  
  acr<-0 #Includede to verify the chain acceptation rate
  alpha_kappa<-3
  beta_kappa<-3
  
  
  

  
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
  SIGMA2_E<-vector("numeric",length = B)
  SIGMA2_M<-vector("numeric",length = B)
  SIGMA2_D<-vector("numeric",length = B)
  
  PHI<-matrix(NA,nrow = B, ncol = m)
  TAU<-vector("numeric",length = B)
  
  KAPPA2_JK<-matrix(NA,nrow = B, ncol = m)
  
  KAPPA2_K<-matrix(NA,nrow = B, ncol = d)
  
  ALPHA_KAPPA<-vector("numeric",length = B)
  BETA_KAPPA<-vector("numeric",length = B)
  
  ACR<-NA
  
  #====================================================================================
  # Functions before the actualization
  #====================================================================================
  
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
  
  
  
  # --- Inside the loop: Actualization ---
  for (b in c(1:tot_it)) {
    
    # 1. Pre-calculate terms needed by C++ samplers
    
    # inv_kappa (for the N students, needed by Betas)
    inv_kappa_val <- 1/kappa2_jk
    inv_kappa_val <- ifelse(is.na(inv_kappa_val) | !is.finite(inv_kappa_val) | inv_kappa_val <= 0, 1e-9, inv_kappa_val)
    inv_kappa <- rep(inv_kappa_val, n_jk)
    
    # r_ijk_no_phi (Residual without the spatial component, needed by phi_sampler)
    r_ijk_no_phi <- y_ijk - rep(beta,times=n) - as.matrix(data_prior)%*%c(beta_E,beta_M,beta_D)
    
    # 2. PHI sampling (Optimized C++)
    phi <- phi_sampler_rcpp(phi, kappa2_jk, d_jk, tau_phi, n_jk, W, r_ijk_no_phi, mun_idx_vec)
    
    # 3. tauPhi_sampler (Simple, kept in R)
    tau_phi<-tauPhi_sampler(phi,d_jk,nu_phi,gamma_phi,m,W) 
    
    # 4.  BETA sampling (Kept in R as it's a simple Normal draw)
    beta <- beta_sampler(kappa2_jk,mu_beta,sigma2_beta,n_jk,y_ijk,data_prior,beta_E,beta_M,beta_D,m,phi)
    
    # 5.  BETAS sampling (Optimized C++)
    
    r_jk_E <- y_ijk - rep(beta,times=n) - as.matrix(cbind(Z_jk,W_k)) %*% c(beta_M,beta_D) - rep(phi,n_jk)
    beta_E <- betaE_sampler_rcpp(X_ijk, r_jk_E, inv_kappa, mu_E, sigma2_E, pe)
    
    r_jk_M <- y_ijk - rep(beta,times=n) - as.matrix(cbind(X_ijk,W_k)) %*% c(beta_E,beta_D) - rep(phi,n_jk)
    beta_M <- betaM_sampler_rcpp(Z_jk, r_jk_M, inv_kappa, mu_M, sigma2_M, pm)
    
    r_jk_D <- y_ijk - rep(beta,times=n) - as.matrix(cbind(X_ijk,Z_jk)) %*% c(beta_E,beta_M) - rep(phi,n_jk)
    beta_D <- betaD_sampler_rcpp(W_k, r_jk_D, inv_kappa, mu_D, sigma2_D, pd)
    
    # 6. Muestreo de KAPPA2_JK (Optimized C++)
    # r_jk needs the full residual
    r_jk <- y_ijk - (rep(beta,times=n)+as.matrix(data_prior)%*%c(beta_E,beta_M,beta_D)+rep(phi,n_jk))
    
    kappa2_jk<-kappa2jk_sampler_rcpp(nu_kappa, n_jk, r_jk, kappa2_k_mun_map, mun_idx_vec, m)
    kappa2_jk<-ifelse(is.na(kappa2_jk) | !is.finite(kappa2_jk) | kappa2_jk <= 0, 1e-9, kappa2_jk)
    
    # 7. Muestreo de KAPPA2_K (Kept in R for now as it involves complex rowsum aggregation)
    kappa2_k<-kappa2k_sampler(alpha_kappa,nu_kappa,beta_kappa,kappa2_jk,d,m_k)
    
    sigma2_beta<-sigma2Beta_sampler(nu_beta,gamma_beta,beta,mu_beta)
    sigma2_E<-sigma2E_sampler(nu_E,gamma_E,beta_E,mu_E,pe)
    sigma2_M<-sigma2M_sampler(nu_M,gamma_M,beta_M,mu_M,pm)
    sigma2_D<-sigma2D_sampler(nu_D,gamma_D,beta_D,mu_D,pd)
    
    
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
      SIGMA2_E[it_stored]<-sigma2_E
      SIGMA2_M[it_stored]<-sigma2_M
      SIGMA2_D[it_stored]<-sigma2_D
      
      
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
  
  THETA<-list(BETA,BETA_E,BETA_M,BETA_D,SIGMA2_BETA,SIGMA2_E,SIGMA2_M,SIGMA2_D,PHI,TAU,KAPPA2_JK,KAPPA2_K,ALPHA_KAPPA,BETA_KAPPA,ACR)
  
  return(THETA)
}        


#================ Running the Gibbs Sampler =================================

tictoc::tic()

THETA<-MCMC(B=15000,nburn=0.1,jumps=5,y_ijk,y_ijk_b,data_prior,X_ijk,Z_jk,W_k,n,n_jk,m_k,m,pe,pm,pd,
            W,d_jk,nu_phi,gamma_phi,mu_beta,nu_beta,gamma_beta,mu_E,nu_E,gamma_E,
            mu_M,nu_M,gamma_M,mu_D,nu_D,gamma_D,nu_kappa,a_alpha_kappa,b_alpha_kappa,
            a_beta_kappa,b_beta_kappa,delta2=1,d, mun_idx_vec, dep_map_k)

tiempo = tictoc::toc()


save(THETA,file = "/Users/macbookpro/Desktop/Tesis/Results/MCMC_model2_standarized.RData")





