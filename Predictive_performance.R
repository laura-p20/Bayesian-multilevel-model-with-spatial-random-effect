
#Set the work directory to get de clean data
setwd("/Users/macbookpro/Desktop/Tesis/Procesed_data")
load("/Users/macbookpro/Desktop/Tesis/Procesed_Data/data_ready.RData")

#Municipios whose info was available in the training process
municipios_used<-data_ready$data_complete

#These files correspond to the depurated data which is going to be used
datos <- read.csv("~/Desktop/Tesis/Procesed_Data/Datos_ProcesadosExamen_Saber_11_2022_2_Prueba.csv")
datos <- datos[!(datos$estu_cod_reside_depto == 99999), ] %>%
  drop_na()

datos <- datos[!(datos$estu_cod_reside_depto == 68 & datos$estu_cod_reside_mcpio == 68264), ]
datos <- datos[!(datos$estu_cod_reside_depto == 88 ), ]
datos <- datos[!(datos$estu_cod_reside_mcpio==27493 ), ]
datos <- datos[, !names(datos) %in% c("homicidios", "codprovincia"
                                      , "ano")]%>%
  arrange(estu_cod_reside_depto,estu_cod_reside_mcpio)%>%
  na.omit(datos) 

datos <- datos %>% 
  semi_join(municipios_used, by = "estu_cod_reside_mcpio")
  
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
Z_jk <- scale(Z_jk)
W_k   <- as.matrix(subset(datos, select = w_k_names))
W_k<-scale(W_k)

#The response vector
y_ijk<-datos$punt_global
y_ijk_b<-datos[,c(3,5,16)]
y_ijk<-as.vector(y_ijk)


municipios_used<-municipios_used%>%
  group_by(estu_cod_reside_mcpio,estu_mcpio_reside)%>%
  summarise(n=n())%>%
  arrange()

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


##===========Extracting the samples of each model----
#Basis model 
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

#=== Prediction=====

y_pred_basic<-rep(0,n)
y_pred_ridge<-rep(0,n)
y_pred_lasso<-rep(0,n)



set.seed(2008)

for(i in c(1:nsams)){
  
  #Generation of the data with basic
  mean<-rep(beta[i],n)+X_ijk%*%beta_E[i,]+Z_jk%*%beta_M[i,]+W_k%*%beta_D[i,]+rep(phi[i,],n_jk)
  v<-rep(sqrt(kappa2_jk[i,]),n_jk)
  y_ppp<-rnorm(n,mean=mean,sd=v)
  
  y_pred_basic<-y_pred_basic+y_ppp*(1/nsams)
  
  
  #------------------------------------------------------------ #
  
  #Generation of the data lasso
  
  mean<-rep(betal[i],n)+X_ijk%*%beta_El[i,]+Z_jk%*%beta_Ml[i,]+W_k%*%beta_Dl[i,]+rep(phil[i,],n_jk)
  v<-rep(sqrt(kappa2_jkl[i,]),n_jk)
  
  y_ppp<-rnorm(n,mean=mean,sd=v)
  
  y_pred_ridge<-y_pred_ridge+y_ppp*(1/nsams)
  
  #------------------------------------------------------------ #
  
  #Generation of the data ridge
  mean<-rep(betar[i],n)+X_ijk%*%beta_Er[i,]+Z_jk%*%beta_Mr[i,]+W_k%*%beta_Dr[i,]+rep(phir[i,],n_jk)
  v<-rep(sqrt(kappa2_jkr[i,]),n_jk)
  
  y_ppp<-rnorm(n,mean=mean,sd=v)
  y_pred_lasso<-y_pred_lasso+y_ppp*(1/nsams)
  
  
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
}



predictive_perfomance<-data.frame(
  model=c("Basic","Ridge","Lasso"),
  MSE=c((1/n)*sum((y_pred_basic-y_ijk)^2),(1/n)*sum((y_pred_ridge-y_ijk)^2),(1/n)*sum((y_pred_lasso-y_ijk)^2)),
  MAE=c((1/n)*sum(abs(y_pred_basic-y_ijk)),(1/n)*sum(abs(y_pred_ridge-y_ijk)),(1/n)*sum(abs(y_pred_lasso-y_ijk))),
  R2=c(1-(sum(abs(y_pred_basic-y_ijk)^2)/(var(y_ijk)*(n-1))),1-(sum(abs(y_pred_ridge-y_ijk)^2)/(var(y_ijk)*(n-1))),1-(sum(abs(y_pred_lasso-y_ijk)^2)/(var(y_ijk)*(n-1))))
  
)


set.seed(230125)
MSE<-0
MAE<-0
R2<-0

for(i in c(1:100)){
  
  y_sample_id<-sample(seq_len(n),n*0.1)
  
  
  
  MSE_0=c((1/n)*sum((y_pred_basic[y_sample_id]-y_ijk[y_sample_id])^2),(1/n)*sum((y_pred_ridge[y_sample_id]-y_ijk[y_sample_id])^2),(1/n)*sum((y_pred_lasso[y_sample_id]-y_ijk[y_sample_id])^2))
  MAE_0=c((1/n)*sum(abs(y_pred_basic[y_sample_id]-y_ijk[y_sample_id])),(1/n)*sum(abs(y_pred_ridge[y_sample_id]-y_ijk[y_sample_id])),(1/n)*sum(abs(y_pred_lasso[y_sample_id]-y_ijk[y_sample_id])))
  R2_0=c(1-(sum(abs(y_pred_basic[y_sample_id]-y_ijk[y_sample_id])^2)/(var(y_ijk[y_sample_id])*(n-1))),1-(sum(abs(y_pred_ridge[y_sample_id]-y_ijk[y_sample_id])^2)/(var(y_ijk[y_sample_id])*(n-1))),1-(sum(abs(y_pred_lasso[y_sample_id]-y_ijk[y_sample_id])^2)/(var(y_ijk[y_sample_id])*(n-1))))
  
  MSE<-MSE+MSE_0/100
  MAE<-MAE+MAE_0/100
  R2<-R2+R2_0/100
  
  # Advance
  if (i %% ceiling(100 / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / 100), "%)\n"))
  }
  
}


predictive_perfomance_boostrap<-data.frame(
  model=c("Basic","Ridge","Lasso"),
  MSE=MSE,
  MAE=MAE,
  R2=R2
)

MAE_MSE_R2<-list(y_pred=y_pred,predictive_perfomance=predictive_perfomance,predictive_perfomance_boostrap=predictive_perfomance_boostrap)
save(MAE_MSE_R2,file = "Predictive_posterior_performance(MSE,MAE,R2).RData")


