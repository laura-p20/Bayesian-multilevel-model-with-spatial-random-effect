#=======================================================================================
#. Description:  
#. These script corresponds to the analysis done after the Gibbs Sampler
#=======================================================================================


#===============Libraries and preparation----

library(dplyr)
library(tidyverse)
library(ggplot2)
library(sf)
library(spdep)
library(dplyr)
library(stringi)
library(stringr)

#Set the work directory to get a previus result of the MCMC using the real data

setwd("/Users/macbookpro/Desktop/Tesis/Results")
load("MCMC_model2_standarized.RData")

#============= Functions -----

CV<-function(x){
  cv<-sd(x)/abs(mean(x))*100
  cat("CV: ",cv,"\n","\n")
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


#===========Extracting the parameter samples----

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

acr<-THETA[[15]]


#=========== Posterior analysis of the covariables----


#Creating a data frame containing the basic posterior inference about the covariables parameters

df <- data.frame(
  variable = c("beta",colnames(beta_E),colnames(beta_M),colnames(beta_D)),
  media = c(mean(beta),apply(beta_E,2,mean),apply(beta_M,2,mean),apply(beta_D,2,mean)),
  lim_inf=c(bayesian_summary(beta)[2],apply(beta_E, 2, bayesian_summary)[2,],apply(beta_M, 2, bayesian_summary)[2,],apply(beta_D, 2, bayesian_summary)[2,]),
  lim_sup=c(bayesian_summary(beta)[3],apply(beta_E, 2, bayesian_summary)[3,],apply(beta_M, 2, bayesian_summary)[3,],apply(beta_D, 2, bayesian_summary)[3,]),
  Nivel=c("Intercepto", rep("Estudiante",ncol(beta_E)),rep("Municipal",ncol(beta_M)),rep("Departamental",ncol(beta_D)))
)


df$variable <- recode(df$variable,
                      "beta"                        = "Intercepto",
                      "fami_educacionmadre_modif"   = "Educación madre",
                      "computador"                  = "Computador en casa",
                      "internet"                    = "Acceso a internet",
                      "etnia"                       = "Pertenencia étnica",
                      "libros_11_25"                = "Libros (11–25)",
                      "libros_26_100"               = "Libros (26–100)",
                      "libros_mas100"               = "Libros (>100)",
                      "estrato_1"                   = "Estrato 1",
                      "estrato_2"                   = "Estrato 2",
                      "estrato_3"                   = "Estrato 3",
                      "estrato_4"                   = "Estrato 4",
                      "estrato_5"                   = "Estrato 5",
                      "estrato_6"                   = "Estrato 6",
                      "Genero_mujer"                = "Género (mujer)",
                      "Calendario_A"                = "Calendario A",
                      "Calendario_B"                = "Calendario B",
                      "cole_privado"                = "Colegio privado",
                      "trabaja_menos_de_10_horas"   = "Trabaja (<10h)",
                      "trabaja__11_a_20_horas"      = "Trabaja (11–20h)",
                      "trabaja__21_a_30_horas"      = "Trabaja (21–30h)",
                      "trabaja_mas_de_30_horas"     = "Trabaja (>30h)",
                      "docenttotal_alumtotal"       = "Tasa docentes/estudiantes",
                      "RISK_VICTIM_2022"            = "Riesgo de victimización",
                      "Homi_x_100k_habitantes"      = "Homicidios (x100k hab.)",
                      "porc_alumn_col_publico"      = "% alumnos colegio público",
                      "terrorismot"                 = "Terrorismo",
                      "hurto"                       = "Hurto",
                      "secuestros"                  = "Secuestros",
                      "discapital"                  = "Distancia a capital",
                      "PIB_percapita_DPTO"          = "PIB per cápita (DPTO)",
                      "proporcio_pob_rural"         = "% población rural",
                      "X..municipios.con.riesgo"    = "% municipios con riesgo",
                      "Homicidios_ponderado_x_100k" = "Homicidios ponderado (x100k)"
)

#A plot in which is possible to visualice posterior means and the crebilibity intervals

ggplot(df, aes(y = reorder(variable, media), x = media, color = Nivel)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lim_inf, xmax = lim_sup), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "Media posterior", y = "Parámetro", color = "Nivel", title ="Ranking de la media posterior de las covariables e intercepto \n" ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right") +
  scale_color_manual(values = c("#D9586C","#3B4CC0", "#788DCE","#BDBDD0"))+
  theme_minimal(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )   # o usa scale_color_manual()


#A more detailed information

bayesian_summary2(beta)
colnames(beta_E)
apply(beta_E,2,bayesian_summary2)
colnames(beta_M)
apply(beta_M,2,bayesian_summary2)
colnames(beta_D)
apply(beta_D,2,bayesian_summary2)

#========= Variance parameters of the covariables----

df <- data.frame(
  variable = c("sigma2_beta","sigma2_E","sigma2_M","sigma2_D"),
  media = c(bayesian_summary(sigma2_beta)[1],bayesian_summary(sigma2_E)[1],bayesian_summary(sigma2_M)[1],bayesian_summary(sigma2_D)[1]),
  li = c(bayesian_summary(sigma2_beta)[2],bayesian_summary(sigma2_E)[2],bayesian_summary(sigma2_M)[2],bayesian_summary(sigma2_D)[2]),
  ls = c(bayesian_summary(sigma2_beta)[3],bayesian_summary(sigma2_E)[3],bayesian_summary(sigma2_M)[3],bayesian_summary(sigma2_D)[3]),
  Nivel=c("Intercepto", "Estudiante","Municipal","Departamental")
)

ggplot(df, aes(y = reorder(variable, media), x = media, color = Nivel)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = li, xmax = ls), height = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "Media posterior", y = "Parámetro", color = "Nivel",title="Medias posteriores de los términos de \n varianza de las covariables \n ") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right") +
  scale_y_discrete(
    labels = c(
      "sigma2_beta" = expression(sigma^2[beta]),
      "sigma2_D"    = expression(sigma^2[D]),
      "sigma2_M"    = expression(sigma^2[M]),
      "sigma2_E"    = expression(sigma^2[E])
    ))+
  scale_color_manual(values = c("#D9586C","#3B4CC0", "#788DCE","#BDBDD0"))+
  theme_minimal(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )    # o usa scale_color_manual()



#=========== Spatial effects----

#Set the work directory to get de clean data
setwd("/Users/macbookpro/Desktop/Tesis/Procesed_data")

#These files correspond to the depurated data which is going to be used
datos <- read.csv('Datos_ProcesadosExamen_Saber_11_2022_2_Entrenamiento.csv')
datos <- datos[!(datos$estu_cod_reside_depto == 99999), ] %>%
  drop_na()

datos <- datos[!(datos$estu_cod_reside_depto == 68 & datos$estu_cod_reside_mcpio == 68264), ]
datos <- datos[!(datos$estu_cod_reside_depto == 88 ), ]
datos <- datos[, !names(datos) %in% c("homicidios", "codprovincia"
                                      , "ano")]%>%
  arrange(estu_cod_reside_depto,estu_cod_reside_mcpio)%>%
  na.omit(datos)


# Here the Municipios whose information is available are extracted 

municipios_data<-datos%>%
  group_by(estu_cod_reside_mcpio)%>%
  summarise(n=n())%>%
  arrange()


#====== Upload Municipios shapefiles

setwd("/Users/macbookpro/Desktop/Tesis/shapes_municipios")
municipios<-st_read("MGN_ADM_MPIO_GRAFICO.shp")

municipios<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

# As its necessary only the information about the Municipios whose information is available
# their shapefiles are selected. 

municipios2<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  inner_join(municipios_data,by=c("mpio_cdpmp"="estu_cod_reside_mcpio"))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

#Obtaining the  phi posterior mean
mean_phi<-apply(phi, 2, mean)
phi_mean<-as.data.frame(cbind(mean_phi,municipios2$mpio_cdpmp))
names(phi_mean)<-c("mean_phi" ,"cod_mpio")

df_mapas<-municipios%>%
  left_join(phi_mean,by=c("mpio_cdpmp"="cod_mpio"))

df_mapas$mean_phi<-ifelse(is.na(df_mapas$mean_phi),0,df_mapas$mean_phi)


#Plotting the map


ggplot(df_mapas) +
  geom_sf(aes(fill = mean_phi), color = "gray70", size = 0.1) +
  scale_fill_gradient(
    low = "#D9586C",   # color para valores bajos
    high = "#3B4CC0",  # color para valores altos
    name = "Media posterior"
  ) +
  labs(
    title = "Medias posteriores del efecto aleatorio \n espacial por municipio \n",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )
#======== Tau_phi----

bayesian_summary2(sqrt(tau_phi))


#======== Municipal variance =======

mean_kappa2_jk<-apply(kappa2_jk, 2, mean)
kappa2_jk_mean<-as.data.frame(cbind(sqrt(mean_kappa2_jk),municipios2$mpio_cdpmp))
names(kappa2_jk_mean)<-c("mean_kappa2_jk" ,"cod_mpio")

df_mapas<-municipios%>%
  left_join(kappa2_jk_mean,by=c("mpio_cdpmp"="cod_mpio"))

df_mapas$mean_kappa2_jk<-ifelse(is.na(df_mapas$mean_kappa2_jk),0,df_mapas$mean_kappa2_jk)

#Plotting the map


ggplot(df_mapas) +
  geom_sf(aes(fill = mean_kappa2_jk), color = "gray70", size = 0.1) +
  scale_fill_gradient(
    low = "#D9586C",   # color para valores bajos
    high = "#3B4CC0",  # color para valores altos
    name = "Media posterior"
  ) +
  labs(
    title = "Medias posteriores de las desviaciónes \n estándar municipales \n",
  ) +
  theme_minimal(base_family = "Helvetica")+
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

#======= Departamental variance ====


setwd("/Users/macbookpro/Desktop/Tesis/shapes_departamentos")

departamentos<-st_read("MGN_DPTO_POLITICO.shp")

departamentos<-departamentos[!(departamentos$DPTO_CCDGO==88),]

mean_kappa2_k<-apply(kappa2_k, 2, mean)
kappa2_k_mean<-as.data.frame(cbind(sqrt(mean_kappa2_k),departamentos$DPTO_CCDGO))
names(kappa2_k_mean)<-c("mean_kappa2_k" ,"cod_depto")

df_mapas<-departamentos%>%
  left_join(kappa2_k_mean,by=c("DPTO_CCDGO"="cod_depto"))

df_mapas$mean_kappa2_k<-ifelse(is.na(df_mapas$mean_kappa2_k),0,df_mapas$mean_kappa2_k)
df_mapas$mean_kappa2_k <- as.numeric(df_mapas$mean_kappa2_k)

ggplot(df_mapas) +
  geom_sf(aes(fill = mean_kappa2_k), color = "gray70", size = 0.1) +
  scale_fill_gradient(
    high = "#D9586C",   # color para valores bajos
    low = "#3B4CC0",    # color para valores altos
    name = "Media posterior"
  ) +
  labs(
    title = "Desviaciones estándar departamentales"
  ) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

#========= Variance parameters =======================

bayesian_summary2(beta_kappa)
bayesian_summary2(alpha_kappa)

