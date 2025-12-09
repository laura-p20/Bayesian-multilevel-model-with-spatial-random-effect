#=======================================================================================
#. Description:  
#. The simulation data set is built in this script 
#=======================================================================================


#===============Libraries and preparation----

library(dplyr)
library(MASS)

#Set the work directory to get a previus result of the MCMC using the real data

setwd("/Users/macbookpro/Desktop/Tesis/Results")
load("MCMC_model2_real_data_standarized_.RData")
setwd("/Users/macbookpro/Desktop/Tesis/Procesed_Data")
load("data_ready.RData")

X_ijk<-data_ready$X_ijk
Z_jk<-scale(data_ready$Z_jk)
W_k<-scale(data_ready$W_k)
n<-data_ready$n
n_jk<-data_ready$n_jk
m<-data_ready$m


#Upload the spatial effect
phi<-THETA[[9]]
phi_post<-apply(phi, 2, mean)

#Upload the municipal variance 
kappa2_jk<-THETA[[11]]
kappa2_jk_post<-apply(kappa2_jk,2,mean)

# Simulation proposals ----
# Crear una lista para almacenar los data.frames de parámetros por nivel
scenarios_parameter <- list()

# --- 1. Parámetros de Intercepto y Nivel Estudiante ---
scenarios_parameter$student <- data.frame(
  Parameter = c("Intercepto_beta", "Educacion_madre", "Computador", "Internet", "Etnia", 
                 "libros_11_25", "libros_26_100", "libros_mas100", 
                "estrato_1", "estrato_2", "estrato_3", "estrato_4", "estrato_5", 
                "estrato_6","Genero", "Calendario_A", "Calendario_B", "cole_privado", 
                "trabaja_menos_de_10_horas", "trabaja_11_a_20_horas", 
                "trabaja_21_a_30_horas", "trabaja_mas_de_30_horas"),
  
  scenario_1 = c(196.16, 18.13, 10.36, 8.93, -14.34, -7.58, 8.48, 20.06, 20.27, 
                            20.36, 18.75, 17.00, 14.45, 11.37, 2.53, 23.81, 19.07, 12.82, 
                            -13.42, -12.85, -16.16, -26.14),
  
  scenario_2 = c(250, 0.50, 1.00, 1.00, -0.5, -0.250, 1.00, 2.00, 2.00, 
                           45, 45, 45, 45, 45, 45, 2.00, 1.00, 1.00, 
                           -0.4250, -0.4250, -0.4250, -0.4250),
  
  scenario_3 = c(353, 12.13, 12.00, 12.00, 0, 0, -9.48, 10.06, 10.27, 
                               0, 0, 0, 0, 0, 0, 10.8, 19.07, 13.82, 
                               -13.42, -13.85, -16.16, -26.14)
)


# --- 2. Parámetros de Nivel Municipal ---
scenarios_parameter$Municipal <- data.frame(
  Parameter = c("Tasa_docentes_por_estudiante", "Riesgo_de_victimizacion", 
                "Homicidios_100mil", "Pct_estudiantes_colegio_publico", 
                "Terrorismo", "hurto", "Secuestros", "Distancia_a_la_capital"),
  
  scenario_1 = c(2.90, 8.7, 0.06, 0.04, 0.34, 0.00, 1.12, 0.05),
  
  scenario_2= c(0, 0.001, 0.00, 0.00, 0.00, 0.00, 0.10, 0.001),
  
  scenario_3 = c(-5.00, -9.54,-0.250, -0.8, -0.034, -0.00025, 0.8, 0.005)

  # Note: Se usan valores base ya que los CV% son bajos, indicando alta precisión.
)


# --- 3. Parámetros de Nivel Departamental ---
scenarios_parameter$Departamental <- data.frame(
  Parameter = c("PIB_per_capita_departamental", "Proporcion_de_poblacion_rural", 
                "Municipios_con_riesgo", "Homicidios_ponderados"),
  
  scenario_1 = c(0.23, 0.06, 0.004, 0.012),
  
  scenario_2= c(0.005, 0.00, 0.00, 0.00),
  
  scenario_3 = c(0.23, -0.06, -0.04, -0.12)
)


# --- Impresión del resultado (Opcional, para verificar) ---
scenarios_parameter$phi_post<-phi_post
scenarios_parameter$kappa2_jk_post<-kappa2_jk_post
print(scenarios_parameter)
setwd("/Users/macbookpro/Desktop/Tesis/Procesed_Data/")
save(scenarios_parameter,file="simulation_values.RData")



#---- Dataa simulation of response variable

responses<-matrix(NA,ncol = 5,nrow = n)
kappa2_jk_post<-rep(sqrt(kappa2_jk_post),n_jk)
set.seed(2301)

for(i in c(2:4)){
  
  responses[,i-1]<-scenarios_parameter$student[[i]][1]+ X_ijk%*%scenarios_parameter$student[[i]][-1]+Z_jk%*%scenarios_parameter$Municipal[[i]]+W_k%*%scenarios_parameter$Departamental[[i]]+rep(phi_post,n_jk)+rnorm(n,rep(0,n),kappa2_jk_post)
  a<-responses[,i-1]
  b<-summary(a)
  cat("Escenario",i-1,"\n","Resúmen",b,"\n","\n")
}

y_ijk_3<-responses[,3]

save(y_ijk_3,file="simulation_response_scene_3.RData")
load("~/Desktop/Tesis/Procesed_Data/simulation_response_scene_1.RData")
load("~/Desktop/Tesis/Procesed_Data/simulation_response_scene_2.RData")
load("~/Desktop/Tesis/Procesed_Data/simulation_response_scene_3.RData")


