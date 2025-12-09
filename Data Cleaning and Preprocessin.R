library(dplyr)
library(readxl)
library(tidyverse)
library(tibble)
library(stringr)
library(stringi)

#Function to verify if there is any missing
missings<-function(x){
  missings<-sum(is.na(x))
  return(missings)
}

#Setting the work directory
setwd("/Users/macbookpro/Desktop/Tesis/Modelos Bayesianos/Bases")


#Student's Variables----

#Uploading the ICFES training database 

icfes<-read_delim("Examen_Saber_11_20222.txt",delim = ";",escape_double = FALSE, trim_ws = TRUE)

#Applying filters to exclude foreign students and students from San Andrés

icfes2<-icfes %>%
  select(everything()) %>%
  filter(estu_nacionalidad=="COLOMBIA" & estu_pais_reside=="COLOMBIA" & estu_cod_reside_depto != 99999 & cole_cod_depto_ubicacion!= 88)%>%
  mutate(fami_educacionmadre_modif=case_when(fami_educacionmadre =="Educación profesional completa" | fami_educacionmadre == "Postgrado"~ 1,
                                           TRUE ~ 0))%>%
  select(c("estu_consecutivo","cole_mcpio_ubicacion","cole_cod_mcpio_ubicacion","cole_depto_ubicacion","cole_cod_depto_ubicacion",
           "estu_genero","estu_tieneetnia","fami_estratovivienda","fami_numlibros","fami_tieneinternet","fami_tienecomputador","fami_educacionmadre_modif","punt_global"))%>%
  mutate(cole_cod_mcpio_ubicacion=as.numeric(cole_cod_mcpio_ubicacion),cole_cod_depto_ubicacion=as.numeric(cole_cod_depto_ubicacion))%>%
  arrange(cole_cod_depto_ubicacion,cole_cod_mcpio_ubicacion) %>%
  drop_na()

#Correct BELÉN DE BAJIRÁ's code
icfes2$cole_cod_mcpio_ubicacion[icfes2$cole_cod_mcpio_ubicacion == 27086] <- 27086

#Select and recode the student`s variable.

X_ijk<-icfes2%>%
  mutate( #Libros
    fami_numlibros = factor(fami_numlibros,
                            levels = c("11 A 25 LIBROS",
                                       "0 A 10 LIBROS",
                                       "26 A 100 LIBROS",
                                       "MÁS DE 100 LIBROS")
    )
  )%>%
  mutate( #CATEGORIA BASE "0 A 10 LIBROS"
    libros_11_25  = as.integer(fami_numlibros == "11 A 25 LIBROS"),
    libros_26_100 = as.integer(fami_numlibros == "26 A 100 LIBROS"),
    libros_mas100 = as.integer(fami_numlibros == "MÁS DE 100 LIBROS")
  )%>%
  mutate(  #CATEGORIA BASE "Sin Estrato"
    fami_estratovivienda = factor(
      fami_estratovivienda,
      levels = c(
        "Estrato 1",   
        "Estrato 2",
        "Estrato 3",
        "Estrato 4",
        "Estrato 5",
        "Estrato 6",
        "Sin Estrato"
      )
    ),
    estrato_1 = as.integer(fami_estratovivienda == "Estrato 1"),
    estrato_2 = as.integer(fami_estratovivienda == "Estrato 2"),
    estrato_3 = as.integer(fami_estratovivienda == "Estrato 3"),
    estrato_4 = as.integer(fami_estratovivienda == "Estrato 4"),
    estrato_5 = as.integer(fami_estratovivienda == "Estrato 5"),
    estrato_6 = as.integer(fami_estratovivienda == "Estrato 6")
  )%>%
  mutate( #Computador e Internet
    genero = as.integer(estu_genero == "F"),
    computador = as.integer(fami_tienecomputador == "Si"),
    internet   = as.integer(fami_tieneinternet   == "Si"),
    etnia      = as.integer(estu_tieneetnia      == "Si"  )
  )%>%  
  select(c("cole_mcpio_ubicacion","cole_cod_mcpio_ubicacion","cole_depto_ubicacion","cole_cod_depto_ubicacion","punt_global","fami_educacionmadre_modif","genero","computador",
           "internet","etnia",paste0(rep("estrato_",times=6),c(1:6)),,paste0(rep("libros_",times=3),c("11_25","26_100","mas100"))))


# Municipal Variables----
##Cede Educación

cede_edu<-read_excel("PANEL_DE_EDUCACION(2022).xlsx")

cede_edu1<-cede_edu%>%
  select(c("codmpio","ano","docen_total", "alumn_total"))%>%
  filter(ano==2021)%>%
  mutate(docen_estu=docen_total/alumn_total)%>%
  arrange(codmpio)%>%
  drop_na()

cede_edu1[cede_edu1$codmpio == 27086,1] <- 27493


##Cede Panel General

cede_general<-read_excel("PANEL_CARACTERISTICAS_GENERALES(2022).xlsx")

cede_general1<-cede_general%>%
  select(c("coddepto","codmpio","ano","pobl_tot","pobl_rur"))%>%
  filter(ano==2021)%>%
  mutate(porcen_rural=pobl_rur/pobl_tot*100)%>%
  arrange(coddepto,codmpio)%>%
  drop_na()


cede_general1$codmpio[cede_general1$codmpio==27086]<- 27493


##Cede Panel Violencia

#Data about homicides is taken with respect to the 2019 due to the missing data from 2021 (2021:297Na 2019:2Na)

cede_violencia<-read_excel("PANEL_CONFLICTO_Y_VIOLENCIA(2022).xlsx")

cede_violencia1<-cede_violencia%>%
  select(ano,codmpio,homicidios)%>%
  filter(ano==2019)%>%
  arrange(codmpio)%>%
  drop_na()


cede_violencia1$codmpio[cede_violencia1$codmpio==27086]<- 27493


#DIVIPOLA: the dataset contains all the Colombia's Municipios and Departamentos

codigos<-read_excel("DIVIPOLA_Municipios.xlsx",range = "A11:D1133")
names(codigos)<-c("cod_dep","departamento","cod_mun","municipio")

#Normalizing the text

codigos<-within(codigos,{
  departamento<-str_squish(departamento)
  departamento<-stri_trans_general(departamento, "Latin-ASCII")
  municipio<-str_squish(municipio)
  municipio<-stri_trans_general(municipio, "Latin-ASCII")
  })%>%
  mutate(cod_dep=as.numeric(cod_dep),cod_mun=as.numeric(cod_mun))

##LAFT 
laft<-read_excel("LAFT.xlsx",sheet = 4)
laft<-laft%>%
      select(c("CODIGO","RISK_VICTIM_2022"))

#Normalizing the text
laft<-within(laft,{ 
  CODIGO<-str_squish(CODIGO)
  CODIGO<-stri_trans_general(CODIGO, "Latin-ASCII")
  Dep<-str_sub(CODIGO,1,str_locate(CODIGO, "-")[,1]-1)
Mun<-str_sub(CODIGO,str_locate(CODIGO, "-")[,1]+1,length(CODIGO))
}
)

laft1<-laft%>%
  left_join(codigos,by=c("Dep"="departamento","Mun"="municipio"))%>%
  mutate(cod_mun=as.numeric(cod_mun))%>%
  drop_na()%>%
  as.data.frame()

#The dataframe below consist of all the Munipical level variables

Z_jk<-cede_violencia1%>%
  left_join(cede_general1,by=c("codmpio"))%>%
  left_join(cede_edu1,by=c("codmpio"))%>%
  left_join(laft1,by=c("codmpio"="cod_mun"))%>%
  mutate(tasa_homicidios=100000*homicidios/pobl_tot)%>%
  arrange(coddepto,codmpio)%>%
  select(c("codmpio","docen_estu","tasa_homicidios","RISK_VICTIM_2022"))%>%
  na.omit()

#Departamental variables----

#Treatment of PIB database
pib<-read_excel("DANE - PIB.xlsx",sheet = 4,range = "A9:U43",col_names = T,col_types = c("numeric","text",rep("numeric",times=19)))%>%
  select("Código Departamento (DIVIPOLA)" ,"DEPARTAMENTOS","2022")%>%
  drop_na()

names(pib)<-c("cod_dep","DEPARTAMENTOS","PIB")

pib<-pib%>%
  filter(cod_dep!=88)%>%
  arrange(cod_dep)

#General CEDE batabase grouped by Departamento
cede_general_dep<-cede_general1%>%
  group_by(coddepto)%>%
  filter(coddepto!=88)%>%
  summarise(tot_pob=sum(pobl_tot),tot_pob_rural=sum(pobl_rur))%>%
  mutate(porcen_rural=100*tot_pob_rural/tot_pob)

MOE<-read_excel("MOE.xlsx")%>%
  drop_na()

MOE<-within(MOE,{
  Depto<-str_squish(Depto)
})

#Violence CEDE batabase grouped by Departamento

cede_violencia_dep<-cede_violencia1%>%
  left_join(codigos,by=c("codmpio"="cod_mun"))%>%
  group_by(cod_dep)%>%
  summarise(sum_hom=sum(homicidios),n_k=n())%>%
  mutate(prom_hom=sum_hom/n_k)

#The matriz containing de Departamento's Information
 
W_k<-pib %>%
  left_join(MOE,by=c("DEPARTAMENTOS"="Depto"))%>%
  left_join(cede_general_dep,by=c("cod_dep"="coddepto"))%>%
  left_join(cede_violencia_dep,by=c("cod_dep"))%>%
  mutate(pib_pc=PIB/1000000)%>%
  select(cod_dep,DEPARTAMENTOS,porcen_rural,pib_pc,"% municipios con riesgo","prom_hom")

names(W_k)<-c("cod_dep","DEPARTAMENTOS","porcen_rural","pib_pc","Porcentaje municipios con riesgo","prom_hom")

W_k$`Porcentaje municipios con riesgo`<-as.numeric(W_k$`Porcentaje municipios con riesgo`)

#Save just the students whose muncipal information is available

Full_dataset_clean<-X_ijk%>%
  inner_join(Z_jk,by=c("cole_cod_mcpio_ubicacion"="codmpio"))%>%
  inner_join(W_k,by=c("cole_cod_depto_ubicacion"="cod_dep"))

#Training and testing division----

set.seed(1278)

id<-sample(c(1:nrow(Full_dataset_clean)),nrow(Full_dataset_clean)*0.5)

Training_dataset_clean<-Full_dataset_clean[id,]%>%
  arrange(cole_cod_depto_ubicacion,cole_cod_mcpio_ubicacion)%>%
  select(-DEPARTAMENTOS)

Test_dataset_clean<-Full_dataset_clean[-id,]%>%
  arrange(cole_cod_depto_ubicacion,cole_cod_mcpio_ubicacion)%>%
  select(-DEPARTAMENTOS)

# Number of students from the j municipio and k departamento for train and test dataset:

# Note: It's mandatory to create n_jk_test and n_jk_train because the number of students from j municipio
# in the trainning dataset is different from the number of students from j municipio in the test dataset.

n_jk_train<-Training_dataset_clean%>%
  select("cole_mcpio_ubicacion","cole_cod_mcpio_ubicacion","cole_depto_ubicacion","cole_cod_depto_ubicacion")%>%
  group_by(cole_cod_mcpio_ubicacion,cole_cod_depto_ubicacion)%>%
  summarise(n_jk=n())%>%
  arrange(cole_cod_depto_ubicacion,cole_cod_mcpio_ubicacion)

n_jk_test<-Test_dataset_clean%>%
  select("cole_mcpio_ubicacion","cole_cod_mcpio_ubicacion","cole_depto_ubicacion","cole_cod_depto_ubicacion")%>%
  group_by(cole_cod_mcpio_ubicacion,cole_cod_depto_ubicacion)%>%
  summarise(n_jk=n())%>%
  arrange(cole_cod_depto_ubicacion,cole_cod_mcpio_ubicacion)

# Number of students from the  k departamento for train and test dataset:

# Note: It's mandatory to create n_jk_test and n_jk_train because the number of students from j municipio
# in the trainning dataset is different from the number of students from j municipio in the test dataset.

n_k_train<-Training_dataset_clean%>%
  select("cole_mcpio_ubicacion","cole_cod_mcpio_ubicacion","cole_depto_ubicacion","cole_cod_depto_ubicacion")%>%
  group_by(cole_cod_depto_ubicacion)%>%
  summarise(n_k=n())%>%
  arrange(cole_cod_depto_ubicacion)

n_k_test<-Test_dataset_clean%>%
  select("cole_mcpio_ubicacion","cole_cod_mcpio_ubicacion","cole_depto_ubicacion","cole_cod_depto_ubicacion")%>%
  group_by(cole_cod_depto_ubicacion)%>%
  summarise(n_k=n())%>%
  arrange(cole_cod_depto_ubicacion)

#Number of municipios inside every departament k:
# The number of municipios inside every departament is constant, and is not afected by the dataset division 
m_k<-n_jk_train%>%
  group_by(cole_cod_depto_ubicacion)%>%
  summarise(n_k=n())


X_ijk_train<-Training_dataset_clean[,c(1:19)]

Z_jk<-Training_dataset_clean[,c(seq(1,4),seq(20,22))]%>%
  select("cole_cod_mcpio_ubicacion","docen_estu","tasa_homicidios","RISK_VICTIM_2022")%>%
  group_by(cole_cod_mcpio_ubicacion,docen_estu,tasa_homicidios,RISK_VICTIM_2022)%>%
  summarise(n_jk=n())%>%
  select(-n_jk)


X_ijk_test<-Test_dataset_clean[,c(1:19)]

Z_jk<-Test_dataset_clean[,c(seq(1,4),seq(20,22))]%>%
  select("cole_cod_mcpio_ubicacion","docen_estu","tasa_homicidios","RISK_VICTIM_2022")%>%
  group_by(cole_cod_mcpio_ubicacion,docen_estu,tasa_homicidios,RISK_VICTIM_2022)%>%
  summarise(n_jk=n())%>%
  select(-n_jk)



#Saving the files in RData format---

Train_dataset_clean<-list("Training_dataset_clean"=Training_dataset_clean,"X_ijk"=X_ijk_train,"Municipal_data"=Z_jk,"Departamental_data"=W_k,"n_jk"=n_jk_train,"n_k"=n_k_train,"m_k"=m_k)
Test_dataset_clean<-list("Test_dataset_clean"=Test_dataset_clean,"X_ijk"=X_ijk_test,"Municipal_data"=Z_jk,"Departamental_data"=W_k,"n_jk"=n_jk_test,"n_k"=n_k_train,"m_k"=m_k)


save(Train_dataset_clean,file="/Users/macbookpro/Desktop/Tesis/Procesed_Data/Training_dataset.RData")
save(Test_dataset_clean,file="/Users/macbookpro/Desktop/Tesis/Procesed_Data/Test_dataset.RData")



