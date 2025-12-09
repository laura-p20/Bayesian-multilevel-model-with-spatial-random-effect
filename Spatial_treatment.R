#==============================================================================
#             Municipal  Neighbors 
#------------------------------------------------------------------------------
#  The adjacency information within Municipios, based on the neighborhood relationship, 
#  is built in these script
#==============================================================================


#========== Libraries

library(sf)
library(spdep)
library(dplyr)
library(stringi)
library(stringr)

#====== Upload the training data set 

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

municpios_data<-municipios_data$estu_cod_reside_mcpio


#====== Upload Municipios shapefiles

setwd("/Users/macbookpro/Desktop/Tesis/shapes_municipios")
municipios<-st_read("MGN_ADM_MPIO_GRAFICO.shp")

# As its necessary only the information about the Municipios whose information is available
# their shapefiles are selected. 

municipios<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  inner_join(municipios_data,by=c("mpio_cdpmp"="estu_cod_reside_mcpio"))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

#===== Built the neighborhood relationship


#Include to valid de shapefile
municipios <- st_make_valid(municipios)

#A list of each municipio containing their neighbors index ()
neighbors_municipios<- st_touches(municipios, sparse = TRUE)


#A vector cotaining each Munipio's number of neighbors 
d_jk<-rep(NA,times=nrow(municipios))

for (i in c(1:nrow(municipios))) {
  d_jk[i]<-length(neighbors_municipios[[i]])
}
#There is only one municipio which has no neighbors, but that's just really because there is
#no information about its neighbors, then in order to do not indeterminate future calculus, the
#count of the number of neighbors for it is set to 1. 

d_jk<-ifelse(d_jk==0,1,d_jk)

#To save the information
neighbors_municipios<-list(neighbors=neighbors_municipios,"d_jk"=d_jk)
save(neighbors_municipios,file="/Users/macbookpro/Desktop/Tesis/Procesed_Data/Spatial_Component.RData")



#============= Analysis of neighbors out of the departament

neighbors<-data.frame(
  Municipio=municipios$mpio_cnmbr,
  Neighbors=d_jk,
  Non_Departamental_neighbors=rep(NA,nrow(municipios))
)

for (i in c(1:nrow(municipios))) {
  dep<-municipios$dpto_ccdgo[i]
  a<-rep(NA,times=neighbors$Neighbors[i])
  a<-ifelse(str_sub(municipios$mpio_cdpmp[neighbors_municipios[[i]]],1,2)==dep,0,1)
  neighbors$Non_Departamental_neighbors[i]<-sum(a)
}

neighbors$Non_Departamental_neighbors_rate<-neighbors$Non_Departamental_neighbors/neighbors$Neighbors

sum(neighbors$Non_Departamental_neighbors_rate>0.4)/nrow(municipios)

