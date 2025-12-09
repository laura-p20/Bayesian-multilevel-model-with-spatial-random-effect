#==================================================#
#. Descriptive maps
#==================================================#


#.  Libraries and settings=========
library(dplyr)
library(tidyverse)
library(cluster)
library(ggplot2)
library(sf)
library(spdep)
library(dplyr)
library(stringi)
library(stringr)
suppressMessages(suppressWarnings(library(coda)))
library(NbClust)
library(mclust)
library(lattice)
library(reshape2)
library(RColorBrewer)
library(colorspace)
library(ggspatial)
library(prettymapr)
library(terra)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)

#Upload the information

load("/Users/macbookpro/Desktop/Tesis/Procesed_Data/data_ready.RData")

# Municipal shapes
municipios<-st_read("/Users/macbookpro/Desktop/Tesis/shapes_municipios/MGN_ADM_MPIO_GRAFICO.shp")
municipios<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

# Departamental shapes
departamentos<-st_read("~/Desktop/Tesis/shapes_departamentos/MGN_DPTO_POLITICO.shp")
departamentos<-departamentos[!(departamentos$DPTO_CCDGO==88),]%>%
  mutate(DPTO_CCDGO=as.numeric(DPTO_CCDGO))



X_ijk<-cbind(data_ready$X_ijk,data_ready$id_dep_mun)
Z_jk<-cbind(data_ready$Z_jk,data_ready$id_dep_mun)
W_k<-cbind(data_ready$W_k,data_ready$id_dep_mun)
y_ijk<-data.frame(punt_global=data_ready$punt_global,data_ready$id_dep_mun)
n<-data_ready$n
n_jk<-data_ready$n_jk
m<-data_ready$m
d<-data_ready$d
index<-data_ready$id_dep_mun

#======= Global score    ========
##========  Municipal Sample global score========

y_ijk_mun<-y_ijk%>%
  group_by(estu_cod_reside_mcpio)%>%
  summarise(mean_mun=mean(punt_global))%>%
  arrange(desc(mean_mun))

global_mun<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  left_join(y_ijk_mun,by=c("mpio_cdpmp"="estu_cod_reside_mcpio"))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

global_mun$mean_mun<-ifelse(is.na(global_mun$mean_mun),mean(global_mun$mean_mun,na.rm=TRUE),global_mun$mean_mun)


# Dowload the other countries silhoutes
mundo <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = global_mun, 
          aes(fill = mean_mun), 
          color = "#BEBEBE", 
          lwd = 0.015) +
  
  # Color
  scale_fill_continuous_sequential(palette="BurgYl") +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "",subtitle = "",fill="Mean") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"  )

##======== Departamental Sample global score========


y_ijk_dep<-y_ijk%>%
  group_by(estu_cod_reside_depto)%>%
  summarise(mean_dep=mean(punt_global))%>%
  arrange(desc(mean_dep))

global_dep<-departamentos%>%
  left_join(y_ijk_dep,by=c("DPTO_CCDGO"="estu_cod_reside_depto"))



ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = global_dep, 
          aes(fill = mean_dep), 
          color = "#BEBEBE", 
          lwd = 0.025) +
  
  # Color
  scale_fill_continuous_sequential(palette="BurgYl") +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "",subtitle = "",fill="Mean") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"  )


#======== Covariates =========
## Departamental level======

W_k_agg<-W_k%>%
  group_by(estu_cod_reside_depto)%>%
  summarise(PIB=mean(PIB_percapita_DPTO),pro_pob_rural=mean(proporcio_pob_rural),
            perce_mun_risk=mean(X..municipios.con.riesgo),homicides=mean(Homicidios_ponderado_x_100k))

covariate_dep<-departamentos%>%
  left_join(W_k_agg,by=c("DPTO_CCDGO"="estu_cod_reside_depto"))



ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = covariate_dep, 
          aes(fill = PIB), 
          color = "#BEBEBE", 
          lwd = 0.025) +
  
  # Color
  scale_fill_continuous_sequential(palette="BurgYl") +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "",subtitle = "",fill="Million \n (COP)") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"  )


ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = covariate_dep, 
          aes(fill = pro_pob_rural), 
          color = "#BEBEBE", 
          lwd = 0.025) +
  
  # Color
  scale_fill_continuous_sequential(palette="BurgYl") +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "",subtitle = "",fill="%") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"  )


ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = covariate_dep, 
          aes(fill = perce_mun_risk), 
          color = "#BEBEBE", 
          lwd = 0.025) +
  
  # Color
  scale_fill_continuous_sequential(palette="BurgYl") +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "",subtitle = "",fill="%") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"  )



ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = covariate_dep, 
          aes(fill = homicides), 
          color = "#BEBEBE", 
          lwd = 0.025) +
  
  # Color
  scale_fill_continuous_sequential(palette="BurgYl") +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "",subtitle = "",fill="%") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"  )


ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = covariate_dep, 
          aes(fill = perce_mun_risk), 
          color = "#BEBEBE", 
          lwd = 0.025) +
  
  # Color
  scale_fill_continuous_sequential(palette="BurgYl") +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "",subtitle = "",fill="Homicides per\n100k inhabitants") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"  )



#======== Covariates =========
## Municipal level======


Z_jk_agg<-Z_jk%>%
  group_by(estu_cod_reside_mcpio)%>%
  summarise(docen=mean(docenttotal_alumtotal),victim=mean(RISK_VICTIM_2022),homi=mean(Homi_x_100k_habitantes),
            porcen_alu=mean(porc_alumn_col_publico),terror=mean(terrorismot),hurto=mean(hurto),
            secuestro=mean(secuestros),distan=mean(discapital))

covariate_mun<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  left_join(Z_jk_agg,by=c("mpio_cdpmp"="estu_cod_reside_mcpio"))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

covariate_mun$victim<-ifelse(is.na(covariate_mun$victim),mean(covariate_mun$victim,na.rm=TRUE),covariate_mun$victim)
covariate_mun$porc_col<-ifelse(is.na(covariate_mun$porc_col),mean(covariate_mun$porc_col,na.rm=TRUE),covariate_mun$porc_col)
covariate_mun$terrorismo<-ifelse(is.na(covariate_mun$terrorismo),mean(covariate_mun$terrorismo,na.rm=TRUE),covariate_mun$terrorismo)
covariate_mun$distan<-ifelse(is.na(covariate_mun$distan),mean(covariate_mun$distan,na.rm=TRUE),covariate_mun$distan)


ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = covariate_mun, 
          aes(fill = hurto), 
          color = "#BEBEBE", 
          lwd = 0.015) +
  
  # Color
  scale_fill_continuous_sequential(palette="BurgYl") +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "",subtitle = "",fill="Victimization \nRisk index") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"  )




