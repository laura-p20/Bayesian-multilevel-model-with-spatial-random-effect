###############################################################################
#. This file contains the posterior inference and results build up by the 
#. MCMC samples of the Ridge models.
#  Here is algo the maps and the posterior departamental and municipal Rankings
################################################################################

#===========Libraries and settings====

setwd("/Users/macbookpro/Desktop/Tesis/Results")

#===============Libraries and preparation----

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
  
#Upload the covariables matrices and some numbers

load("/Users/macbookpro/Desktop/Tesis/Procesed_Data/data_ready.RData")

X_ijk<-data_ready$X_ijk
Z_jk<-data_ready$Z_jk
Z_jk<-scale(Z_jk)
W_k<-data_ready$W_k
W_k<-scale(W_k)
n<-data_ready$n
n_jk<-data_ready$n_jk
m<-data_ready$m
d<-data_ready$d
index<-data_ready$id_dep_mun

y_ijk<-data_ready$punt_global


#Municipios whose info is available
municipios_used<-data_ready$data_complete
municipios_used<-municipios_used%>%
  group_by(estu_cod_reside_mcpio,estu_mcpio_reside)%>%
  summarise(n=n())%>%
  arrange()

municipios_used<-municipios_used[,-3]

departamentos_used<-data_ready$data_complete%>%
  group_by(estu_cod_reside_depto,estu_depto_reside)%>%
  summarise(n=n())%>%
  arrange()
departamentos_used<-departamentos_used[,-3]


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

IC_0.025<-function(x){
  return(quantile(x,c(0.025)))
}
 IC_0.975<-function(x){
   return(quantile(x,c(0.975)))
 }

#=========== Extracting the MCMC samples of the Ridge Model ----

load("~/Desktop/Tesis/Results/MCMC_ridge_real_data.RData")
#load("~/Desktop/Tesis/Results/MCMC_model2_real_data_standarized.RData")
beta<-THETA[[1]]
beta_E<-THETA[[2]]
beta_M<-THETA[[3]]
beta_D<-THETA[[4]]

sigma2_beta<-THETA[[5]]
lambda2_E<-THETA[[6]]
lambda2_M<-THETA[[7]]
lambda2_D<-THETA[[8]]

phi<-THETA[[9]]
tau_phi<-THETA[[10]]

kappa2_jk<-THETA[[11]]
kappa2_k<-THETA[[12]]

alpha_kappa<-THETA[[13]]
beta_kappa<-THETA[[14]]

acr<-THETA[[15]]

nsams<-length(beta)
#================== Global score posterior mean =======

# Bayesian stimation of the parameters

beta_post<-mean(beta)
beta_E_post<-apply(beta_E, 2, stat_mode)
beta_M_post<-apply(beta_M, 2, stat_mode)
beta_D_post<-apply(beta_D, 2, stat_mode)

phi_post<-apply(phi,2,mean)

#The Global score posterior mean disgregated

mun_post_mean<-matrix(NA,nrow = nsams,ncol=m)
dep_post_mean<-matrix(NA,nrow = nsams,ncol=d)


for (i in c(1:nsams)) {
  global_mean<-rep(beta[i],n)+X_ijk%*%beta_E[i,]+Z_jk%*%beta_M[i,]+W_k%*%beta_D[i,]+rep(phi[i,],n_jk)
  mun_post_mean[i,]<-tapply(global_mean,index$estu_cod_reside_mcpio,mean)
  dep_post_mean[i,]<-tapply(global_mean,index$estu_cod_reside_depto,mean)
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
}



df_mun<-data.frame(
  media_post=apply(mun_post_mean, 2,mean),
  lim_inf=apply(mun_post_mean, 2,IC_0.025),
  lim_sup=apply(mun_post_mean, 2,IC_0.975)
)

df_mun<-cbind(municipios_used,df_mun)
df_mun$signif_dif<-ifelse(df_mun$lim_inf>250,"1",ifelse(df_mun$lim_sup<250,"2","3"))

df_mun<-df_mun%>%
  arrange(desc(media_post))
  


df_dep<-data.frame(
  media_post=apply(dep_post_mean, 2,mean),
  lim_inf=apply(dep_post_mean, 2,IC_0.025),
  lim_sup=apply(dep_post_mean, 2,IC_0.975)
)

df_dep<-cbind(departamentos_used,df_dep)
df_dep$signif_dif<-ifelse(df_dep$lim_inf>250,"1",ifelse(df_dep$lim_sup<250,"2","3"))

#A plot in which is possible to visualice posterior means and the crebilibity intervals
# of the municipal posterior mean
ggplot(df_mun[c(1:20,1091:1111),], aes(y = reorder(estu_mcpio_reside, media_post), 
                   x = media_post, 
                   color = signif_dif)) +
  geom_point(size = 1.5) +
  geom_errorbarh(aes(xmin = lim_inf, xmax = lim_sup), height = 0.3) +
  geom_vline(xintercept = 250, linetype = "dashed", color = "gray50") +
  labs(x = "Posterior Mean", y = "Departament", color = "", title = "Global Score Posterior Mean",subtitle = "Top & Bottom 20") +
  scale_color_manual(values = c( "#70284a", "#dc7176","#f2a28a")) +
  theme_minimal()+
  theme(
    plot.title = element_text(color = "#7C3B5E",face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40",size = 12, hjust = 0.5),
    axis.title = element_text(color = "grey40",face = "plain"),
    legend.position = "none"
  )



#A plot in which is possible to visualice posterior means and the crebilibity intervals

ggplot(df_mun, aes(y = reorder(estu_depto_reside, media_post), 
                   x = media_post, 
                   color = signif_dif)) +
  geom_point(size = 1.5) +
  geom_errorbarh(aes(xmin = lim_inf, xmax = lim_sup), height = 0.3) +
  geom_vline(xintercept = 250, linetype = "dashed", color = "gray50") +
  labs(x = "Posterior Mean", y = "Departament", color = "", title = "Global Score Posterior Mean",subtitle = "Departamental Ranking") +
  scale_color_manual(values = c( "#70284a", "#dc7176","#f2a28a")) +
  theme_minimal(base_size = 13, base_family = "Helvetica") +
  theme(
    plot.title = element_text(color = "#7C3B5E",face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40",size = 12, hjust = 0.5),
    axis.title = element_text(color = "grey40",face = "plain"),
    legend.position = "none"
  )



#================ Clustering according to the posterior mean ========

##===== Departamental Cluster====
n_clusters<-rep(NA,nsams)
adjacency_matrix_dep<-matrix(0,nrow = d,ncol = d)

for (i in c(1:nsams)) {
  data<-cbind(dep_post_mean[i,],sqrt(kappa2_k[i,]))
  a<-NbClust(data,method = "kmeans",index="silhouette")
  n_clusters[i]<- a$Best.nc[1]
  k_means<-kmeans(data,centers = n_clusters[[i]])
  
  adjacency_matrix_iter<-matrix(NA,nrow = d,ncol = d)
  
  for (j in c(1:d)) {
    adjacency_matrix_iter[j,]<-c(rep(0,j),as.numeric(k_means$cluster[j]==k_means$cluster[-c(1:j)])*1/nsams)
  }
  
  adjacency_matrix_dep<-adjacency_matrix_dep+adjacency_matrix_iter
  
  
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }

}

adjacency_matrix_dep<-adjacency_matrix_dep+diag(1,nrow = d)+t(adjacency_matrix_dep)

colnames(adjacency_matrix_dep)<-departamentos$estu_depto_reside
rownames(adjacency_matrix_dep)<-departamentos$estu_depto_reside


idx_order<-data.frame(
  idx=c(1:d),
  post_mean=colMeans(dep_post_mean)
)%>%
arrange(desc(post_mean))

adjacency_matrix_dep<-adjacency_matrix_dep[idx_order$idx,idx_order$idx]


adjacency_dep_df <- melt(adjacency_matrix_dep)
names(adjacency_dep_df ) <- c("x", "y", "p")

adjacency_dep_df$y <- factor(adjacency_dep_df$y, 
                             levels = rev(unique(adjacency_dep_df$y)))

base_plot <- ggplot(adjacency_dep_df, aes(x = x, y = y, fill = p)) +
  geom_tile() + 
  coord_fixed() + 
  theme_minimal() +
  labs(title = "Heat Map- Adjacency Matrix",subtitle = "Departamental relationship", x = NULL, y = NULL,fill="Probability")+theme(
    plot.title = element_text(
      color = "#7C3B5E",   
      size  = 18,          
      face  = "bold",     
      hjust = 0.5          
    ),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
  )
base_plot<-base_plot +
scale_fill_continuous_sequential(palette = "BurgYl")+
scale_x_discrete(labels = NULL) + 
  theme(axis.text.x = element_blank())



#Cluster estimation based on adjacency matrix 

clust_label_dep<-Mclust(adjacency_matrix_dep,verbose = FALSE)
clusters_df_dep<-data.frame(
  cluster=clust_label_dep$classification,
  post_mean=df_dep$media_post[idx_order$idx],
  dep_cod=df_dep$estu_cod_reside_depto[idx_order$idx]
    )

departamentos<-st_read("~/Desktop/Tesis/shapes_departamentos/MGN_DPTO_POLITICO.shp")
departamentos<-departamentos[!(departamentos$DPTO_CCDGO==88),]%>%
  mutate(DPTO_CCDGO=as.numeric(DPTO_CCDGO))

df_mapas<-departamentos%>%
  left_join(clusters_df_dep,by=c("DPTO_CCDGO"="dep_cod"))



# Dowload the other countries silhoutes
mundo <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = df_mapas, 
          aes(fill = as.factor(cluster)), 
          color = "#BEBEBE", 
          lwd = 0.025) +
  
  # Color
  scale_fill_discrete_sequential(
    palette = "BurgYl", 
    rev = FALSE,
    name = "Cluster"
  ) +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "Departamental Segmentation",subtitle = "Clustering based on posterior mean estimates") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# Cluster analysis

ggplot(clusters_df_dep, aes(x = as.factor(cluster), y = post_mean, fill = as.factor(cluster), colour = as.factor(cluster))) +
  geom_boxplot(
    alpha = 0.7,          # Transparencia media
    outlier.shape = NA   # Ocultamos outliers  # Borde gris oscuro uniforme
  ) +
  geom_jitter(
    aes(color = as.factor(cluster)), 
    width = 0.1, 
    alpha = 0.25
  ) +
  scale_color_discrete_sequential(
    palette = "BurgYl",
    rev = FALSE,
    name = "Cluster"
  )+
  scale_fill_discrete_sequential(
    palette="BurgYl",
    rev=FALSE,
    name="Cluster"
  )+
  theme_minimal()+
  labs(title = "Departamental posterior mean",subtitle = "Summary across clusters", x = NULL, y = NULL,fill="Probability")+theme(
    plot.title = element_text(
      color = "#7C3B5E",   
      size  = 18,          
      face  = "bold",     
      hjust = 0.5          
    ),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
  )



tapply(clusters_df_dep$post_mean, clusters_df_dep$cluster, summary)



##===== Municipal Cluster====


adjacency_matrix_mun<-matrix(0,nrow = m,ncol = m)

for (i in c(1:nsams)) {
  
  #Data prep for each iter
  data<-cbind(mun_post_mean[i,],sqrt(kappa2_jk[i,]))
  
  #Real K-means
  k_means<-kmeans(data,centers = 3)
  
  adjacency_matrix_iter<-matrix(NA,nrow = m,ncol =m )

  for (j in c(1:m)) {
    adjacency_matrix_iter[j,]<-c(rep(0,j),as.numeric(k_means$cluster[j]==k_means$cluster[-c(1:j)])*1/nsams)
  }

  adjacency_matrix_mun<-adjacency_matrix_mun+adjacency_matrix_iter
  
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
  
}

adjacency_matrix_mun<-adjacency_matrix_mun+diag(1,nrow = m)+t(adjacency_matrix_mun)


colnames(adjacency_matrix_mun)<-municipios_used$estu_mcpio_reside
rownames(adjacency_matrix_mun)<-municipios_used$estu_mcpio_reside


idx_order_mun<-data.frame(
  idx=c(1:m),
  post_mean=colMeans(mun_post_mean)
)%>%
  arrange(desc(post_mean))


adjacency_matrix_mun3<-adjacency_matrix_mun[idx_order_mun$idx,idx_order_mun$idx]


adjacency_mun_df <- melt(adjacency_matrix_mun3)
names(adjacency_mun_df) <- c("x", "y", "p")

adjacency_mun_df$y <- factor(adjacency_mun_df$y, 
                             levels = rev(unique(adjacency_mun_df$y)))

base_plot2 <- ggplot(adjacency_mun_df, aes(x = x, y = y, fill = p)) +
  geom_raster() + 
  coord_fixed() + 
  theme_minimal() +
  labs(title = "Heat Map- Adjacency Matrix",subtitle = "Municipal relationship", x = NULL, y = NULL,fill="Probability")+theme(
    plot.title = element_text(
      color = "#7C3B5E",   
      size  = 18,          
      face  = "bold",      
      hjust = 0.5          
    ),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
  )
base_plot2<-base_plot2 +
  scale_fill_continuous_sequential(palette = "BurgYl")+
  scale_x_discrete(labels = NULL) + 
  theme(axis.text.x = element_blank())+
  scale_y_discrete(labels = NULL) + 
  theme(axis.text.y = element_blank())



#Cluster estimation based on adjacency matrix 

clust_label_mun<-Mclust(adjacency_matrix_mun3,verbose = FALSE)
clusters_df_mun<-data.frame(
  cluster=clust_label_mun$classification,
  post_mean=df_mun$media_post[idx_order_mun$idx],
  mun_cod=df_mun$estu_cod_reside_mcpio[idx_order_mun$idx]
)

municipios<-st_read("/Users/macbookpro/Desktop/Tesis/shapes_municipios/MGN_ADM_MPIO_GRAFICO.shp")

municipios<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

# As its necessary only the information about the Municipios whose information is available
# their shapefiles are selected. 

municipios2<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  left_join(clusters_df_mun,by=c("mpio_cdpmp"="mun_cod"))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

municipios2$cluster<-ifelse(is.na(municipios2$cluster),9,municipios2$cluster)


# Dowload the other countries silhoutes
mundo <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = municipios2, 
          aes(fill = as.factor(cluster)), 
          color = "#BEBEBE", 
          lwd = 0.025) +
  
  # Color
  scale_fill_discrete_sequential(
    palette = "BurgYl", 
    rev = FALSE,
    name = "Cluster"
  ) +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "Municipal Segmentation",subtitle = "Clustering based on posterior mean estimates") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# Cluster analysis

tapply(clusters_df_mun$post_mean, clusters_df_mun$cluster, summary)

ggplot(clusters_df_mun, aes(x = as.factor(cluster), y = post_mean, fill = as.factor(cluster), colour = as.factor(cluster))) +
  geom_boxplot(
    alpha = 0.7,          # Transparencia media
    outlier.shape = NA   # Ocultamos outliers  # Borde gris oscuro uniforme
  ) +
  geom_jitter(
    aes(color = as.factor(cluster)), 
    width = 0.1, 
    alpha = 0.25
  ) +
  scale_color_discrete_sequential(
    palette = "BurgYl",
    rev = FALSE,
    name = "Cluster"
  )+
  scale_fill_discrete_sequential(
    palette="BurgYl",
    rev=FALSE,
    name="Cluster"
  )+
  theme_minimal()+
  labs(title = "Municipal posterior mean",subtitle = "Summary across clusters", x = NULL, y = NULL,fill="Probability")+theme(
    plot.title = element_text(
      color = "#7C3B5E",   
      size  = 18,          
      face  = "bold",     
      hjust = 0.5          
    ),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
  )


# =======  Spatial Cluster ===========

load("~/Desktop/Tesis/Results/MCMC_model2_real_data_standarized.RData")
beta<-THETA[[1]]
beta_E<-THETA[[2]]
beta_M<-THETA[[3]]
beta_D<-THETA[[4]]

sigma2_beta<-THETA[[5]]
lambda2_E<-THETA[[6]]
lambda2_M<-THETA[[7]]
lambda2_D<-THETA[[8]]

phi<-THETA[[9]]
tau_phi<-THETA[[10]]

kappa2_jk<-THETA[[11]]
kappa2_k<-THETA[[12]]

alpha_kappa<-THETA[[13]]
beta_kappa<-THETA[[14]]

acr<-THETA[[15]]

nsams<-length(beta)
n_clusters<-rep(NA,nsams)
adjacency_matrix_spatial<-matrix(0,nrow = m,ncol = m)
k_range <- 2:10

for (i in c(1:nsams)) {
  data<-cbind(phi[i,])
  
  dist_matrix <- dist(data) 
  
  # 2. Calcular Silhouette (Bucle interno optimizado)
  avg_sil_width <- numeric(length(k_range))
  
  for (k_idx in seq_along(k_range)) {
    k <- k_range[k_idx]
    # nstart bajo aquí para velocidad de decisión
    km_temp <- kmeans(data, centers = k, nstart = 10) 
    sil <- silhouette(km_temp$cluster, dist_matrix)
    avg_sil_width[k_idx] <- mean(sil[, 3])
  }
  
  # 3. Seleccionar ganador
  best_k_idx <- which.max(avg_sil_width)
  n_clusters[i] <- k_range[best_k_idx]
  
  k_means<-kmeans(data,centers = n_clusters[i])
  
  adjacency_matrix_iter<-matrix(NA,nrow = m,ncol = m)
  
  for (j in c(1:m)) {
    adjacency_matrix_iter[j,]<-c(rep(0,j),as.numeric(k_means$cluster[j]==k_means$cluster[-c(1:j)])*1/nsams)
  }
  
  adjacency_matrix_spatial<-adjacency_matrix_spatial+adjacency_matrix_iter
  
  
  # Advance
  if (i %% ceiling(nsams / 10) == 0) {
    cat(paste0("Iteración ", i, " de ", nsams, " (", round(100 * i / nsams), "%)\n"))
  }
  
}


adjacency_matrix_spatial_1<-adjacency_matrix_spatial+diag(1,nrow = m)+t(adjacency_matrix_spatial)

spatial_adj_matrix<-list(adjacency_matrix_spatial_1)
save(spatial_adj_matrix,file="Spatial_adj_matrix.RData")

colnames(adjacency_matrix_spatial_1)<-municipios_used$estu_mcpio_reside
rownames(adjacency_matrix_spatial_1)<-municipios_used$estu_mcpio_reside


idx_order_mun<-data.frame(
  idx=c(1:m),
  post_mean=colMeans(phi)
)%>%
  arrange(desc(post_mean))


adjacency_matrix_spatial_1<-adjacency_matrix_spatial_1[idx_order_mun$idx,idx_order_mun$idx]


adjacency_mun_spatial_df <- melt(adjacency_matrix_spatial_1)
names(adjacency_mun_spatial_df) <- c("x", "y", "p")

adjacency_mun_spatial_df$y <- factor(adjacency_mun_spatial_df$y, 
                                     levels = rev(unique(adjacency_mun_spatial_df$y)))

base_plot2 <- ggplot(adjacency_mun_spatial_df, aes(x = x, y = y, fill = p)) +
  geom_raster() + 
  coord_fixed() + 
  theme_minimal() +
  labs(title = "Heat Map – Adjacency Matrix",subtitle = "Municipal Spatial Relationship", x = NULL, y = NULL,fill="Probability")+theme(
    plot.title = element_text(
      color = "#7C3B5E",   
      size  = 18,          
      face  = "bold",      
      hjust = 0.5          
    ),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
  )
base_plot2<-base_plot2 +
  scale_fill_continuous_sequential(palette = "BurgYl")+
  scale_x_discrete(labels = NULL) + 
  theme(axis.text.x = element_blank())+
  scale_y_discrete(labels = NULL) + 
  theme(axis.text.y = element_blank())





#Cluster estimation based on adjacency matrix 

clust_label_mun<-Mclust(adjacency_matrix_spatial_1,verbose = FALSE)
clusters_df_mun<-data.frame(
  cluster=clust_label_mun$classification,
  post_mean=idx_order_mun$post_mean,
  mun_cod=municipios_used$estu_cod_reside_mcpio[idx_order_mun$idx]
)

municipios<-st_read("/Users/macbookpro/Desktop/Tesis/shapes_municipios/MGN_ADM_MPIO_GRAFICO.shp")

municipios<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

# As its necessary only the information about the Municipios whose information is available
# their shapefiles are selected. 

municipios2<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  left_join(clusters_df_mun,by=c("mpio_cdpmp"="mun_cod"))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

municipios2$cluster<-ifelse(is.na(municipios2$cluster),3,municipios2$cluster)


# Dowload the other countries silhoutes
mundo <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = municipios2, 
          aes(fill = as.factor(cluster)), 
          color = "#BEBEBE", 
          lwd = 0.025) +
  
  # Color
  scale_fill_discrete_sequential(
    palette = "BurgYl", 
    rev = FALSE,
    name = "Cluster"
  ) +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "Municipal Segmentation",subtitle = "Clustering based on spatial random effects") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )



ggplot(clusters_df_mun, aes(x = as.factor(cluster), y = post_mean, fill = as.factor(cluster), colour = as.factor(cluster))) +
  geom_boxplot(
    alpha = 0.7,          # Transparencia media
    outlier.shape = NA   # Ocultamos outliers  # Borde gris oscuro uniforme
  ) +
  geom_jitter(
    aes(color = as.factor(cluster)), 
    width = 0.1, 
    alpha = 0.25
  ) +
  scale_color_discrete_sequential(
    palette = "BurgYl",
    rev = FALSE,
    name = "Cluster"
  )+
  scale_fill_discrete_sequential(
    palette="BurgYl",
    rev=FALSE,
    name="Cluster"
  )+
  theme_minimal()+
  labs(title = "Municipal Spatial Random Effect",subtitle = "Summary across clusters", x = NULL, y = NULL,fill="Probability")+theme(
    plot.title = element_text(
      color = "#7C3B5E",   
      size  = 18,          
      face  = "bold",     
      hjust = 0.5          
    ),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
  )

tapply(clusters_df_mun$post_mean, clusters_df_mun$cluster, summary)
# ===== Random spatial effect- posterior mean=======

phi_mean<-data.frame(
  mean_phi=colMeans(phi),
  cod_mun=municipios_used$estu_cod_reside_mcpio
)


municipios<-st_read("/Users/macbookpro/Desktop/Tesis/shapes_municipios/MGN_ADM_MPIO_GRAFICO.shp")

municipios<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

# As its necessary only the information about the Municipios whose information is available
# their shapefiles are selected. 

municipios2<-municipios%>%
  filter(dpto_ccdgo!=88)%>%
  mutate(mpio_cdpmp=as.numeric(mpio_cdpmp))%>%
  left_join(phi_mean,by=c("mpio_cdpmp"="cod_mun"))%>%
  arrange(dpto_ccdgo,mpio_ccdgo)

municipios2$mean_phi<-ifelse(is.na(municipios2$mean_phi),mean(phi_mean$mean_phi)+0.05,municipios2$mean_phi)


# Dowload the other countries silhoutes
mundo <- ne_countries(scale = "medium", returnclass = "sf")
mi_paleta_personalizada <- c(
  "#70284a",
  "#823451",
  "#944059",# Extremo Burdeos oscuro
  "#c37177",
  "#ce8886",
  "#e5b7a6",
  "#fbe6c5", # Centro neutro (gris muy claro)
  "#dfcbc9",
  "#bcb0ce",
  "#9d96d2",
  "#7e7bd6",
  "#5e60db",
  "#3f45df"

  # Extremo Salmón/Melocotón
)

ggplot() +
  #Background layer
  geom_sf(data = mundo, 
          fill = "#E2E2E2",  # Land color
          color = NA) +      # No borders
  
  #Data layer
  geom_sf(data = municipios2, 
          aes(fill = mean_phi), 
          color = "#BEBEBE", 
          lwd = 0.025) +
  
  # Color
  scale_fill_gradientn(
    colors = mi_paleta_personalizada,
    name = "Posterior \n mean",
    # Opcional: si quieres que el color central (#eeeeee) coincida 
    # exactamente con el 0 de tus datos, descomenta la línea de abajo:
     values = scales::rescale(c(min(municipios2$mean_phi), 0, max(municipios2$mean_phi)))
  ) +
  
  # --- ZOOM  ---
  # Doing the zoom over Colombia
  coord_sf(
    xlim = c(-79.2, -66.8), 
    ylim = c(-4.3, 12.6),   
    expand = FALSE          
  ) +
  
  # --- AESTHETIC ---
  labs(title = "Spatial Random Effects",subtitle = "Municipal level posterior means") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color del Mar
    plot.title = element_text(color = "#7C3B5E", face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 12,hjust = 0.5 ),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )



m<-list(adjacency_matrix_mun=adjacency_matrix_mun3,adjacency_matrix_dep=adjacency_matrix_dep,df_mun=df_mun,df_dep=df_dep,dep_post_mean=dep_post_mean,mun_post_mean=mun_post_mean)
save(m,file="/Users/macbookpro/Desktop/Tesis/Results/Matriz_adjacencia_departamental.RData")
load("/Users/macbookpro/Desktop/Tesis/Results/Matriz_adjacencia_departamental.RData")

