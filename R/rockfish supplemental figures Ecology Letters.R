library(dbplyr)
library(bcdata)
library(bcmaps)
library(sf)
library(sp)
library(ggspatial) # this is for adding the arrows
# library(rgdal) #use this for data conversion
library(ggrepel) # to offset labels using geom_sf_label_repel  --> Not done here
library(riverdist) # to snap points to River --> Not done here
# library(bcmapsdata)
library(viridis)
library(ggnewscale)
library(tidyr)
library(dplyr)
library(sdmTMB)
library(sdmTMBextra)
library(ggplot2)
#remotes::install_github("strengejacke/ggeffects")
library(ggeffects)
library(cowplot)
library(emmeans)

rockfish <- readRDS(file="R/sdm_model_data.rds")
rockfish$rca_distkm

ggplot(data=rockfish,aes(x=rca_distkm)) +
  geom_histogram(binwidth=2.5,colour="black",fill="grey50")+
  ylab("Frequency") + xlab("Distance to nearest RCA (km)") +
  ylim(c(0,20000))+
  theme_classic()
ggsave("Figures/figure s10 boxplot distance.jpeg",dpi=600,units="in",height=5,width=6.5)

quantile(rockfish$cumulative_catch_before_rca)
ggplot(data=rockfish,aes(x=cumulative_catch_before_rca)) +
  geom_histogram(colour="black",fill="grey50",boundary = 0,bins=25)+
  facet_wrap(~species,nrow=5,scale="free")+
  ylim(c(0, NA))+
  ylab("Frequency") + xlab("Cumulative harvest before RCA establishment (kg)") +
  theme_classic()
ggsave("Figures/figure s11 boxplot harvest.jpeg",dpi=600,units="in",height=7,width=10)
