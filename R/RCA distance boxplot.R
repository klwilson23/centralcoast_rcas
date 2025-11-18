library(ggplot2)
library(dplyr)
library(tidyr)

rockfish <- readRDS(file="R/sdm_model_data.rds")
rockfish$species <- gsub("vermillion","vermilion",rockfish$species)
rockfish$species <- gsub("sebastolobus","shortspine_thornyhead",rockfish$species)
range(rockfish$rca_distkm)
ggplot(data=rockfish,aes(x=species,y=rca_distkm)) +
  # geom_violin(trim = FALSE,fill="dodgerblue") + 
  geom_boxplot(fill="dodgerblue",outliers=FALSE,width=0.5) +
  ggbeeswarm::geom_quasirandom(shape = 21,fill="grey10",colour="white",size=0.7,alpha=0.5)+
  ylab("Distance from nearest RCA (km)") + xlab("Species") +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(filename="Figures/RCA distance boxplot.jpeg",dpi=600,units="in",width=8,height=6.5)
