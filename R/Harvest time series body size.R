rockfish <- readRDS(file="R/sdm_model_data_size.rds")
rockfish$cumul_bin_harvest
rockfish$rca_cutoff <- ifelse(rockfish$rca_distkm>=5,"Far","Close")
rockfish$rca_cutoff <- as.factor(rockfish$rca_cutoff)

rockfish$species <- gsub("vermillion","vermilion",rockfish$species)
rockfish$species <- gsub("sebastolobus","shortspine_thornyhead",rockfish$species)

harvest_timeseries <- aggregate(cbind(cumul_bin_harvest,cumul_bin_harvest_pre_rca,cumul_bin_harvest_post_rca)~year+species+rca_cutoff,data=rockfish,FUN=sum,na.rm=TRUE)
print(harvest_timeseries)
table(harvest_timeseries$rca_cutoff)
library(ggplot2)
library(dplyr)
library(tidyr)

ggplot(data=harvest_timeseries,aes(x=year,y=cumul_bin_harvest,colour=rca_cutoff)) +
  geom_line(data=harvest_timeseries,aes(x=year,y=cumul_bin_harvest,colour=rca_cutoff)) +
  facet_wrap(~species,scales="free_y") +
  ylab("Total Commercial Harvest (kg)") + xlab("Year") +
  scale_color_discrete(name="RCA Cutoff") +
  scale_y_continuous(limits = c(0,NA)) +
  theme_minimal() +
  theme(legend.position = "top")
ggsave(filename="Figures/Harvest timeseries body size.jpeg",dpi=600,units="in",width=12,height=7)


ggplot(data=harvest_timeseries,aes(x=year,y=cumul_bin_harvest_pre_rca,colour=rca_cutoff)) +
  geom_line(data=harvest_timeseries,aes(x=year,y=cumul_bin_harvest_pre_rca,colour=rca_cutoff)) +
  facet_wrap(~species,scales="free_y") +
  ylab("Total Commercial Harvest before RCA establishment (kg)") + xlab("Year") +
  scale_color_discrete(name="RCA Cutoff") +
  scale_y_continuous(limits = c(0,NA)) +
  theme_minimal() +
  theme(legend.position = "top")
ggsave(filename="Figures/Harvest timeseries pre RCA body size.jpeg",dpi=600,units="in",width=12,height=5)

rockfish$rca_cutoff <- ifelse(rockfish$rca_distkm>0,">1km away","Inside RCA or within 1km")
rockfish$rca_cutoff <- factor(rockfish$rca_cutoff,levels=c("Inside RCA or within 1km",">1km away"))

harvest_snapshot <- aggregate(cbind(cumulative_catch_before_rca,cumulative_catch_after_rca)~species+rca_name+rca_cutoff,data=rockfish,FUN=sum,na.rm=TRUE)
# harvest_snapshot <- pivot_longer(harvest_snapshot,cols="")
ggplot(data=rockfish,aes(x=cumulative_catch_before_rca,y=cumulative_catch_after_rca,colour=rca_cutoff)) +
  geom_point() +
  facet_wrap(~species,ncol=4,scales="free") +
  xlab("Total Commercial Harvest before RCA establishment (kg/km2)") + ylab("Total Commercial Harvest after RCA establishment (kg/km2)") +
  geom_abline(slope=1,intercept=0,lwd=1,linetype=3)+
  scale_color_discrete(name="Proximity to RCA") +
  scale_y_continuous(limits = c(0,NA),n.breaks = 4) + scale_x_continuous(limits = c(0,NA),n.breaks = 4) +
  theme_minimal() +
  theme(legend.position = "top",strip.text.x=element_text(hjust=0))
ggsave(filename="Figures/Harvest before and after RCA by species body size.jpeg",dpi=600,units="in",width=8,height=8)
