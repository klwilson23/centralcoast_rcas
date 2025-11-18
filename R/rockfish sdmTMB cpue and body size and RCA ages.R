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

quadroot <- function(x) {return((x)^(1/4))}
wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")


rcas <- st_read("Data/RCAs_2019_BCAlbers/RCAs_BCAlbers.shp")
coast_line <- readRDS("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_boundary.rds")
fit_sdm <- readRDS(file="R/sdm_groundfish_re6.rds")
mesh <- readRDS(file="R/sdm_mesh.rds")
pfma <- st_read("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Spatial data/PFMA/DFO_PFMA_SUBAREAS_SP/DFO_SBAREA_polygon.shp")
rockfish <- readRDS(file="R/sdm_model_data.rds")

# make plot for RCA ages
rca_distkm <- c(0,25)
year_seq <- c(2007,2010,2013,2016,2019,2022)-2004
rmax_seq_fact <- c(0.043,0.096,0.142)

# plot 1:
rca_dist_seq <- (-1*(rca_distkm-mean(rockfish$rca_distkm))/sd(rockfish$rca_distkm))
year_protected_seq <- (year_seq-mean(rockfish$year_protected))/sd(rockfish$year_protected)
rmax_seq <- (rmax_seq_fact-mean(rockfish$rmax))/sd(rockfish$rmax)


cpue_trends_ages <- emmeans(fit_sdm, ~ year_protected_scaled+rmax_scale+protection_dist_scale,at = list(year_protected_scaled = year_protected_seq,rmax_scale=rmax_seq,protection_dist_scale=rca_dist_seq))
# emm_df <- as.data.frame(summary(size_trends_ages))
cpue_trends <- broom::tidy(cpue_trends_ages,conf.int = TRUE)

df <- data.frame("RCA_age"=rep(rep(year_seq,times=length(rmax_seq_fact)),times=length(rca_dist_seq)),
                 "rmax"=rep(rep(rmax_seq_fact,each=length(year_seq)),times=length(rca_dist_seq)),
                 "RCA_distance"=rep(rep(rca_distkm,each=length(year_seq)),each=length(rmax_seq_fact)),
                 broom::tidy(cpue_trends_ages,conf.int = TRUE))

df$rmax <- as.factor(df$rmax)
df$RCA_distance <- factor(df$RCA_distance,levels=c(0,25),labels=c("Inside","25 km away"))

p1 <- ggplot(data=df,aes(x=RCA_age,y=estimate,colour=RCA_distance,fill=RCA_distance,group=RCA_distance)) +
  geom_ribbon(data=df,alpha=0.5,aes(x=RCA_age,ymin=conf.low,ymax=conf.high,group=RCA_distance)) +
  geom_line(data=df,aes(x=RCA_age,y=estimate,colour=RCA_distance,group=RCA_distance))+
  facet_grid(cols=vars(rmax),labeller=label_bquote(cols = r[max]==.(as.character(rmax))))+
  xlab("") + ylab("log(CPUE)") +
  theme_minimal() +
  scale_fill_brewer(name="Proximity to nearest RCA",palette=2,type="qual")+
  scale_color_brewer(name="Proximity to nearest RCA",palette=2,type="qual")+
  theme(legend.position="top",text = element_text(size=11),legend.text = element_text(size=10),strip.text=element_text(size=10,hjust=0))
print(p1)

# now body size
fit_sdm <- readRDS(file="R/sdm_groundfish_size_re.rds")
rockfish <- readRDS(file="R/sdm_model_data_size.rds")

rca_dist_seq <- (-1*(rca_distkm-mean(rockfish$rca_distkm))/sd(rockfish$rca_distkm))
year_protected_seq <- (year_seq-mean(rockfish$year_protected))/sd(rockfish$year_protected)
rmax_seq <- (rmax_seq_fact-mean(rockfish$rmax))/sd(rockfish$rmax)


size_trends_ages <- emmeans(fit_sdm, ~ year_protected_scaled+rmax_scale+protection_dist_scale,at = list(year_protected_scaled = year_protected_seq,rmax_scale=rmax_seq,protection_dist_scale=rca_dist_seq))
# emm_df <- as.data.frame(summary(size_trends_ages))
size_trends <- broom::tidy(size_trends_ages,conf.int = TRUE)

df <- data.frame("RCA_age"=rep(rep(year_seq,times=length(rmax_seq_fact)),times=length(rca_dist_seq)),
                 "rmax"=rep(rep(rmax_seq_fact,each=length(year_seq)),times=length(rca_dist_seq)),
                 "RCA_distance"=rep(rep(rca_distkm,each=length(year_seq)),each=length(rmax_seq_fact)),
                 broom::tidy(size_trends_ages,conf.int = TRUE))

df$rmax <- as.factor(df$rmax)
df$RCA_distance <- factor(df$RCA_distance,levels=c(0,25),labels=c("Inside","25 km away"))

p2 <- ggplot(data=df,aes(x=RCA_age,y=estimate,colour=RCA_distance,fill=RCA_distance,group=RCA_distance)) +
  geom_ribbon(data=df,alpha=0.5,aes(x=RCA_age,ymin=conf.low,ymax=conf.high,group=RCA_distance)) +
  geom_line(data=df,aes(x=RCA_age,y=estimate,colour=RCA_distance,group=RCA_distance))+
  facet_grid(cols=vars(rmax),labeller=label_bquote(cols = r[max]==.(as.character(rmax))))+
  xlab("Years since RCA establishment") + ylab("log(total length)") +
  theme_minimal() +
  scale_fill_brewer(name="Proximity to nearest RCA",palette=2,type="qual")+
  scale_color_brewer(name="Proximity to nearest RCA",palette=2,type="qual")+
  theme(legend.position="top",text = element_text(size=11),legend.text = element_text(size=10),strip.text=element_text(size=10,hjust=0))
print(p2)


megaP <- ggpubr::ggarrange(p1,p2,nrow=2,heights=c(1,1),labels=c("(a)","(b)"),common.legend = TRUE)

ggsave(filename = "Figures/RCA ages and cpue and body size.jpeg",plot=megaP,dpi=600,units="in",height=5.5,width=6)
