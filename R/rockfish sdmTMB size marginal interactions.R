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

quadroot <- function(x) {return((x)^(1/4))}

wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")

rcas <- st_read("Data/RCAs_2019_BCAlbers/RCAs_BCAlbers.shp")
coast_line <- readRDS("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_boundary.rds")
fit_sdm <- readRDS(file="R/sdm_groundfish_size_re.rds")
fit_sdm <- sdmTMB:::reload_model(fit_sdm)
mesh <- readRDS(file="R/sdm_mesh_size.rds")
pfma <- st_read("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Spatial data/PFMA/DFO_PFMA_SUBAREAS_SP/DFO_SBAREA_polygon.shp")
rockfish <- readRDS(file="R/sdm_model_data_size.rds")
grid <-readRDS("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_rockfish_grid_4km.rds")
cell_size <- 4000
grid$protected <- ifelse(grid$rca_distkm==0,1,0)
grid$year_protected <- ifelse(grid$rca_name%in%c("Fish Egg Inlet","Goose Island","McMullin Group","Smith Sound"),2004,2005)
grid$upperoceansr <- relevel(grid$upperoceansr,ref="(11) Mainland Fjords")
range(rockfish$cumulative_catch_before_rca)
quantile(rockfish$cumulative_catch_before_rca,probs=c(0,0.025,0.25,0.5,0.75,0.975,1))

range(rockfish$year_protected_scaled)
range(rockfish$protection_dist_scale)
range(rockfish$rmax_scale)
min_kg <- min(rockfish$cumulative_catch[rockfish$cumulative_catch!=0])
min_before_kg <- min(rockfish$cumulative_catch_before_rca[rockfish$cumulative_catch_before_rca!=0])
min_after_kg <- min(rockfish$cumulative_catch_after_rca[rockfish$cumulative_catch_after_rca!=0])


new_grid <- expand.grid("puid"="12411",
                        "species"="quillback",
                        "gear"="dive",
                        "effort"=120,
                        "chase_prob"=0,
                        "depth_seq"=50,
                        "cumulative_catch_before_rca" = c(0,10,100),
                        "year" = c(2010,2015,2023),
                        "rca_distkm" = seq(from=0,to=20,length=20),
                        # "protection_dist_scale" = seq(-2,2,by=0.5),
                        "rmax" = c(0.043,0.096,0.142), # yelloweye, widow, lingcod
                        "week"=32)
new_grid <- merge(new_grid,grid,by="puid",all.x=TRUE)
new_grid$rca_distkm <- new_grid$rca_distkm.x
new_grid$log_depth <- log(new_grid$depth_seq)
new_grid$species <- factor(new_grid$species,levels=levels(rockfish$species))
new_grid$biweekly <- ceiling(new_grid$week / 2)
new_grid$biweekly_scale <- (new_grid$biweekly-mean(rockfish$biweekly))/sd(rockfish$biweekly)
new_grid$depth_cat <- rockfish$depth_cat[match(new_grid$species,rockfish$species)]
new_grid$depth_cat <- factor(new_grid$depth_cat,levels=levels(rockfish$depth_cat))
new_grid$complexity_scaled <- 0
new_grid$rca_name <- factor(new_grid$rca_name,levels=levels(rockfish$rca_name))
new_grid$protection_dist_scale <- (-1*(new_grid$rca_distkm-mean(rockfish$rca_distkm))/sd(rockfish$rca_distkm))
new_grid$year_protected <- ifelse(new_grid$rca_name%in%c("Fish Egg Inlet","Goose Island","McMullin Group","Smith Sound"),2004,2005)
new_grid$year_protected <- new_grid$year-new_grid$year_protected
new_grid$year_protected_scaled <- (new_grid$year_protected-mean(rockfish$year_protected))/sd(rockfish$year_protected)
new_grid$rmax_scale <- (new_grid$rmax-mean(rockfish$rmax))/sd(rockfish$rmax)
new_grid$cumul_catch_before_rca_scaled <- (quadroot(new_grid$cumulative_catch_before_rca) - mean(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%new_grid$species])))/sd(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%new_grid$species]))
new_grid <- st_as_sf(new_grid,coords=c("X","Y"),crs=sf::st_crs(coast_line))
new_grid$body_size <- 1
p <- predict(fit_sdm, newdata = as.data.frame(new_grid),se_fit=TRUE,type="link",re_form=~0,re_form_iid=~0)

p$rmax_cat <- factor(p$rmax,levels=unique(p$rmax),labels=c("Low","Mid","High"))
p$harv_cat <- factor(p$cumulative_catch_before_rca,levels=sort(unique(p$cumulative_catch_before_rca)),labels=c("Unexploited","10 kg","100 kg"))
p$year_cat <- factor(p$year,levels=sort(unique(p$year)),labels=c("5-years","10-years","18-years"))

p_sub <- p#p[p$year_cat!="Mid",]
ggplot(data=p_sub,aes(x=-1*rca_distkm,y=est,group=interaction(year_cat),colour = year_cat)) +
  geom_ribbon(data=p_sub,alpha=0.2,aes(x=-1*rca_distkm,ymin=est-1.0*est_se,ymax=est+1.0*est_se,group=interaction(year_cat),fill = year_cat)) +
  geom_line(data=p_sub,aes(x=-1*rca_distkm,y=est,group=interaction(year_cat),colour = year_cat))+
  facet_grid(rows=vars(harv_cat),cols=vars(rmax_cat),labeller = label_bquote(cols = .(as.character(rmax_cat))~r[max],rows=.(harv_cat)))+
  xlab("Proximity to nearest RCA (km)") + ylab("ln(total length)") +
  theme_minimal() + theme(legend.position="top") +
  theme(text = element_text(size=11),legend.text = element_text(size=10),strip.text=element_text(size=9)) +
  scale_colour_brewer(name="RCA Age",type="qual",palette=2) +  scale_fill_brewer(name="RCA Age",type="qual",palette=2)

ggsave(filename = "Figures/effect of RCA interactions size.jpeg",dpi=600,units="in",height=5,width=6)


new_grid <- expand.grid("puid"="12411",
                        "species"="quillback",
                        "gear"="dive",
                        "effort"=120,
                        "chase_prob"=0,
                        "depth_seq"=50,
                        "cumulative_catch_before_rca" = seq(0,100,length=20),
                        "year" = c(2010,2015,2023),
                        "rca_distkm" = c(0,5,20),
                        "rmax" = c(0.043,0.096,0.142), # yelloweye, lingcod, brown
                        "week"=32)
new_grid <- merge(new_grid,grid,by="puid",all.x=TRUE)
new_grid$rca_distkm <- new_grid$rca_distkm.x
new_grid$log_depth <- log(new_grid$depth_seq)
new_grid$species <- factor(new_grid$species,levels=levels(rockfish$species))
new_grid$biweekly <- ceiling(new_grid$week / 2)
new_grid$biweekly_scale <- (new_grid$biweekly-mean(rockfish$biweekly))/sd(rockfish$biweekly)
new_grid$depth_cat <- rockfish$depth_cat[match(new_grid$species,rockfish$species)]
new_grid$depth_cat <- factor(new_grid$depth_cat,levels=levels(rockfish$depth_cat))
new_grid$complexity_scaled <- 0
new_grid$rca_name <- factor(new_grid$rca_name,levels=levels(rockfish$rca_name))
new_grid$protection_dist_scale <- (-1*(new_grid$rca_distkm-mean(rockfish$rca_distkm))/sd(rockfish$rca_distkm))
new_grid$year_protected <- ifelse(new_grid$rca_name%in%c("Fish Egg Inlet","Goose Island","McMullin Group","Smith Sound"),2004,2005)
new_grid$year_protected <- new_grid$year-new_grid$year_protected
new_grid$year_protected_scaled <- (new_grid$year_protected-mean(rockfish$year_protected))/sd(rockfish$year_protected)
new_grid$rmax_scale <- (new_grid$rmax-mean(rockfish$rmax))/sd(rockfish$rmax)
new_grid$cumul_catch_before_rca_scaled <- (quadroot(new_grid$cumulative_catch_before_rca) - mean(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%new_grid$species])))/sd(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%new_grid$species]))
new_grid <- st_as_sf(new_grid,coords=c("X","Y"),crs=sf::st_crs(coast_line))
new_grid$body_size <- 1
p <- predict(fit_sdm, newdata = as.data.frame(new_grid),se_fit=TRUE,type="link",re_form=~0,re_form_iid=~0)

p$rmax_cat <- factor(p$rmax,levels=unique(p$rmax),labels=c("Low","Mid","High"))
p$year_cat <- factor(p$year,levels=sort(unique(p$year)),labels=c("5-years","10-years","18-years"))
p$rca_cat <- factor(p$rca_distkm,levels=sort(unique(p$rca_distkm)),labels=c("Inside RCA","5km from RCA","20km from RCA"))

p_sub <- p#p[p$year_cat!="Mid",]
ggplot(data=p_sub,aes(x=cumulative_catch_before_rca,y=est,group=interaction(year_cat),colour = year_cat)) +
  geom_ribbon(data=p_sub,alpha=0.2,aes(x=cumulative_catch_before_rca,ymin=est-1.0*est_se,ymax=est+1.0*est_se,group=interaction(year_cat),fill = year_cat)) +
  geom_line(data=p_sub,aes(x=cumulative_catch_before_rca,y=est,group=interaction(year_cat),colour = year_cat))+
  facet_grid(rows=vars(rca_cat),cols=vars(rmax_cat),labeller = label_bquote(cols = .(as.character(rmax_cat))~r[max],rows=.(rca_cat)))+
  xlab("Cumulative commercial harvest (kg)") + ylab("ln(total length)") +
  theme_minimal() + theme(legend.position="top") +
  theme(text = element_text(size=11),legend.text = element_text(size=10),strip.text=element_text(size=9)) +
  scale_colour_brewer(name="RCA Age",type="qual",palette=2) + scale_fill_brewer(name="RCA Age",type="qual",palette=2)

ggsave(filename = "Figures/effect of harvest interactions size.jpeg",dpi=600,units="in",height=5,width=6)


new_grid <- expand.grid("puid"="12411",
                        "species"="quillback",
                        "gear"="dive",
                        "effort"=120,
                        "chase_prob"=0,
                        "depth_seq"=50,
                        "cumulative_catch_before_rca" = c(0,10,100),
                        "year" = seq(2006,2023,length=20),
                        "rca_distkm" = c(0,5,20),
                        "rmax" = c(0.043,0.096,0.142), # yelloweye, lingcod, brown
                        "week"=32)
new_grid <- merge(new_grid,grid,by="puid",all.x=TRUE)
new_grid$rca_distkm <- new_grid$rca_distkm.x
new_grid$log_depth <- log(new_grid$depth_seq)
new_grid$species <- factor(new_grid$species,levels=levels(rockfish$species))
new_grid$biweekly <- ceiling(new_grid$week / 2)
new_grid$biweekly_scale <- (new_grid$biweekly-mean(rockfish$biweekly))/sd(rockfish$biweekly)
new_grid$depth_cat <- rockfish$depth_cat[match(new_grid$species,rockfish$species)]
new_grid$depth_cat <- factor(new_grid$depth_cat,levels=levels(rockfish$depth_cat))
new_grid$complexity_scaled <- 0
new_grid$rca_name <- factor(new_grid$rca_name,levels=levels(rockfish$rca_name))
new_grid$protection_dist_scale <- (-1*(new_grid$rca_distkm-mean(rockfish$rca_distkm))/sd(rockfish$rca_distkm))
new_grid$year_protected <- ifelse(new_grid$rca_name%in%c("Fish Egg Inlet","Goose Island","McMullin Group","Smith Sound"),2004,2005)
new_grid$year_protected <- new_grid$year-new_grid$year_protected
new_grid$year_protected_scaled <- (new_grid$year_protected-mean(rockfish$year_protected))/sd(rockfish$year_protected)
new_grid$rmax_scale <- (new_grid$rmax-mean(rockfish$rmax))/sd(rockfish$rmax)
new_grid$cumul_catch_before_rca_scaled <- (quadroot(new_grid$cumulative_catch_before_rca) - mean(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%new_grid$species])))/sd(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%new_grid$species]))
new_grid <- st_as_sf(new_grid,coords=c("X","Y"),crs=sf::st_crs(coast_line))
new_grid$body_size <- 1

p <- predict(fit_sdm, newdata = as.data.frame(new_grid),se_fit=TRUE,type="link",re_form=~0,re_form_iid=~0)

p$rmax_cat <- factor(p$rmax,levels=unique(p$rmax),labels=c("Low","Mid","High"))
p$rca_cat <- factor(p$rca_distkm,levels=sort(unique(p$rca_distkm)),labels=c("Inside RCA","5km from RCA","20km from RCA"))
p$harv_cat <- factor(p$cumul_catch_before_rca,levels=sort(unique(p$cumul_catch_before_rca)),labels=c("Unexploited","10 kg","100 kg"))

p_sub <- p#p[p$year_cat!="Mid",]
ggplot(data=p_sub,aes(x=year_protected,y=est,group=interaction(rca_cat),colour = rca_cat)) +
  geom_ribbon(data=p_sub,alpha=0.2,aes(x=year_protected,ymin=est-1.0*est_se,ymax=est+1.0*est_se,group=interaction(rca_cat),fill = rca_cat)) +
  geom_line(data=p_sub,aes(x=year_protected,y=est,group=interaction(rca_cat),colour = rca_cat))+
  facet_grid(rows=vars(harv_cat),cols=vars(rmax_cat),labeller = label_bquote(cols = .(as.character(rmax_cat))~r[max],rows=.(harv_cat)))+
  xlab("RCA Age (years since protection)") + ylab("ln(total length)") +
  theme_minimal() + theme(legend.position="top") +
  theme(text = element_text(size=11),legend.text = element_text(size=10),strip.text=element_text(size=9)) +
  scale_colour_brewer(name="Proximity to RCA",type="qual",palette=2) + scale_fill_brewer(name="Proximity to RCA",type="qual",palette=2)

ggsave(filename = "Figures/effect of annual trend interactions size.jpeg",dpi=600,units="in",height=5,width=6)

new_grid <- expand.grid("puid"="12411",
                        "species"="quillback",
                        "gear"="dive",
                        "effort"=120,
                        "chase_prob"=0,
                        "depth_seq"=50,
                        "cumulative_catch_before_rca" = c(0,10,100),
                        "year" = c(2010,2015,2023),
                        "rca_distkm" = c(0,5,20),
                        "rmax" = seq(min(rockfish$rmax),max(rockfish$rmax),length=20), # yelloweye, lingcod, brown
                        "week"=32)
new_grid <- merge(new_grid,grid,by="puid",all.x=TRUE)
new_grid$rca_distkm <- new_grid$rca_distkm.x
new_grid$log_depth <- log(new_grid$depth_seq)
new_grid$species <- factor(new_grid$species,levels=levels(rockfish$species))
new_grid$biweekly <- ceiling(new_grid$week / 2)
new_grid$biweekly_scale <- (new_grid$biweekly-mean(rockfish$biweekly))/sd(rockfish$biweekly)
new_grid$depth_cat <- rockfish$depth_cat[match(new_grid$species,rockfish$species)]
new_grid$depth_cat <- factor(new_grid$depth_cat,levels=levels(rockfish$depth_cat))
new_grid$complexity_scaled <- 0
new_grid$rca_name <- factor(new_grid$rca_name,levels=levels(rockfish$rca_name))
new_grid$protection_dist_scale <- (-1*(new_grid$rca_distkm-mean(rockfish$rca_distkm))/sd(rockfish$rca_distkm))
new_grid$year_protected <- ifelse(new_grid$rca_name%in%c("Fish Egg Inlet","Goose Island","McMullin Group","Smith Sound"),2004,2005)
new_grid$year_protected <- new_grid$year-new_grid$year_protected
new_grid$year_protected_scaled <- (new_grid$year_protected-mean(rockfish$year_protected))/sd(rockfish$year_protected)
new_grid$rmax_scale <- (new_grid$rmax-mean(rockfish$rmax))/sd(rockfish$rmax)
new_grid$cumul_catch_before_rca_scaled <- (quadroot(new_grid$cumulative_catch_before_rca) - mean(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%new_grid$species])))/sd(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%new_grid$species]))
new_grid <- st_as_sf(new_grid,coords=c("X","Y"),crs=sf::st_crs(coast_line))
new_grid$body_size <- 1

p <- predict(fit_sdm, newdata = as.data.frame(new_grid),se_fit=TRUE,type="link",re_form=~0,re_form_iid=~0)

p$rca_cat <- factor(p$rca_distkm,levels=sort(unique(p$rca_distkm)),labels=c("Inside RCA","5km from RCA","20km from RCA"))
p$harv_cat <- factor(p$cumulative_catch_before_rca,levels=sort(unique(p$cumulative_catch_before_rca)),labels=c("Unexploited","10 kg","100 kg"))
p$year_cat <- factor(p$year,levels=sort(unique(p$year)),labels=c("5-years","10-years","18-years"))

p_sub <- p#p[p$year_cat!="Mid",]
ggplot(data=p_sub,aes(x=rmax,y=est,group=interaction(year_cat),colour = year_cat)) +
  geom_ribbon(data=p_sub,alpha=0.2,aes(x=rmax,ymin=est-1.0*est_se,ymax=est+1.0*est_se,group=interaction(year_cat),fill = year_cat)) +
  geom_line(data=p_sub,aes(x=rmax,y=est,group=interaction(year_cat),colour = year_cat))+
  facet_grid(rows=vars(harv_cat),cols=vars(rca_cat))+
  xlab(expression(r[max])) + ylab("ln(total length)") +
  theme_minimal() + theme(legend.position="top") +
  theme(text = element_text(size=11),legend.text = element_text(size=10),strip.text=element_text(size=9)) +
  scale_colour_brewer(name="RCA Age",type="qual",palette=2) + scale_fill_brewer(name="RCA Age",type="qual",palette=2)
ggsave(filename = "Figures/effect of rmax interactions size.jpeg",dpi=600,units="in",height=5,width=6)

