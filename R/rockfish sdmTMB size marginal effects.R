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


wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")

rcas <- st_read("Data/RCAs_2019_BCAlbers/RCAs_BCAlbers.shp")
coast_line <- readRDS("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_boundary.rds")
fit_sdm <- readRDS(file="R/sdm_groundfish_size_re.rds")
# fit_sdm <- sdmTMB:::reload_model(fit_sdm)
mesh <- readRDS(file="R/sdm_mesh_size.rds")
pfma <- st_read("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Spatial data/PFMA/DFO_PFMA_SUBAREAS_SP/DFO_SBAREA_polygon.shp")
rockfish <- readRDS(file="R/sdm_model_data_size.rds")
grid <-readRDS("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_rockfish_grid_4km.rds")
cell_size <- 4000
grid$protected <- ifelse(grid$rca_distkm==0,1,0)
grid$protection_dist <- grid$protected*grid$rca_distkm
grid$protection_dist_scale <- (-1*(grid$rca_distkm-mean(rockfish$rca_distkm))/sd(rockfish$rca_distkm))
grid$year_protected <- ifelse(grid$rca_name%in%c("Fish Egg Inlet","Goose Island","McMullin Group","Smith Sound"),2004,2005)
grid$upperoceansr <- relevel(grid$upperoceansr,ref="(11) Mainland Fjords")

new_grid <- expand.grid("puid"="12411",
                        "year"=2006,
                        "species"="quillback",
                        "gear"="dive",
                        "effort"=120,
                        "chase_prob"=0,
                        "depth_seq"=seq(5,305,by=25),
                        "cumul_catch_before_rca_scaled" = 0,
                        "week"=32,
                        "depth_cat"=levels(rockfish$depth_cat))
new_grid <- merge(new_grid,grid,by="puid",all.x=TRUE)
new_grid$log_depth <- log(new_grid$depth_seq)
new_grid$species <- factor(new_grid$species,levels=levels(rockfish$species))
new_grid$year_protected <- 2005
new_grid$year_protected <- new_grid$year-new_grid$year_protected
new_grid$year_protected_scaled <- (new_grid$year_protected-mean(rockfish$year_protected))/(sd(rockfish$year_protected))
new_grid$biweekly <- ceiling(new_grid$week / 2)
new_grid$biweekly_scale <- (new_grid$biweekly-mean(rockfish$biweekly))/sd(rockfish$biweekly)
# new_grid$depth_cat <- rockfish$depth_cat[match(new_grid$species,rockfish$species)]
new_grid$depth_cat <- factor(new_grid$depth_cat,levels=levels(rockfish$depth_cat))
new_grid$complexity_scaled <- 0
new_grid$rmax_scale <- rockfish$rmax_scale[match(new_grid$species,rockfish$species)]
new_grid$rca_name <- factor(new_grid$rca_name,levels=levels(rockfish$rca_name))
new_grid <- st_as_sf(new_grid,coords=c("X","Y"),crs=sf::st_crs(coast_line))
new_grid$body_size <- 1
p <- predict(fit_sdm, newdata = as.data.frame(new_grid),se_fit=TRUE,type="link",re_form=~0,re_form_iid=~0)

p$depth_cat <- factor(p$depth_cat,levels=levels(p$depth_cat),labels=c("quillback","china","copper/black/deacon/lingcod","yelloweye","vermillion/widow","canary","puget_sound","yellowtail","tiger","dusky-dark"))
table(rockfish$species)
ggplot(data=p,aes(x=depth_seq,y=est,colour=depth_cat)) +
  geom_ribbon(data=p,alpha=0.5,aes(x=depth_seq,ymin=est-1.96*p$est_se,ymax=est+1.96*p$est_se,fill=depth_cat)) +
  geom_line(data=p,aes(x=depth_seq,y=est,colour=depth_cat))+
  facet_wrap(~depth_cat,scales="free_y") +
  xlab("Depth (m)") + ylab("ln(total length)") +
  scale_fill_viridis_d(name="Group") +
  scale_colour_viridis_d(name="Group") +
  theme_minimal() + theme(legend.position="top") +
  theme(text = element_text(size=10),legend.text = element_text(size=8),strip.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave(filename = "Figures/Depth size smoothers.jpeg",dpi=600,units="in",height=5,width=7)
