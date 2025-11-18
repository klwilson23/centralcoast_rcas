# See: https://github.com/bcgov/bcmaps
# see: https://cran.r-project.org/web/packages/bcmaps/bcmaps.pdf
# availble layers in bcmaps:
# https://gist.github.com/ateucher/86674e12ffba66ecce87cceb7bffbf41
# https://github.com/poissonconsulting/fwabc#<https://github.com/poissonconsulting/fwabc>
# https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
mode <- function(x) {
  ux <- unique(na.omit(x))
  tx <- tabulate(match(x, ux))
  if(length(ux) != 1 & sum(max(tx) == tx) > 1) {
    if (is.character(ux)) return(NA_character_) else return(NA_real_)
  }
  max_tx <- tx == max(tx)
  return(ux[max_tx])
}


buffer <- 25000

library(dbplyr)
library(bcdata)
library(dplyr)
library(bcmaps)
library(sf)
library(sp)
library(ggplot2)
library(ggspatial) # this is for adding the arrows
library(rgdal) #use this for data conversion
library(ggrepel) # to offset labels using geom_sf_label_repel  --> Not done here
library(riverdist) # to snap points to River --> Not done here
library(bcmapsdata)
library(viridis)
library(ggnewscale)
library(tidyr)
ocean <- bc_neighbours(ask=FALSE) %>% 
  filter(type=="Ocean")

full_dat <- read.csv("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/CCIRA Tow Video/Full Transcription/R script merging all data/rockfish body sizes.csv",header=TRUE)
lh_traits <- read.csv("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/CCIRA Tow Video/Full Transcription/R script merging all data/species life histories.csv")
full_dat$max_age <- lh_traits$max..age..yrs.[match(full_dat$species,lh_traits$name.in.data.file)]

full_dat$lon <- full_dat$longitude
full_dat$lat <- full_dat$latitude

rcas <- st_read("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Analyses/Groundfish/Spatial GLMM/Data/RCAs_2019_BCAlbers/RCAs_BCAlbers.shp")
data2 <- full_dat[,c("lon","lat")] #%>% mutate(UTM_W=as.numeric(UTM_W),UTM_N=as.numeric(UTM_N)) #
sputm <- SpatialPoints(data2, proj4string=CRS("+proj=longlat"))
spgeo <- spTransform(sputm, CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
data_points <- as_tibble(full_dat) %>% mutate(lat=spgeo@coords[,2],lon=spgeo@coords[,1]) %>%
  st_as_sf(coords=c("lon","lat"),crs=3005)
central_coast <- st_read("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Analyses/Groundfish/Spatial GLMM/Data/MaPP_CC_SpatialPlanning_Area.shp")
## Subset the points data to be a single point for each location (for labelling Locations)
data_point_labels <- data_points %>% group_by(survey_id) %>% mutate("n_species"=length(unique(species))) %>% slice(1)

ncc_extent <- raster::extent(data_points)
ncc_extent@xmin <- ncc_extent@xmin - 3*buffer
ncc_extent@ymin <- ncc_extent@ymin - buffer 
ncc_extent@xmax <- ncc_extent@xmax + buffer 
ncc_extent@ymax <- ncc_extent@ymax + buffer 
plot_area_ncc <- ncc_extent %>%
  st_bbox() %>%                 # Turn into a square
  st_as_sfc(crs=st_crs(data_points))
st_crs(plot_area_ncc) <- st_crs(data_points)

bchres <- bc_bound_hres()
coast_line <- bchres %>% st_intersection(plot_area_ncc)

#Below is to get a dataset of lakes in plot box
lakes_in_plot_area <- bcdc_query_geodata("freshwater-atlas-lakes") %>%
  filter(INTERSECTS(plot_area_ncc)) %>%
  collect() %>%
  st_intersection(plot_area_ncc)

# Ocean colouring - this doesn't work well  because the resolution isn't the same
ocean_colour <- bc_neighbours(ask=FALSE) %>% 
  filter(type=="Ocean") %>% 
  st_intersection(plot_area_ncc)

#alaska <- bcmaps::bc_neighbours() %>% st_intersection(plot_area_ncc)
#alaska_full <- bcmaps::bc_neighbours()
pfma <- st_read("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Spatial data/PFMA/DFO_PFMA_SUBAREAS_SP/DFO_SBAREA_polygon.shp")
pfma <- pfma[pfma$MGMT_AREA %in% c(5,6,7,8,9,10,106,107,108,109,110),]
catch_area <- pfma %>% st_intersection(plot_area_ncc) %>% st_difference(bchres)
unique(catch_area$MGMT_AREA)

data_points$mgmt_area <- pfma$MGMT_AREA[unlist(st_intersects(data_points,pfma))]
data_points$mgmt_subarea <- pfma$LABEL[unlist(st_intersects(data_points,pfma))]

habitat <- raster::raster("Data/rfsubstrate_100m/NCC_substrate_20m.tif")
habitat <- raster::crop(habitat,raster::extent(central_coast))
habitat_qcs <- raster::raster("Data/rfsubstrate_100m/QCS_substrate_20m.tif")
habitat_qcs <- raster::crop(habitat_qcs,raster::extent(central_coast))
habitat_deep <- raster::raster("Data/rfsubstrate_100m/Ranger_RF_Substrate_100m.tif")
habitat_deep <- raster::crop(habitat_deep,raster::extent(central_coast))
ncc_raster <- terra::rast(habitat)
#plot(habitat_deep)
qcs_raster <- terra::rast(habitat_qcs)
deep_raster <- terra::rast(habitat_deep)
obj_raster <- merge(ncc_raster,qcs_raster)
#plot(obj_raster)
levels(obj_raster) <- data.frame(value=c(4,3,2,1),id=c("Mud","Sand","Mixed","Rock"),SUBSTRATE=c("Rock","Mixed","Sand","Mud")) # first column sets the intended order, 2nd column sets the original order, third column sets the categorical order
levels(deep_raster) <- data.frame(value=c(4,3,2,1),id=c("Mud","Sand","Mixed","Rock"),SUBSTRATE=c("Rock","Mixed","Sand","Mud")) # first column sets the intended order, 2nd column sets the original order, third column sets the categorical order

#plot(obj_raster)
#obj_raster <- terra::focal(obj_raster, 3, fun="modal", na.policy="only")
data_points$avg_habitat_dfo <- sapply(data_points$geometry, function(x) {mode(terra::extract(obj_raster, terra::vect(x))[[2]])})
data_points$avg_habitat_dfo[is.na(data_points$avg_habitat_dfo)] <- sapply(data_points$geometry[is.na(data_points$avg_habitat_dfo)], function(x) {mode(terra::extract(deep_raster, terra::vect(x))[[2]])})
table(data_points$avg_habitat_dfo,useNA="ifany")
#sapply(data_points$geometry[1:2], function(x) {mean(as.numeric(factor((terra::nearby(obj_raster, terra::vect(x))[[2]]),levels=c("Mud","Sand","Mixed","Rock"))))})
data_points$avg_habitat_dfo[which(is.na(data_points$avg_habitat_dfo))] <- sapply(which(is.na(data_points$avg_habitat_dfo)), FUN = function(x) {mode(terra::extract(obj_raster, terra::vect(data_points$geometry[!is.na(data_points$avg_habitat_dfo)][which.min(sf::st_distance(data_points$geometry[!is.na(data_points$avg_habitat_dfo)], data_points$geometry[x]))]))[[2]])})
table(data_points$avg_habitat_dfo,useNA="ifany")

saveRDS(data_points,"~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_rockfish_bodysize_spatial.rds")

bc_sub <- bchres %>% st_intersection(plot_area_ncc)
bc_fish <- bchres %>% st_intersection(catch_area)

#saveRDS(coast_line,"~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_boundary.rds")
alaska <- st_transform(USAboundaries::us_states(states = "Alaska", resolution = "high"),st_crs(plot_area_ncc))
alaska_full <- alaska %>% st_intersection(bchres)
alaska <- alaska  %>% st_intersection(plot_area_ncc)
#sf::st_drop_geometry()
sub_data_point_labels <- data_point_labels[sample(1:nrow(data_point_labels),8,replace=FALSE),]
## The fun part! we get to make a map
## Plot map -- Order of plotting MATTERS
cc_rca <- sf::st_filter(rcas,data_point_labels)
cc_rca <- cc_rca %>% st_intersection(catch_area)
cc_rca_labels <- cc_rca[!duplicated(st_drop_geometry(cc_rca[,1:2])),]

ncc_coast <- ggplot() +
  #geom_sf(data = catch_area, aes(fill=as.factor(MGMT_AREA))) +
  #scale_fill_manual(name="Management areas",values=c("#7570b3","#d95f02","#e7298a","#a6761d","#6a3d9a","#fdbf6f","#bd0026","#e6ab02","#cab2d6","grey25","orange")) +
  new_scale("fill") +
  geom_sf(data=alaska, fill='grey85') +
  geom_sf(data = lakes_in_plot_area, fill = "#1f78b4",colour=NA) +     #Plot Lakes
  geom_sf(data = coast_line, fill = "grey60") +                 #Plot coastline
  geom_sf(data=cc_rca,fill="yellow",alpha=0.5,col="black",lwd=0.2)+
  # geom_sf(data=data_point_labels,pch=21,fill="white") +
  geom_sf(data=sf::st_jitter(data_point_labels,1000),aes(fill=n_species),pch=21) +
  geom_sf(data=cc_rca,fill="yellow",alpha=0.4,col=NA,lwd=0.2)+
  #ggsflabel::geom_sf_label_repel(data=sub_data_point_labels,aes(label=survey_id),size=1.5, force = 1, nudge_x = -2, seed = 10,max.overlaps=20)+ #Add labels
  ggsflabel::geom_sf_label_repel(data=cc_rca_labels,aes(label=NAME),size=1.5, force = 1, nudge_x = -2, seed = 10,max.overlaps=20)+ #Add labels
  geom_sf(data = plot_area_ncc, alpha = 0,colour='black') +        #Plot area box
  coord_sf(expand = FALSE) +                                    #Expands box to axes
  xlab('Longitude') + ylab('Latitude') +                        #Axis titles
  annotation_scale(location = "bl", width_hint = 0.5) +         #Rose Compass
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.5, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering,
                         height = unit(1,"cm"), width = unit(1, "cm"))+
  scale_fill_viridis_c(name="Species richness",limits=c(0,12),option="F") +
  theme(panel.background = element_rect('lightblue1'), panel.grid.major = element_line('lightblue1'),legend.position="top",legend.box.just="center",legend.box="horizontal",legend.justification = "center",legend.key.size=unit(1, "lines"),legend.margin = margin(c(0,0,0,0),unit="lines"),legend.title=element_text(size=7),legend.text = element_text(size=5))

bc_neigh <- bc_neighbours(ask=FALSE)
bc_neigh <- bc_neigh[bc_neigh$name%in%c("Alaska"),]
bc_map <- ggplot() +
  geom_sf(data = bc_bound(ask=FALSE), fill = "grey10",colour=NA) +
  geom_sf(data=bc_neigh, fill='grey50',colour=NA) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(fill=adjustcolor("white",0.90), colour="black")) +
  theme(legend.justification = c(0, 1),legend.position = c(0, .95)) +
  theme(text = element_text(family = "Futura-Medium"),legend.title = element_text(family = "Futura-Bold", size = 10),legend.text = element_text(family = "Futura-Medium", size = 10))
bc <- bc_map + geom_rect(aes(xmin=ncc_extent@xmin,ymin=ncc_extent@ymin,xmax=ncc_extent@xmax,ymax=ncc_extent@ymax),fill = "tomato", colour = "grey70",alpha=0.5)

library(cowplot)
ncc_inset <- ggdraw(ncc_coast) +
  draw_plot({
    bc},
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = 0.12, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.735,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.15, 
    height = 0.15)

#Save the plot
ggsave('~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Analyses/Groundfish/Spatial GLMM/Figures/ncc_rockfish_bodysize_map.pdf',plot=ncc_inset,width = 6, height = 7,units='in',dpi=800)
ggsave('~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Analyses/Groundfish/Spatial GLMM/Figures/ncc_rockfish_bodysize_map.jpeg',plot=ncc_inset,width = 6, height = 7,units='in',dpi=800)

