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
library(cowplot)

data_points <- readRDS("~/Dropbox (CCIRA)/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_rockfish_spatial.rds")
pfma <- st_read("~/Dropbox (CCIRA)/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Spatial data/PFMA/DFO_PFMA_SUBAREAS_SP/DFO_SBAREA_polygon.shp")
rcas <- st_read("~/Dropbox (CCIRA)/Ecological Research - test CCIRA business acct/Analyses/Groundfish/Spatial GLMM/Data/RCAs_2019_BCAlbers/RCAs_BCAlbers.shp")
central_coast <- st_read("~/Dropbox (CCIRA)/Ecological Research - test CCIRA business acct/Analyses/Groundfish/Spatial GLMM/Data/MaPP_CC_SpatialPlanning_Area.shp")
grid <-readRDS("~/Dropbox (CCIRA)/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_rockfish_grid_4km.rds")
cell_size <- 4000
grid$protected <- ifelse(grid$rca_distkm==0,1,0)
grid$year_protected <- ifelse(grid$rca_name%in%c("Fish Egg Inlet","Goose Island","McMullin Group","Smith Sound"),2004,2005)
grid$upperoceansr <- relevel(grid$upperoceansr,ref="(11) Mainland Fjords")

buffer <- 25000

## Subset the points data to be a single point for each location (for labelling Locations)
data_point_labels <- data_points %>% group_by(survey_id) %>% mutate("n_species"=length(unique(species))) %>% slice(1)

pfma <- pfma[pfma$MGMT_AREA %in% c(5,6,7,8,9,10,106,107,108,109,110),]

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

catch_area <- pfma %>% st_intersection(plot_area_ncc) %>% st_difference(bchres)
catch_union <- catch_area %>% st_union()
bc_sub <- bchres %>% st_intersection(plot_area_ncc)
bc_fish <- bchres %>% st_intersection(catch_area)


grid_cc <- grid %>% st_intersection(catch_area)

grid_spacing <- 4000  # size of squares, in units of the CRS (i.e. meters for 5514)

polygony <- st_make_grid(plot_area_ncc, square = T, cellsize = c(grid_spacing, grid_spacing)) %>% # the grid, covering bounding box
  st_sf() # not really required, but makes the grid nicer to work with later

grid_cc <- polygony %>% sf::st_crop(plot_area_ncc)
grid_cc <- st_join(grid_cc,grid,join=st_contains)

grid_cc$sample_size <- lengths(st_intersects(grid_cc,data_point_labels))
grid_cc$sample_size[grid_cc$sample_size==0] <- NA
grid_cc <- grid_cc[!is.na(grid_cc$sample_size),]
alaska <- st_transform(USAboundaries::us_states(states = "Alaska", resolution = "high"),st_crs(plot_area_ncc))
alaska_full <- alaska %>% st_intersection(bchres)
alaska <- alaska  %>% st_intersection(plot_area_ncc)


sub_data_point_labels <- data_point_labels[sample(1:nrow(data_point_labels),8,replace=FALSE),]
## The fun part! we get to make a map
## Plot map -- Order of plotting MATTERS
cc_rca <- sf::st_filter(rcas,data_point_labels)
cc_rca <- cc_rca %>% st_intersection(catch_area)
cc_rca_labels <- cc_rca[!duplicated(st_drop_geometry(cc_rca[,1:2])),]
cc_rca_labels$NAME <- c("Kitasu Bay","West Aristazabal I.","Goose Island","McMullin Group","West Calvert","Fish Egg Inlet","Smith Sound")

ncc_coast <- ggplot() +
  geom_sf(data=alaska, fill='grey85') +
  geom_sf(data=grid_cc,aes(fill=sample_size),colour="black") +
  geom_sf(data = coast_line, fill = "grey60") +                 #Plot coastline
  geom_sf(data=cc_rca,fill="yellow",alpha=0.6,col="black",lwd=0.3)+
  # geom_sf(data=cc_rca,fill="yellow",alpha=0.4,col=NA,lwd=0.2)+
  # geom_sf_text(data=sf::st_jitter(cc_rca_labels,amount=2000),aes(label=NAME),size=3,colour="white", nudge_x = 10000,nudge_y = 1000,check_overlap = FALSE,hjust=0,fontface = "bold",bg.color = "black",)+ #Add labels
  geom_text_repel(data=sf::st_jitter(cc_rca_labels,amount=0),aes(label=NAME,, geometry = geometry),nudge_x = -8000,nudge_y = -2000,size=3.2,colour="white",hjust=1,fontface = "plain",bg.color = "black",stat = "sf_coordinates",segment.colour="grey15",min.segment.length=0.25)+ #Add labels
  
  geom_sf(data = plot_area_ncc, alpha = 0,colour='black') +        #Plot area box
  coord_sf(expand = FALSE) +                                    #Expands box to axes
  xlab('Longitude') + ylab('Latitude') +                        #Axis titles
  scale_fill_viridis_c(name="No. of surveys",trans="log10",option="F") +
  annotation_scale(location = "bl", width_hint = 0.5) +         #Rose Compass
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.5, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering,
                         height = unit(1,"cm"), width = unit(1, "cm"))+
  theme(panel.background = element_rect('lightblue1'), panel.grid.major = element_line('lightblue1'),legend.position="top",legend.box.just="center",legend.box="horizontal",legend.justification = "center",legend.key.size=unit(1, "lines"),legend.margin = margin(c(0,0,0,0),unit="lines"),legend.title=element_text(size=12),legend.text = element_text(size=5)) +
  theme(legend.location = "top")

bc_neigh <- bc_neighbours(ask=FALSE)
bc_neigh <- bc_neigh[bc_neigh$name%in%c("Alaska"),]
bc_map <- ggplot() +
  geom_sf(data = bc_bound(ask=FALSE), fill = "grey10",colour=NA) +
  geom_sf(data=bc_neigh, fill='grey50',colour=NA) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(fill=adjustcolor("white",0.90), colour="black")) +
  theme(legend.justification = c(0, 1),legend.position = c(0, .95)) +
  theme(text = element_text(family = "Futura-Medium"),legend.title = element_text(family = "Futura-Bold", size = 12),legend.text = element_text(family = "Futura-Medium", size = 12))
bc <- bc_map + geom_rect(aes(xmin=ncc_extent@xmin,ymin=ncc_extent@ymin,xmax=ncc_extent@xmax,ymax=ncc_extent@ymax),fill = "tomato", colour = "grey70",alpha=0.5)

ncc_inset <- ggdraw(ncc_coast) +
  draw_plot({
    bc},
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = 0.118, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.698,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.15, 
    height = 0.15)
ncc_inset

#Save the plot
ggsave('~/Dropbox (CCIRA)/Ecological Research - test CCIRA business acct/Analyses/Groundfish/Spatial GLMM/Figures/ncc_rockfish_map_update.jpeg',plot=ncc_inset,width = 6, height = 7,units='in',dpi=800)
knitr::plot_crop('~/Dropbox (CCIRA)/Ecological Research - test CCIRA business acct/Analyses/Groundfish/Spatial GLMM/Figures/ncc_rockfish_map_update.jpeg')
