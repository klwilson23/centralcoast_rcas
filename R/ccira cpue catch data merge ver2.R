library(tidyr)
library(dplyr)
library(sf)
# merge back with Sean Anderson's data pull
rockfish <- readRDS("Data/ccira_sdmTMB_data.rds")
rockfish$ccira_row_id <- 1:nrow(rockfish)
rockfish$date <- lubridate::as_date(rockfish$date)
rockfish <- sf::st_as_sf(rockfish)

rockfish_harvest <- readRDS("Data/cumulative-catch-all-1.0km.rds")
rockfish_harvest$date <- lubridate::as_date(rockfish_harvest$date)
rockfish_harvest <- st_transform(rockfish_harvest, st_crs(rockfish))

rockfish_full <- dplyr::left_join(rockfish,rockfish_harvest[,c("ccira_row_id","before_rca","cumulative_catch")] %>% st_drop_geometry(),by="ccira_row_id")

rockfish_wide <- tidyr::pivot_wider(rockfish_full,names_from='before_rca',names_prefix = "before_rca",values_from="cumulative_catch")
rockfish_wide <- rockfish_wide[,!colnames(rockfish_wide)=="before_rcaNA"]
colnames(rockfish_wide)[colnames(rockfish_wide)%in%c("before_rcaTRUE","before_rcaFALSE")] <- c("cumulative_catch_after_rca","cumulative_catch_before_rca")

aggregate(cbind(cumulative_catch_before_rca,cumulative_catch_after_rca,cumul_bin_harvest,harvest)~species,data=rockfish_wide,sum,na.rm=TRUE)

rockfish_wide$cumulative_catch_before_rca[is.na(rockfish_wide$cumulative_catch_before_rca)] <- 0
rockfish_wide$cumulative_catch_after_rca[is.na(rockfish_wide$cumulative_catch_after_rca)] <- 0
rockfish_wide$cumulative_catch <- rockfish_wide$cumulative_catch_before_rca + rockfish_wide$cumulative_catch_after_rca

min_kg <- min(rockfish_wide$cumulative_catch[rockfish_wide$cumulative_catch!=0])
min_before_kg <- min(rockfish_wide$cumulative_catch_before_rca[rockfish_wide$cumulative_catch_before_rca!=0])
min_after_kg <- min(rockfish_wide$cumulative_catch_after_rca[rockfish_wide$cumulative_catch_after_rca!=0])

rockfish_wide <- rockfish_wide %>%
  group_by(species) %>%
  mutate("cumul_catch_scaled"=as.numeric(scale((cumulative_catch)^(1/4),center=TRUE,scale=TRUE)),
         "cumul_catch_before_rca_scaled"=as.numeric(scale((cumulative_catch_before_rca)^(1/4),center=TRUE,scale=TRUE)),
         "cumul_catch_after_rca_scaled"=as.numeric(scale((cumulative_catch_after_rca)^(1/4),center=TRUE,scale=TRUE))) %>%
  ungroup()
rockfish_wide$cumul_catch_scaled[is.na(rockfish_wide$cumul_catch_scaled)] <- 0
rockfish_wide$cumul_catch_before_rca_scaled[is.na(rockfish_wide$cumul_catch_before_rca_scaled)] <- 0
rockfish_wide$cumul_catch_after_rca_scaled[is.na(rockfish_wide$cumul_catch_after_rca_scaled)] <- 0

plot(log(rockfish_wide$cumulative_catch_after_rca),log(rockfish_wide$cumul_bin_harvest_post_rca))
plot(rockfish_wide$cumul_catch_before_rca_scaled,rockfish_wide$cumul_catch_after_rca_scaled)

ggplot(data=rockfish_wide,aes(x=log(cumulative_catch_before_rca+0.5*min_before_kg),y=log(counts/effort))) +
  geom_point()+
  geom_smooth()+
  facet_wrap(~species,scales="free")

ggplot(data=rockfish_wide,aes(x=(cumulative_catch_before_rca)^(1/4),y=log(counts/effort))) +
  geom_point()+
  geom_smooth()+
  facet_wrap(~species,scales="free")

ggplot(data=rockfish_wide,aes(x=(cumulative_catch_before_rca),y=log(counts/effort))) +
  geom_point()+
  geom_smooth()+
  facet_wrap(~species,scales="free")

ggplot(data=rockfish_wide,aes(x=sqrt(cumulative_catch_before_rca),y=log(counts/effort))) +
  geom_point()+
  geom_smooth()+
  facet_wrap(~species,scales="free")

ggplot(data=rockfish_wide,aes(x=sqrt(cumulative_catch_before_rca),y=(cumulative_catch_before_rca)^(1/4))) +
  geom_point()+
  geom_smooth()+
  facet_wrap(~species,scales="free")

saveRDS(rockfish_wide,"Data/ccira_sdmTMB_data_catch.rds")
