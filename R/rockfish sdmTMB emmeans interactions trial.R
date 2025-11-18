# Example: Getting coefficient values for continuous interactions using emmeans
# This demonstrates how to extract coefficients at various values of
# continuous predictors in an interaction
library(sdmTMB)
library(sf)
library(dplyr)
library(emmeans)

quadroot <- function(x) {return((x)^(1/4))}
wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")

fit_sdm <- readRDS(file="R/sdm_groundfish_re6.rds")
fit_sdm <- sdmTMB:::reload_model(fit_sdm)
mesh <- readRDS(file="R/sdm_mesh.rds")
rockfish <- readRDS(file="R/sdm_model_data.rds")
grid <-readRDS("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_rockfish_grid_4km.rds")

cell_size <- 4000
grid$protected <- ifelse(grid$rca_distkm==0,1,0)
grid$year_protected <- ifelse(grid$rca_name%in%c("Fish Egg Inlet","Goose Island","McMullin Group","Smith Sound"),2004,2005)
grid$upperoceansr <- relevel(grid$upperoceansr,ref="(11) Mainland Fjords")
new_grid <- expand.grid("puid"="12411",
                        "species"="quillback",
                        "gear"="dive",
                        "effort"=120,
                        "chase_prob"=0,
                        "depth_seq"=50,
                        "cumulative_catch_before_rca" = seq(0,100,length=20),
                        "year" = c(2006,2013,2023),
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
new_grid <- st_as_sf(new_grid,coords=c("X","Y"),crs=sf::st_crs(rockfish))

# Method 3: Contrast between slopes
# Compare temperature effect at shallow vs deep water
ref_grid(fit_sdm)

fit_sdm$formula

rca_dist_seq <- (-1*(c(0,1,5,15,25)-mean(rockfish$rca_distkm))/sd(rockfish$rca_distkm))
year_protected_seq <- ((c(2005,2010,2015,2020,2025)-2004)-mean(rockfish$year_protected))/sd(rockfish$year_protected)
rmax_seq <- (c(0.043,0.096,0.142)-mean(rockfish$rmax))/sd(rockfish$rmax)
cumul_catch_before_rca_seq <- (quadroot(c(0,10,100,1000,10000)) - mean(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%"quillback"])))/sd(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%"quillback"]))

rca_slopes_over_time <- emtrends(fit_sdm, ~ year_protected_scaled, var = "protection_dist_scale",
                                 at = list(year_protected_scaled = year_protected_seq,cumul_catch_before_rca_scaled=cumul_catch_before_rca_seq,rmax_scale =rmax_seq))
print(list("RCA distance effect along recent (+2 std. dev) and early years (-2 std. devs):",
           rca_slopes_over_time))

# Compare temperature effect at shallow vs deep water
temp_contrast <- pairs(emtrends(fit_sdm, ~ year_protected_scaled, var = "protection_dist_scale",
                                at = list(year_protected_scaled = rev(range(year_protected_seq)),
                                          cumul_catch_before_rca_scaled=cumul_catch_before_rca_seq,
                                          rmax_scale =rmax_seq)))
print("Contrast in RCA distance effect between recent (+2 std. dev) and early years (-2 std. devs):")
print(temp_contrast)


# Method 4: Plot the interaction using emmip
emmip(fit_sdm, year_protected_scaled ~ protection_dist_scale, var = "protection_dist_scale",
      at = list(year_protected_scaled = year_protected_seq,protection_dist_scale=rca_dist_seq,cumul_catch_before_rca_scaled=cumul_catch_before_rca_seq,rmax_scale =rmax_seq),CIs = TRUE)




rca_dist_seq <- (-1*(c(0,1,5,15,25)-mean(rockfish$rca_distkm))/sd(rockfish$rca_distkm))
year_protected_seq <- ((c(2005,2010,2015,2020,2025)-2004)-mean(rockfish$year_protected))/sd(rockfish$year_protected)
rmax_seq <- (c(0.043,0.096,0.142)-mean(rockfish$rmax))/sd(rockfish$rmax)
cumul_catch_before_rca_seq <- (quadroot(c(0,10,100,1000,10000)) - mean(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%"quillback"])))/sd(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%"quillback"]))

rca_slopes_over_time <- emtrends(fit_sdm, ~ year_protected_scaled | rmax_scale, var = "protection_dist_scale",
                                 at = list(year_protected_scaled = year_protected_seq,rmax_scale=rmax_seq))
df <- data.frame("RCA_age"=rep(c(2005,2010,2015,2020,2025)-2004,times=3),"rmax"=rep(c(0.043,0.096,0.142),each=5),tidy(rca_slopes_over_time))
df$ui <- df$protection_dist_scale.trend+1.96*df$std.error
df$li <- df$protection_dist_scale.trend-1.96*df$std.error

set.seed(2023)
df <- df %>%
  group_by(rmax) %>%
  mutate(jit=jitter(x=0,factor=50)) %>%
  ungroup()
df$RCA_jitter <- df$RCA_age + (df$jit-mean(df$jit))
print(list("RCA distance effect along recent (+2 std. dev) and early years (-2 std. devs):",
           rca_slopes_over_time))
p2 <- ggplot(data=df,aes(x=protection_dist_scale.trend,y=RCA_jitter,fill=as.factor(rmax),colour=as.factor(rmax))) +
  geom_point()+
  geom_segment(aes(x=li,xend=ui,y=RCA_jitter)) +
  ylab("RCA age") + xlab("Effect of spatial protection (distance from nearest RCA)") +
  theme_minimal() + theme(legend.position="top") +
  scale_y_continuous(trans = "reverse") + scale_fill_viridis_d(name=expression(r[max]),option="D") + scale_colour_viridis_d(name=expression(r[max]),option="D")+
  theme(text = element_text(size=11),legend.text = element_text(size=10),strip.text=element_text(size=9))
p2
mega_plot <- ggpubr::ggarrange(p1,p2,ncol=2,widths=c(3,2),heights=c(1,1))

ggsave(filename = "Figures/effect of RCA interactions emmeans.jpeg",plot=p2,dpi=600,units="in",height=5,width=5.5)
