library(dbplyr)
library(bcdata)
library(bcmaps)
library(sf)
library(sp)
library(ggplot2)
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
library(glmmTMB)
library(lme4)
# remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE)
# remotes::install_github("pbs-assess/sdmTMBextra", dependencies = TRUE)
# remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,ref="dispformula2")
# remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,ref="refactor-re")
# devtools::install_github("pbs-assess/sdmTMB", dependencies = TRUE,ref="re-slope")
# pak::pak("pbs-assess/sdmTMB", dependencies = TRUE)
# pak::pak("pbs-assess/sdmTMBextra", dependencies = TRUE)
trial <- FALSE
rockfish <- readRDS("Data/ccira_sdmTMB_data_catch.rds")
rockfish <- as.data.frame(rockfish)
rockfish$species <- factor(rockfish$species,levels=sort(unique(rockfish$species)))
rockfish$rmax_cat <- ifelse(rockfish$rmax_scale>0.5,"fast",ifelse(rockfish$rmax_scale < -0.5,"slow","average"))
rockfish$rca_rmax <- factor(paste(rockfish$rca_name,rockfish$rmax_cat,sep="-"))
rockfish$rca_species <- factor(paste(rockfish$rca_name,rockfish$species,sep="-"))
rockfish$log_depth <- log(rockfish$depth)
rca_group_names <- unique(rockfish$rca_species)
rca_predictions <- data.frame("rca_group"=rca_group_names,"species"=stringr::str_split_i(rca_group_names,"-",2),"RCA"=stringr::str_split_i(rca_group_names,"-",1))


ggplot(data=rockfish) +
  geom_boxplot(aes(x=avg_habitat_dfo,y=complexity)) +
  facet_wrap(~gear,scales="free_y")
ggplot(data=rockfish[rockfish$gear=="dive",]) +
  geom_boxplot(aes(x=ifelse(rca_distkm>0,"outside","inside"),y=complexity)) +
  #geom_smooth(aes(x=ifelse(rca_distkm>0,"outside","inside"),y=complexity))+
  facet_wrap(~rca_name)
ggplot(data=rockfish[rockfish$species%in%c("quillback"),]) +
  geom_boxplot(aes(x=as.factor(month),y=log(counts/effort))) +
  #geom_smooth(aes(x=month,y=cpue))+
  facet_wrap(~gear,scales="free_y")

bc_ocean <- bcmaps::bc_neighbours(ask=FALSE) %>% filter(type=="Ocean")
bc_land <- bcmaps::bc_neighbours(ask=FALSE) %>% filter(type=="Province")
bc_coast <- bcmaps::bc_bound()

purrr::map(rockfish, ~sum(is.na(.)))
sf_rockfish <- st_as_sf(rockfish)

catch_area <- sf_rockfish %>% st_buffer(dist = 200) %>% st_bbox() %>% st_as_sfc()
bc_fish <- bc_coast %>% st_intersection(catch_area)

coast_line <- readRDS("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_boundary.rds")

#mesh <- make_mesh(rockfish, xy_cols = c("X_km", "Y_km"), type="kmeans",n_knots = 300)
mesh <- make_mesh(rockfish, xy_cols = c("X_km", "Y_km"), cutoff=3)

barrier_mesh <- sdmTMBextra::add_barrier_mesh(mesh, coast_line, range_fraction = 0.1,proj_scaling = 1000, plot = TRUE)
ggplot(coast_line) +
  geom_sf() +
  geom_sf(data = barrier_mesh$mesh_sf[barrier_mesh$normal_triangles, ], size = 1, colour = "blue") +
  geom_sf(data = barrier_mesh$mesh_sf[barrier_mesh$barrier_triangles, ], size = 1, colour = "green")
# maybe drop 'habitat' as variable; should we add 'protection_dist_scale:rmax_scale'

# trial lme4 v. sdmTMB
if(trial)
{
  fit_spp_re_nonsf <- sdmTMB(
    counts ~ year_protected_scaled + gear + complexity_scaled + protection_dist_scale + rmax_scale + cumul_catch_before_rca_scaled + year_protected_scaled:rmax_scale + protection_dist_scale:cumul_catch_before_rca_scaled + protection_dist_scale:year_protected_scaled + protection_dist_scale:rmax_scale + (1 + protection_dist_scale | rca_name) + (1 + protection_dist_scale | species),
    offset = log(rockfish$effort/(1-rockfish$chase_prob)),
    data = rockfish,
    mesh = barrier_mesh,
    family = nbinom2(link="log"),
    # control = sdmTMBcontrol(profile=TRUE),
    spatial="off",
    reml = FALSE,
    silent=FALSE)
  tidy_fit <- tidy(fit_spp_re_nonsf,conf.int=TRUE,effects=c("ran_vals"))
    
  fit_glmm <- glmmTMB::glmmTMB(counts ~ year_protected_scaled + gear + complexity_scaled + protection_dist_scale + rmax_scale + cumul_catch_before_rca_scaled + year_protected_scaled:rmax_scale + protection_dist_scale:cumul_catch_before_rca_scaled + protection_dist_scale:year_protected_scaled + protection_dist_scale:rmax_scale + (1 + protection_dist_scale | rca_name) + (1 + protection_dist_scale | species),
  offset = log(rockfish$effort/(1-rockfish$chase_prob)),
  data = rockfish,  family = nbinom2(link="log"))

  tinytest::expect_equal(mean(diag(cor(ranef(fit_spp_re_nonsf)[[1]]$rca_species, ranef(fit_glmm)[[1]]$rca_species))), 1)
  
  tinytest::expect_equal(ranef(fit_spp_re_nonsf)[[1]]$rca_species[,1], ranef(fit_glmm)[[1]]$rca_species[,1],1)
  }

fit_spp_re_6 <- sdmTMB(
  counts ~ year_protected_scaled + gear + s(log_depth,by=depth_cat) + complexity_scaled + s(biweekly_scale) + protection_dist_scale + rmax_scale + cumul_catch_before_rca_scaled + year_protected_scaled:rmax_scale + protection_dist_scale:cumul_catch_before_rca_scaled + protection_dist_scale:year_protected_scaled + protection_dist_scale:rmax_scale + (1 + protection_dist_scale | rca_name) + (1 + protection_dist_scale | species),
  offset = log(rockfish$effort/(1-rockfish$chase_prob)),
  # experimental = list(slope_group = "rca_species", slope_covariate = "protection_dist_scale"),
  #dispformula = ~gear,
  #spatial_varying = ~ 0 + protection_dist_scale,
  data = rockfish,
  mesh = barrier_mesh,
  time = "species",
  #extra_time = c(2008, 2009, 2010, 2011, 2012),
  #family = delta_truncated_nbinom2(),
  family = nbinom2(link="log"),
  spatiotemporal = "iid",
  # control = sdmTMBcontrol(profile=TRUE),
  spatial="on",
  reml = FALSE,
  silent=FALSE)
edf6 <- sdmTMB::cAIC(fit_spp_re_6, what = "EDF")
summary(fit_spp_re_6)
sanity(fit_spp_re_6)
round(summary(fit_spp_re_6$sd_report, select = c("fixed"), p.value = TRUE),4)
tidy_fit6 <- tidy(fit_spp_re_6,conf.int=TRUE,effects=c("fixed"))
tidy_fit6$group_name <- "global"
tidy_fit6$par_name <- tidy_fit6$term
tidy_fit6$level_ids <- "global"

ran_pars <- tidy(fit_spp_re_6,conf.int=TRUE,effects=c("ran_pars"))
ran_pars <- ran_pars[,!colnames(ran_pars)=="model"]
ran_pars$group_name <- "global"
ran_pars$par_name <- ran_pars$term
ran_pars$level_ids <- "global"
tidy_fit6 <- rbind(tidy_fit6,ran_pars)

re_tidy <-tidy(fit_spp_re_6,conf.int=TRUE,effects=c("ran_vals"))

re_tidy$par_name <- paste(re_tidy$term,re_tidy$level_ids,sep="-")
colnames(re_tidy)[colnames(re_tidy)=="conf.hi"] <- "conf.high"
re_tidy <- re_tidy[,!colnames(re_tidy)%in%c("model")]
tidy_fit6 <- rbind(tidy_fit6,re_tidy)
print(tidy_fit6,n=199)
print(tidy(fit_spp_re_6,conf.int=TRUE,effects=c("ran_vals")),n=70)

tidy(fit_spp_re_6,conf.int=TRUE,effects=c("ran_vcov"))

n_sims <- 10000
# tmb_sd <- fit_spp_re_6$sd_report
samps <- sdmTMB:::rmvnorm_prec(fit_spp_re_6$tmb_obj$env$last.par.best, fit_spp_re_6$sd_report, n_sims)
pars <- c(fit_spp_re_6$sd_report$par.fixed, fit_spp_re_6$sd_report$par.random)
pn <- names(pars)
b_j <- which(row.names(samps) == "b_j")[which(tidy_fit6$par_name=="protection_dist_scale")]
re_b_pars <- which(row.names(samps) == "re_b_pars") # intercepts followed by random slopes here...
tidy_re_b_pars <- which(tidy_fit6$group_name!="global")
re <- samps[re_b_pars, ]
plot(rowMeans(re),tidy_fit6$estimate[tidy_re_b_pars])
abline(a=0,b=1)
rca_names <- sort(unique(rockfish$rca_name))
species <- sort(unique(rockfish$species))
fixed_rca_intercepts <- FALSE
if(fixed_rca_intercepts)
{
  rca_intercepts_idx <- 0#seq(1,length(rca_names),by=1)
  rca_slopes_idx <- seq(1,length(rca_names),by=1)
  rca_intercepts <- re[rca_intercepts_idx,]
  rca_slopes <- re[rca_slopes_idx,]
  
  species_intercepts_idx <- seq(length(rca_names)+1,nrow(re),by=2)
  species_slopes_idx <- seq(length(rca_names)+2,nrow(re),by=2)
  species_intercepts <- re[species_intercepts_idx,]
  species_slopes <- re[species_slopes_idx,] 
}else
{
  rca_intercepts_idx <- seq(1,2*length(rca_names),by=2)
  rca_slopes_idx <- seq(2,2*length(rca_names),by=2)
  rca_intercepts <- re[rca_intercepts_idx,]
  rca_slopes <- re[rca_slopes_idx,]
  
  species_intercepts_idx <- seq(2*length(rca_names)+1,nrow(re),by=2)
  species_slopes_idx <- seq(2*length(rca_names)+2,nrow(re),by=2)
  species_intercepts <- re[species_intercepts_idx,]
  species_slopes <- re[species_slopes_idx,] 
}

b <- samps[b_j,,drop=FALSE]
rca_names_idx <- as.numeric(rep(rca_names,each=length(species)))
species_idx <- as.numeric(rep(species,time=length(rca_names)))

plot(rca_intercepts,rca_slopes)
plot(species_intercepts,species_slopes)
out <- matrix(nrow = nrow(rca_slopes)*nrow(species_slopes), ncol = ncol(rca_slopes))
for (i in 1:ncol(rca_slopes)) {
  for(j in 1:nrow(out)) {
    out[j,i] <- b[i] + rca_slopes[rca_names_idx[j],i] + species_slopes[species_idx[j],i]
  }
}
rca_names_df <- rca_names[rca_names_idx]
species_df <- species[species_idx]

rca_species <- data.frame(rca_names_df,species_df)
rca_species <- apply( rca_species , 1 , paste , collapse = "-" )

out_subset <- out[match(rca_predictions$rca_group,rca_species),]
rca_predictions$mean <- apply(out_subset,1,mean)
rca_predictions$ui <- apply(out_subset,1,quantile,probs=0.975)
rca_predictions$li <- apply(out_subset,1,quantile,probs=0.025)
rca_predictions$p_above <- apply(out_subset,1,function(x){sum(x>0)/length(x)})
rca_predictions$p_below <- apply(out_subset,1,function(x){sum(x<0)/length(x)})

rca_predictions$rca_effect <- ifelse(rca_predictions$p_above>=0.95,"rca benefits",ifelse(rca_predictions$p_below>0.95,"rca detriments","not significant"))

# https://waves-vagues.dfo-mpo.gc.ca/Library/40876159.pdf
# Mastrandrea et al 2010
rca_predictions$rca_benefits_ipcc <- ifelse(rca_predictions$p_above>=0.99,"very likely",
                                            ifelse(rca_predictions$p_above>=0.90,"very likely",
                                                   ifelse(rca_predictions$p_above>=0.66,"likely",
                                                          ifelse(rca_predictions$p_above>=0.33,"as likely as not",
                                                                 ifelse(rca_predictions$p_above>=0.10,"unlikely",
                                                                        ifelse(rca_predictions$p_above>=0.01,"very unlikely","very unlikely"))))))

rca_predictions$rca_benefits_ipcc <- factor(rca_predictions$rca_benefits_ipcc,levels=c("very likely","likely","as likely as not","unlikely","very unlikely"))

table(rca_predictions$rca_benefits_ipcc,rca_predictions$species)
table(rca_predictions$rca_benefits_ipcc,rca_predictions$RCA)
table(rca_predictions$rca_benefits_ipcc)

table(rca_predictions$rca_effect,rca_predictions$species)
table(rca_predictions$rca_effect,rca_predictions$RCA)
table(rca_predictions$rca_effect)
rca_predictions$sample_size <- sapply(unique(rca_predictions$rca_group),function(x){sum(paste(rockfish$rca_name,rockfish$species,sep="-") %in% x)})
rca_predictions$totals <- sapply(unique(rca_predictions$rca_group),function(x){length(rockfish$counts[paste(rockfish$rca_name,rockfish$species,sep="-") %in% x])})
rca_predictions[rca_predictions$species=="yelloweye",]
rca_predictions[rca_predictions$species=="black",]
rca_predictions[rca_predictions$species=="quillback",]
rca_predictions[rca_predictions$species=="copper",]
rca_predictions[rca_predictions$species=="brown",]
rca_predictions[rca_predictions$species=="lingcod",]


rca_predictions[order(rca_predictions$mean,decreasing =TRUE),]

s_spp_re6 <- simulate(fit_spp_re_6, nsim = 1000,type="mle-mvn")
r_spp_re6 <- sdmTMB::dharma_residuals(s_spp_re6, fit_spp_re_6, plot = TRUE)
pred_fixed6 <- fit_spp_re_6$family$linkinv(predict(fit_spp_re_6)$est_non_rf)

r_nb6 <- DHARMa::createDHARMa(simulatedResponse = s_spp_re6,observedResponse = (fit_spp_re_6$data$counts),fittedPredictedResponse = pred_fixed6)
DHARMa::plotResiduals(r_nb6) 
DHARMa::plotQQunif(r_nb6)
DHARMa::plotResiduals(r_nb6, rank = TRUE, quantreg = FALSE, smoothScatter = TRUE)
jpeg(filename = "Figures/count_cpue_residuals.jpeg",width=8,height=6,units="in",res=600)
plot(r_nb6)
dev.off()
dharma_df <- data.frame("y"=r_nb6$observedResponse/sd(r_nb6$observedResponse),"x"=rowMeans(r_nb6$simulatedResponse),"resids"=r_nb6$scaledResiduals)

DHARMa::plotResiduals(r_nb6) 
DHARMa::plotResiduals(r_nb6,form = rockfish$log_depth)
DHARMa::plotQQunif(r_nb6,form = rockfish$log_depth)

saveRDS(fit_spp_re_6,file="R/sdm_groundfish_re6.rds")

saveRDS(rockfish,file="R/sdm_model_data.rds")

saveRDS(barrier_mesh,file="R/sdm_mesh.rds")
saveRDS(edf6,file="R/sdm_edf6.rds")
saveRDS(s_spp_re6,file="R/sdm_dharma_re6.rds")

saveRDS(rca_predictions,file="R/sdm_rca_predictions.rds")

# edf6 <- readRDS(file="~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Analyses/Groundfish/Spatial GLMM/R/sdm_edf6.rds")

