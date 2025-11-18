library(sdmTMB)
library(emmeans)
library(sf)
library(dplyr)

quadroot <- function(x) {return((x)^(1/4))}
wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")

fit_sdm <- readRDS(file="R/sdm_groundfish_re6.rds")
fit_sdm <- sdmTMB:::reload_model(fit_sdm)
mesh <- readRDS(file="R/sdm_mesh.rds")
rockfish <- readRDS(file="R/sdm_model_data.rds")
coast_line <- readRDS("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_boundary.rds")
grid <-readRDS("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_rockfish_grid_4km.rds")
cell_size <- 4000
grid$protected <- ifelse(grid$rca_distkm==0,1,0)
grid$year_protected <- ifelse(grid$rca_name%in%c("Fish Egg Inlet","Goose Island","McMullin Group","Smith Sound"),2004,2005)
grid$upperoceansr <- relevel(grid$upperoceansr,ref="(11) Mainland Fjords")
pfma <- st_read("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Spatial data/PFMA/DFO_PFMA_SUBAREAS_SP/DFO_SBAREA_polygon.shp")

new_grid <- expand.grid("puid"="12411",
                        "species"="quillback",
                        "gear"="dive",
                        "effort"=120,
                        "chase_prob"=0,
                        "depth_seq"=50,
                        "cumulative_catch_before_rca" = 10,
                        "year" = c(2015),
                        "rca_distkm" = 5,
                        # "protection_dist_scale" = seq(-2,2,by=0.5),
                        "rmax" = c(0.096), # yelloweye, lingcod, brown
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

# Method 3: Contrast between slopes
# Compare temperature effect at shallow vs deep water

library(sdmTMB)
library(dplyr)

emmeans_sdmTMB <- function(model, specs, at = NULL, 
                           data = NULL, level = 0.95,new_df=TRUE,new_grid=NULL) {
  if (is.null(data)) data <- model$data
  data <- data[, !sapply(data, is.character)]
  # Build reference grid
  grid <- data %>% 
    mutate(across(where(is.numeric), ~ mean(., na.rm = TRUE))) %>%
    mutate(across(where(is.factor), ~ first(.)))
  max_l <- max(unlist(lapply(at,length)))
  grid <- grid[1,] %>% slice(rep(row_number(), max_l))
  
  if(new_df){grid <- new_grid}
  if (!is.null(at)) {
    for (nm in names(at)) {
      grid[[nm]] <- at[[nm]]
    }
  }
  
  reference_grid <- expand.grid(grid[unique(c(specs, names(at)))], KEEP.OUT.ATTRS = FALSE)
  reference_grid <- full_join(reference_grid,grid)
  # Predict
  preds <- predict(model, newdata = as.data.frame(reference_grid),offset=log(reference_grid$effort/(1-reference_grid$chase_prob)),se_fit=TRUE,type="link",re_form=~0,re_form_iid=~0)
  reference_grid$fit <- preds$est
  reference_grid$se_fit <- preds$est_se
  
  # Summarize by specs
  out <- reference_grid %>%
    group_by(across(names(at))) %>%   # no all_of() needed anymore
    summarise(
      emmean = mean(fit),
      SE = sqrt(mean(se_fit^2)),
      lower.CL = emmean - qnorm(1 - (1 - level)/2) * SE,
      upper.CL = emmean + qnorm(1 - (1 - level)/2) * SE,
      .groups = "drop"
    )
  
  contr <- out %>%
    rename(level = !!specs) %>%
    arrange(level) %>%
    mutate(dummy = 1) %>%
    inner_join(., ., by = "dummy", suffix = c(".1", ".2")) %>%
    filter(level.1 < level.2) %>%
    transmute(
      contrast = paste(level.2, "-", level.1),
      estimate = emmean.2 - emmean.1,
      SE = sqrt(SE.1^2 + SE.2^2),   # assumes independence of estimates
      lower.CL = estimate - qnorm(1 - (1 - level)/2) * SE,
      upper.CL = estimate + qnorm(1 - (1 - level)/2) * SE,
      z.ratio = estimate / SE,
      p.value = 2 * (1 - pnorm(abs(z.ratio)))
    )
  
  return(list(emmeans = out, contrasts = contr))
}



emmeans_sdmTMB(model=fit_sdm,specs = "cumul_catch_before_rca_scaled", at = list(protection_dist_scale = seq(-2, 2, by = 0.5)),new_grid=new_grid,new_df=TRUE)
