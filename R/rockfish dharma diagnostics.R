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

fit_sdm <- readRDS(file="R/sdm_groundfish_re6.rds")
mesh <- readRDS(file="R/sdm_mesh.rds")
rockfish <- readRDS(file="R/sdm_model_data.rds")

s_spp_re <- simulate(fit_sdm, nsim = 1000,type="mle-mvn")
r_spp_re <- sdmTMB::dharma_residuals(s_spp_re, fit_sdm, plot = TRUE)
pred_fixed <- fit_sdm$family$linkinv(predict(fit_sdm)$est_non_rf)

r_nb <- DHARMa::createDHARMa(simulatedResponse = s_spp_re,observedResponse = (fit_sdm$data$counts),fittedPredictedResponse = pred_fixed)
DHARMa::plotResiduals(r_nb) 
DHARMa::plotQQunif(r_nb)
DHARMa::plotResiduals(r_nb, rank = TRUE, quantreg = FALSE, smoothScatter = TRUE)

dharma_df <- data.frame("scaled_residuals"=r_nb$scaledResiduals,"x"=rowMeans(r_nb$simulatedResponse),"predicted_values"=fitted(fit_sdm),"y"=r_nb$observedResponse/sd(r_nb$observedResponse),"resids"=r_nb$scaledResiduals,"sim_quantiles"=r_spp_re$expected,"obs_quantiles"=r_spp_re$observed)


ggplot(dharma_df, aes(x = sim_quantiles, y = obs_quantiles)) +
  geom_point() +
  # geom_smooth(method = "loess", se = FALSE) + # Add a LOESS smoother
  geom_abline(slope=1,intercept=0, linetype = "dashed", color = "red") + # Reference line for scaled residuals
  geom_vline(xintercept=0.5, linetype = "dashed", color = "grey20") + # Reference line for scaled residuals
  geom_hline(yintercept=0.5, linetype = "dashed", color = "grey20") + # Reference line for scaled residuals
  labs(x = "Theoretical quantiles",
       y = "Sample quantiles") +
  theme_classic()

ggsave("Figures/sdm_tmb_dharma_cpue.jpeg",dpi=600,width=7,height=6,units="in")

fit_sdm <- readRDS(file="R/sdm_groundfish_size_re.rds")
mesh <- readRDS(file="R/sdm_mesh_size.rds")
rockfish <- readRDS(file="R/sdm_model_data_size.rds")

s_spp_re <- simulate(fit_sdm, nsim = 1000,type="mle-mvn")
r_spp_re <- sdmTMB::dharma_residuals(s_spp_re, fit_sdm, plot = TRUE)
pred_fixed <- fit_sdm$family$linkinv(predict(fit_sdm)$est_non_rf)

r_nb <- DHARMa::createDHARMa(simulatedResponse = s_spp_re,observedResponse = (fit_sdm$data$counts),fittedPredictedResponse = pred_fixed)
DHARMa::plotResiduals(r_nb) 
DHARMa::plotQQunif(r_nb)
DHARMa::plotResiduals(r_nb, rank = TRUE, quantreg = FALSE, smoothScatter = TRUE)

dharma_df <- data.frame("scaled_residuals"=r_nb$scaledResiduals,"x"=rowMeans(r_nb$simulatedResponse),"predicted_values"=fitted(fit_sdm),"y"=r_nb$observedResponse/sd(r_nb$observedResponse),"resids"=r_nb$scaledResiduals,"sim_quantiles"=r_spp_re$expected,"obs_quantiles"=r_spp_re$observed)


ggplot(dharma_df, aes(x = sim_quantiles, y = obs_quantiles)) +
  geom_point() +
  # geom_smooth(method = "loess", se = FALSE) + # Add a LOESS smoother
  geom_abline(slope=1,intercept=0, linetype = "dashed", color = "red") + # Reference line for scaled residuals
  geom_vline(xintercept=0.5, linetype = "dashed", color = "grey20") + # Reference line for scaled residuals
  geom_hline(yintercept=0.5, linetype = "dashed", color = "grey20") + # Reference line for scaled residuals
  
  labs(x = "Theoretical quantiles",
       y = "Sample quantiles") +
  theme_classic()

ggsave("Figures/sdm_tmb_dharma_size.jpeg",dpi=600,width=7,height=6,units="in")
