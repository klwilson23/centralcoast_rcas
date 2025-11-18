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
fit_sdm <- readRDS(file="R/sdm_groundfish_size_re.rds")
mesh <- readRDS(file="R/sdm_mesh.rds")
pfma <- st_read("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Spatial data/PFMA/DFO_PFMA_SUBAREAS_SP/DFO_SBAREA_polygon.shp")
rockfish <- readRDS(file="R/sdm_model_data_size.rds")
grid <-readRDS("~/CCIRA Dropbox/Kyle Wilson/Ecological Research - test CCIRA business acct/Master Data/Groundfish, sponge, coral/Data Merger 2022/ccira_rockfish_grid_4km.rds")

cell_size <- 4000
grid$protected <- ifelse(grid$rca_distkm==0,1,0)
grid$year_protected <- ifelse(grid$rca_name%in%c("Fish Egg Inlet","Goose Island","McMullin Group","Smith Sound"),2004,2005)
grid$upperoceansr <- relevel(grid$upperoceansr,ref="(11) Mainland Fjords")
range(rockfish$cumulative_catch_before_rca)
quantile(quadroot(rockfish$cumulative_catch_before_rca),probs=c(0.25,0.75,0.95))^4
range(rockfish$year_protected_scaled)
range(rockfish$protection_dist_scale)
range(rockfish$rmax_scale)
min_kg <- min(rockfish$cumulative_catch[rockfish$cumulative_catch!=0])
min_before_kg <- min(rockfish$cumulative_catch_before_rca[rockfish$cumulative_catch_before_rca!=0])
min_after_kg <- min(rockfish$cumulative_catch_after_rca[rockfish$cumulative_catch_after_rca!=0])

rca_distkm <- c(0,1,5,10,15,25,40,50)
rca_dist_seq <- (-1*(rca_distkm-mean(rockfish$rca_distkm))/sd(rockfish$rca_distkm))
year_seq <- c(2007,2010,2013,2016,2019,2022)-2004
year_protected_seq <- (year_seq-mean(rockfish$year_protected))/sd(rockfish$year_protected)
rmax_seq_fact <- c(0.043,0.096,0.142)
rmax_seq <- (rmax_seq_fact-mean(rockfish$rmax))/sd(rockfish$rmax)
# catch_seq <- c(0,5000)
# cumul_catch_before_rca_seq <- (quadroot(catch_seq) - mean(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%"quillback"])))/sd(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%"quillback"]))

cumul_catch_before_rca_seq <- c(-1,1)

catch_seq <- round((cumul_catch_before_rca_seq*sd(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%"quillback"])) + mean(quadroot(rockfish$cumulative_catch_before_rca[rockfish$species%in%"quillback"])))^4,0)

size_trends <- emmeans(fit_sdm, ~ protection_dist_scale,at = list(protection_dist_scale = rca_dist_seq),type = "response")

rca_distkm_2 <- seq(from=0,to=50,by=5)
rca_dist_seq2 <- (-1*(rca_distkm_2-mean(rockfish$rca_distkm))/sd(rockfish$rca_distkm))
size_trends2 <- emmeans(fit_sdm, ~ protection_dist_scale,at = list(protection_dist_scale = rca_dist_seq2),type = "response")
emm_df2 <- as.data.frame(summary(size_trends2))
emm_df2$rca_distkm <- rca_distkm_2

openxlsx::write.xlsx(emm_df2,file="Results/body sizes at distances.xlsx")

emm_df <- as.data.frame(summary(size_trends))
# size_trends <- broom::tidy(size_trends,conf.int = TRUE)

p1 <- ggplot(data=emm_df,aes(x=-1*rca_distkm,y=response),fill="#7fcdbb",colour="#7fcdbb") +
  geom_ribbon(data=emm_df,alpha=0.5,aes(x=-1*rca_distkm,ymin=lower.CL,ymax=upper.CL),fill="#7fcdbb",colour="#7fcdbb") +
  geom_line(data=emm_df,aes(x=-1*rca_distkm,y=response),colour="#7fcdbb")+
  xlab("Proximity to nearest RCA (km)") + ylab("Total length (cm)") +
  theme_minimal() +
  theme(legend.position="none",text = element_text(size=11),legend.text = element_text(size=10),strip.text=element_text(size=9))
print(p1)
rca_slopes_over_time <- emtrends(fit_sdm, ~ year_protected_scaled | rmax_scale | cumul_catch_before_rca_scaled, var = "protection_dist_scale",
                                 at = list(year_protected_scaled = year_protected_seq,rmax_scale=rmax_seq,
                                           cumul_catch_before_rca_scaled = cumul_catch_before_rca_seq))
df <- data.frame("RCA_age"=rep(rep(year_seq,times=length(rmax_seq_fact)),times=length(cumul_catch_before_rca_seq)),
                 "rmax"=rep(rep(rmax_seq_fact,each=length(year_seq)),times=length(cumul_catch_before_rca_seq)),
                 "cumul_catch"=rep(rep(catch_seq,each=length(year_seq)),each=length(rmax_seq_fact)),
                 broom::tidy(rca_slopes_over_time))

quant_interval <- qnorm(0.95) #inner 90% interval
likely_interval <- qnorm(0.83) #inner 90% interval

df$ui <- df$protection_dist_scale.trend+quant_interval*df$std.error
df$li <- df$protection_dist_scale.trend-quant_interval*df$std.error

df$ui_likely <- df$protection_dist_scale.trend+likely_interval*df$std.error
df$li_likely <- df$protection_dist_scale.trend-likely_interval*df$std.error

set.seed(2023)
df <- df %>%
  group_by(rmax) %>%
  mutate(jit=jitter(x=0,factor=50),RCA_jitter=RCA_age + jit) %>%
  ungroup()
df$RCA_jitter <- df$RCA_age + (df$jit-mean(df$jit))
df$ipcc_uncertainty <- ifelse(df$li>0,"Very likely",ifelse(df$li_likely>0,"Likely","Not very likely"))
df$ipcc_uncertainty <- factor(df$ipcc_uncertainty,levels=c("Very likely","Likely","Not very likely"))
print(list("RCA distance effect along recent (+2 std. dev) and early years (-2 std. devs):",rca_slopes_over_time))
# df$RCA_age <- factor(df$RCA_age,levels=unique(df$RCA_age),labels=c(paste("RCA age",unique(df$RCA_age),"year(s)",sep=" ")))
df$rmax <- factor(df$rmax,levels=(unique(df$rmax)))
df$RCA_age <- factor(df$RCA_age,levels=rev(unique(df$RCA_age)))
df$cumul_catch <- factor(df$cumul_catch,levels=catch_seq,labels=c("-1\u03C3 of harvest","+1\u03C3 of harvest"))
p2 <- ggplot(data=df,aes(x=protection_dist_scale.trend,y=RCA_age,group=as.factor(rmax),fill=as.factor(rmax),colour=as.factor(rmax),shape=ipcc_uncertainty)) +
  geom_segment(aes(x=li,xend=ui,y=RCA_age,colour=as.factor(rmax)),position = position_dodge(width = 0.5)) +
  geom_segment(aes(x=li_likely,xend=ui_likely,y=RCA_age,colour=as.factor(rmax)),position = position_dodge(width = 0.5),lwd=1.25) +
  # geom_pointrange(aes(xmin=li,xmax=ui,y=as.factor(RCA_age)),position=position_jitter(seed=2023)) +
  # geom_jitter(colour="black",position=position_jitter(seed=2023))+
  facet_wrap(~cumul_catch,nrow=1)+
  geom_point(colour="black",position = position_dodge(width = 0.5)) +
  ylab("RCA Age") + xlab("Effect of proximity to nearest RCA (km)") +
  theme_minimal() + theme(legend.position="top") +
  scale_y_discrete(position = "left") +
  geom_vline(xintercept=0,lty=3)+
  annotate("text",label="More effective RCAs",x=0.05,y=6.7,size=2.75) +
  annotate("text",label="Less effective RCAs",x=-0.05,y=6.7,size=2.75) +
  coord_cartesian(clip = "off")+
  geom_segment(aes(x = 0.01, y = 6.4, xend = 0.05, yend = 6.4), 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),colour="black") +
  geom_segment(aes(x = -0.01, y = 6.4, xend = -0.05, yend = 6.4), 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),colour="black") +
  # scale_y_discrete(limits = rev) +
  # scale_fill_viridis_d(name=expression(r[max]),option="D") + 
  # scale_fill_manual(name="Likelihood of impact",values=c("dodgerblue","orange")) +
  # scale_colour_manual(name="Likelihood of impact",values=c("dodgerblue","orange")) + 
  scale_colour_viridis_d(name=expression(r[max]),option="D") +
  scale_fill_viridis_d(name=expression(r[max]),option="D") +
  scale_shape_manual(name="Likelihood of impact",values=c(24,23,25)) +
  guides(shape="none")+
  theme(text = element_text(size=11),legend.text = element_text(size=10),strip.text=element_text(size=9,hjust = 0.4))
p2

megaP <- cowplot::align_plots(p1, p2,align="v",axis="tblr")

cowplot::ggdraw(megaP[[1]]) + cowplot::draw_plot(megaP[[2]])

megaP1_align <- cowplot::ggdraw(megaP[[1]]) 
megaP2_align <- cowplot::ggdraw(megaP[[2]])

megaP <- ggpubr::ggarrange(megaP1_align,megaP2_align,nrow=2,heights=c(0.5,1),labels=c("(a)","(b)"))

ggsave(filename = "Figures/marginal effects of RCA interactions emmeans body size.jpeg",plot=megaP,dpi=600,units="in",height=6,width=7)

