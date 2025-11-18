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
count_rca <- readRDS("R/sdm_rca_predictions.rds")
size_rca <- readRDS("R/sdm_size_rca_predictions.rds")
rockfish <- readRDS(file="R/sdm_model_data.rds")
rockfish_size <- readRDS(file="R/sdm_model_data_size.rds")
species <- sort(unique(rockfish$species))

count_rca$metric <- "CPUE"
size_rca$metric <- "Body Size"
rca_effects <- rbind(count_rca,size_rca)
rca_effects$p <- sqrt((pmax(0,pmin(1,(rca_effects$ui-0)/(rca_effects$ui-rca_effects$li))) - 0.5)^2) / 0.5
rca_effects$metric <- factor(rca_effects$metric,levels=c("CPUE","Body Size"))
rca_effects$species[rca_effects$species=="dusky"] <- "dusky-dark"
rca_effects$species <- factor(rca_effects$species,levels=rev(levels(rockfish$species)))

rca_effects$rca_benefits_ipcc <- ifelse(rca_effects$p_above>=0.99,"Very likely",
                                            ifelse(rca_effects$p_above>=0.90,"Very likely",
                                                   ifelse(rca_effects$p_above>=0.66,"Likely",
                                                          ifelse(rca_effects$p_above>=0.33,"As likely as not",
                                                                 ifelse(rca_effects$p_above>=0.10,"Unlikely",
                                                                        ifelse(rca_effects$p_above>=0.01,"Very unlikely","Very unlikely"))))))

rca_effects$rca_benefits_ipcc <- factor(rca_effects$rca_benefits_ipcc,levels=c("Very likely","Likely","As likely as not","Unlikely","Very unlikely"))

round(table(rca_effects$rca_benefits_ipcc,rca_effects$metric)/colSums(table(rca_effects$rca_benefits_ipcc,rca_effects$metric)),2)
table(rca_effects$rca_benefits_ipcc,rca_effects$RCA)
table(rca_effects$rca_benefits_ipcc,rca_effects$species)
table(rca_effects$rca_benefits_ipcc,rca_effects$species,rca_effects$metric)

rca_effects$species <- gsub("vermillion","vermilion",rca_effects$species)
rca_effects$species <- gsub("sebastolobus","shortspine_thornyhead",rca_effects$species)

rockfish$species <- gsub("vermillion","vermilion",rockfish$species)
rockfish$species <- gsub("sebastolobus","shortspine_thornyhead",rockfish$species)

rockfish_size$species <- gsub("vermillion","vermilion",rockfish_size$species)
rockfish_size$species <- gsub("sebastolobus","shortspine_thornyhead",rockfish_size$species)

species <- sort(unique(rockfish$species))


openxlsx::write.xlsx(rca_effects,file="Data/rca effects by species and metric.xlsx")
rca_effects$RCA <- gsub("West Aristazabal Island","West Aristazabal",rca_effects$RCA)
rca_effects$metric <- factor(rca_effects$metric,levels=c("CPUE","Body Size"),labels=c("log(CPUE)","log(total length)"))
p1 <- ggplot(data=rca_effects[rca_effects$metric=="log(CPUE)",],aes(x=mean,y=species,fill=rca_benefits_ipcc))+
  geom_vline(xintercept=0,colour="red",linetype=2)+
  geom_segment(aes(x = li,xend=ui,y=species),colour="black")+
  geom_point(pch=21,size=2)+
  facet_grid(rows=vars(metric),cols=vars(RCA)) +
  scale_fill_brewer(name="Pr(\u03B2)>0",type="div",palette="PRGn",direction = 1,,limits=c("Very likely","Likely","As likely as not","Unlikely","Very unlikely")) +
  xlab("") + ylab("Species") +
  theme_minimal() +
  theme(legend.position="top") +
  scale_x_continuous(n.breaks = 4) +
  scale_y_discrete(limits=rev)+
  theme(text = element_text(size=12),,axis.text.x = element_text(size=10),axis.text.y = element_text(size=8.5),strip.text.x=element_text(size=7,hjust = -0),panel.grid.minor.x = element_blank(),strip.clip = "off")

p2 <- ggplot(data=rca_effects[rca_effects$metric=="log(total length)",],aes(x=mean,y=species,fill=rca_benefits_ipcc))+
  geom_vline(xintercept=0,colour="red",linetype=2)+
  geom_segment(aes(x = li,xend=ui,y=species),colour="black")+
  geom_point(pch=21,size=2)+
  facet_grid(rows=vars(metric),cols=vars(RCA)) +
  scale_fill_brewer(name="Pr(\u03B2)>0",type="div",palette="PRGn",direction = 1,limits=c("Very likely","Likely","As likely as not","Unlikely","Very unlikely")) +
  xlab("Effect of proximity to nearest RCA") + ylab("Species") +
  theme_minimal() +
  theme(legend.position="top") +
  scale_x_continuous(n.breaks = 3) +
  scale_y_discrete(limits=rev) +
  theme(text = element_text(size=12),axis.text.x = element_text(size=10),axis.text.y = element_text(size=8.5),strip.text.x=element_text(size=7,hjust = -0),panel.grid.minor.x = element_blank(),strip.clip = "off")

ggpubr::ggarrange(p1,p2,nrow=2,ncol=1,common.legend = TRUE,labels=c("(a)","(b)"),heights=c(nrow(rca_effects[rca_effects$metric=="log(CPUE)",]),nrow(rca_effects[rca_effects$metric=="log(total length)",])))
ggsave(filename="Figures/rca_effects_scatter.jpeg",dpi=600,height=8,width=7)


rca_table <- rca_effects %>%
  group_by(species,metric) %>%
  summarise("Positive_effects"=sum(if_else(rca_benefits_ipcc=="Very likely" | rca_benefits_ipcc=="Likely",1,0)),"No. sites"=n())

cpue_dist <- rockfish %>%
  group_by(species) %>%
  summarise("Mean distance"=mean(rca_distkm),"Median distance"=median(rca_distkm),"No. inside RCA"=sum(rca_distkm==0))

size_dist <- rockfish_size %>%
  group_by(species) %>%
  summarise("Mean distance"=mean(rca_distkm),"Median distance"=median(rca_distkm),"No. inside RCA"=sum(rca_distkm==0))
cpue_dist$metric <- "log(CPUE)"
size_dist$metric <- "log(total length)"
trait_dist <- rbind(cpue_dist,size_dist)
rca_table <- left_join(rca_table,trait_dist,by=c("species","metric"))
print(rca_table,n=1)

fit_sdm <- readRDS(file="R/sdm_groundfish_size_re.rds")
# fit_sdm <- sdmTMB:::reload_model(fit_sdm)
tidy_fit <- tidy(fit_sdm,conf.int=TRUE,effects=c("fixed"))
tidy_fit$group_name <- "global"
tidy_fit$par_name <- tidy_fit$term
tidy_fit$level_ids <- "global"

ran_pars <- tidy(fit_sdm,conf.int=TRUE,effects=c("ran_pars"))
ran_pars <- ran_pars[,!colnames(ran_pars)=="model"]
ran_pars$group_name <- "global"
ran_pars$par_name <- ran_pars$term
ran_pars$level_ids <- "global"
tidy_fit <- rbind(tidy_fit,ran_pars)

re_tidy <-tidy(fit_sdm,conf.int=TRUE,effects=c("ran_vals"))
re_tidy$par_name <- paste(re_tidy$term,re_tidy$level_ids,sep="-")
colnames(re_tidy)[colnames(re_tidy)=="conf.hi"] <- "conf.high"
re_tidy <- re_tidy[,!colnames(re_tidy)%in%c("model")]
tidy_fit <- rbind(tidy_fit,re_tidy)
print(tidy_fit,n=199)
n_sims <- 10000
samps <- sdmTMB:::rmvnorm_prec(fit_sdm$tmb_obj$env$last.par.best, fit_sdm$sd_report, n_sims)
b_j <- which(row.names(samps) == "b_j")[which(tidy_fit$par_name=="protection_dist_scale")]
re_b_pars <- which(row.names(samps) == "re_b_pars") # intercepts followed by random slopes here...
re <- samps[re_b_pars, ]
rca_names <- sort(unique(rockfish$rca_name))
fixed_rca_intercepts <- FALSE
if(fixed_rca_intercepts)
{
  rca_intercepts_idx <- 0#seq(1,length(rca_names),by=1)
  rca_slopes_idx <- seq(1,length(rca_names),by=1)
  rca_intercepts <- re[rca_intercepts_idx,]
  rca_slopes <- re[rca_slopes_idx,]
}else
{
  rca_intercepts_idx <- seq(1,2*length(rca_names),by=2)
  rca_slopes_idx <- seq(2,2*length(rca_names),by=2)
  rca_intercepts <- re[rca_intercepts_idx,]
  rca_slopes <- re[rca_slopes_idx,]
}

b <- samps[b_j,,drop=FALSE]
rca_names_idx <- as.numeric(rep(rca_names,each=length(species)))
out <- matrix(nrow = nrow(rca_slopes), ncol = ncol(rca_slopes))
for (i in 1:ncol(rca_slopes)) {
  for(j in 1:nrow(out)) {
    out[j,i] <- b[i] + rca_slopes[j,i]
  }
}
size_RCA <- apply(out,1,function(x){sum(x>0)/length(x)})


fit_sdm <- readRDS(file="R/sdm_groundfish_re6.rds")
# fit_sdm <- sdmTMB:::reload_model(fit_sdm)
tidy_fit <- tidy(fit_sdm,conf.int=TRUE,effects=c("fixed"))
tidy_fit$group_name <- "global"
tidy_fit$par_name <- tidy_fit$term
tidy_fit$level_ids <- "global"

ran_pars <- tidy(fit_sdm,conf.int=TRUE,effects=c("ran_pars"))
ran_pars <- ran_pars[,!colnames(ran_pars)=="model"]
ran_pars$group_name <- "global"
ran_pars$par_name <- ran_pars$term
ran_pars$level_ids <- "global"
tidy_fit <- rbind(tidy_fit,ran_pars)

re_tidy <-tidy(fit_sdm,conf.int=TRUE,effects=c("ran_vals"))
re_tidy$par_name <- paste(re_tidy$term,re_tidy$level_ids,sep="-")
colnames(re_tidy)[colnames(re_tidy)=="conf.hi"] <- "conf.high"
re_tidy <- re_tidy[,!colnames(re_tidy)%in%c("model")]
tidy_fit <- rbind(tidy_fit,re_tidy)
print(tidy_fit,n=75)
samps <- sdmTMB:::rmvnorm_prec(fit_sdm$tmb_obj$env$last.par.best, fit_sdm$sd_report, n_sims)
b_j <- which(row.names(samps) == "b_j")[which(tidy_fit$par_name=="protection_dist_scale")]
re_b_pars <- which(row.names(samps) == "re_b_pars") # intercepts followed by random slopes here...
re <- samps[re_b_pars, ]
rca_names <- sort(unique(rockfish$rca_name))
fixed_rca_intercepts <- FALSE
if(fixed_rca_intercepts)
{
  rca_intercepts_idx <- 0#seq(1,length(rca_names),by=1)
  rca_slopes_idx <- seq(1,length(rca_names),by=1)
  rca_intercepts <- re[rca_intercepts_idx,]
  rca_slopes <- re[rca_slopes_idx,]
}else
{
  rca_intercepts_idx <- seq(1,2*length(rca_names),by=2)
  rca_slopes_idx <- seq(2,2*length(rca_names),by=2)
  rca_intercepts <- re[rca_intercepts_idx,]
  rca_slopes <- re[rca_slopes_idx,]
}

b <- samps[b_j,,drop=FALSE]
rca_names_idx <- as.numeric(rep(rca_names,each=length(species)))
out <- matrix(nrow = nrow(rca_slopes), ncol = ncol(rca_slopes))
for (i in 1:ncol(rca_slopes)) {
  for(j in 1:nrow(out)) {
    out[j,i] <- b[i] + rca_slopes[j,i]
  }
}
cpue_RCA <- apply(out,1,function(x){sum(x>0)/length(x)})
names(cpue_RCA) <- rca_names
names(size_RCA) <- rca_names
cpue_RCA;size_RCA
sebastes <- data.frame("species"="Sebastes spp.","metric"=c("log(CPUE)","log(total length)"),"Positive_effects"=c(5,7),"No. sites"=c(7,7))
colnames(sebastes) <- c("species","metric","Positive_effects","No. sites")
any_cpue_dist <- rockfish %>%
  summarise("Mean distance"=mean(rca_distkm),"Median distance"=median(rca_distkm),"No. inside RCA"=sum(rca_distkm==0,na.rm=TRUE))

any_size_dist <- rockfish_size %>%
  summarise("Mean distance"=mean(rca_distkm),"Median distance"=median(rca_distkm),"No. inside RCA"=sum(rca_distkm==0,na.rm=TRUE))
any_cpue_dist$species <- "Sebastes spp."
any_cpue_dist$metric <- "log(CPUE)"
any_size_dist$species <- "Sebastes spp."
any_size_dist$metric <- "log(total length)"
any_dist <- rbind(any_cpue_dist,any_size_dist)
sebastes_table <- left_join(sebastes,any_dist,by=c("species","metric"))
mean_rmax <- rockfish %>%
  group_by(species) %>% summarise("rmax"=mean(rmax))
mean_rmax <- rbind(mean_rmax,data.frame("species"="Sebastes spp.","rmax"=mean(rockfish$rmax)))

rca_table_effects <- rbind(rca_table,sebastes_table)
rca_table_effects <- left_join(rca_table_effects,mean_rmax)
rca_table_effects <- rca_table_effects[order(rca_table_effects$rmax),]
print(rca_table_effects,n=50)
openxlsx::write.xlsx(rca_table_effects,file="Data/species RCA and distances.xlsx")
