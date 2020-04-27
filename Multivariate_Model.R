# #################################################################### ####
# Title: Comparing the effects of seeded diversity on species richness ####
#        and above ground biomass , and community assembly             #### 
# Authors: Emma R Ladouceur & Shane A. Blowes                           ####
# #################################################################### ####

rm(list=ls())
detach("package:ggplot2", unload=TRUE)
detach("package:plyr", unload=TRUE)
library(tidyverse)
library(brms)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(bayesplot)

#setwd('~/Desktop/Academic/R code/SeedAdditionSynthesis/')
setwd('~/Dropbox/Projects/SeedAddDraft/')
plot<-read.csv("./Data/SeedAdd_Plot_Level.csv", header=TRUE) %>%
  as_tibble()

plot$fyr.trt<-as.factor(plot$yr.trt)
plot$seed.rich<-as.numeric(as.character(plot$seed.rich))
plot$site<-as.factor(plot$site)
plot$block<-as.factor(plot$block)
# Centered seed richness
plot$seed.rich.m<-plot$seed.rich-mean(plot$seed.rich)
# log-biomass
plot$l.biomass <- log(plot$biomass.plot)
plot$Experiment<-plot$Experiment_
head(plot)

summary(plot)
plot2<-plot %>% drop_na(live_mass)

# load model object
load("./Model Fits/multi.Rdata") # object name: seedadd.multi

 seedadd.multi <- brm(mvbind(rich.plot, l.biomass) ~  seed.rich.m + (seed.rich.m | p | Experiment), 
               data = plot, cores = 4, chains = 4)


 setwd('~/Dropbox/Projects/SeedAdd/Model_fits/')
 save(seedadd.multi, file = './multi.Rdata')
 
 
# check  correlation coefficients between variables
summary(seedadd.multi)

plot(seedadd.multi)

pmb2<-pp_check(seedadd.multi, resp = 'lbiomass')+ theme_classic()
pmr2<-pp_check(seedadd.multi, resp = 'richplot')+ theme_classic()
# Figure S1 f
(pmr2 | pmb2)


mm1<-residuals(seedadd.multi)
mm1
plot <- cbind(plot,
              residual_mm1_rich = mm1[,,'richplot'][,'Estimate'],
              residual_mm1_biomass = mm1[,,'lbiomass'][,'Estimate'])

par(mfrow=c(2,2))
with(plot, plot(Experiment, residual_mm1_rich))
with(plot, plot(Experiment, residual_mm1_biomass))

with(plot, plot(site, residual_mm1_rich))
with(plot, plot(site, residual_mm1_biomass))

par(mfrow=c(2,2))
with(plot, plot(block, residual_mm1_rich))
with(plot, plot(block, residual_mm1_biomass))

with(plot, plot(fyr.trt, residual_mm1_rich))
with(plot, plot(fyr.trt, residual_mm1_biomass))
par(mfrow=c(1,2))
with(plot, plot(seed.rich, residual_mm1_rich))
with(plot, plot(seed.rich, residual_mm1_biomass))


# make sure to detach plyr at top

mm_fitted <- cbind(seedadd.multi$data,
                   fitted(seedadd.multi, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(plot %>% distinct(Experiment, seed.rich, seed.rich.m, biomass.plot,rich.plot),
             by = c('Experiment', 'seed.rich.m','rich.plot')) 


mm_fixef <- fixef(seedadd.multi)


mm_coef <- coef(seedadd.multi)
mm_coef
mm_coefr <-  bind_cols(mm_coef$Experiment[,,'richplot_Intercept'] %>% 
                         as_tibble() %>% 
                         mutate(Intercept = Estimate,
                                Intercept_lower = Q2.5,
                                Intercept_upper = Q97.5,
                                Experiment = rownames(mm_coef$Experiment[,,'richplot_Intercept'])) %>% 
                         select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                       mm_coef$Experiment[,,'richplot_seed.rich.m'] %>% 
                         as_tibble() %>% 
                         mutate(Slope = Estimate,
                                Slope_lower = Q2.5,
                                Slope_upper = Q97.5) %>% 
                         select(-Estimate, -Est.Error, -Q2.5, -Q97.5))  %>% 
  inner_join(plot %>% 
               group_by(Experiment) %>% 
               summarise(xmin = min(seed.rich),
                         xmax = max(seed.rich),
                         cxmin = min(seed.rich.m),
                         cxmax = max(seed.rich.m)),
             by = 'Experiment') 



mm_coefb <-  bind_cols(mm_coef$Experiment[,,'lbiomass_Intercept'] %>% 
                         as_tibble() %>% 
                         mutate(Intercept = Estimate,
                                Intercept_lower = Q2.5,
                                Intercept_upper = Q97.5,
                                Experiment = rownames(mm_coef$Experiment[,,'lbiomass_Intercept'])) %>% 
                         select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                       mm_coef$Experiment[,,'lbiomass_seed.rich.m'] %>% 
                         as_tibble() %>% 
                         mutate(Slope = Estimate,
                                Slope_lower = Q2.5,
                                Slope_upper = Q97.5) %>% 
                         select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  inner_join(plot %>% 
               group_by(Experiment) %>% 
               summarise(xmin = min(seed.rich),
                         xmax = max(seed.rich),
                         cxmin = min(seed.rich.m),
                         cxmax = max(seed.rich.m)),
             by = 'Experiment')

# use grid_arrange_shared_legend function at the beginning of Main_Analysis.R
# to create figures

#Delta plot
mm_coefr2<-mm_coefr[,c(-1,-2,-3,-8,-9,-10,-11)]
mm_coefb2<-mm_coefb[,c(-1,-2,-3,-8,-9,-10,-11)]
colnames(mm_coefr2)
colnames(mm_coefb2)
names(mm_coefr2) <- c("Experiment","R.Slope","R.Slope_lower","R.Slope_upper")
names(mm_coefb2) <- c("Experiment","B.Slope","B.Slope_lower","B.Slope_upper")
m.delta.coefs<-bind_cols(mm_coefr2,mm_coefb2)



library(plyr)
m.delta.coefs$Study<-revalue(m.delta.coefs$Experiment, c("ASGA_Michigan"="Michigan.us", "California_Invade"="California.I.us","California_Prop_Limi"="California.P.L.us","CCR_04"="CedarCreek4.us","CCR_093"="CedarCreek93.us","Germany_Montane"="Montane.de","Halle"="Halle.de","Jena"="Jena.de","Jena2"="Jena2.de","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow.us","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field.us","Texas_Temple_Prarie"="Texas.Templ.Prairie.us"))
m.delta.coefs$Study<-factor(as.character(m.delta.coefs$Study))

ggplot(data=m.delta.coefs, aes(x=R.Slope, y=B.Slope,color=Study)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = B.Slope_lower, ymax = B.Slope_upper,colour = Study), 
                width = 0, size = 0.75,alpha=0.5) +
  geom_errorbarh(aes(xmin = R.Slope_lower, xmax = R.Slope_upper,colour = Study), 
                 width = 0, size = 0.75,alpha=0.5) +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  labs(x = 'Species Richness Estimate',
       y = expression(paste('Biomass Estimate [log(g/',m^2, ')]'))) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"),
        legend.position="bottom")

# Results are qualitatively consistent with univaraite models
# Estimates vary slightly, but uncertainty has the same implications
