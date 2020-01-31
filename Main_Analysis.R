# #################################################################### ####
# Title: Comparing the effects of seeded diversity on species richness ####
#        and above ground biomass , and community assembly             #### 
# Authors: Emma R Ladouceur & Shane A. Blowes                          ####
# Details: Figures 1-3                                                 ####
# #################################################################### ####

# Libraries 
rm(list=ls())
detach("package:ggplot2", unload=TRUE)
detach("package:plyr", unload=TRUE)
library(tidyverse)
library(brms)
library(ggplot2)
library(gridExtra)
library(grid)
library(bayesplot)

# TIDYVERSE SHARED LEGEND (to make plots)
# https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

# Data
setwd('~/Desktop/Academic/R code/SeedAdditionSynthesis/')
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


# load model objects
load("./Model Fits/rich.Rdata") # object name: seedadd.rich
load("./Model Fits/biomass.Rdata") # object name: seedadd.biomass

# seedadd.rich<- brm(rich.plot ~  seed.rich.m + (seed.rich.m | Experiment/site/block/fyr.trt), 
#                  data = plot, cores = 4, chains = 4)

# seedadd.biomass <- brm(l.biomass ~  seed.rich.m + (seed.rich.m | Experiment/site/block/fyr.trt),
#                    data = plot, 
#                    chains = 4, cores = 4)


# richness model summary
summary(seedadd.rich)

# Figure S1a
color_scheme_set("darkgray")
pp_check(seedadd.rich)+ theme_classic() # predicted vs. observed values

# models residuals
m1<-residuals(seedadd.rich)
m1<-as.data.frame(m1)
rr.plot<-cbind(plot,m1$Estimate)

par(mfrow=c(3,2))
with(rr.plot, plot(Experiment, m1$Estimate))
with(rr.plot, plot(site, m1$Estimate))
with(rr.plot, plot(block, m1$Estimate))
with(rr.plot, plot(fyr.trt, m1$Estimate))
with(rr.plot, plot(seed.rich, m1$Estimate))

# make sure to detach plyr at top

# for plotting fixed effects
S_seed_fitted <- cbind(seedadd.rich$data,
                       fitted(seedadd.rich, re_formula = NA)) %>% 
  as_tibble() %>% 
  # join with plot data for figures
  inner_join(plot %>% distinct(Experiment, site,block, fyr.trt, seed.rich, seed.rich.m, rich.plot),
             by = c('Experiment','site','block','fyr.trt', 'seed.rich.m', 'rich.plot'))


# fixed effect coefficients
S_seed_fixef <- fixef(seedadd.rich)


S_seed_exp_coef <- coef(seedadd.rich)
S_seed_exp_coef 
S_seed_exp_coef2 <-  bind_cols(S_seed_exp_coef$Experiment[,,'Intercept'] %>% 
                                 as_tibble() %>% 
                                 mutate(Intercept = Estimate,
                                        Intercept_lower = Q2.5,
                                        Intercept_upper = Q97.5,
                                        Experiment = rownames(S_seed_exp_coef$Experiment[,,'Intercept'])) %>% 
                                 select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                               S_seed_exp_coef$Experiment[,,'seed.rich.m'] %>% 
                                 as_tibble() %>% 
                                 mutate(Slope = Estimate,
                                        Slope_lower = Q2.5,
                                        Slope_upper = Q97.5) %>% 
                                 select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  # join with min and max of the x-values
  inner_join(plot %>% 
               group_by(Experiment) %>% 
               summarise(xmin = min(seed.rich),
                         xmax = max(seed.rich),
                         cxmin = min(seed.rich.m),
                         cxmax = max(seed.rich.m)),
             by = 'Experiment')


# Biomass model
summary(seedadd.biomass)


# inspection of chain diagnostics
plot(seedadd.biomass) 

# Figure S1 b
pp_check(seedadd.biomass)+ theme_classic()

# residuals
m2<-residuals(seedadd.biomass)
m2<-as.data.frame(m2)
rb.plot<-cbind(plot,m2$Estimate)

par(mfrow=c(2,3))
with(rb.plot, plot(Experiment, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(site, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(block, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(fyr.trt, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(seed.rich, m2$Estimate));abline(h=0, lty=2)


biomass_fitted <- cbind(seedadd.biomass$data,
                        fitted(seedadd.biomass, re_formula = NA)) %>% 
  as_tibble() %>% 
  inner_join(plot %>% distinct(Experiment, site,block, fyr.trt, seed.rich, seed.rich.m, biomass.plot,l.biomass),
             by = c('Experiment','site','block','fyr.trt', 'seed.rich.m', 'l.biomass'))


biomass_fixef <- fixef(seedadd.biomass)

biomass_exp_coef <- coef(seedadd.biomass)

biomass_exp_coef2 <-  bind_cols(biomass_exp_coef$Experiment[,,'Intercept'] %>% 
                                  as_tibble() %>% 
                                  mutate(Intercept = Estimate,
                                         Intercept_lower = Q2.5,
                                         Intercept_upper = Q97.5,
                                         Experiment = rownames(biomass_exp_coef$Experiment[,,'Intercept'])) %>% 
                                  select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                                biomass_exp_coef$Experiment[,,'seed.rich.m'] %>% 
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


# use grid_arrange_shared_legend function at the beginning of  this script
# to create figures

library(plyr)
S_seed_fitted$Study<-revalue(S_seed_fitted$Experiment, c("ASGA_Michigan"="Michigan.us", "California_Invade"="California.I.us","California_Prop_Limi"="California.P.L.us","CCR_04"="CedarCreek4.us","CCR_093"="CedarCreek93.us","Germany_Montane"="Montane.de","Halle"="Halle.de","Jena"="Jena.de","Jena2"="Jena2.de","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow.us","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field.us","Texas_Temple_Prarie"="Texas.Templ.Prairie.us"))
S_seed_exp_coef2$Study<-revalue(S_seed_exp_coef2$Experiment, c("ASGA_Michigan"="Michigan.us", "California_Invade"="California.I.us","California_Prop_Limi"="California.P.L.us","CCR_04"="CedarCreek4.us","CCR_093"="CedarCreek93.us","Germany_Montane"="Montane.de","Halle"="Halle.de","Jena"="Jena.de","Jena2"="Jena2.de","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow.us","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field.us","Texas_Temple_Prarie"="Texas.Templ.Prairie.us"))

S_seed_fitted$Study<-factor(as.character(S_seed_fitted$Study))
S_seed_exp_coef2$Study<-factor(as.character(S_seed_exp_coef2$Study))


biomass_fitted$Study<-revalue(biomass_fitted$Experiment, c("ASGA_Michigan"="Michigan.us", "California_Invade"="California.I.us","California_Prop_Limi"="California.P.L.us","CCR_04"="CedarCreek4.us","CCR_093"="CedarCreek93.us","Germany_Montane"="Montane.de","Halle"="Halle.de","Jena"="Jena.de","Jena2"="Jena2.de","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow.us","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field.us","Texas_Temple_Prarie"="Texas.Templ.Prairie.us"))
biomass_exp_coef2$Study<-revalue(biomass_exp_coef2$Experiment, c("ASGA_Michigan"="Michigan.us", "California_Invade"="California.I.us","California_Prop_Limi"="California.P.L.us","CCR_04"="CedarCreek4.us","CCR_093"="CedarCreek93.us","Germany_Montane"="Montane.de","Halle"="Halle.de","Jena"="Jena.de","Jena2"="Jena2.de","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow.us","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field.us","Texas_Temple_Prarie"="Texas.Templ.Prairie.us"))

biomass_fitted$Study<-factor(as.character(biomass_fitted$Study))
biomass_exp_coef2$Study<-factor(as.character(biomass_exp_coef2$Study))


# Figure 1 prep

S_seed_fixef_df<-as.data.frame(S_seed_fixef)
biomass_fixef_df<-as.data.frame((biomass_fixef))
S_seed_fixef_df$Model<-'Richness'
biomass_fixef_df$Model<-'Biomass'
fixef.all<-bind_rows(S_seed_fixef_df,biomass_fixef_df)


rc<-ggplot() + 
  geom_point(data = S_seed_exp_coef2, aes(x = Study, y = Slope,colour = Study),size = 4) +
  geom_errorbar(data = S_seed_exp_coef2, aes(x = Study,ymin = Slope_lower,
                                             ymax = Slope_upper,colour = Study),
                width = 0, size = 1.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_hline(data = filter(fixef.all, Model=='Richness'),
             aes(yintercept = Estimate[2]), size = 1.2) +
  geom_rect(data = filter(fixef.all, Model=='Richness'),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(x = 'Study',
       y = 'Slope', title= "a) Richness") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  scale_x_discrete(limits = rev(levels(S_seed_exp_coef2$Study)))+coord_flip() + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),legend.position="bottom",axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())

bc<-ggplot() + 
  geom_point(data = biomass_exp_coef2, aes(x = Study, y = Slope,colour = Study),size = 4) +
  geom_errorbar(data = biomass_exp_coef2, aes(x = Study,ymin = Slope_lower,
                                              ymax = Slope_upper,colour = Study),
                width = 0, size = 1.5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_hline(data = filter(fixef.all, Model=='Biomass'),
             aes(yintercept = Estimate[2]), size = 1.2) +
  geom_rect(data = filter(fixef.all, Model=='Biomass'),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(x = 'Study',
       y = 'Slope', title = expression(paste('b) Biomass [log(g/',m^2, ')]'))) +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  scale_x_discrete(limits = rev(levels(biomass_exp_coef2$Study)))+coord_flip() + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),legend.position="bottom",axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())

# Figure 1
grid_arrange_shared_legend(rc,bc,nrow=1)

# Figure 2 a) Richness

r1<-ggplot() +
  geom_point(data = S_seed_fitted,
             aes(x = seed.rich, y = rich.plot,
                 colour = Study),
             size = 1.2, position = position_jitter(width = 0.95, height = 0.95)) +
  geom_segment(data = S_seed_exp_coef2,
               aes(x = xmin, 
                   xend = xmax,
                   y = (Intercept + Slope * cxmin),
                   yend = (Intercept + Slope * cxmax),
                   group = Study,
                   colour = Study),
               size = 1.2) +
  # uncertainy in fixed effect
  geom_ribbon(data = S_seed_fitted,
              aes(x = seed.rich, ymin = Q2.5, ymax = Q97.5),
              alpha = 0.3) +
  # fixed effect
  geom_line(data = S_seed_fitted,
            aes(x = seed.rich, y = Estimate),
            size = 1.5) +
  #scale_y_continuous(trans = 'log', breaks = c(2,4,8,16.24,36,48,64)) +
  labs(x = 'Number of species added ',
       y = 'Species richness', title= 'a) Richness') +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"))


# Figure 2 b) Biomass

b1<-ggplot() +
  geom_point(data = biomass_fitted,
             aes(x = seed.rich, y = biomass.plot,
                 colour = Study),
             size = 1.2, position = position_jitter(width = 0.95, height = 0.95)) +
  geom_segment(data = biomass_exp_coef2,
               aes(x = xmin, 
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax),
                   group = Study,
                   colour = Study),
               size = 1.2) +
  geom_ribbon(data = biomass_fitted,
              aes(x = seed.rich, 
                  ymin = exp(Q2.5), 
                  ymax = exp(Q97.5)),
              alpha = 0.3) +
  # fixed effect
  geom_line(data = biomass_fitted,
            aes(x = seed.rich, y = exp(Estimate)),
            size = 1.5) +
  scale_y_continuous(trans = 'log', breaks = c(8, 64, 512, 1024, 2048)) +
  labs(x = 'Number of species added ',
       y = expression(paste('Biomass [log(g/',m^2, ')]')), title='b) Biomass') +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"))

# Figure 2
grid_arrange_shared_legend(r1,b1,nrow=1)



# Figure 3 prep
S_seed_exp_coef3<-S_seed_exp_coef2[,c(-1,-2,-3,-8,-9,-10,-11,-13)]
biomass_exp_coef3<-biomass_exp_coef2[,c(-1,-2,-3,-8,-9,-10,-11,-13)]
colnames(S_seed_exp_coef2)
colnames(biomass_exp_coef3)
names(S_seed_exp_coef3) <- c("Experiment","R.Slope","R.Slope_lower","R.Slope_upper","Study")
names(biomass_exp_coef3) <- c("Experiment","B.Slope","B.Slope_lower","B.Slope_upper","Study")
delta.coefs<-bind_cols(S_seed_exp_coef3,biomass_exp_coef3)


delta.coefs$Study<-revalue(delta.coefs$Experiment, c("ASGA_Michigan"="Michigan.us", "California_Invade"="California.I.us","California_Prop_Limi"="California.P.L.us","CCR_04"="CedarCreek4.us","CCR_093"="CedarCreek93.us","Germany_Montane"="Montane.de","Halle"="Halle.de","Jena"="Jena.de","Jena2"="Jena2.de","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow.us","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field.us","Texas_Temple_Prarie"="Texas.Templ.Prairie.us"))

delta.coefs$Study<-factor(as.character(delta.coefs$Study))


# Figure 3
ggplot(data=delta.coefs, aes(x=R.Slope, y=B.Slope,color=Study)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = B.Slope_lower, ymax = B.Slope_upper,colour = Study), width = 0, size = 0.75,alpha=0.5) +
  geom_errorbarh(aes(xmin = R.Slope_lower, xmax = R.Slope_upper,colour = Study), width = 0, size = 0.75,alpha=0.5) +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  labs(x = 'Species Richness Slope',
       y = expression(paste('Biomass Slope [log(g/',m^2, ')]'))) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),legend.position="bottom")




