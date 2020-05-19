# #################################################################### ####
# Title: Reducing dispersal limitation via seed addition leads to      ####
#        increased species richness, but not aboveground biomass       #### 
# Authors: Emma R Ladouceur & Shane A. Blowes                          ####
# Details: Main Figures 1-3                                            ####
# #################################################################### ####

# Libraries 
library(tidyverse)
library(brms)
library(ggplot2)
library(bayesplot)
library(patchwork)


# Data
#setwd('~/Desktop/Academic/R code/SeedAdditionSynthesis/')
setwd('~/Data/')
plot<-read.csv("./SeedAdd_Plot_Level.csv", header=TRUE) %>%
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
# 
# seedadd.biomass <- brm(l.biomass ~  seed.rich.m + (seed.rich.m | Experiment/site/block/fyr.trt),
#                    data = plot,
#                    chains = 4, cores = 4)

setwd('~/Model_fits/')
save(seedadd.rich, file = './rich.Rdata')
save(seedadd.biomass, file = './biomass.Rdata')


# richness model summary
summary(seedadd.rich)
# inspection of chain diagnostics
plot(seedadd.biomass) 

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



library(plyr)
S_seed_fitted$Experiment<-revalue(S_seed_fitted$Experiment, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))
S_seed_exp_coef2$Experiment<-revalue(S_seed_exp_coef2$Experiment, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))

S_seed_fitted$Experiment<-factor(as.character(S_seed_fitted$Experiment))
S_seed_exp_coef2$Experiment<-factor(as.character(S_seed_exp_coef2$Experiment))


biomass_fitted$Experiment<-revalue(biomass_fitted$Experiment, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))
biomass_exp_coef2$Experiment<-revalue(biomass_exp_coef2$Experiment, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))

biomass_fitted$Experiment<-factor(as.character(biomass_fitted$Experiment))
biomass_exp_coef2$Experiment<-factor(as.character(biomass_exp_coef2$Experiment))


# Figure 1 & 2 prep

S_seed_fixef_df<-as.data.frame(S_seed_fixef)
biomass_fixef_df<-as.data.frame((biomass_fixef))
S_seed_fixef_df$Model<-'Richness'
biomass_fixef_df$Model<-'Biomass'
fixef.all<-bind_rows(S_seed_fixef_df,biomass_fixef_df)
fixef.all

# Figure 1 b) Richness
fig1b<-ggplot() + 
  geom_point(data = S_seed_exp_coef2, aes(x = Experiment, y = Slope,colour = Experiment),size = 2) +
  geom_errorbar(data = S_seed_exp_coef2, aes(x = Experiment,ymin = Slope_lower,
                                             ymax = Slope_upper,colour = Experiment),
                width = 0, size = 1) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_hline(data = filter(fixef.all, Model=='Richness'),
             aes(yintercept = Estimate[2]), size = 1.2) +
  geom_rect(data = filter(fixef.all, Model=='Richness'),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(x = 'Experiment',
       y = 'Species Richness / species of seed added', title= "",
       subtitle="b)") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  scale_x_discrete(limits = rev(levels(S_seed_exp_coef2$Experiment)))+coord_flip() + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),legend.position="none",
                   #axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())

# Figure 2 b) Biomass
fig2b<-ggplot() + 
  geom_point(data = biomass_exp_coef2, aes(x = Experiment, y = Slope,colour = Experiment),size = 2) +
  geom_errorbar(data = biomass_exp_coef2, aes(x = Experiment,ymin = Slope_lower,
                                              ymax = Slope_upper,colour = Experiment),
                width = 0, size = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_hline(data = filter(fixef.all, Model=='Biomass'),
             aes(yintercept = Estimate[2]), size = 1.2) +
  geom_rect(data = filter(fixef.all, Model=='Biomass'),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(x = 'Experiment',
       y = expression(paste( 'Biomass [log(g/',m^2,')] / species of seed added')), title = '',
       subtitle= "b)") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  scale_x_discrete(limits = rev(levels(biomass_exp_coef2$Experiment)))+coord_flip() + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),legend.position="none",
                   #axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())


# Figure 1 a) Richness

fig1a<-ggplot() +
  geom_point(data = S_seed_fitted,
             aes(x = seed.rich, y = rich.plot,
                 colour = Experiment),
             size = 1.2, position = position_jitter(width = 0.95, height = 0.95)) +
  geom_segment(data = S_seed_exp_coef2,
               aes(x = xmin, 
                   xend = xmax,
                   y = (Intercept + Slope * cxmin),
                   yend = (Intercept + Slope * cxmax),
                   group = Experiment,
                   colour = Experiment),
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
  labs(x = 'Number of species of seed added ',
       y = 'Species richness', title= 'Species Richness', subtitle= "a)") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),
                                                                                                                             legend.position="bottom")


# Figure 2 a) Biomass

fig2a<-ggplot() +
  geom_point(data = biomass_fitted,
             aes(x = seed.rich, y = biomass.plot,
                 colour = Experiment),
             size = 1.2, position = position_jitter(width = 0.95, height = 0.95)) +
  geom_segment(data = biomass_exp_coef2,
               aes(x = xmin, 
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax),
                   group = Experiment,
                   colour = Experiment),
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
  labs(x = 'Number of species of seed added ',
       y = expression(paste('Biomass [log(g/',m^2, ')]')), title='Biomass', subtitle="a)") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),
                     legend.position="bottom")

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

m.legend<-g_legend(ba)


(fig1a+ theme(legend.position="none") | fig1b)/(m.legend)  +
  plot_layout(heights = c(10,1))

(fig2a+ theme(legend.position="none") | fig2b)/(m.legend)  +
  plot_layout(heights = c(10,1))



# Figure 3 prep
S_seed_exp_coef3<-S_seed_exp_coef2[,c(-1,-2,-3,-8,-9,-10,-11)]
biomass_exp_coef3<-biomass_exp_coef2[,c(-1,-2,-3,-8,-9,-10,-11)]
colnames(S_seed_exp_coef3)
colnames(biomass_exp_coef3)
names(S_seed_exp_coef3) <- c("Experiment","R.Slope","R.Slope_lower","R.Slope_upper")
names(biomass_exp_coef3) <- c("Experiment","B.Slope","B.Slope_lower","B.Slope_upper")
delta.coefs<-bind_cols(S_seed_exp_coef3,biomass_exp_coef3)


delta.coefs$Experiment<-revalue(delta.coefs$Experiment, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))

delta.coefs$Experiment<-factor(as.character(delta.coefs$Experiment))


# Figure 3

ggplot(data=delta.coefs, aes(x=R.Slope, y=B.Slope,color=Experiment)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin = B.Slope_lower, ymax = B.Slope_upper,colour = Experiment), width = 0, size = 0.75,alpha=0.5) +
  geom_errorbarh(aes(xmin = R.Slope_lower, xmax = R.Slope_upper,colour = Experiment), width = 0, size = 0.75,alpha=0.5) +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  labs(x = 'Species Richness Slope / number of species of seed added ',
       y = expression(paste('Biomass Slope [log(g/',m^2, ')] / number of species of seed added'))) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),legend.position="bottom")




