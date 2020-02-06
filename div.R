# #################################################################### ####
# Title: Comparing the effects of seeded diversity on species richness ####
#        and above ground biomass , and community assembly             #### 
# Authors: Emma R Ladouceur & Shane A. Blowes                           ####
# Details: Figures S3                                                  ####
# #################################################################### ####

# Libraries
rm(list=ls())
detach("package:ggplot2", unload=TRUE)
detach("package:plyr", unload=TRUE)
library(tidyverse)
library(tibble)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(brms)
library(vegan)
library(bayesplot)


# Data
#setwd('~/Desktop/Academic/R code/SeedAdditionSynthesis/'))
setwd('~/Dropbox/Projects/SeedAddDraft/')
spdat<-read.csv("./Data/SeedAdd_Sp_level.csv", header=TRUE) 

spdat.bm<- spdat %>% select(unique.id_,species,biomass.sp)

# wide format
plot.wide.bm<-spread(spdat.bm,species,biomass.sp)

# replace NA with 0
plot.wide.bm[is.na(plot.wide.bm)] <- 0

# make column 1 row names
rownames(plot.wide.bm) <- plot.wide.bm[,1]
#remove old column
plot.wide.bm <- plot.wide.bm[,-1]

bm.seed.pie<-diversity(plot.wide.bm, index = 'invsimpson')
bm.seed.pie<-as.data.frame(bm.seed.pie)

bm.seed.pie<-setDT(bm.seed.pie, keep.rownames = TRUE)[]

bm.seed.pie$biomass.pie<-bm.seed.pie$bm.seed.pie
bm.seed.pie$unique.id_<-bm.seed.pie$rn
bm.seed.pie <- bm.seed.pie[,-1]
bm.seed.pie <- bm.seed.pie[,-1]

sp.pie<-left_join(spdat,bm.seed.pie)

colnames(sp.pie)

seed.pie1<- sp.pie %>% select(unique.id_,seed.rich,Experiment_,site,block,yr.trt,biomass.pie, rich.plot,biomass.plot,biomass.sp)

# plot level
seed.pie<-distinct(seed.pie1,unique.id_, .keep_all= TRUE)

seed.pie$fyr.trt<-as.factor(seed.pie$yr.trt)
seed.pie$seed.rich<-as.numeric(as.character(seed.pie$seed.rich))
seed.pie$site<-as.factor(seed.pie$site)
seed.pie$block<-as.factor(seed.pie$block)
# Centered seed richness
seed.pie$seed.rich.m<-seed.pie$seed.rich-mean(seed.pie$seed.rich)


seed.pie$Experiment<-seed.pie$Experiment_
seed.pie$l.b.pie <- log(seed.pie$biomass.pie)

summary(seed.pie)

#write.csv(seed.pie,"~/Dropbox/Projects/SeedAdd/Data/seed.pie.csv")


# Diversity Models

#setwd('~/Desktop/Academic/R code/SeedAdditionSynthesis/')
setwd('~/Dropbox/Projects/SeedAddDraft/')
spdat<-read.csv("./Data/seed.pie.csv", header=TRUE) %>%
  as_tibble()


levels(seed.pie$Experiment)
seed.pie$fyr.trt<-as.factor(seed.pie$yr.trt)
seed.pie$seed.rich<-as.numeric(as.character(seed.pie$seed.rich))
seed.pie$site<-as.factor(seed.pie$site)
seed.pie$block<-as.factor(seed.pie$block)
seed.pie$seed.rich.m<-seed.pie$seed.rich-mean(seed.pie$seed.rich)
seed.pie$l.b.pie <- log(seed.pie$biomass.pie)

hist(seed.pie$biomass.pie)
hist(seed.pie$l.b.pie)

# load model object
load("./Model Fits/seedpie.Rdata") # object name: m.l.seed.pie

# m.l.seed.pie <- brm(l.b.pie ~  seed.rich.m + (seed.rich.m | Experiment/site/block/fyr.trt),
#                     data = seed.pie, cores = 4, chains = 4)


summary(m.l.seed.pie)

plot(m.l.seed.pie)

# Figure S1 c
color_scheme_set("darkgray")
pp_check(m.l.seed.pie) + theme_classic()

# residuals from log model
m2<-residuals(m.l.seed.pie)
m2<-as.data.frame(m2)
s.p.<-seed.pie %>% drop_na(l.b.pie)

rb.plot<-cbind(s.p.,m2$Estimate)
head(rb.plot)
par(mfrow=c(3,2))
with(rb.plot, plot(Experiment, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(site, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(block, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(fyr.trt, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(seed.rich, m2$Estimate));abline(h=0, lty=2)

# make sure to detach plyr at top

l.pie_fitted <- cbind(m.l.seed.pie$data,
                    fitted(m.l.seed.pie, re_formula = NA))  %>%
  as_tibble() %>%
   inner_join(seed.pie %>% distinct(Experiment, site,block, fyr.trt, seed.rich, seed.rich.m, biomass.pie,l.b.pie),
            by = c('Experiment','site','block','fyr.trt', 'seed.rich.m', 'l.b.pie'))


l.pie_fixef <- fixef(m.l.seed.pie)

l.pie_exp_coef <- coef(m.l.seed.pie)

l.pie_exp_coef2 <-  bind_cols(l.pie_exp_coef$Experiment[,,'Intercept'] %>% 
                              as_tibble() %>% 
                              mutate(Intercept = Estimate,
                                     Intercept_lower = Q2.5,
                                     Intercept_upper = Q97.5,
                                     Experiment = rownames(l.pie_exp_coef$Experiment[,,'Intercept'])) %>% 
                              select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                            l.pie_exp_coef$Experiment[,,'seed.rich.m'] %>% 
                              as_tibble() %>% 
                              mutate(Slope = Estimate,
                                     Slope_lower = Q2.5,
                                     Slope_upper = Q97.5) %>% 
                              select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) 

pie_exp_coef2$Experiment<-as.factor(as.character(pie_exp_coef2$Experiment))
seed.pie$seed.rich<-as.numeric(as.character(seed.pie$seed.rich))

l.pie_exp_coef3 <- l.pie_exp_coef2 %>%  
  inner_join(seed.pie %>% 
               group_by(Experiment) %>% 
               summarise(xmin = min(seed.rich),
                         xmax = max(seed.rich),
                         cxmin = min(seed.rich.m),
                         cxmax = max(seed.rich.m)),
             by = 'Experiment')



# use grid_arrange_shared_legend function at the beginning of Main_Analysis.R
# to create figures

library(plyr)
l.pie_fitted$Experiment<-as.factor(as.character(l.pie_fitted$Experiment))
l.pie_fitted$Study<-revalue(l.pie_fitted$Experiment, c("ASGA_Michigan"="Michigan.us", "California_Invade"="California.I.us","California_Prop_Limi"="California.P.L.us","CCR_04"="CedarCreek4.us","CCR_093"="CedarCreek93.us","Germany_Montane"="Montane.de","Halle"="Halle.de","Jena"="Jena.de","Jena2"="Jena2.de","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow.us","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field.us","Texas_Temple_Prarie"="Texas.Templ.Prairie.us"))
l.pie_exp_coef3$Study<-revalue(l.pie_exp_coef3$Experiment, c("ASGA_Michigan"="Michigan.us", "California_Invade"="California.I.us","California_Prop_Limi"="California.P.L.us","CCR_04"="CedarCreek4.us","CCR_093"="CedarCreek93.us","Germany_Montane"="Montane.de","Halle"="Halle.de","Jena"="Jena.de","Jena2"="Jena2.de","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow.us","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field.us","Texas_Temple_Prarie"="Texas.Templ.Prairie.us"))

l.pie_fitted$Study<-factor(as.character(l.pie_fitted$Study))
l.pie_exp_coef3$Study<-factor(as.character(l.pie_exp_coef3$Study))

# regression
dr <- ggplot() +
  geom_point(data = l.pie_fitted,
             aes(x = seed.rich, y = biomass.pie,
                 colour = Study),
             size = 1.2, position = position_jitter(width = 0.95, height = 0.95)) +
  geom_segment(data = l.pie_exp_coef3,
               aes(x = xmin,
                   xend = xmax,
                   y = exp(Intercept + Slope * cxmin),
                   yend = exp(Intercept + Slope * cxmax),
                   group = Study,
                   colour = Study),
               size = 1.2) +
  geom_ribbon(data = l.pie_fitted,
              aes(x = seed.rich, ymin = exp(Q2.5), ymax = exp(Q97.5)),
              alpha = 0.3) +
  geom_line(data = l.pie_fitted,
            aes(x = seed.rich, y = exp(Estimate)),
            size = 1.5) +
  scale_y_continuous(trans = 'log', breaks = c(2,4,8, 16,24)) +
  labs(x = 'Number of species added',
       y = 'Diversity', title= 'b) Evenness Regression') +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF","#F9B90AFF" , "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"))

dr

#coefs 

#fixed
l.pie_fixef_df<-as.data.frame((l.pie_fixef))
l.pie_fixef_df

l.pie_exp_coef3$Study<-revalue(l.pie_exp_coef3$Experiment, c("ASGA_Michigan"="Michigan.us", "California_Invade"="California.I.us","California_Prop_Limi"="California.P.L.us","CCR_04"="CedarCreek4.us","CCR_093"="CedarCreek93.us","Germany_Montane"="Montane.de","Halle"="Halle.de","Jena"="Jena.de","Jena2"="Jena2.de","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow.us","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field.us","Texas_Temple_Prarie"="Texas.Templ.Prairie.us"))

l.pie_exp_coef3$Study<-as.factor(as.character(l.pie_exp_coef3$Study))

dc<-ggplot() + 
  geom_point(data = l.pie_exp_coef3, aes(x = Study, y = Slope,colour = Study),size = 4) +
  geom_errorbar(data = l.pie_exp_coef3, aes(x = Study,ymin = Slope_lower,
                                          ymax = Slope_upper,colour = Study),
                width = 0, size = 1.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_hline(data = filter(l.pie_fixef_df),
             aes(yintercept = Estimate[2]), size = 1.2) +
  geom_rect(data = filter(l.pie_fixef_df),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(x = 'Experiment',
       y = 'Slope',title= 'a) Evenness Coefficients') +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF","#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  scale_x_discrete(limits = rev(levels(l.pie_exp_coef3$Study)))+
  coord_flip() + 
  theme_bw()+theme(axis.text.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),legend.position="bottom")

# Figure S3
grid_arrange_shared_legend(dc,dr,nrow=1)


