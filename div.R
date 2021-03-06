# #################################################################### ####
# Title: Reducing dispersal limitation via seed addition leads to      ####
#        increased species richness, but not aboveground biomass       #### 
# Authors: Emma R Ladouceur & Shane A. Blowes                          ####
# Details: Figures S3 : SPIE                                           ####
# #################################################################### ####


# Libraries
library(tidyverse)
library(tibble)
library(data.table)
library(ggplot2)
library(brms)
library(vegan)
library(bayesplot)
library(patchwork)


# Data

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

summary(seed.pie)

#write.csv(seed.pie,"~/Data/seed.pie.csv")


# Diversity Models

seed.pie.d<-read.csv("./Data/seed.pie.csv", header=TRUE) %>%
  as_tibble()


levels(seed.pie.d$Experiment)
seed.pie.d$Experiment<-seed.pie.d$Experiment_
seed.pie.d$fyr.trt<-as.factor(seed.pie.d$yr.trt)
seed.pie.d$seed.rich<-as.numeric(as.character(seed.pie.d$seed.rich))
seed.pie.d$site<-as.factor(seed.pie.d$site)
seed.pie.d$block<-as.factor(seed.pie.d$block)
seed.pie.d$seed.rich.m<-seed.pie.d$seed.rich-mean(seed.pie.d$seed.rich)
seed.pie.d$l.b.pie <- log(seed.pie.d$biomass.pie)

hist(seed.pie.d$biomass.pie)
hist(seed.pie.d$l.b.pie)

# load model object
load("./Model Fits/seed.pie.Rdata") # object name: m.seed.pie

# seed.pie <- brm(biomass.pie ~  seed.rich.m + (seed.rich.m | Experiment/site/block/fyr.trt),
#                      data = seed.pie, cores = 4, chains = 4)
# 
# save(seed.pie, file = './Model Fits/seed.pie.Rdata')

summary(seed.pie)

plot(m.seed.pie)

# Figure S1 c
color_scheme_set("darkgray")
pp_check(seed.pie) + theme_classic()

# residuals from model
m2<-residuals(m.seed.pie)
m2<-as.data.frame(m2)
s.p.<-seed.pie %>% drop_na(biomass.pie)

rb.plot<-cbind(s.p.,m2$Estimate)
head(rb.plot)
par(mfrow=c(3,2))
with(rb.plot, plot(Experiment, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(site, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(block, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(fyr.trt, m2$Estimate));abline(h=0, lty=2)
with(rb.plot, plot(seed.rich, m2$Estimate));abline(h=0, lty=2)

# make sure to detach plyr at top

pie_fitted <- cbind(seed.pie$data,
                    fitted(seed.pie, re_formula = NA))  %>%
  as_tibble() %>%
   inner_join(seed.pie.d %>% distinct(Experiment, site,block, fyr.trt, seed.rich, seed.rich.m, biomass.pie),
            by = c('Experiment','site','block','fyr.trt', 'seed.rich.m','biomass.pie'))


pie_fixef <- fixef(seed.pie)

pie_exp_coef <- coef(seed.pie)

pie_exp_coef2 <-  bind_cols(pie_exp_coef$Experiment[,,'Intercept'] %>% 
                              as_tibble() %>% 
                              mutate(Intercept = Estimate,
                                     Intercept_lower = Q2.5,
                                     Intercept_upper = Q97.5,
                                     Experiment = rownames(pie_exp_coef$Experiment[,,'Intercept'])) %>% 
                              select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                            pie_exp_coef$Experiment[,,'seed.rich.m'] %>% 
                              as_tibble() %>% 
                              mutate(Slope = Estimate,
                                     Slope_lower = Q2.5,
                                     Slope_upper = Q97.5) %>% 
                              select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) 

pie_exp_coef2$Experiment<-as.factor(as.character(pie_exp_coef2$Experiment))
seed.pie$seed.rich<-as.numeric(as.character(seed.pie$seed.rich))


pie_exp_coef3 <- pie_exp_coef2 %>%  
  left_join(seed.pie.d %>% 
               group_by(Experiment) %>% 
               summarise(xmin = min(seed.rich),
                         xmax = max(seed.rich),
                         cxmin = min(seed.rich.m),
                         cxmax = max(seed.rich.m)),
             by = 'Experiment')



library(plyr)
pie_fitted$Experiment<-as.factor(as.character(pie_fitted$Experiment))
pie_fitted$Experiment<-plyr::revalue(pie_fitted$Experiment, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))
pie_exp_coef3$Experiment<-plyr::revalue(pie_exp_coef3$Experiment, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))

pie_fitted$Experiment<-factor(as.character(pie_fitted$Experiment))
pie_exp_coef3$Experiment<-factor(as.character(pie_exp_coef3$Experiment))



# Figure S3 a) : Evenness Regression
dr <- ggplot() +
  geom_point(data = pie_fitted,
             aes(x = seed.rich, y = biomass.pie,
                 colour = Experiment),
             size = 1.2, position = position_jitter(width = 0.95, height = 0.95)) +
  geom_segment(data = pie_exp_coef3,
               aes(x = xmin,
                   xend = xmax,
                   y = (Intercept + Slope * cxmin),
                   yend = (Intercept + Slope * cxmax),
                   group = Experiment,
                   colour = Experiment),
               size = 1.2) +
  geom_ribbon(data = pie_fitted,
              aes(x = seed.rich, ymin = Q2.5, ymax = Q97.5),
              alpha = 0.3) +
  geom_line(data = pie_fitted,
            aes(x = seed.rich, y = Estimate),
            size = 1.5) +
  labs(x = 'Number of species of seed added',
       y = expression(S[PIE]), title= ' Evenness', subtitle="a)") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF","#F9B90AFF" , "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),
                                                                                                                             legend.position="bottom")

dr

# coefs 

#fixed
pie_fixef_df<-as.data.frame(pie_fixef)
pie_fixef_df

pie_exp_coef3$Experiment<-plyr::revalue(pie_exp_coef3$Experiment, c("ASGA_Michigan"="Michigan.us", "California_Invade"="California.I.us","California_Prop_Limi"="California.P.L.us","CCR_04"="CedarCreek4.us","CCR_093"="CedarCreek93.us","Germany_Montane"="Montane.de","Halle"="Halle.de","Jena"="Jena.de","Jena2"="Jena2.de","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow.us","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field.us","Texas_Temple_Prarie"="Texas.Temple.Prairie.us"))

pie_exp_coef3$Experiment<-as.factor(as.character(pie_exp_coef3$Experiment))

# Figure S3 b) : Evenness coefficients
dc<-ggplot() + 
  geom_point(data = pie_exp_coef3, aes(x = Experiment, y = Slope,colour = Experiment),size = 2) +
  geom_errorbar(data = pie_exp_coef3, aes(x = Experiment,ymin = Slope_lower,
                                          ymax = Slope_upper,colour = Experiment),
                width = 0, size = 1) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_hline(data = filter(pie_fixef_df),
             aes(yintercept = Estimate[2]), size = 1.2) +
  geom_rect(data = filter(pie_fixef_df),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(x = 'Experiment',
       y=expression(S[PIE] * "/ number of species of seed added") , title= ' ', subtitle= "b)") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF","#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  scale_x_discrete(limits = rev(levels(pie_exp_coef3$Experiment)))+
  coord_flip() + 
  theme_bw()+theme(axis.text.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),legend.position="bottom",
                   axis.ticks.y=element_blank())

dc


# extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

d.legend<-g_legend(dr)

# use patchwork for legend
(dr+ theme(legend.position="none") | dc+ theme(legend.position="none"))/(d.legend)  +
  plot_layout(heights = c(10,1))





