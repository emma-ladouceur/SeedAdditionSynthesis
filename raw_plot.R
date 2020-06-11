# #################################################################### ####
# Title: Reducing dispersal limitation via seed addition leads to      ####
#        increased species richness, but not aboveground biomass       #### 
# Authors: Emma R Ladouceur & Shane A. Blowes                          ####
# Details: Figures S1 : Raw data                                       ####
# #################################################################### ####


# Libraries
library(ggplot2)
library(ggforce)
library(patchwork)
library(tidyverse)

setwd('~/Data/')
plot<-read.csv("./Data/SeedAdd_Plot_Level.csv", header=TRUE) %>%
  as_tibble()


colnames(plot)

plot$Treatment<-plot$trt

plot$Experiment<-plyr::revalue(plot$Experiment_, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))
plot$Experiment<-factor(as.character(plot$Experiment))

 # Figure S1 a) Richness
rp<-ggplot(plot,aes(x=Experiment, y=rich.plot)) +
geom_jitter( aes(shape = Treatment,color = Treatment),
             position = position_jitterdodge(jitter.width = 0.4,jitter.height = 0.4, dodge.width = 0.8), alpha=0.4,size=1.2) +
  stat_summary(aes(shape=Treatment),
               fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom = "pointrange",  size = 0.4,
               position = position_dodge(0.8),color="black") +
  scale_color_manual(values =  c("#A1C720FF","#15983DFF"))  + 
  theme_classic()+theme(axis.text.x = element_text(size=9, angle=7), plot.margin=margin(t=4,1,1,1, "lines"),
                        legend.direction = "horizontal", legend.position = c(0.5,1.2) )+
  labs(title = "a) Plot Species Richness") + ylab("Species Richness") 


# Figure S1 b) Biomass
bp<-ggplot(plot,aes(x= Experiment, y= biomass.plot)) +
  geom_jitter( aes(shape = Treatment,color = Treatment),
               position = position_jitterdodge(jitter.width = 0.4,jitter.height = 0.4, dodge.width = 0.8), alpha=0.4,size=1.2) +
  stat_summary(aes(shape= Treatment),
               fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom = "pointrange",  size = 0.4,
               position = position_dodge(0.8)) + ylim(0,1100) +
  scale_color_manual(values =  c("#A1C720FF","#15983DFF")) + 
  theme_classic()+theme(axis.text.x = element_text(size=9, angle=7),legend.position="none") +

  labs(title = 'b) Plot Biomass')  + ylab(expression(paste('Biomass (g/',m^2, ')'))) 



(rp)/(bp)

