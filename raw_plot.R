# #################################################################### ####
# Title: Comparing the effects of seeded diversity on species richness ####
#        and above ground biomass , and community assembly             #### 
# Authors: Emma R Ladouceur & Shane A. Blowes                           ####
# Details: Figure S2                                                   ####
# #################################################################### ####

rm(list=ls())
detach("package:ggplot2", unload=TRUE)
detach("package:plyr", unload=TRUE)
library(ggplot2)
library(ggforce)
library(grid)
library(gridExtra)

setwd('~/Desktop/Academic/R code/SeedAdditionSynthesis/')
plot<-read.csv("./Data/SeedAdd_Plot_Level.csv", header=TRUE) %>%
  as_tibble()


colnames(plot)

plot$Treatment<-plot$trt

# use grid_arrange_shared_legend function at the beginning of Main_Analysis.R
# to create figures

library(plyr)
plot$Study<-revalue(plot$Experiment_, c("ASGA_Michigan"="Michigan.us", "California_Invade"="California.I.us","California_Prop_Limi"="California.P.L.us","CCR_04"="CedarCreek4.us","CCR_093"="CedarCreek93.us","Germany_Montane"="Montane.de","Halle"="Halle.de","Jena"="Jena.de","Jena2"="Jena2.de","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow.us","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field.us","Texas_Temple_Prarie"="Texas.Templ.Prairie.us"))
plot$Study<-factor(as.character(plot$Study))


rp<-ggplot(plot,aes(x=Study, y=rich.plot)) +
geom_jitter( aes(shape = Treatment,color = Treatment),
             position = position_jitterdodge(jitter.width = 0.4,jitter.height = 0.4, dodge.width = 0.8), alpha=0.2,size=1.2) +
   geom_sina(aes(color = Treatment), size = 0.7 , alpha=0.2) +
  stat_summary(aes(shape=Treatment),
               fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom = "pointrange",  size = 0.4,
               position = position_dodge(0.8)) +
  scale_color_manual(values =  c("#A1C720FF","#15983DFF"))  + 
  theme_classic()+theme(axis.text.x = element_text(size=9, angle=7))+
  labs(title = "a) Plot Species Richness") + ylab("Species Richness") 


bp<-ggplot(plot,aes(x= Study, y= biomass.plot)) +
  geom_jitter( aes(shape = Treatment,color = Treatment),
               position = position_jitterdodge(jitter.width = 0.4,jitter.height = 0.4, dodge.width = 0.8), alpha=0.2,size=1.2) +
    geom_sina(aes(color = Treatment), size = 0.7 , alpha=0.2)+
  stat_summary(aes(shape= Treatment),
               fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom = "pointrange",  size = 0.4,
               position = position_dodge(0.8)) + ylim(0,1100) +
  scale_color_manual(values =  c("#A1C720FF","#15983DFF")) + 
  theme_classic()+theme(axis.text.x = element_text(size=9, angle=7)) +
  labs(title = 'b) Plot Biomass')  + ylab(expression(paste('Biomass (g/',m^2, ')'))) 

# Figure S2
grid_arrange_shared_legend(rp,bp,ncol=1,nrow=2)

