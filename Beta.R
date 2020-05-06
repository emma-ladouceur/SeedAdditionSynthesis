# #################################################################### ####
# Title: Comparing the effects of seeded diversity on species richness ####
#        and above ground biomass , and community assembly             #### 
# Authors: Emma R Ladouceur & Shane A. Blowes                          ####
# Details: Figure S4                                                   ####
# #################################################################### ####


rm(list=ls())
detach("package:ggplot2", unload=TRUE)
detach("package:plyr", unload=TRUE)


library(tidyr)
library(dplyr)
library(betapart)
library(tibble) 
library(tidyverse)
library(bayesplot)
library(patchwork)
library(ggplot2)



setwd('~/Desktop/Academic/R code/SeedAdditionSynthesis/')
plot<-read.csv("./Data/SeedAdd_Plot_Level.csv", header=TRUE) %>%
  as_tibble()

plot$unique.id<-plot$unique.id_

# select certain columns from plot
plot2 <- plot %>% 
  select(unique.id, Experiment_, site, block,yr.trt, trt, seed.rich) 


plot2 %>% 
  distinct(trt)

beta_pairs <- function(x){
  # function to return the dissimilarities (turnover and nestedness component)
  # for each control treatment comparison in x
  
  # return dataframe with dissimilarities and the treatment magnitude (seed.rich)

  # separate out the control and treatment plots
  contr.plots = x %>% 
    filter(trt=='Control')
  
  # fix for treatment labels
  trt.plots = x %>% 
    filter(trt=='Seeds' | trt=='Seed')
  
  out <- tibble()
  if(nrow(contr.plots)>0){
    for(i in 1:nrow(contr.plots)){
      beta = beta.pair(bind_rows(contr.plots %>% 
                                   slice(i) %>% 
                                   select(-seed.rich, -trt), 
                                 trt.plots %>% 
                                   select(-seed.rich, -trt)),
                       index.family = 'jaccard')
      # buid the data we want for analysis
      out <- bind_rows(out,
                       tibble(
                         seed.rich = trt.plots$seed.rich,
                         jtu = as.matrix(beta$beta.jtu)[-1,1],
                         jne = as.matrix(beta$beta.jne)[-1,1],
                         group = i)
      )
    }
  }  
  # escape for no controls
  else{
    out = tibble(
      seed.rich = NA,
      jtu = NA,
      jne = NA,
      group = NA)
  }
  return(out)
}


wide.df <- bind_rows(
  # 1) ASGA_Michigan.w
  left_join(ASGA_Michigan.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich),
  # 2) California_Invade.w
  left_join(California_Invade.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich),
  # 3) California_Prop_Limiti.w
  left_join(California_Prop_Limi.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich),
  # 4) CCR_04.w
  left_join(CCR_04.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich),
  # 5) CCR_93.w
  left_join(CCR_093.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich),
  # 6) Germany_Montane.w
  left_join(Germany_Montane.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich),
  # 7) Halle.w
  left_join(Halle.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich),
  # 8) Jena.w
  left_join(Jena.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich),
  # 9) Jena2.w
  left_join(Jena2.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich),
  # 10) Kansas_KUFS_LTER_Hay_Meadow_Exp_2.w
  left_join(Kansas_KUFS_LTER_Hay_Meadow_Exp_2.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich),
  # 11) Kansas_KUFS_LTER_Old_Field_Exp_1.w
  left_join(Kansas_KUFS_LTER_Old_Field_Exp_1.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich),
  # 12) Texas_Temple_Prarie.w
  left_join(Texas_Temple_Prarie.w, plot2, 
            by = 'unique.id') %>% 
    group_by(Experiment_, site, block, yr.trt) %>% 
    nest_legacy(starts_with('sp_'), trt, seed.rich)
)


# calculate the beta components
wide.df <- wide.df %>% 
  mutate(beta = purrr::map(data, ~ beta_pairs(.x)))


beta.df = wide.df %>% 
  unnest_legacy(beta) %>%
  unite(col = pw_beta_group,
        c(Experiment_, site, block, yr.trt, group), sep = '_', remove = F) %>% 
  select(-group)


# write.csv(beta.df,"./Data/beta.df.csv")
setwd('~/Dropbox/Projects/SeedAddDraft/')
beta<-read.csv("./Data/beta.df.csv", header=TRUE) %>%
  as_tibble()

# sb
beta<-read.csv("~/Dropbox/SeedAdd/Data/beta.df.csv", header=TRUE) %>%
  as_tibble()

beta$Experiment<-beta$Experiment_
beta$fyr.trt<-as.factor(beta$yr.trt)
beta$seed.rich<-as.numeric(as.character(beta$seed.rich))
beta$site<-as.factor(beta$site)
beta$block<-as.factor(beta$block)
beta<-beta %>%
  drop_na() 

nrow(beta)
nrow(beta2)
View(beta2)


# Load model objects
load("./Model Fits/betat.Rdata") # object name: turnover.zoib
load("~/Dropbox/SeedAdd/Model_fits/betat.Rdata") # object name: turnover.zoib
load("./Model Fits/betan.Rdata") # object name: nested.zib
load("~/Dropbox/SeedAdd/Model_fits/betan.Rdata") # object name: nested.zib

# turnover.zoib <- brm(jtu ~  seed.rich + (seed.rich | Experiment/site/block/fyr.trt),
#                          family=zero_one_inflated_beta(),
#                          data = beta,
#                          inits = '0',
#                          cores = 4, chains = 4)
# 
# setwd('~/Dropbox/Projects/SeedAdd/Model_fits/')
# save(turnover.zoib, file = './betat.Rdata')

# Turnover model
summary(turnover.zoib)

plot(turnover.zoib) 

# Figure S1 d
pp_check(turnover.zoib) + theme_classic()

betat_fitted <- cbind(turnover.zoib$data,
                      fitted(turnover.zoib, re_formula = NA)) %>% 
  as_tibble() 


betat_fixef <- fixef(turnover.zoib)

betat_exp_coef <- coef(turnover.zoib)

betat_exp_coef2 <-  bind_cols(betat_exp_coef$Experiment[,,'Intercept'] %>% 
                                as_tibble() %>% 
                                mutate(Intercept = Estimate,
                                       Intercept_lower = Q2.5,
                                       Intercept_upper = Q97.5,
                                       Experiment = rownames(betat_exp_coef$Experiment[,,'Intercept'])) %>% 
                                select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                              betat_exp_coef$Experiment[,,'seed.rich'] %>% 
                                as_tibble() %>% 
                                mutate(Slope = Estimate,
                                       Slope_lower = Q2.5,
                                       Slope_upper = Q97.5) %>% 
                                select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  inner_join(beta %>% 
               group_by(Experiment) %>% 
               summarise(xmin = min(seed.rich),
                         xmax = max(seed.rich)),
             by = 'Experiment')


# Nestedness model

# nested.zib <- brm(jne ~  seed.rich + (seed.rich | Experiment/site/block/fyr.trt),
#                       family=zero_inflated_beta(),
#                       data = beta,
#                       inits = '0',
#                       cores = 4, chains = 4)
# 
# setwd('~/Dropbox/Projects/SeedAdd/Model_fits/')
# save(nested.zib, file = './betan.Rdata')

summary(nested.zib)
plot(nested.zib) 

# Figure S1 e
pp_check(nested.zib)+ theme_classic()

# residuals
n1<-residuals(nested.zib)
n1<-as.data.frame(n1)

rn.plot<-cbind(beta,n1$Estimate)

par(mfrow=c(2,3))
with(rn.plot, plot(Experiment, n1$Estimate))
with(rn.plot, plot(site, n1$Estimate))
with(rn.plot, plot(block, n1$Estimate))
with(rn.plot, plot(fyr.trt, n1$Estimate))
with(rn.plot, plot(pw_beta_group, n1$Estimate))


# make sure to detach plyr at top

# fixed effects
betan_fitted <- cbind(nested.zib$data,
                      fitted(nested.zib, re_formula = NA)) %>% 
  as_tibble() 

betan_fixef <- fixef(nested.zib)

betan_exp_coef <- coef(nested.zib)
betan_exp_coef 

betad<-beta%>%distinct(Experiment,seed.rich)
View(betad)
# this gets us the coefficients for the varying intercepts and slopes
betan_exp_coef2 <-  bind_cols(betan_exp_coef$Experiment[,,'Intercept'] %>% 
                                as_tibble() %>% 
                                mutate(Intercept = Estimate,
                                       Intercept_lower = Q2.5,
                                       Intercept_upper = Q97.5,
                                       Experiment = rownames(betan_exp_coef$Experiment[,,'Intercept'])) %>% 
                                select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                              betan_exp_coef$Experiment[,,'seed.rich'] %>% 
                                as_tibble() %>% 
                                mutate(Slope = Estimate,
                                       Slope_lower = Q2.5,
                                       Slope_upper = Q97.5) %>% 
                                select(-Estimate, -Est.Error, -Q2.5, -Q97.5)) %>% 
  inner_join(beta %>% 
               group_by(Experiment) %>% 
               summarise(xmin = min(seed.rich),
                         xmax = max(seed.rich)),
             by = 'Experiment')


# use grid_arrange_shared_legend function at the beginning of Main_Analysis.R
# to create figures


# beta plots

colnames(betat_exp_coef2)
colnames(betan_exp_coef2)
betat_exp_coef2$Model<-'Turnover'
betan_exp_coef2$Model<-'Nestedness'

#fixed
betat_fixef_df<-as.data.frame(betat_fixef)
betan_fixef_df<-as.data.frame((betan_fixef))
betat_fixef_df$Model<-'Turnover'
betan_fixef_df$Model<-'Nestedness'
fixef.all<-bind_rows(betat_fixef_df,betan_fixef_df)
fixef.all

library(plyr)
betat_exp_coef2$Experiment<-revalue(betat_exp_coef2$Experiment, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))
betan_exp_coef2$Experiment<-revalue(betan_exp_coef2$Experiment, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))

betat_exp_coef2$Experiment<-factor(as.character(betat_exp_coef2$Experiment))
betan_exp_coef2$Experiment<-factor(as.character(betan_exp_coef2$Experiment))

betat_fitted$Experiment<-revalue(betat_fitted$Experiment, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))
betat_fitted$Experiment<-factor(as.character(betat_fitted$Experiment))
betan_fitted$Experiment<-revalue(betan_fitted$Experiment, c("ASGA_Michigan"="Michigan", "California_Invade"="California.1","California_Prop_Limi"="California.2","CCR_04"="Cedar.Creek.4","CCR_093"="Cedar.Creek.93","Germany_Montane"="Montane","Halle"="Halle","Jena"="Jena","Jena2"="Jena.2","Kansas_KUFS_LTER_Hay_Meadow_Exp_2"="Kansas.Hay.Meadow","Kansas_KUFS_LTER_Old_Field_Exp_1"="Kansas.Old.Field","Texas_Temple_Prarie"="Texas.Temple.Prairie"))
betan_fitted$Experiment<-factor(as.character(betan_fitted$Experiment))


tc<-ggplot() + 
  geom_point(data = betat_exp_coef2, aes(x = Experiment, y = Slope,colour = Experiment),size = 2) +
  geom_errorbar(data = betat_exp_coef2, aes(x = Experiment,ymin = Slope_lower,
                                            ymax = Slope_upper,colour = Experiment),
                width = 0, size = 1) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_hline(data = filter(fixef.all, Model=='Turnover'),
             aes(yintercept = Estimate[2]), size = 1.2) +
  geom_rect(data = filter(fixef.all, Model=='Turnover'),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(x = 'Experiment',
       y = 'Change in Turnover / species of seed added', subtitle= "b) ") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  scale_x_discrete(limits = rev(levels(betat_exp_coef2$Experiment)))+
  coord_flip() + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),legend.position="bottom",axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())


nc<-ggplot() + 
  geom_point(data = betan_exp_coef2, aes(x = Experiment, y = Slope,colour = Experiment),size = 2) +
  geom_errorbar(data = betan_exp_coef2, aes(x = Experiment,ymin = Slope_lower,
                                            ymax = Slope_upper,colour = Experiment),
                width = 0, size = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_hline(data = filter(fixef.all, Model=='Nestedness'),
             aes(yintercept = Estimate[2]), size = 1.2) +
  geom_rect(data = filter(fixef.all, Model=='Nestedness'),
            aes(xmin = -Inf, xmax = Inf,
                ymin = Q2.5[2], ymax = Q97.5[2]),
            alpha = 0.3) +
  labs(x = 'Experiment',
       y = 'Change in Nestedness / species of seed added', subtitle = "d) ") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF", "#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  scale_x_discrete(limits = rev(levels(betan_exp_coef2$Experiment)))+coord_flip() + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),legend.position="bottom",#axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())

nc
colnames(betat_fitted)

betat_exp_coef2
betat.line<-betat_exp_coef2 %>% filter(xmin < xmax)
betat.point<-betat_exp_coef2 %>% filter(xmin==xmax)

betan.line<-betan_exp_coef2 %>% filter(xmin < xmax)
betan.point<-betan_exp_coef2 %>% filter(xmin==xmax)

#"#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF","#F9B90AFF" , "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF"

color_line <- c("#EC579AFF", "#149BEDFF","#8F2F8BFF")
 color_point <- c("#FA6B09FF","#EE0011FF" , "#15983DFF","#A1C720FF","#0C5BB0FF","#F9B90AFF", "#16A08CFF" ,"#6A7F93FF","#9A703EFF" )


betat.point

btr <-
ggplot() +
  geom_point(data = betat_fitted,
             aes(x = seed.rich, y = jtu,
                 colour = Experiment), #alpha=0.4,
             size = 1, position = position_jitter(width = 0.1)) +
  geom_segment(data = betat.line,
               aes(x = xmin,
                   xend = xmax,
                   y = plogis(Intercept + Slope * xmin),
                   yend = plogis(Intercept + Slope * xmax),
                   colour = Experiment),
                   #group = Experiment,
                    colour = color_line,
               size = 1.2) +
  geom_point(data = betat.point,
             aes(x = xmax, y = plogis(Intercept + Slope)), 
             fill=color_point,shape=21, size=3.5,stroke=1,
             color="black") +
  geom_ribbon(data = betat_fitted,
              aes(x = seed.rich, ymin = Q2.5, ymax = Q97.5),
              alpha = 0.3) +
  geom_line(data = betat_fitted,
            aes(x = seed.rich, y = Estimate),
            size = 1.5) +
  labs(x = 'Number of species of seed added',
       y = 'Turnover', title= 'Turnover', subtitle="a)") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF","#F9B90AFF" , "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),
                     legend.position="bottom")
btr

summary(turnover.zoib)
summary(nested.zib)

View(betat_exp_coef2)

bnr <- ggplot() +
  geom_point(data = betan_fitted,
             aes(x = seed.rich, y = jne,
                 colour = Experiment),
             size = 1.2, position = position_jitter(width = 0.1)) +
  geom_segment(data = betan.line,
               aes(x = xmin,
                   xend = xmax,
                   y = plogis(Intercept + Slope * xmin),
                   yend = plogis(Intercept + Slope * xmax)),
               #group = Experiment,
               colour = color_line,
               size = 1.2) +
  geom_point(data = betan.point,
             aes(x = xmax, y = plogis(Intercept + Slope)), 
             fill=color_point,shape=21, size=3.5,stroke=1,
             color="black") +
  geom_ribbon(data = betan_fitted,
              aes(x = seed.rich, ymin = Q2.5, ymax = Q97.5),
              alpha = 0.3) +
  geom_line(data = betan_fitted,
            aes(x = seed.rich, y = Estimate),
            size = 1.5) +
  labs(x = 'Number of species of seed added',
       y = 'Nestedness', title= ' Nestedness', subtitle="c)") +
  scale_colour_manual(values = c( "#EE0011FF" , "#EC579AFF", "#15983DFF", "#149BEDFF", "#0C5BB0FF", "#8F2F8BFF","#F9B90AFF" , "#16A08CFF" ,"#6A7F93FF","#FA6B09FF","#A1C720FF","#9A703EFF" ))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),
                     legend.position="bottom")
bnr



# extract legend
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
# make legend an object
b.legend<-g_legend(nc)

# use patchwork to arrange figures with single legend
( btr+ theme(legend.position="none") | tc+ theme(legend.position="none")  ) / (bnr+ theme(legend.position="none")| nc+ theme(legend.position="none") ) /(b.legend)  +
  plot_layout(heights = c(10,10,2.3))



