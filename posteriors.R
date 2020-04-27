# #################################################################### ####
# Title: Comparing the effects of seeded diversity on species richness ####
#        and above ground biomass , and community assembly             #### 
# Authors: Emma R Ladouceur & Shane A. Blowes                          ####
# Details: Supplementary Information Section 4                         ####
# #################################################################### ####


library(tidyverse)
library(brms)
library(ggridges)
library(gridExtra)
library(grid)


#setwd('~/Desktop/Academic/R code/SeedAdditionSynthesis/')
setwd('~/Dropbox/Projects/SeedAddDraft/')
load("./Model Fits/rich.Rdata") # object name: seedadd.rich
load("./Model Fits/biomass.Rdata") # object name: seedadd.biomass



# study-levels (use model with fewest missing values)
study_levels <- seedadd.rich$data %>% 
  as_tibble() %>% 
  distinct(Experiment) %>% 
  mutate(level1 =  Experiment,
         level = gsub(' ', '.', level1)) %>%
  nest_legacy(level)



study_sample_posterior <- study_levels %>%
  mutate(rich = purrr::map(data, ~posterior_samples(seedadd.rich, 
                                                    pars = paste('r_Experiment[', as.character(.x$level), ',seed.rich.m]', sep=''),
                                                    exact = TRUE,
                                                    subset = floor(runif(n = 1000,
                                                                         min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         rich.overall = purrr::map(data, ~posterior_samples(seedadd.rich, 
                                                            pars = 'b_seed.rich.m',
                                                            exact = TRUE,
                                                            subset = floor(runif(n = 1000,
                                                                                 min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         bm = purrr::map(data, ~posterior_samples(seedadd.biomass, 
                                                  pars = paste('r_Experiment[', as.character(.x$level), ',seed.rich.m]', sep=''),
                                                  exact = TRUE,
                                                  subset = floor(runif(n = 1000,
                                                                       min = 1, max = 2000))) %>% unlist() %>% as.numeric()),
         bm.overall = purrr::map(data, ~posterior_samples(seedadd.biomass, 
                                                          pars = 'b_seed.rich.m',
                                                          exact = TRUE,
                                                          subset = floor(runif(n = 1000,
                                                                               min = 1, max = 2000))) %>% unlist() %>% as.numeric()))



rich_fixef <- fixef(seedadd.rich)
bm_fixef <- fixef(seedadd.biomass)


rich_posterior <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest_legacy(rich,rich.overall) %>% 
  mutate(response = 'rich.plot',
         rich_global_slope = rich.overall) 

bm_posterior <- study_sample_posterior  %>% 
  select(-data) %>% 
  unnest_legacy(bm,bm.overall) %>% 
  mutate(response = 'l.biomass',
         bm_global_slope = bm.overall) 


bm_posterior$c.v.b<-bm_posterior$bm+bm_posterior$bm_global_slope
rich_posterior$c.v.r<-rich_posterior$rich+rich_posterior$rich_global_slope


p_corr <- cor.test(x=bm_posterior$c.v.b, y=rich_posterior$c.v.r, method = 'spearman')
p_corr

# Spearman's rank correlation rho
# 
# data:  bm_posterior$c.v.b and rich_posterior$c.v.r
# S = 2.3412e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1870921 


