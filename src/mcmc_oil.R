library(tidyverse);library (rstatix);library (stringr);library(viridis)
library(ggpattern); library (vegan) ;library (ggrepel)
library(lme4); library(glmmTMB); library (DHARMa) 

# Data preparation ####
read.csv("data/species_oil.csv", sep =",")%>%
  dplyr::select(Taxon, species, family, community, ecology, PERoil, seedmass, ratio, FDD, GDD, under_snow)%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, species, community, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(seedmass = log(seedmass), 
         PERoil = log(PERoil))%>%
  na.omit()-> oil_df

hist(oil_df$seedmass) # normal distribution after log transformation
hist(oil_df$PERoil) # almost normal distribution after log transformation
hist(oil_df$ratio) # normal distributed
hist(oil_df$GDD)
hist(oil_df$under_snow) # not really normally distributed
# GLMM with dharma ####
# escoger entre family = Gamma(link="log"),  gaussian
# ratio with gaussian family
# PERoil (log transform) with Gamma family
# GDD gamma family
a <- glmmTMB(PERoil ~ under_snow + (1| familia) ,  family = Gamma(link="log"), data= oil_df) 
#b <- glmmTMB(ratio  ~ seedmass + (1| familia) , family = gaussian, data= glm) 
summary(a)
Anova(a)
residuals <- simulateResiduals (a) ; plot(residuals)
# to compare models
AIC(a, b)

# simple lm
lm(PERoil~under_snow, data=oil_df)-> c
summary(c)
# take into account phylogeny! 
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree.tree")), 
                    ape::read.tree("results/tree.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(PERoil ~ under_snow, # PERoil
                   random = ~ animal + ID + community,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_df,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g3
x11()
plot(g3)
summary(g3) 
