library(tidyverse);library (rstatix);library (stringr);library(viridis)
library(ggpattern); library (vegan) ;library (ggrepel)
library(lme4); library(glmmTMB); library (DHARMa) 

############################################ BIOLOGICAL TRADE-OFFS #################################################
# seed mass (log transformed) #####
oil_data%>%
  merge(seedmass)%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, species, community, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Lseedmass = log(mean), 
         LPERoil = log(PERoil),
         Lratio=log(ratio))%>%
  na.omit()-> oil_seedmass # 
# data distribution
hist(oil_seedmass$seedmass) 
hist(oil_seedmass$Lseedmass)# normal distribution after log transformation
hist(oil_seedmass$PERoil)
hist(oil_seedmass$LPERoil) # almost normal distribution after log transformation
hist(oil_seedmass$ratio) # almost normal distributed
hist(oil_seedmass$Lratio)# almost normal distributed

##### GLMM with dharma 
a <- glmmTMB(LPERoil~ Lseedmass + (1| familia) ,  family = Gamma(link="log"), data= oil_seedmass) 
a <- glmmTMB(Lratio  ~ Lseedmass + (1| familia) , family = Gamma(link="log"), data= oil_seedmass) 
summary(a)
residuals <- simulateResiduals (a) ; plot(residuals)
# to compare modelsAIC(a, b)
# Gamma(link = "log"), model plots GOOD in oil content still not good with ratio
# oil content (log transformed) and ratio 

##### MCMC 
# take into account phylogeny! 
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_oil.tree")), 
                    ape::read.tree("results/tree_oil.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(Lratio ~ Lseedmass, # ratio 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_seedmass,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g3
x11()
plot(g3)
summary(g3) 

# longevity (raw germination curves and p50) in 2024 longevity analysis script more options####
# right now 23 species with longevity data from Pavia and oil data content 
### Read data raw germination for MCMC#
read.csv("data/2022/germination22.csv", sep =";") %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/2022/species22.csv")) %>% 
  merge(oil_data, by=c("Taxon", "species", "community", "family", "ecology"))%>% # from data handling script 
  rename(familia = family)%>%
  convert_as_factor(Taxon, species, community, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID) %>%
  mutate(LPERoil = log(PERoil),
         Lratio=log(ratio))%>%
  na.omit()-> df22_new # remove species without sp_pref (Minuartia CF) and without oil content (C.ramosissimum, G.verna, G.campestris and P. pyrenaica)

unique(df22_new$familia)
### Read tree

phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree22.tree")), 
                    ape::read.tree("results/tree22.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

# shorter iterations
#nite = 10000
#nthi = 10
#nbur = 100

### Set priors for germination models

priors <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        #G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        #G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

### All species model

MCMCglmm::MCMCglmm(cbind(germinated, seeds - germinated) ~
                     scale(ageing) * scale(Lratio), #scale(ageing) * scale(LPERoil) + scale(ageing)*scale(ratio)
                   random = ~ animal + ID ,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = df22_new,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1

# save(m1, file = "results/mcmc.Rdata")
x11()
plot(m1)

# load("results/mcmc.Rdata")
summary(m1)

### Random and phylo

# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"]) 

mean(m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"])) %>% round(2)
coda::HPDinterval(m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"]))[, 2] %>% round(2)

# Random effects animal
summary(m1)$Gcovariances[1, 1] %>% round(2) 
summary(m1)$Gcovariances[1, 2] %>% round(2) 
summary(m1)$Gcovariances[1, 3] %>% round(2)

# Random effects species ID
summary(m1)$Gcovariances[2, 1] %>% round(2)
summary(m1)$Gcovariances[2, 2] %>% round(2) 
summary(m1)$Gcovariances[2, 3] %>% round(2) 

# Random effects community
summary(m1)$Gcovariances[3, 1] %>% round(2)
summary(m1)$Gcovariances[3, 2] %>% round(2) 
summary(m1)$Gcovariances[3, 3] %>% round(2)


#### P50 MCMC
# read data
oil_data %>%
  merge(p50)%>% #23 species
  rename(familia = family)%>%
  convert_as_factor(Taxon, species, community, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID) %>%
  mutate(LPERoil = log(PERoil),
         Lratio=log(ratio),
         Lp50 = log(p50))%>%
  na.omit()-> genstat22_new
hist(genstat22_new$p50) # normally distributed??
hist(genstat22_new$Lp50)
hist(genstat22_new$PERoil) #not normally distributed
hist(genstat22_new$LPERoil)# almost normally distributed
hist(genstat22_new$ratio) # normally distributed
hist(genstat22_new$Lratio) # normally distributed

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(p50 ~  scale(Lratio),
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = genstat22_new,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g3
x11()
plot(g3)
summary(g3)

#### P50 GLMM
# Gamma(link="log")
a <- glmmTMB(p50 ~ PERoil + (1| familia) ,  family = gaussian, data= genstat22_new) #LPERoil
a <- glmmTMB(p50  ~ Lratio + (1| familia) , family = gaussian, data= genstat22_new)  #Lratio
summary(a)
residuals <- simulateResiduals (a) ; plot(residuals)

# T50 ()####
oil_data%>%
  merge(t50)%>% # 29 species
  rename(familia = family)%>%
  convert_as_factor(Taxon, species, community, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(LPERoil = log(PERoil),
         Lratio=log(ratio),
         LT50 = log(T50))%>%
  na.omit()-> oil_t50 # 
hist(oil_t50$T50) # not normally dsitributed
hist(oil_t50$LT50)
##GLMM with dharma 
a <- glmmTMB(T50 ~ LPERoil + (1| familia) ,  family = Gamma(link="log"), data= oil_t50) 
a <- glmmTMB(T50 ~ Lratio + (1| familia) , family = Gamma(link="log"), data= oil_t50) 
summary(a)
residuals <- simulateResiduals (a) ; plot(residuals)
# to compare modelsAIC(a, b)
##### MCMC 
# take into account phylogeny! 
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_oil.tree")), 
                    ape::read.tree("results/tree_oil.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(T50 ~ Lratio , # ratio 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_t50,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g3
x11()
plot(g3)
summary(g3) 

################ ECOLOGICAL TRADE_OFFS #######################
# Ecology/distribution  ####
str(oil_data)
oil_data%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, species, community, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(LPERoil = log(PERoil),
         Lratio=log(ratio))%>%
  na.omit()-> oil_ecology

# Gamma(link="log")
a <- glmmTMB(LPERoil ~ ecology + (1| familia) ,  family = Gamma(link="log"), data= oil_ecology) 
a <- glmmTMB(Lratio  ~ ecology  + (1| familia) , family = Gamma(link="log"), data= oil_ecology) 
summary(a)
residuals <- simulateResiduals (a) ; plot(residuals)
# to compare modelsAIC(a, b)

##### MCMC 
# take into account phylogeny! 
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_oil.tree")), 
                    ape::read.tree("results/tree_oil.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(LPERoil ~ ecology, # ratio 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_ecology,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g3
x11()
plot(g3)
summary(g3) 
# GDD ####
oil_data%>%
  merge(sp_pref)%>% # 29 species
  rename(familia = family)%>%
  convert_as_factor(Taxon, species, community, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(LPERoil = log(PERoil),
         Lratio=log(ratio))%>%
  na.omit()-> oil_sp_pref #

str(oil_df)
# Gamma(link="log")
a <- glmmTMB(LPERoil ~ GDD*community + (1| familia) ,  family = gaussian, data= oil_sp_pref) 
a <- glmmTMB(Lratio  ~ GDD*community  + (1| familia) , family = Gamma(link="log"), data= oil_sp_pref) 
summary(a)
residuals <- simulateResiduals (a) ; plot(residuals)
# to compare modelsAIC(a, b)

##### MCMC 
# take into account phylogeny! 
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_oil.tree")), 
                    ape::read.tree("results/tree_oil.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(Lratio ~ GDD*community, # ratio 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_sp_pref,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g3
x11()
plot(g3)
summary(g3) 
# FDD ####
str(oil_sp_pref)
# Gamma(link="log")
a <- glmmTMB(LPERoil ~ FDD*community + (1| familia) ,  family = Gamma(link="log"), data= oil_sp_pref) 
a <- glmmTMB(Lratio  ~ FDD*community  + (1| familia) , family = Gamma(link="log"), data= oil_sp_pref) 
summary(a)
residuals <- simulateResiduals (a) ; plot(residuals)
# to compare modelsAIC(a, b)

##### MCMC 
# take into account phylogeny! 
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_oil.tree")), 
                    ape::read.tree("results/tree_oil.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(LPERoil ~ FDD*community, # ratio 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_sp_pref,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g3
x11()
plot(g3)
summary(g3) 
