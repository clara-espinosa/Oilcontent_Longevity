library(tidyverse);library (rstatix);library (stringr);library(viridis)
library(ggpattern); library (vegan) ;library (ggrepel)
library(lme4); library(glmmTMB); library (DHARMa) 

read.csv("data/species_traits_summary.csv")%>%# from script header data handling
  
############################################ BIOLOGICAL TRADE-OFFS #################################################
# seed mass (log transformed) #####
read.csv("data/species_traits_summary.csv")%>%# from script header data handling
  dplyr::select(Taxon, community, family, oil.content, ratio, mass_50)%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Lseedmass = log(mass_50), 
         Loil.content = log(oil.content),
         Lratio=log(ratio))%>%
  as.data.frame()-> oil_seedmass # 47 species 49 accessions

# data distribution
hist(oil_seedmass$mass_50) 
hist(oil_seedmass$Lseedmass)# normal distribution after log transformation
hist(oil_seedmass$oil.content)
hist(oil_seedmass$Loil.content) # almost normal distribution after log transformation
hist(oil_seedmass$ratio) # almost  normal distributed
hist(oil_seedmass$Lratio)# normal distributed

### MCMC 
# take into account phylogeny! 
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_oil.tree")), 
                    ape::read.tree("results/tree_oil.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

plot(ape::read.tree("results/tree_oil.tree"))
### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))

# mcmc glmm
MCMCglmm::MCMCglmm(Lratio  ~ Lseedmass, #  Loil.content   Lratio 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_seedmass,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g4
g3 # Loil.content ~ Lseedmass
g4 # Lratio ~ Lseedmass

# None show a significant relationship although negative trends

x11()
plot(g4)   
summary(g4)  

# longevity (raw germination curves and p50 (square root transformed)) in 2024 longevity analysis script more options####
read.csv("data/longevity/species.csv", sep =",")-> longevity_sp
setdiff(oil_data$Taxon,longevity_sp$Taxon)
setdiff(longevity_sp$Taxon,oil_data$Taxon)
### Read and build data raw germination for MCMC#
read.csv("data/longevity/germination.csv", sep =",") %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/longevity/species.csv", sep =",")) %>% 
  merge(oil_data, by=c("Taxon", "community", "family", "ecology"))%>% # from data handling script 
  group_by(Taxon, family,ageing) %>%
  summarise(seeds = sum (seeds, na.rm = T),germinated = sum(germinated, na.rm = TRUE), 
            oil.content = mean(oil.content), ratio = mean(ratio)) %>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID) %>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio))%>%
  na.omit()%>%
  as.data.frame()-> df24 

unique(df24$familia) # 16 families
unique(df24$Taxon) #31 species

# data distribution
hist(df24$oil.content)
hist(df24$Loil.content) # almost normal distribution after log transformation
hist(df24$ratio) # almost  normal distributed
hist(df24$Lratio)# normal distributed


### Read tree

phangorn::nnls.tree(cophenetic(ape::read.tree("results/treelongevity.tree")), 
                    ape::read.tree("results/treelongevity.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL
plot(ape::read.tree("results/treelongevity.tree"))
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
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

### All species model

MCMCglmm::MCMCglmm(cbind(germinated, seeds - germinated) ~
                     scale(ageing) * scale(Lratio), #scale(ageing) * scale(Loil.content) + scale(ageing)*scale(Lratio)
                   random = ~ animal + ID ,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = df24,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m2

# save(m1, file = "results/mcmc.Rdata")
m1 # germ ~ ageing * Loil.content
m2 # germ ~ ageing * Lratio

x11()
plot(m2)

# load("results/mcmc.Rdata")
summary(m2)

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
read.csv("data/species_traits_summary.csv", sep =",")%>% 
  dplyr::select(Taxon,  family, community, p50, oil.content, ratio)%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID) %>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio),
         Lp50 = log(p50), 
         SRp50 = sqrt(p50), )%>%
  na.omit()%>%
  as.data.frame()-> oil_p50 #n=33

unique(oil_p50$familia) # 16 families
unique(oil_p50$Taxon) #33 levels (silene ciliata and thymus in both communities)
setdiff(oil_data$Taxon,oil_p50$Taxon)
setdiff(longevity_sp$Taxon,oil_p50$Taxon)
setdiff(oil_p50$Taxon,oil_data$Taxon)

hist(oil_p50$p50) # normally distributed?
hist(oil_p50$Lp50)
hist(oil_p50$SRp50) # normally distributed'
hist(oil_p50$oil.content) #not normally distributed
hist(oil_p50$Loil.content)# almost normally distributed
hist(oil_p50$ratio) # not normally distributed
hist(oil_p50$Lratio) # almost normally distributed

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# model
MCMCglmm::MCMCglmm( SRp50 ~ Lratio , # Loil.content  Lratio
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_p50,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g6
g5 # square root p50 ~ Loil.content
g6 # square root p50 ~ Lratio

x11()
plot(g6)
summary(g6)


# T50/EHS (from germination phenology experiment) ####
read.csv("data/species_traits_summary.csv", sep =",")%>% 
  dplyr::select(Taxon, family, community, T50_mean, 
                EHS_mean ,oil.content, ratio)%>%
  na.omit() %>%# 
  group_by(Taxon, community, family)%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio),
         LEHS_mean = log (EHS_mean),
         LT50_mean = log(T50_mean))%>%
  as.data.frame()-> oil_t50 # 36 species
unique(oil_t50$Taxon)
hist(oil_t50$T50_mean) # not normally dsitributed
hist(oil_t50$LT50_mean)
hist(oil_t50$EHS_mean) # not normally dsitributed
hist(oil_t50$LEHS_mean)  # normally dsitributed
hist(oil_p50$oil.content) #not normally distributed
hist(oil_p50$Loil.content)# almost normally distributed
hist(oil_p50$ratio) # almost normally distributed
hist(oil_p50$Lratio) # normally distributed

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
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# model
MCMCglmm::MCMCglmm(LEHS_mean~ Lratio , # Lratio  Loil.content
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_t50,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g10

g7 # T50_mean ~ Loil.content
g8 # T50_mean ~ Lratio
g9 # LEHS_mean ~ Loil.content 
g10 # LEHS_mean~ Lratio  


x11()
plot(g10)
summary(g10) 

########################### ECOLOGICAL TRADE_OFFS KEEP working on!!!####################### 
# Ecology####
read.csv("data/species_traits_summary.csv")%>% #n= 49
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio))%>%
  na.omit()-> oil_ecology # 48 species
setdiff(oil_ecology$Taxon, oil_data$Taxon)
setdiff(oil_data$Taxon,oil_ecology$Taxon)
unique(oil_ecology$Taxon)

hist(oil_ecology$oil.content) # not normally dsitributed
hist(oil_ecology$Loil.content)
hist(oil_ecology$ratio)
hist(oil_ecology$Lratio)

##### MCMC 
# take into account phylogeny! 
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_oil.tree")), 
                    ape::read.tree("results/tree_oil.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL
plot(ape::read.tree("results/tree_oil.tree"))
### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(Lratio ~ ecology, # Lratio Loil.content
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_ecology,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g12

g11 #Loil.content ~ ecology
g12 # Lratio ~ ecology
x11()
plot(g12)
summary(g12) 
# GDD ####
read.csv("data/species_traits_summary.csv")%>%
  dplyr::select(Taxon, community, family, oil.content, ratio, GDD, FDD, Snw)%>%
  na.omit () %>% # 36 species
  rename(familia = family)%>%
  convert_as_factor(Taxon, community, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio),
         LFDD=log(FDD),
         LGDD=log(GDD),
         LSnw=log(Snw))%>%
  as.data.frame()-> oil_sp_pref #46
unique(oil_sp_pref$Taxon)

hist(oil_sp_pref$oil.content) # not normally dsitributed
hist(oil_sp_pref$Loil.content) # normally distributed
hist(oil_sp_pref$ratio)
hist(oil_sp_pref$Lratio) # clearly normally distributed
hist(oil_sp_pref$FDD)
hist(oil_sp_pref$LFDD)# somewhat normally distributed
hist(oil_sp_pref$GDD) # bimodal (one for each community most likely)
hist(oil_sp_pref$LGDD) # still not normally distributed
hist(oil_sp_pref$Snw)
hist(oil_sp_pref$LSnw) # somewhat normally distributed?

##### MCMC 
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
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))

# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(Lratio~ Snw, # Lratio Loil.content
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_sp_pref,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g18

g13 # Loil.content ~ GDD
g14 # Lratio ~ GDD
g15 # Loil.content ~ FDD
g16 # Lratio ~ FDD
g17 # Loil.content ~ Snow
g18 # Lratio ~ Snow

x11()
plot(g18)
summary(g18) 

