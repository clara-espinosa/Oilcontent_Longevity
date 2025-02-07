library(tidyverse);library (rstatix);library (stringr)
library (vegan) ;library(MCMCglmm)

### MCMC take into account phylogeny ####
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_alpinedata.tree")), 
                    ape::read.tree("results/tree_alpinedata.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

plot(ape::read.tree("results/tree_alpinedata.tree"))
### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
############################ BIOLOGICAL CORRELATES ######################
# seed mass (log transformed) #####
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
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

# mcmc glmm
MCMCglmm::MCMCglmm(Lseedmass ~ ratio , #  
                   random = ~ animal,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_seedmass,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g4
g3 # Lseedmass ~ oil.content
g4 # Lseedmass ~ ratio

# None show a significant relationship although negative trends

x11()
plot(g4)   
summary(g4)  

### Random and phylo

# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- g4$VCV[,"animal"]/(g4$VCV[,"animal"] + g4$VCV[,"units"]) 

mean(g4$VCV[,"animal"]/(g4$VCV[,"animal"] + g4$VCV[,"units"])) %>% round(2)
coda::HPDinterval(g4$VCV[,"animal"]/(g4$VCV[,"animal"] + g4$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(g4$VCV[,"animal"]/(g4$VCV[,"animal"] + g4$VCV[,"units"]))[, 2] %>% round(2)
lambda <- g4$VCV[,"animal"]/(g4$VCV[,"animal"] + g4$VCV[,"units"]) 


# Random effects animal
summary(g4)$Gcovariances[1, 1] %>% round(2) 
summary(g4)$Gcovariances[1, 2] %>% round(2) 
summary(g4)$Gcovariances[1, 3] %>% round(2)

# longevity (p50)####
# read data
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  dplyr::select(Taxon,  family, community, p50, oil.content, ratio)%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID) %>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio),
         Lp50 = log(p50), 
         SRp50 = sqrt(p50), )%>%
  na.omit()%>%
  as.data.frame()-> oil_p50 #n=37

unique(oil_p50$familia) # 17 families
unique(oil_p50$Taxon) #35 levels (silene ciliata and thymus in both communities

hist(oil_p50$p50) # quite normally distributed

# model
MCMCglmm::MCMCglmm(p50 ~ ratio , # 
                   random = ~ animal,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_p50,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g6
g5 # p50 ~ oil.content
g6 # p50 ~ ratio

x11()
plot(g6)
summary(g6)

### Random and phylo

# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- g6$VCV[,"animal"]/(g6$VCV[,"animal"] + g6$VCV[,"units"]) 

mean(g6$VCV[,"animal"]/(g6$VCV[,"animal"] + g6$VCV[,"units"])) %>% round(2)
coda::HPDinterval(g6$VCV[,"animal"]/(g6$VCV[,"animal"] + g6$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(g6$VCV[,"animal"]/(g6$VCV[,"animal"] + g6$VCV[,"units"]))[, 2] %>% round(2)
lambda <- g6$VCV[,"animal"]/(g6$VCV[,"animal"] + g6$VCV[,"units"]) 


# Random effects animal
summary(g6)$Gcovariances[1, 1] %>% round(2) 
summary(g6)$Gcovariances[1, 2] %>% round(2) 
summary(g6)$Gcovariances[1, 3] %>% round(2)


# Environmental Heat Sum (EHS, log transformed) ####
read.csv("data/species_traits.csv", sep =",")%>% 
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  dplyr::select(Taxon, family, community,EHS_mean ,oil.content, ratio)%>%
  na.omit() %>%# 
  group_by(Taxon, community, family)%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio),
         LEHS_mean = log (EHS_mean))%>%
  as.data.frame()-> oil_EHS # 34 species 36 accessions
unique(oil_EHS$Taxon)
hist(oil_EHS$EHS_mean) # not normally dsitributed
hist(oil_EHS$LEHS_mean)  # normally dsitributed

# model
MCMCglmm::MCMCglmm(LEHS_mean~ ratio , # Lratio  Loil.content
                   random = ~ animal,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_EHS,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g8


g7 # LEHS_mean ~ oil.content 
g8 # LEHS_mean~ ratio  

x11()
plot(g8)
summary(g8) 

### Random and phylo

# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- g8$VCV[,"animal"]/(g8$VCV[,"animal"] + g8$VCV[,"units"]) 

mean(g8$VCV[,"animal"]/(g8$VCV[,"animal"] + g8$VCV[,"units"])) %>% round(2)
coda::HPDinterval(g8$VCV[,"animal"]/(g8$VCV[,"animal"] + g8$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(g8$VCV[,"animal"]/(g8$VCV[,"animal"] + g8$VCV[,"units"]))[, 2] %>% round(2)
lambda <- g8$VCV[,"animal"]/(g8$VCV[,"animal"] + g8$VCV[,"units"]) 


# Random effects animal
summary(g8)$Gcovariances[1, 1] %>% round(2) 
summary(g8)$Gcovariances[1, 2] %>% round(2) 
summary(g8)$Gcovariances[1, 3] %>% round(2)

########################### ECOLOGICAL DRIVERS ####################### 
# GDD, FDD and SNOW  ####
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  dplyr::select(Taxon, community, family, oil.content, ratio, GDD, FDD, Snw)%>%
  na.omit () %>% 
  rename(familia = family)%>%
  convert_as_factor(Taxon, community, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio))%>%
  as.data.frame()-> oil_eco #46
unique(oil_eco$Taxon)

hist(oil_eco$oil.content) # not normally dsitributed
hist(oil_eco$Loil.content) # normally distributed
hist(oil_eco$ratio)
hist(oil_eco$Lratio) # clearly normally distributed


##### MCMC 
MCMCglmm::MCMCglmm(Lratio~ Snw, # Lratio Loil.content
                   random = ~ animal,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_eco,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g14

g9 # Loil.content ~ GDD
g10 # Lratio ~ GDD
g11 # Loil.content ~ FDD
g12 # Lratio ~ FDD
g13 # Loil.content ~ Snow
g14 # Lratio ~ Snow

x11()
plot(g14)
summary(g14) 

### Random and phylo

# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- g14$VCV[,"animal"]/(g14$VCV[,"animal"] + g14$VCV[,"units"]) 

mean(g14$VCV[,"animal"]/(g14$VCV[,"animal"] + g14$VCV[,"units"])) %>% round(2)
coda::HPDinterval(g14$VCV[,"animal"]/(g14$VCV[,"animal"] + g14$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(g14$VCV[,"animal"]/(g14$VCV[,"animal"] + g14$VCV[,"units"]))[, 2] %>% round(2)
lambda <- g14$VCV[,"animal"]/(g14$VCV[,"animal"] + g14$VCV[,"units"]) 


# Random effects animal
summary(g14)$Gcovariances[1, 1] %>% round(2) 
summary(g14)$Gcovariances[1, 2] %>% round(2) 
summary(g14)$Gcovariances[1, 3] %>% round(2)
