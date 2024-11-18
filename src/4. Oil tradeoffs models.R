library(tidyverse);library (rstatix);library (stringr);library(viridis)
library(ggpattern); library (vegan) ;library (ggrepel)
library(lme4); library(glmmTMB); library (DHARMa) 


############################################ BIOLOGICAL TRADE-OFFS #################################################
# seed mass (log transformed) #####
read.csv("data/seed_mass.csv")%>%
  merge(oil_data, by= c("Taxon", "community"))%>%
  dplyr::select(Taxon, community, family, oil.content, ratio, mass_50)%>%
  na.omit()%>%
  group_by(Taxon, family)%>%
  summarise(mass_50 = mean(mass_50), oil.content= mean(oil.content), ratio=mean(ratio))%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Lseedmass = log(mass_50), 
         Loil.content = log(oil.content),
         Lratio=log(ratio))%>%
  na.omit()%>%
  as.data.frame()-> oil_seedmass # 47 species 
unique(oil_seedmass$animal)
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
                        #G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# mcmc glmm
MCMCglmm::MCMCglmm(Loil.content  ~ Lseedmass, #  Loil.content   Lratio 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_seedmass,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g1
g1 # Loil.content ~ Lseedmass
g2 # Lratio ~ Lseedmass

# None show a significant relationship although negative trends

x11()
plot(g1)   
summary(g1)  

# longevity (raw germination curves and p50) in 2024 longevity analysis script more options####
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
  na.omit()%>%# remove species without sp_pref (Minuartia CF) and without oil content (C.ramosissimum, G.verna, G.campestris and P. pyrenaica)
  as.data.frame()-> df24 

unique(df24$familia)
unique(df24$Taxon) #31 species
setdiff(oil_data$Taxon,df24$Taxon)
setdiff(longevity_sp$Taxon,df24$Taxon)
setdiff(df22_new$Taxon,oil_data$Taxon)
# data distribution
hist(df24 $oil.content)
hist(df24 $Loil.content) # almost normal distribution after log transformation
hist(df24 $ratio) # almost  normal distributed
hist(df24 $Lratio)# normal distributed


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
                        #G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        #G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

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
plot(m1)

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
  group_by(Taxon, community, family)%>%
  summarise(p50 = mean(p50), oil.content= mean (oil.content), ratio = mean(ratio))%>%
  na.omit() %>% # 
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID) %>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio),
         Lp50 = log(p50))%>%
  na.omit()%>%
  as.data.frame()-> oil_p50 #n=29

unique(oil_p50$familia)
unique(oil_p50$Taxon) #27 levels (silene ciliata and thymus in both communities)
setdiff(oil_data$Taxon,oil_p50$Taxon)
setdiff(longevity_sp$Taxon,oil_p50$Taxon)
setdiff(oil_p50$Taxon,oil_data$Taxon)

hist(oil_p50$p50) # normally distributed
hist(oil_p50$Lp50)
hist(oil_p50$oil.content) #not normally distributed
hist(oil_p50$Loil.content)# almost normally distributed
hist(oil_p50$ratio) # normally distributed
hist(oil_p50$Lratio) 

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# model
MCMCglmm::MCMCglmm( p50 ~ Loil.content , # Loil.content  Lraio
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_p50,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g4
g3 # p50 ~ Loil.content
g4 # p50 ~ Lratio

x11()
plot(g4)
summary(g4)


# T50 ()####
read.csv("data/species_traits_summary.csv", sep =",")%>% 
  dplyr::select(Taxon, family, community, mean_T50, F_T50, S_T50, oil.content, ratio)%>%
  na.omit() %>%# 47 species
  group_by(Taxon, family)%>%
  summarise(mean_T50 = mean(mean_T50), F_T50 = mean(F_T50), S_T50 = mean(S_T50), 
            oil.content= mean (oil.content), ratio = mean(ratio))%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio),
         LS_T50 = log(S_T50),
         LF_T50 = log(F_T50),
         Lmean_T50 = log(mean_T50))%>%
  na.omit()%>%
  as.data.frame()-> oil_t50 # 33 species
hist(oil_t50$mean_T50) # not normally dsitributed
hist(oil_t50$Lmean_T50)
hist(oil_t50$F_T50)
hist(oil_t50$LF_T50)
hist(oil_t50$S_T50)
hist(oil_t50$ LS_T50)
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
                        #G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# model
MCMCglmm::MCMCglmm(mean_T50 ~ Lratio , # Lratio  Loil.content
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_t50,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g6

g5 # mean_T50 ~ Loil.content
g6 # mean_T50 ~ Lratio
g7 # F_T50 ~ Loil.content not updated
g8 # F_T50 ~ Lratio  not updated
g9 # S_T50 ~ Loil.content not updated
g10 # S_T50 ~ Lratio not updated

x11()
plot(g6)
summary(g6) 

# 1 model with all biological correlates too many explanatory variables for n= 25 remove?? #########
read.csv("data/species_traits_summary.csv", sep =",")%>% 
  dplyr::select(Taxon,community , family, ecology, oil.content, ratio, mass_50, p50, mean_T50, F_T50, S_T50)%>% 
  na.omit()%>%
  group_by(Taxon, family)%>%
  summarise(mean_T50 = mean(mean_T50), F_T50 = mean(F_T50), S_T50 = mean(S_T50), 
            oil.content= mean (oil.content), ratio = mean(ratio),
            p50 = mean(p50), mass_50 = mean(mass_50))%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Lmass_50 = log(mass_50),
         Lp50 = log(p50),
         Lmean_T50 = log(mean_T50),
         LF_T50 = log(F_T50),
         LS_T50 = log(S_T50),
         Loil.content = log(oil.content),
         Lratio=log(ratio))%>%
  na.omit()-> oil_bio # N = 25

oil_bio%>%
  dplyr::select(oil.content, ratio, mass_50, p50, mean_T50, F_T50, S_T50, Lmass_50:Lratio)%>%
  gather(trait, values, oil.content:Lratio)%>%
  ggplot()+
  geom_histogram(aes(values, fill = trait))+
  facet_wrap(~trait, scales = "free")

### Read tree

phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_oil.tree")), 
                    ape::read.tree("results/tree_oil.tree"), method = "ultrametric") -> 
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

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# correct glmm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(ratio ~  scale(Lmass_50)+scale(Lp50)+scale(LS_T50),
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_bio,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> gbio_ratio
gbio_oil # Loil.content ~  scale(Lmass_50)+scale(Lp50)+scale(LS_T50)
gbio_ratio # ratio ~  scale(Lmass_50)+scale(Lp50)+scale(LS_T50)
x11()
plot(gbio_oil)
summary(gbio_ratio) # nada significativo si ponemos todo junto en Loil.content ni en Lratio!!
rm(ratio)
rm(gbioratio)

########################### ECOLOGICAL TRADE_OFFS KEEP working on!!!####################### 
# Ecology/distribution ####
read.csv()
str(oil_data)
oil_data%>% #n= 50
  group_by(Taxon, family, ecology)%>%
  summarise(oil.content = mean(oil.content), ratio = mean(ratio))%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio))%>%
  na.omit()-> oil_ecology # 48 species
setdiff(oil_ecology$Taxon, oil_seedmass$Taxon)
setdiff(oil_seedmass$Taxon,oil_ecology$Taxon)
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
                        #G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(Loil.content ~ ecology, # Lratio 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_ecology,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g11

g11 #Loil.content ~ ecology
g12 # Lratio ~ ecology
x11()
plot(g3)
summary(g3) 
# GDD ####
# try quadratic explanatory variables
read.csv("data/species_traits_summary.csv")%>%
  na.omit () %>% # 36 species
  rename(familia = family)%>%
  convert_as_factor(Taxon, community, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio))%>%
  na.omit()-> oil_sp_pref #


# Gamma(link="log")
a <- glmmTMB(Loil.content ~ GDD^2 + (1| familia) ,  family = gaussian, data= oil_sp_pref) 
a <- glmmTMB(Lratio  ~ GDD^2  + (1| familia) , family = Gamma(link="log"), data= oil_sp_pref) 
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
MCMCglmm::MCMCglmm(Lratio ~ GDD, # ratio Loil.content
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_sp_pref,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g14
g13 #Lratio ~ GDD
g14 #Loil.content ~ GDD
x11()
plot(g13)
summary(g14) 
# FDD ####
str(oil_sp_pref)
# Gamma(link="log")
a <- glmmTMB(Loil.content ~ FDD*community + (1| familia) ,  family = Gamma(link="log"), data= oil_sp_pref) 
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
MCMCglmm::MCMCglmm(Loil.content ~ FDD, # Lratio 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_sp_pref,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g15

g15 #Loil.content ~ FDD
g16 #Lratio ~ FDD
x11()
plot(g15)
summary(g15) 

# 1 model with all ecological correlates ####
oil_data%>% #n= 36
  merge(sp_pref, by=c("Taxon", "community"))%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, species, community, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio))%>%
  na.omit()-> oil_eco #
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
MCMCglmm::MCMCglmm(Lratio ~ ecology*FDD + ecology*GDD, # ratio 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_eco,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> geco
x11()
plot(geco)
summary(geco)
