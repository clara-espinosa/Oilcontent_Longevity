library(tidyverse);library(rstatix)
library(vegan);library(MCMCglmm);library(ape)


# Data format preparation for mcmc###
read.csv("data/oil_regionaldata.csv", sep = ",")%>%
  rename(seed_mass= X50seed_mass_mg)%>%
  group_by(Taxon, family, ecology)%>%
  summarise(oil_content = mean(oil_content), seed_mass = mean(seed_mass))%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia, ecology) %>%
  mutate(ecology= fct_relevel(ecology, "Strict lowland", "Generalist", "Strict alpine"))%>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil_content = log(oil_content),
         Lseed_mass=log(seed_mass))%>%
  na.omit()%>%
  as.data.frame()-> regional_oil_analysis

# data exploration and visualization
hist(regional_oil_analysis$oil_content) # clearly not normal
hist(regional_oil_analysis$seed_mass)# clearly not normal
hist(regional_oil_analysis$Loil_content) # normally distributed
hist(regional_oil_analysis$Lseed_mass)# normally distributed
qqnorm(regional_oil_analysis$oil_content)

# descriptive statistics
regional_oil_analysis%>% 
  dplyr::select(Taxon, ecology, oil_content, seed_mass)%>%
  get_summary_stats(seed_mass)

##### USE MCMC to take into account phylogeny! 
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_regionaldata.tree")), 
                    ape::read.tree("results/tree_regionaldata.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL
plot(ape::read.tree("results/tree_regionaldata.tree"))
### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))

# 1. COMPARISON between altitudinal preferences#####
str(regional_oil_analysis)
MCMCglmm::MCMCglmm(Loil_content ~ ecology, 
                   random = ~ animal,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = regional_oil_analysis,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g1

plot(g1)
summary(g1) 
  
### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- g1$VCV[,"animal"]/(g1$VCV[,"animal"] + g1$VCV[,"units"]) 

mean(g1$VCV[,"animal"]/(g1$VCV[,"animal"] + g1$VCV[,"units"])) %>% round(2)
coda::HPDinterval(g1$VCV[,"animal"]/(g1$VCV[,"animal"] + g1$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(g1$VCV[,"animal"]/(g1$VCV[,"animal"] + g1$VCV[,"units"]))[, 2] %>% round(2)

# Random effects animal
summary(g1)$Gcovariances[1, 1] %>% round(2) 
summary(g1)$Gcovariances[1, 2] %>% round(2) 
summary(g1)$Gcovariances[1, 3] %>% round(2)


# 2. Relationship between oil content and seed mass (trend but not significant, too much variation) ####
str(regional_oil_analysis)
MCMCglmm::MCMCglmm(Lseed_mass ~oil_content, 
                   random = ~ animal,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = regional_oil_analysis,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g2

plot(g2)
summary(g2) # no significant effect of oil content in seed mass (or viceversa)

### Random and phylo

# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- g2$VCV[,"animal"]/(g2$VCV[,"animal"] + g2$VCV[,"units"]) 

mean(g2$VCV[,"animal"]/(g2$VCV[,"animal"] + g2$VCV[,"units"])) %>% round(2)
coda::HPDinterval(g2$VCV[,"animal"]/(g2$VCV[,"animal"] + g2$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(g2$VCV[,"animal"]/(g2$VCV[,"animal"] + g2$VCV[,"units"]))[, 2] %>% round(2)

# Random effects animal
summary(g2)$Gcovariances[1, 1] %>% round(2) 
summary(g2)$Gcovariances[1, 2] %>% round(2) 
summary(g2)$Gcovariances[1, 3] %>% round(2)
