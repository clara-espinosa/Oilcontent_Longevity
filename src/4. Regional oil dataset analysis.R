library(tidyverse);library(readxl);library(rstatix)
library(vegan);library(glmmTMB);library(DHARMa)
library(phylosignal);library(phylobase);library(ape);library(tidytree)

# oil content data (own + literature)
read.csv("data/oil_fulldataset.csv", sep = ";")%>% rename(seed_mass= X50seed_mass_mg)-> oil_fulldataset # 80 species
unique(oil_fulldataset$family) # 19 families

# data exploration and visualization
hist(oil_fulldataset$oil_content) # clearly not normal
hist(oil_fulldataset$seed_mass)# clearly not normal

# Data format preparation for mcmc###

read.csv("data/oil_fulldataset.csv", sep = ";")%>%
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
  as.data.frame()-> full_oil_analysis

# descriptive statistics
full_oil_analysis%>% 
  dplyr::select(Taxon, ecology, oil_content, seed_mass)%>%
  #group_by(ecology)%>%
  get_summary_stats(seed_mass)

hist(full_oil_analysis$Loil_content) # normally distributed
hist(full_oil_analysis$Lseed_mass)# normally distributed

##### USE MCMC to take into account phylogeny! 
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_oil_fulldataset.tree")), 
                    ape::read.tree("results/tree_oil_fulldataset.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL
plot(ape::read.tree("results/tree_oil_fulldataset.tree"))
### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)

# 1. COMPARISON between altitudinal preferences, NO DIFFERENCES RELATED TO altitude #####
# correct glm? ASK EDUARDO!!!
str(full_oil_analysis)
MCMCglmm::MCMCglmm(Loil_content ~ ecology, 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = full_oil_analysis,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g1

plot(g1)
summary(g1) 
  
# 2. Relationship between oil content and seed mass (trend but not significant, too much variation) ####
str(full_oil_analysis)
MCMCglmm::MCMCglmm(Loil_content ~ Lseed_mass, 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = full_oil_analysis,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g2

plot(g2)
summary(g2) # no significant effect of seed mass in oil content (or viceversa)
