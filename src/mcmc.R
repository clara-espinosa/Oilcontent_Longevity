library(tidyverse)

### Read data

raw_df <-read.csv("data/germination.csv", sep =";") 
df <- raw_df %>%
  gather(raw_df, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv", sep =";")) %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>%
  na.omit %>% 
  unite("ecology",habitat:micro, sep=" ", remove = FALSE) %>%
  mutate(micro=factor(micro)) %>% 
  mutate(micro=fct_relevel(micro,c("neutral","snowbed","fellfield"))) %>%
  arrange(micro)

### calculate germination indices
library (GerminaR)

#change data structure
dat <- raw_df %>% 
  mutate(across(c(code, ageing), as.factor))
str(dat)
gerind <- ger_summary(SeedN = "seeds", evalName = "D", data=dat[3:7])

# choose germination indices 
#              GRP - germination percentage (from 0 to 100%)
#              MGR - mean germination rate (time units)
#              SYN - syncronization index (from 0 to 1)
gerind <- gerind %>% 
  select(grp, mgr, syn)

# join germination indices with dat info (code and ageing days)
gerind <- bind_cols(dat[,1:2], gerind)
str(gerind)
# write.csv (gerind, file = "data/indices.csv")
# add variable to the dataframe
df <- df %>% 
  merge(read.csv("data/indices.csv", sep =";"))
  

### Read tree

phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree.tree")), 
                    ape::read.tree("results/tree.tree"), rooted = TRUE) -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations

nite = 1000000
nthi = 100
nbur = 100000

# shorter iterations
nite = 10000
nthi = 10
nbur = 100

### Set priors for germination models

priors <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

### All species model

MCMCglmm::MCMCglmm(cbind(germinated, seeds - germinated) ~
                    scale(ageing)*micro + scale(ageing)*habitat,
                   random = ~ animal + ID + bedrock + site:bedrock,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = df,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1

# save(m1, file = "results/mcmc.Rdata")
x11()
plot(m1)

# load("results/mcmc.Rdata")
summary(m1)

### GLM without phylogeny

glm(cbind(germinated, seeds - germinated) ~
      ageing * micro + ageing *habitat, family = "binomial", data = df) -> m2
summary(m2)

### Overall percentages

df %>%
  group_by(ageing, micro) %>%
  summarise(p = sum(germinated) / sum(seeds))

df %>%
  group_by(ageing, habitat) %>%
  summarise(p = sum(germinated) / sum(seeds))


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

# Random effects bedrock

summary(m1)$Gcovariances[3, 1] %>% round(2)
summary(m1)$Gcovariances[3, 2] %>% round(2) 
summary(m1)$Gcovariances[3, 3] %>% round(2)

# Random effects site:bedrock

summary(m1)$Gcovariances[4, 1] %>% round(2)
summary(m1)$Gcovariances[4, 2] %>% round(2) 
summary(m1)$Gcovariances[4, 3] %>% round(2)

### Gaussian priors

priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))

# Gaussian model
MCMCglmm::MCMCglmm(syn ~ scale(ageing) * micro+ scale(ageing)*habitat,
                   random = ~ animal + ID + bedrock + site:bedrock,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = df,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g1
# save(m1, file = "results/mcmc.Rdata")
x11()
plot(g1)

# load("results/mcmc.Rdata")
summary(g1)
