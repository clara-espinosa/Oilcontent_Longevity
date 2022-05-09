library(tidyverse)

### Read data

read.csv("data/germination.csv") %>%
  gather(days, germinated, D7:D28) %>%
  group_by(code, ageing, sown) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv"), by = "code") %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>%
  na.omit -> df

### Read tree

phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree.tree")), 
                    ape::read.tree("results/tree.tree"), rooted = TRUE) -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations

nite = 1000000
nthi = 100
nbur = 100000

nite = 10000
nthi = 10
nbur = 100

### Set priors for germination models

priors <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

### All species model

MCMCglmm::MCMCglmm(cbind(germinated, sown - germinated) ~
                     scale(ageing) * snow + scale(ageing) * habitat + scale(ageing) * bedrock,
                   random = ~ animal + ID + site,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = df,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1

# save(m1, file = "results/mcmc.Rdata")
plot(m1)

# load("results/mcmc.Rdata")
summary(m1)

### GLM without phylogeny

glm(cbind(germinated, sown - germinated) ~
      ageing * snow, family = "binomial", data = df) -> m2
summary(m2)

### Overall percentages

df %>%
  group_by(ageing, snow) %>%
  summarise(p = sum(germinated) / sum(sown))

df %>%
  group_by(ageing, bedrock) %>%
  summarise(p = sum(germinated) / sum(sown))

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

# Random effects site

summary(m1)$Gcovariances[3, 1] %>% round(2)
summary(m1)$Gcovariances[3, 2] %>% round(2) 
summary(m1)$Gcovariances[3, 3] %>% round(2)