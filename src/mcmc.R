library(tidyverse);library (rstatix);library (stringr)

library(ggpattern)
#### SPECIES PREFERENCES ####
read.csv("data/sp_pref_picos.csv", sep =";")%>%
  select(species:Snw) -> pref_picos
read.csv("data/sp_pref_villa.csv", sep =";")%>%
  select(species:Snw) -> pref_villa

rbind(pref_picos, pref_villa) -> sp_pref
##################  LONGEVITY 2022  #########################
### Read data
raw22 <-read.csv("data/2022/germination22.csv", sep =";") 
raw22 %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/2022/species22.csv", sep =";")) %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>%
  na.omit %>% 
  unite("ecology",distribution:microhabitat, sep=" ", remove = FALSE) %>%
  mutate(microhabitat=factor(microhabitat)) %>% 
  mutate(microhabitat=fct_relevel(microhabitat,c("neutral","snowbed","fellfield"))) %>%
  arrange(microhabitat) -> df22

### calculate germination indices
library (GerminaR)

#change data structure
dat22 <- raw22 %>% 
  mutate(across(c(code, ageing), as.factor))
str(dat22)
gerind22 <- ger_summary(SeedN = "seeds", evalName = "D", data=dat22[3:7])

# choose germination indices 
#              GRP - germination percentage (from 0 to 100%)
#              MGR - mean germination rate (time units)
#              SYN - syncronization index (from 0 to 1)
gerind22 <- gerind22 %>% 
  select(grp, mgr, syn)

# join germination indices with dat info (code and ageing days)
gerind22 <- bind_cols(dat22[,1:2], gerind22)
str(gerind22)
# write.csv (gerind22, file = "data/2022/indices22.csv")
# add variable to the dataframe
df22 <- df22 %>% 
  merge(read.csv("data/2022/indices22.csv", sep =";"))
  

### Read tree

phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree22.tree")), 
                    ape::read.tree("results/tree22.tree"), rooted = TRUE) -> 
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
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

### All species model

 MCMCglmm::MCMCglmm(cbind(germinated, seeds - germinated) ~
                    scale(ageing)*microhabitat + scale(ageing) * distribution,
                   random = ~ animal + ID + bedrock + site:bedrock,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = df22,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1

# save(m1, file = "results/mcmc.Rdata")
x11()
plot(m1)

# load("results/mcmc.Rdata")
summary(m1)

### GLM without phylogeny

glm(cbind(germinated, seeds - germinated) ~
      ageing * microhabitat + ageing *distribution, family = "binomial", data = df22) -> m2
summary(m2)

### Overall percentages

df22 %>%
  group_by(ageing, microhabitat) %>%
  summarise(p = sum(germinated) / sum(seeds))

df22 %>%
  group_by(ageing, distribution) %>%
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
MCMCglmm::MCMCglmm(grp ~ scale(ageing),
                   random = ~ animal + ID + bedrock + site:bedrock,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = df22,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g1
# save(m1, file = "results/mcmc.Rdata")
x11()
plot(g1)
glm (grp ~ scale(ageing)*microhabitat*distribution, family = "gaussian", data = df22) -> grp
summary(grp)

# load("results/mcmc.Rdata")
summary(g1)

#### GENSTAT DATA 2022 ANALISYS ####
genstat22 <-read.csv("data/2022/genstat22.csv", sep =";") 
str (genstat22)
genstat22$code <- as.factor(genstat22$code)
genstat22$species <- as.factor(genstat22$species)
genstat22$familia <- as.factor(genstat22$familia)
genstat22$site <- as.factor(genstat22$site)
genstat22$bedrock <- as.factor(genstat22$bedrock)
genstat22$micro <- as.factor(genstat22$micro)
genstat22$distribution <- as.factor(genstat22$distribution)
View (genstat22)

genstat22 <- genstat22 %>%
mutate(ID = gsub(" ", "_", species), animal = ID)

# compare p50, Ki, Slope!!
glm (slope ~ distribution*micro, family = "gaussian", data = genstat22) -> m3
summary(m3)
# Normal glm show no effect of microhabitat preference and distribution on p50, Ki and slope values

# take into account phylogeny! 
# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(slope ~ micro*distribution,
                   random = ~ animal + ID + bedrock + site:bedrock,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = genstat22,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g3
summary(g3)
#p50, Ki and slope no differences found according to micro and distribution

### MCMC with species preferences as explanatory ####
raw22 %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/2022/species22.csv", sep =";")) %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>%
  na.omit %>% 
 

gerind22%>%
  merge(read.csv("data/2022/species22.csv", sep =";")) %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>%
  merge(sp_pref)%>% 
  na.omit ()-> test1
str(test1)

# try out normal glm (very mixd not clear results)
glm (grp ~ ageing*Snw, family = "gaussian", data = test1) -> test    #*FDD *GDD*Snw
summary(test)

ggplot (test1, aes(x=Snw, y= syn, color= ageing)) +
  geom_point()+
  geom_smooth(method = "lm", se=FALSE, level = 0.9)+
  labs( title= "")
#visualization show some differential patterns in grp (germination percentage) from day 30 ageing responses according to FDD/GDD/Snow
# with mgr (mean germination rate) da0 species from more snow places go slower than other, with FDD and GDD contradictor responses
##################  LONGEVITY 2023 ################

### Read data
raw23 <-read.csv("data/2023/longevity23.csv", sep =";") 
raw23 %>%
  mutate(viable = seeds-empty) %>%
  group_by(species, code, ageing) %>%
  mutate(viable = sum(viable)) %>%
  gather(scores, germinated, D0:D62) %>%
  group_by(species,code,  ageing) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE), viable = first(viable)) %>%
  merge(read.csv("data/2023/species23.csv", sep =";")) %>%
  rename(seeds = viable)%>% 
  mutate(ID = gsub(" ", "_", species), animal = ID)%>%
  select(!(family))-> df23

str(df23)

### calculate germination indices
library (GerminaR)

#change data structure
raw23 %>%
  mutate(viable = seeds-empty) %>%
  group_by(species, code, ageing) %>%
  mutate(seeds = sum(viable)) %>%
  mutate(ageing = factor(ageing),
         code = factor(code),
         species = factor(species))%>%
  select(species, code, ageing, seeds, D0:D62) %>%
  as.data.frame()-> dat23
str(dat23)
gerind23 <- ger_summary(SeedN = "seeds", evalName = "D", data=dat23[4:13])

# choose germination indices 
#              GRP - germination percentage (from 0 to 100%)
#              MGR - mean germination rate (time units)
#              SYN - syncronization index (from 0 to 1)
gerind23 <- gerind23 %>% 
  select(grp, mgr, syn)

# join germination indices with dat info (code and ageing days)
gerind23 <- bind_cols(dat23[,1:3], gerind23)
str(gerind23)
# write.csv (gerind23, file = "data/2023/indices23.csv")
# add variable to the dataframe
df23 %>% 
  merge(read.csv("data/2023/indices23.csv", sep =";")) %>% 
  rename(seeds = viable) ->df23
  


### Read tree

phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree23.tree")), 
                    ape::read.tree("results/tree23.tree"), rooted = TRUE) -> 
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
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

### All species model

MCMCglmm::MCMCglmm(cbind(germinated, seeds - germinated) ~
                     scale(ageing)*microhabitat_preference + scale(ageing) * community,
                   random = ~ animal + ID + collection_site,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = df23,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1

# save(m1, file = "results/mcmc.Rdata")
x11()
plot(m1)

# load("results/mcmc.Rdata")
summary(m1)

### GLM without phylogeny

glm(cbind(germinated, seeds - germinated) ~
      ageing * microhabitat + ageing *distribution, family = "binomial", data = df23) -> m2
summary(m2)

### Overall percentages

df23 %>%
  group_by(ageing, microhabitat) %>%
  summarise(p = sum(germinated) / sum(seeds))

df23 %>%
  group_by(ageing, distribution) %>%
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
MCMCglmm::MCMCglmm(grp ~ scale(ageing),
                   random = ~ animal + ID + bedrock + site:bedrock,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = df,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g1
# save(m1, file = "results/mcmc.Rdata")
x11()
plot(g1)
glm (grp ~ scale(ageing)*microhabitat*distribution, family = "gaussian", data = df) -> grp
summary(grp)

# load("results/mcmc.Rdata")
summary(g1)

##### PERSISTENCE 23 #####
### Read data
library(viridis)
str(persistence)
x11()
persistence %>%
  filter(retrieval_season == "Spring_23" |  retrieval_season == "Autumn_23")%>%
  select(retrieval_season, community, microhabitat_buried, species, seeds_initial, bag1:bag3)%>%
  convert_as_factor(retrieval_season, community, microhabitat_buried, species) %>%
  mutate(species = fct_relevel (species, "Armeria duriaei", "Dianthus langeanus", "Plantago holosteum",
                                "Luzula caespitosa", "Phyteuma hemisphaericum", "Silene ciliata", 
                                "Androsace villosa",  "Carex sempervirens","Gypsophila repens", 
                                "Armeria cantabrica","Festuca glacialis","Jasione cavanillesii"))%>%
  mutate(retrieval_season = fct_relevel (retrieval_season, "Spring_23","Autumn_23"))%>%
  mutate(microhabitat_buried = recode_factor(microhabitat_buried, "Warm" = "Fellfield", 
                               "Crio"="Fellfield", 
                               "Cold" = "Snowbed", 
                               "Snow" = "Snowbed"))%>%
  gather(bag, germ, bag1:bag3)%>%
  mutate(germ_per=germ/10)%>%
  filter(species=="Armeria duriaei")%>%
  ggplot(aes(microhabitat_buried, germ_per, pattern=retrieval_season, fill=microhabitat_buried))+
  geom_boxplot(color="black")+
  geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", 
                       pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_key_scale_factor = 0.6) +
  scale_fill_manual (values = c ("chocolate2", "deepskyblue3"), guide = "none") + #
  labs (y= "Germination proportion", x= "Microhabitat buried")+ #title= "Field germination",
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  facet_wrap(~species, ncol=3)+
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (size = 22),
        strip.text.x = element_text( size = 14, face = "italic"),# face = "bold",
        strip.text.y = element_text(size = 14, angle = 360),
        legend.position = "right", #bottom
        plot.tag.position = c(0,1),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size = 13, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text (size=15), 
        axis.title.x = element_text (size=15))


persistence %>%
  select(species, community, retrieval_season, site_buried, microhabitat_buried, 
         seeds_initial, seeds_identifiables, D0)%>%
  filter (retrieval_season == "Spring_23")%>%
  na.omit()%>%
  group_by(community, microhabitat_buried, species)%>%
  summarise(seeds_identifiables = sum(seeds_identifiables), 
            germinated = sum(D0))%>%
  mutate(germPER = (germinated/seeds_identifiables))%>%
  filter(community =="temperate") %>%
  ggplot(aes(microhabitat_buried, germPER, color= species)) +
  geom_point(size =2, position = "jitter") +
  scale_color_viridis_d() +
  ylim (c(0,1))+
  labs (title = "Temperate")+
  theme_classic()
  
