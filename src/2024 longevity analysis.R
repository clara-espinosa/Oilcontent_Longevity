# ANALISI 2024 LONGEVITY DATA (PAVIA) WITH SP PREFRENCES, OIL CONTENT AND SEED MASS AS EXPLANATORY VARIABLES

library(tidyverse);library (rstatix);library (stringr);library(viridis)
library(ggpattern); library (vegan) ;library (ggrepel)

### analisys using raw germination data, new explicatives variables and phylogeny via MCMC GLMM ####
### Read data #
read.csv("data/2022/germination22.csv", sep =";") %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(header) %>% # from data handling script 
  mutate(ID = gsub(" ", "_", species), animal = ID) %>%
  na.omit %>% # remove species without sp_pref (Minuartia CF) and without oil content (C.ramosissimum, G.verna, G.campestris and P. pyrenaica)
  merge(read.csv("data/2022/indices22.csv", sep =","), by = c("code", "ageing"))-> df22_new

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
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        #G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

### All species model

MCMCglmm::MCMCglmm(cbind(germinated, seeds - germinated) ~
                     scale(ageing)*scale(GDD) + scale(ageing) * scale(PERoil) + scale(ageing)*scale(ratio) , #+ scale(ageing) * scale(meanseedmass) 
                   random = ~ animal + ID + community,
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

# Gaussian model for grp, mgr and syn ####
### Gaussian priors 
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))


MCMCglmm::MCMCglmm(grp ~ scale(ageing)*scale(GDD) + scale(ageing) * scale(meanseedmass) + scale(ageing) * scale(PERoil),
                   random = ~ animal + ID + community,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = df22_new,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g1
# save(m1, file = "results/mcmc.Rdata")
x11()
plot(g1)
glm (grp ~ scale(ageing)*microhabitat*distribution, family = "gaussian", data = df) -> grp
summary(grp)

# load("results/mcmc.Rdata")
summary(g1)

#### Genstat as response variables with new explicative variables ####

read.csv("data/2022/genstat22.csv", sep =",")%>%
  left_join(header, by = c("species", "community", "code")) %>%
  dplyr::select(species, code, community, site, familia, distribution, Ki:upper95, bio1:PERoil, ratio)%>%
  convert_as_factor(species, code, community, familia,  distribution) %>%
  mutate(ID = gsub(" ", "_", species), animal = ID)%>%
  na.omit()-> genstat22_new

# compare p50, Ki, Slope!!
glm (p50 ~ GDD+meanseedmass+PERoil+ratio, family = "gaussian", data = genstat22_new) -> m3
summary(m3)
# slope affected only by PERoil
# Ki nothing significant
# p50 affected only by PERoil

# take into account phylogeny! 
### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# correct glm? ASK EDUARDO!!!
MCMCglmm::MCMCglmm(p50 ~ scale(GDD) + scale(PERoil),
                   random = ~ animal + ID + community,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = genstat22_new,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g3
x11()
plot(g3)
summary(g3)
#p50, Ki and slope no differences found according to micro and distribution

#PCA taking into account p50 ki and slope?
read.csv("data/2022/genstat22.csv", sep =",")%>%
  left_join(header, by = c("species", "community", "code")) %>%
  dplyr::select(species, code, community, site, familia, microhabitat, distribution, Ki:upper95, bio1:PERoil)%>%
  convert_as_factor(species, code, community, familia, microhabitat, distribution)%>%
  mutate(species = make.cepnames(species))%>%
  na.omit()%>%
  dplyr:: select(species, code, community, site, familia, microhabitat, distribution, # after checkin correlation remove>0.7
                 Ki, slope, p50, bio7:meanseedmass, PERoil)-> genstatPCA

genstatPCA[,8:16] %>%
  FactoMineR::PCA() -> pca_genstat

pca_genstat$var$contrib
pca_genstat$eig
genstatPCA[ ,8:16] %>%  cor()

cbind((genstatPCA %>%  dplyr::select(community, species, familia, microhabitat, distribution)), data.frame(pca_genstat$ind$coord[, 1:4])) %>%
  mutate(species = factor(species)) -> pcaInds

pca_genstat$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") %>%
  mutate(Variable = fct_recode(Variable, "Snow" = "Snw"))-> pcaVars

### Plot PCA
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = community), color = "black", show.legend = T, size = 4, shape = 21) + # family
  geom_label(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 4) +
  #geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 5, segment.size= 1,
  # point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  geom_text_repel (data = pcaInds, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(title= "PCA genstat results", tag = "")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        plot.tag.position = c(0.02,1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca_genstat$eig[1, 2], 0),
                                  "% variance explained)", sep = "")) + #, limits = c(-5, 5)
  scale_y_continuous(name = paste("Axis 2 (", round(pca_genstat$eig[2, 2], 0), 
                                  "% variance explained)", sep = "")) #, limits = c(-4, 4)

### exploration visualization ####
# ageing x final germ
read.csv("data/2022/germination22.csv", sep =";") %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE))%>%
  merge(header) %>%
  mutate(germPER = germinated/seeds)%>%
  group_by(community, ageing) %>%
  summarise(germPER)%>%
  mutate(ageing= as.factor(ageing))%>%
  ggplot(aes(x=ageing, y = germPER, fill= ageing), color="black")+
  geom_boxplot(show.legend = F)+
  scale_fill_viridis_d()+
  labs(x= "Ageing days", y = "Final germination")+
  #facet_wrap(~community)+
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position = c(0,1),
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "right", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))

# germination curves facet x community 
read.csv("data/2022/germination22.csv", sep =";") %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE))%>%
  merge(header) %>%
  mutate(germPER = germinated/seeds)%>%
  ggplot(aes(x=ageing, y = germPER, group = species, color= GDD))+ #GDD
  geom_smooth(method = "loess", se = FALSE, linewidth= 1.2)+
  scale_color_viridis ()+
  facet_wrap(~community)+
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position = c(0,1),
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "right", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))  


# p50 scatterplot (genstat) 
read.csv("data/2022/genstat22.csv", sep =",")%>%
  left_join(header, by = c("species", "community", "code")) %>%
  dplyr::select(species, code, community, site, familia,  distribution, Ki, slope, p50, GDD, PERoil)%>%
  convert_as_factor(species, code, community, familia,  distribution) %>%
  ggplot(aes(x=GDD, y=p50), color="black")+ #, fill=species
  geom_point(size= 4, shape=21)+
  labs(x= "Growing degree days (GDD)", y = "p50 (days)")+
  #geom_text_repel(aes(x=PERoil, y=p50,label=species))+
  geom_smooth(method="lm")+
  scale_fill_viridis_d ()+ #direction =-1
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position = c(0,1),
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_text(size = 13), 
         axis.text.x= element_text (size = 12, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "right", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))

# slope scatterplot (genstat) 
read.csv("data/2022/genstat22.csv", sep =",")%>%
  left_join(header, by = c("species", "community", "code")) %>%
  dplyr::select(species, code, community, site, familia, distribution, Ki, slope, p50, GDD, PERoil)%>%
  convert_as_factor(species, code, community, familia,  distribution) %>%
  ggplot(aes(x=PERoil, y=slope, fill=GDD), color="black")+
  geom_point(size= 4, shape=21)+
  geom_smooth(method="glm")+
  geom_text(aes(x=PERoil, y=slope,label=species))+
  scale_fill_viridis ()+
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position = c(0,1),
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "right", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))
