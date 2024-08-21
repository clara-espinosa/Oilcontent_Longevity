library(tidyverse);library (rstatix);library (stringr)
library(ggpattern); library (vegan) ;library (ggrepel)

### read species oil content data ####
read.csv("data/species_oil.csv")%>%
  full_join (oil_data_full, by=c("Taxon", "community", "family"))%>% # from format_oil_data script
  select(community, Taxon,family, ecology, oil.content, ratio)-> oil_data # 50 species
unique(oil_data$Taxon)
###################################### BIOLOGICAL TRADE-OFFS ################################################
## read species with seed mass data ####
read.csv("data/seed_mass.csv", sep = ",")%>%
  group_by(Taxon, community)%>%
  get_summary_stats(mass_50)%>%
  dplyr::select(Taxon, community, n, mean)%>%
  rename(mass_50 = mean)%>%
  as.data.frame()-> seedmass # 66 species
### read longevity data (p50) ####
read.csv("data/longevity/genstat.csv", sep= ",")%>%
  dplyr::select(Taxon, community, p50)%>%
  group_by(Taxon, community)%>%
  summarise(p50 = mean(p50))-> p50 # 42 species

### read germ data (t50/cold_strat_germ) ####
read.csv("data/germ_t50_data.csv")%>% # so far only t50 = !!
  dplyr:: select(Taxon, community, mean_T50, F_T50, S_T50)->t50  #45 species with t50 trait
  #full_join(read.csv("data/germ_drivers_data.csv")) -> germ_traits # %>%  
  # na.omit() # 
###################################### ECOLOGICAL TRADE-OFFS ################################################
#### Load SPECIES PREFERENCES data####
# based from 80 iButtons floristic relevés
# original scripts and data in Omaña80 and Picos Github repository (respectively)
#script name: indices and glms

# FROM ibuttons data: in picos from 02/10/2018 to 07/08/2019; in villabandin from 12/07/2021 to 29/05/2022
# FDD: sum of degrees (daily mean) when daily Tmean is below 0 ºC (in absolute values)
# GDD: sum of degrees (daily mean) when daily Tmean is above 5 ºC 
# Snowdays = sum of days with snow cover, calculated as days when Tmax <0.5 and Tmin > -0.5 ºC
# These variables (altogether with bio1, bio2 and bio7) calculated mean x day, then month and finally year (whole datset)

# From inventory data: 80 plots per each community, at each plot species coverage estimated in %
# To calculate species preferences several filters applied:
# consider only plots where each species has at least 10% coverage
# climatic variables weighted per coverage 

read.csv("data/sp_pref.csv", sep =",") -> sp_pref # 127 species

# merge all info from species and create 1 header dataset for plotting 
oil_data %>%
  full_join(p50,by= c("Taxon", "community") )%>%
  left_join(seedmass, by= c("Taxon", "community"))%>%
  left_join(t50, by= c("Taxon", "community"))%>%
  left_join(sp_pref, by = c("Taxon", "community"))%>%
  dplyr::select(Taxon_full, Taxon, community, family, ecology, 
                oil.content, ratio, mass_50, p50, mean_T50, F_T50, S_T50, GDD, FDD, Snw)%>%
  write.csv("data/species_traits_summary.csv")


##### EXTRA ###################
#### PCA and correlations (check our previous categorical classification based on DCA (fell vs Snow vs neutral)
# when considered both communities together previous categorical classification DO NOT match
# try dividing communities?
header_clim[,6:11] %>%
  FactoMineR::PCA() -> pca_sp_pref
  
pca_sp_pref$var$contrib
pca_sp_pref$eig
header_clim[,6:11] %>%  cor()

cbind((header_clim %>%  dplyr::select(community, species, familia, microhabitat, distribution)), data.frame(pca_sp_pref$ind$coord[, 1:4])) %>%
  mutate(species = factor(species)) -> pcaInds

pca_sp_pref$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") %>%
  mutate(Variable = fct_recode(Variable, "Snow" = "Snw"))-> pcaVars

### Plot PCA
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = microhabitat), color = "black", show.legend = T, size = 3, shape = 21) + # family
  geom_label(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 4) +
  #geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 5, segment.size= 1,
                  # point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  geom_text_repel (data = pcaInds, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(title= "Species climatic preferences", tag = "")+
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
  scale_x_continuous(name = paste("Axis 1 (", round(pca_sp_pref$eig[1, 2], 0),
                                  "% variance explained)", sep = "")) + #, limits = c(-5, 5)
  scale_y_continuous(name = paste("Axis 2 (", round(pca_sp_pref$eig[2, 2], 0), 
                                  "% variance explained)", sep = "")) #, limits = c(-4, 4)
# only temperate
# categorical classification doesn't match either, but seems like there is a pattern
header_clim%>%
  filter(community == "Temperate") -> header_clim_tem

header_clim_tem[,6:11] %>%
  FactoMineR::PCA() -> pca_sp_pref

pca_sp_pref$var$contrib
pca_sp_pref$eig
pca_sp_pref[,6:11] %>%  cor()

cbind((header_clim_tem %>%  dplyr::select(community, species, familia, microhabitat, distribution)), data.frame(pca_sp_pref$ind$coord[, 1:4])) %>%
  mutate(species = factor(species)) -> pcaInds

pca_sp_pref$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") %>%
  mutate(Variable = fct_recode(Variable, "Snow" = "Snw"))-> pcaVars

### Plot PCA
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = distribution), color = "black", show.legend = T, size = 4, shape = 21) + # family
  geom_label(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 4) +
  #geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 5, segment.size= 1,
  # point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  geom_text_repel (data = pcaInds, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(title= "Species climatic preferences", tag = "")+
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
  scale_x_continuous(name = paste("Axis 1 (", round(pca_sp_pref$eig[1, 2], 0),
                                  "% variance explained)", sep = "")) + #, limits = c(-5, 5)
  scale_y_continuous(name = paste("Axis 2 (", round(pca_sp_pref$eig[2, 2], 0), 
                                  "% variance explained)", sep = "")) #, limits = c(-4, 4)
# only mediterranean
# categorical classification doesn't match either, but seems like there is a pattern
# seem to be a king of gradient between strict alpine and generalist
header_clim%>%
  filter(community == "Mediterranean") -> header_clim_med

header_clim_med[,6:11] %>%
  FactoMineR::PCA() -> pca_sp_pref

pca_sp_pref$var$contrib
pca_sp_pref$eig
pca_sp_pref[,6:11] %>%  cor()

cbind((header_clim_med %>%  dplyr::select(community, species, familia, microhabitat, distribution)), data.frame(pca_sp_pref$ind$coord[, 1:4])) %>%
  mutate(species = factor(species)) -> pcaInds

pca_sp_pref$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") %>%
  mutate(Variable = fct_recode(Variable, "Snow" = "Snw"))-> pcaVars

### Plot PCA
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = distribution), color = "black", show.legend = T, size = 4, shape = 21) + # family
  geom_label(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 4) +
  #geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 5, segment.size= 1,
  # point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  geom_text_repel (data = pcaInds, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(title= "Species climatic preferences", tag = "")+
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
  scale_x_continuous(name = paste("Axis 1 (", round(pca_sp_pref$eig[1, 2], 0),
                                  "% variance explained)", sep = "")) + #, limits = c(-5, 5)
  scale_y_continuous(name = paste("Axis 2 (", round(pca_sp_pref$eig[2, 2], 0), 
                                  "% variance explained)", sep = "")) #, limits = c(-4, 4)
