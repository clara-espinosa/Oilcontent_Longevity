library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan);library

#### oil data transformation #####
# oil data mg oil/g sample
read.csv("data/oil_mg_wide.csv")%>%
  dplyr:: select(Taxon:C24.0)%>%
  mutate(species=make.cepnames(species))%>%# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
  gather(oil_type, oil_type_mg, 7:26)-> oil_mg # transform to long format

read.csv("data/oil_mg_wide.csv")%>%
  dplyr:: select(Taxon:C24.0)%>%
  mutate(species=make.cepnames(species))%>%
  group_by (family, community)%>%
  get_summary_stats()
  
# use total amount of oil as explanatory variable in PCA
oil_mg%>%
  convert_as_factor(Taxon, species, family, community)%>%
  group_by(species, community, family) %>%
  summarise(PERoil= (sum(oil_type_mg)/10))%>%
  data.frame()-> PERoil_sp

# percentage of dfferent types of oil
read.csv("data/oil_per_wide.csv")%>%
  dplyr:: select(Taxon:C24.0)%>%
  mutate(species=make.cepnames(species))%>%
  merge (PERoil_sp, by = c("species", "family", "community"))-> oil_data_pca #%>%
  #dplyr::select(species:year_analisis, C14.0:C20.0, C22.0, C22.1n9, C24.0, PERoil)-> oil_data_pca # remove correlated >0.7

# trans form into long format
read.csv("oil_per_wide.csv")%>%
  dplyr:: select(Taxon:C24.0)%>%
  mutate(species=make.cepnames(species))%>%
  gather(oil_type, oil_type_PER, 7:26) -> oil_per

read.csv("oil_per_wide.csv")%>%
  dplyr:: select(Taxon:C24.0)%>%
  mutate(species=make.cepnames(species))%>%
  group_by (family, community)%>%
  get_summary_stats()

# merge mg oil and % oil into long format
oil_mg %>%
  merge(oil_per)%>%
  merge(PERoil_sp, by= c("community", "species", "family"))%>%
  write.csv("oil_data_long.csv")

## PCA without Teesdalia/PERoil ####
str(oil_data_pca)
oil_data_pca%>%
  filter(!species=="Teesconf")-> oil_data_pca
oil_data_pca[, 7:27]%>%
  FactoMineR::PCA() -> pca_oil

pca_oil$var$contrib
pca_oil$eig
oil_data_pca[, 7:27]%>%  cor() # remove C12:0, C20:1n9, C20:2n6 (hihgly correlated with other variables with higher contributions to axes)

cbind((oil_data_pca[, 7:27] %>% dplyr::select(family, species, community)), data.frame(pca_oil$ind$coord[, 1:2])) %>%
  mutate(species = factor(species)) -> pcaInds

pca_oil$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") -> pcaVars

### Plot PCA 
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = family), color = "black", show.legend = T, size = 5, shape = 21) + # family
  #geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 5) +
  #geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 5, segment.size= 1,
  # point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  geom_text_repel (data = pcaInds, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =5, max.overlaps = 15) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(title= "FAME PCA species", tag = "")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, color = "black"),
        plot.tag.position = c(0.02,1),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca_oil$eig[1, 2], 0),
                                  "% variance explained)", sep = "") ) + #,limits = c(-5, 5)
  scale_y_continuous(name = paste("Axis 2 (", round(pca_oil$eig[2, 2], 0), 
                                  "% variance explained)", sep = "")) #, limits = c(-4, 4)


# PCA with header data ####
oil_data_pca %>%
  merge(read.csv("data/species_oil.csv"))%>%
  dplyr::select(community, Taxon, species, family, ecology,PERoil:seedmass, FDD: under_snow)%>%
  convert_as_factor(community, Taxon, species, family, ecology)%>%
  as.data.frame()-> oil_data_pca2

str(oil_data_pca2)
oil_data_pca2[, 6:30]%>%  cor() # remove C12:0, C20:1n9, C20:2n6 (hihgly correlated with other variables with higher contributions to axes)

oil_data_pca2[, 6:30]%>%
  FactoMineR::PCA() -> pca_oil

pca_oil$var$contrib
pca_oil$eig
oil_data_pca2$family
cbind((oil_data_pca2[, 1:5]), data.frame(pca_oil$ind$coord[, 1:2])) %>%
  mutate(species = factor(species)) -> pcaInds

pca_oil$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") -> pcaVars

### Plot PCA
x11()
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = family), color = "black", show.legend = T, size = 5, shape = 21) + # family
  #geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 5) +
  #geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 5, segment.size= 1,
  # point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  geom_text_repel (data = pcaInds, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =5, max.overlaps = 15) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(title= "FAME PCA species", tag = "")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, color = "black"),
        plot.tag.position = c(0.02,1),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca_oil$eig[1, 2], 0),
                                  "% variance explained)", sep = "") ) + #,limits = c(-5, 5)
  scale_y_continuous(name = paste("Axis 2 (", round(pca_oil$eig[2, 2], 0), 
                                  "% variance explained)", sep = "")) #, limits = c(-4, 4)

