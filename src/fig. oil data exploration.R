# Figure 1 oil data exploration
library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan);library(ggpubr)

# Panel A and B : PCA variables + species 
# JOIN FAME composition (percentages relative to the total content of oil) merge with total oil for PCA
read.csv("data/oil_per_wide.csv")%>%
  dplyr:: select(Taxon:C24.0)%>%
  mutate(species=make.cepnames(species))%>%
  merge (PERoil_sp, by = c("species"))-> oil_data_pca #%>%
#dplyr::select(species:year_analisis, C14.0:C20.0, C22.0, C22.1n9, C24.0, PERoil)-> oil_data_pca # remove correlated >0.7

str(oil_data_pca)
oil_data_pca[, 7:27]%>%
  FactoMineR::PCA() -> pca_oil

pca_oil$var$contrib
pca_oil$eig
oil_data_pca[, 7:27]%>%  cor() # remove C12:0, C20:1n9, C20:2n6 (hihgly correlated with other variables with higher contributions to axes)

# PCa with only the FA with more than 3% relative proportion within 
oil_data_pca%>%
  dplyr::select(C16.0, C18.0, C18.1n9c, C18.2n6c, C18.3n3, C18.3n6, C20.1n9, C22.1n9, PERoil)%>% #FAME >3% relative abundace explica bastante mas
  FactoMineR::PCA() -> pca_oil


oil_data_pca%>%
  dplyr::select(C16.0, C18.0, C18.1n9c, C18.2n6c, C18.3n3, C18.3n6, C20.1n9, C22.1n9, PERoil)%>%
  cor() 

cbind(oil_data_pca, data.frame(pca_oil$ind$coord[, 1:2])) %>%
  mutate(species = factor(species))%>%
  mutate(family = factor(family))-> pcaInds

pca_oil$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") -> pcaVars


### PANEL A) Plot PCA species#########
x11()
pcaInds%>%
  as.data.frame()%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Apiaceae", "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", 
                              "Fabaceae", "Salicaceae","Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))->pcaInds 
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  #coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(fill = family), color = "black", show.legend = T, size = 4, shape = 21) + # family
  #geom_text_repel (data = pcaInds, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4, max.overlaps = 5) +
  ggthemes::theme_tufte(base_size=12) + 
  scale_fill_manual(values=rainbow(17))+
  guides(fill=guide_legend(ncol=1, keywidth=0.1,keyheight=0.1,default.unit="cm")) +
  labs( tag = "A)")+#title= "PCA species",
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        legend.position = "left", 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10, color = "black"),
        plot.margin = unit(c(0, 0.2,0,0), "cm")) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca_oil$eig[1, 2], 0),
                                  "% variance explained)", sep = "") ) + #,limits = c(-5, 5)
  scale_y_continuous(name = paste("Axis 2 (", round(pca_oil$eig[2, 2], 0), 
                                  "% variance explained)", sep = "")) -> fig3A;fig3A#, limits = c(-4, 4)
### PANEL B) Plot PCA variables#########
x11()
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  #coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  #geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 3) +
  geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 4, segment.size= 1,
   point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  ggthemes::theme_tufte(base_size=12) + 
  labs( tag = "B)")+#title= "PCA variables",
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        plot.margin = unit(c(0, 0,0.2,0), "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10, color = "black"))+
  scale_x_continuous(name = paste("Axis 1 (", round(pca_oil$eig[1, 2], 0),
                                  "% variance explained)", sep = "") ) + #,limits = c(-5, 5)
  scale_y_continuous(name = paste("Axis 2 (", round(pca_oil$eig[2, 2], 0), 
                                  "% variance explained)", sep = ""))-> fig3B;fig3B #, limits = c(-4, 4)

# combine A + B panels ###
library(patchwork)
fig3A + fig3B + 
  plot_layout(guides = 'auto')-> fig3AB;fig3AB

#### PANEL C) Total oil content (%) per species ####
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/species_oil.csv"))%>%
  mutate(Taxon = as.factor(Taxon))%>%
  mutate(Taxon = fct_relevel(Taxon,"Luzula caespitosa" , "Kobresia myosuroides", "Carex sempervirens",
                             "Neoschischkinia truncatula", "Koeleria vallesiana", "Helictochloa marginata", 
                             "Festuca glacialis", "Festuca summilusitana", "Saxifraga paniculata", 
                             "Saxifraga conifera", "Sedum brevifolium", "Sedum anglicum", "Salix breviserrata",
                             "Anthyllis vulneraria", "Teesdalia conferta", "Helianthemum canum", "Helianthemum urrielense",
                             "Androsace villosa", "Thymus praecox", "Plantago alpina", "Plantago holosteum",
                             "Conopodium majus", "Phyteuma hemisphaericum", "Jurinea humilis", "Solidago virgaurea", 
                             "Armeria cantabrica", "Armeria duriaei", "Spergula morisonii", "Minuartia recurva", 
                             "Arenaria erinacea", "Gypsophila repens", "Dianthus langeanus", "Silene ciliata", 
                             "Silene acaulis"))%>% 
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Apiaceae", "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", 
                              "Fabaceae", "Salicaceae","Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>% 
  group_by(Taxon,family, order) %>%
  summarise(PERoil = mean(PERoil))%>%
  ggplot(aes(x=Taxon, y= PERoil, fill=family ) )+ #family
  #geom_point(size = 4, shape=21, color = "black", show.legend = F)+
  geom_bar(stat = "identity", show.legend = F, color="black")+ #problem with Thymus and S. ciliata (2 accesions summ % oil)
  coord_flip()+
  scale_fill_manual (values=rainbow(17))+
  ggthemes::theme_tufte(base_size=12) + 
  labs( tag = "C)", y= "Oil content (%)")+ #title= "Species oil content (%)",
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        legend.position = "", 
        plot.margin = unit(c(0, 0.2,0,0), "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))->fig3C;fig3C

#### PANEL D) oil types percentatges per sp (% specific oil /total oil ID) ####
library(viridis)
read.csv("data/oil_data_long.csv")%>%
  mutate(ID = paste(species, community))%>%
  group_by(Taxon, family, oil_type)%>%
  summarise(oil_type_PER = mean(oil_type_PER)) %>%
  mutate(Taxon = as.factor(Taxon))%>%
  mutate(Taxon = fct_relevel(Taxon,"Luzula caespitosa" , "Kobresia myosuroides", "Carex sempervirens",
                             "Neoschischkinia truncatula", "Koeleria vallesiana", "Helictochloa marginata", 
                             "Festuca glacialis", "Festuca summilusitana", "Saxifraga paniculata", 
                             "Saxifraga conifera", "Sedum brevifolium", "Sedum anglicum", "Salix breviserrata",
                             "Anthyllis vulneraria", "Teesdalia conferta", "Helianthemum canum", "Helianthemum urrielense",
                             "Androsace villosa", "Thymus praecox", "Plantago alpina", "Plantago holosteum",
                             "Conopodium majus", "Phyteuma hemisphaericum", "Jurinea humilis", "Solidago virgaurea", 
                             "Armeria cantabrica", "Armeria duriaei", "Spergula morisonii", "Minuartia recurva", 
                             "Arenaria erinacea", "Gypsophila repens", "Dianthus langeanus", "Silene ciliata", 
                             "Silene acaulis"))%>% 
  filter (oil_type_PER >3) %>% # if filter above 3% change geom_bar position = fill
  ggplot(aes(x=Taxon, y= oil_type_PER, fill= oil_type))+
  geom_bar(position = "fill", stat = "identity",  color = "black", show.legend = T)+ #position = "stack", 
  coord_flip()+
  scale_fill_viridis_d (name= "", guide = guide_legend (title.position = "top",direction = "horizontal")) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "D)", y= "Oil content (relative %)")+# title= "FA types  (>3%)", 
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        legend.position = "bottom", 
        legend.title = element_text(size = 10),
        legend.key.size = unit(1,"line"),
        legend.text = element_text(size = 10, color = "black"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = unit(c(0, 0.2,0,0), "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black")) -> fig3D;fig3D

### PANEL E) unsaturated vs saturated fatty acids composition per species ####
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/oil_type_specification.csv"), by= "oil_type")%>%
  dplyr::select(oil_type, community, Taxon, family, oil_type_PER, simple_name, saturation_type)%>%
  group_by(Taxon, family, saturation_type)%>%
  summarise(oil_type_PER = sum(oil_type_PER))%>%
  mutate(Taxon = as.factor(Taxon))%>%
  mutate(Taxon = fct_relevel(Taxon,"Luzula caespitosa" , "Kobresia myosuroides", "Carex sempervirens",
                             "Neoschischkinia truncatula", "Koeleria vallesiana", "Helictochloa marginata", 
                             "Festuca glacialis", "Festuca summilusitana", "Saxifraga paniculata", 
                             "Saxifraga conifera", "Sedum brevifolium", "Sedum anglicum", "Salix breviserrata",
                             "Anthyllis vulneraria", "Teesdalia conferta", "Helianthemum canum", "Helianthemum urrielense",
                             "Androsace villosa", "Thymus praecox", "Plantago alpina", "Plantago holosteum",
                             "Conopodium majus", "Phyteuma hemisphaericum", "Jurinea humilis", "Solidago virgaurea", 
                             "Armeria cantabrica", "Armeria duriaei", "Spergula morisonii", "Minuartia recurva", 
                             "Arenaria erinacea", "Gypsophila repens", "Dianthus langeanus", "Silene ciliata", 
                             "Silene acaulis"))%>%
  ggplot(aes(x=Taxon, y= oil_type_PER, fill= saturation_type))+
  geom_bar(position = "fill", stat = "identity",  color = "black", show.legend = T)+
  scale_fill_manual (name= "FA types", values = c("brown", "chocolate1"),
                        guide = guide_legend (title.position = "top",direction = "horizontal", reverse = T)) +
  coord_flip()+
  ggthemes::theme_tufte(base_size=12) + 
  labs( tag = "E)", y= "Oil content (relative %)")+ #title= "UFA vs SFA",
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0, 0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"))-> fig3E;fig3E

# combine panels ####

ggpubr::ggarrange(fig3C, fig3D, fig3E,ncol =3, nrow= 1,common.legend = FALSE, widths = c(1.5,1,1),align = "h")->fig3CDE;fig3CDE

ggpubr::ggarrange(fig3AB, fig3CDE,ncol =1, nrow= 2,common.legend = FALSE,heights = c(1,1.5), align = "h")
