# Figure 1 oil data exploration
library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan);library(ggpubr)

# descriptive statistics
View (oil_data_full)
oil_data_full%>%

#### PANEL A) Total oil content (%) per species ####
#order like family fct relevel
poales <- c("#78E208", "#45A747", "#396F3E")
rosids <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220","#fcb622", "darkgoldenrod3", "#f68b08","#ff7125" ,"#c87107",  "#8E5005") #, 
asterids <- c("#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f")
col <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220", "#fcb622","darkgoldenrod3", "#f68b08", "#ff7125" ,"#c87107",  "#8E5005",
         "#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f",
         "#78E208", "#45A747", "#396F3E")
read.csv("data/species_oil.csv")%>%
  merge (oil_data_full)%>%
  mutate(Taxon = as.factor(Taxon))%>%
  mutate(Taxon = fct_relevel(Taxon,"Luzula caespitosa" , "Kobresia myosuroides", "Carex sempervirens",
                             "Koeleria vallesiana","Sesleria caerulea","Helictochloa marginata", "Festuca glacialis",  
                             "Festuca summilusitana", "Avenella flexuosa","Neoschischkinia truncatula", 
                             "Saxifraga paniculata", "Saxifraga oppositifolia" ,"Saxifraga conifera",
                             "Sempervivum arachnoideum","Sedum brevifolium", "Sedum anglicum", 
                             "Helianthemum canum", "Helianthemum urrielense",
                             "Teesdalia conferta","Iberis carnosa", "Salix breviserrata","Anthyllis vulneraria",  
                             "Armeria cantabrica", "Armeria duriaei", "Spergula morisonii",
                             "Minuartia recurva", "Minuartia verna", "Cerastium ramosissimum",
                             "Arenaria erinacea", "Gypsophila repens", "Dianthus langeanus", 
                             "Silene ciliata", "Silene acaulis","Androsace villosa",
                             "Gentiana verna","Plantago alpina", "Plantago holosteum",
                             "Veronica nummularia", "Pedicularis pyrenaica", "Thymus praecox", 
                             "Conopodium majus", "Dethawia splendens", "Jasione cavanillesii", 
                             "Phyteuma hemisphaericum", "Jurinea humilis", "Solidago virgaurea", 
                             "Phalacrocarpum oppositifolium"))%>% 
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
  "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
  "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
    "Crassulaceae", "Saxifragaceae", "Poaceae", 
   "Cyperaceae", "Juncaceae"))%>% 
  group_by(Taxon, family)%>%
  summarise(oil.content = mean(oil.content))%>%
  ggplot(aes(x=Taxon, y= oil.content, fill=family ) )+ #family
  #geom_point(size = 4, shape=21, color = "black", show.legend = F)+
  geom_bar(stat = "identity", show.legend = F, color="black")+ #problem with Thymus and S. ciliata (2 accesions summ % oil)
  coord_flip()+
  scale_fill_manual (values=col)+
  ggthemes::theme_tufte(base_size=12) + 
  labs( tag = "A)", y= "Oil content (%)")+ #title= "Species oil content (%)",
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        legend.position = "", 
        plot.background = element_blank(),
        plot.margin = unit(c(0, 0.2,0,0), "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))->fig3A;fig3A

#### PANEL B) oil types percentatges per sp (% specific oil /total oil ID) ####
library(viridis)
library(scales)

show_col(viridis_pal()(6))

read.csv("data/species_oil.csv")%>%
  merge(oil_data_full)%>%
  dplyr::select(Taxon, family, C16.0, C18.2n6c, C18.1n9c, C18.3n3, C22.1n9, C18.3n6)%>%
  gather(oil_type, oil_PER, C16.0:C18.3n6)%>%
  group_by(Taxon, family, oil_type)%>%
  summarise(oil_PER = mean(oil_PER)) %>%
  mutate(Taxon = as.factor(Taxon))%>%
  mutate(Taxon = fct_relevel(Taxon,"Luzula caespitosa" , "Kobresia myosuroides", "Carex sempervirens",
                             "Koeleria vallesiana","Sesleria caerulea","Helictochloa marginata", "Festuca glacialis",  
                             "Festuca summilusitana", "Avenella flexuosa","Neoschischkinia truncatula", 
                             "Saxifraga paniculata", "Saxifraga oppositifolia" ,"Saxifraga conifera",
                             "Sempervivum arachnoideum","Sedum brevifolium", "Sedum anglicum", 
                             "Helianthemum canum", "Helianthemum urrielense",
                             "Teesdalia conferta","Iberis carnosa", "Salix breviserrata","Anthyllis vulneraria",  
                             "Armeria cantabrica", "Armeria duriaei", "Spergula morisonii",
                             "Minuartia recurva", "Minuartia verna", "Cerastium ramosissimum",
                             "Arenaria erinacea", "Gypsophila repens", "Dianthus langeanus", 
                             "Silene ciliata", "Silene acaulis","Androsace villosa",
                             "Gentiana verna","Plantago alpina", "Plantago holosteum",
                             "Veronica nummularia", "Pedicularis pyrenaica", "Thymus praecox", 
                             "Conopodium majus", "Dethawia splendens", "Jasione cavanillesii", 
                             "Phyteuma hemisphaericum", "Jurinea humilis", "Solidago virgaurea", 
                             "Phalacrocarpum oppositifolium"))%>% 
  group_by(Taxon)%>%
  mutate(other= (100-sum(oil_PER)))%>%
  spread(oil_type, oil_PER)%>%
  gather(oil_type, oil_PER, other:C22.1n9)%>%
  #filter (oil_PER>10) %>% # if filter above 3% change geom_bar position = fill
  #mutate(oil_type= ifelse(oil_mean<10, "other",oil_type ))%>%
  #group_by(Taxon, oil_type)%>%
  #summarise(oil_PER= sum(oil_PER))%>%
  ggplot(aes(x=Taxon, y= oil_PER, fill= oil_type))+
  geom_bar(position = position_stack(reverse=T), stat = "identity",  color = "black", show.legend = T)+ #position = "stack", 
  coord_flip()+
  scale_fill_manual(name= "", values = c("#440154FF", "#414497FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF", "lightgrey"), 
                    guide = guide_legend (title.position = "top",direction = "horizontal")) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "B)", y= "Oil content (relative %)")+# title= "FA types  (>3%)", 
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.key.size = unit(1,"line"),
        legend.text = element_text(size = 10, color = "black"),
        legend.margin = margin(t = 0, unit='cm'),
        legend.box.margin = unit(c(0, 0,0,0), "cm"),
        plot.margin = unit(c(0, 0.2,0,0), "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black")) -> fig3B;fig3B

### PANEL C) unsaturated vs saturated fatty acids composition per species ####
read.csv("data/species_oil.csv")%>%
  merge(oil_data_full)%>%
  dplyr::select(Taxon, family, UFA, SFA)%>%
  gather(saturation_type, oil_PER, UFA:SFA)%>%
  group_by(Taxon, family, saturation_type)%>%
  summarise(oil_PER = mean(oil_PER)) %>%
  mutate(Taxon = as.factor(Taxon))%>%
  mutate(Taxon = fct_relevel(Taxon,"Luzula caespitosa" , "Kobresia myosuroides", "Carex sempervirens",
                             "Koeleria vallesiana","Sesleria caerulea","Helictochloa marginata", "Festuca glacialis",  
                             "Festuca summilusitana", "Avenella flexuosa","Neoschischkinia truncatula", 
                             "Saxifraga paniculata", "Saxifraga oppositifolia" ,"Saxifraga conifera",
                             "Sempervivum arachnoideum","Sedum brevifolium", "Sedum anglicum", 
                             "Helianthemum canum", "Helianthemum urrielense",
                             "Teesdalia conferta","Iberis carnosa", "Salix breviserrata","Anthyllis vulneraria",  
                             "Armeria cantabrica", "Armeria duriaei", "Spergula morisonii",
                             "Minuartia recurva", "Minuartia verna", "Cerastium ramosissimum",
                             "Arenaria erinacea", "Gypsophila repens", "Dianthus langeanus", 
                             "Silene ciliata", "Silene acaulis","Androsace villosa",
                             "Gentiana verna","Plantago alpina", "Plantago holosteum",
                             "Veronica nummularia", "Pedicularis pyrenaica", "Thymus praecox", 
                             "Conopodium majus", "Dethawia splendens", "Jasione cavanillesii", 
                             "Phyteuma hemisphaericum", "Jurinea humilis", "Solidago virgaurea", 
                             "Phalacrocarpum oppositifolium"))%>% 
  ggplot(aes(x=Taxon, y= oil_PER, fill= saturation_type))+
  geom_bar(position = "fill", stat = "identity",  color = "black", show.legend = T)+
  scale_fill_manual (name= "FA types", values = c("brown1", "brown4"),
                        guide = guide_legend (title.position = "top",direction = "horizontal", reverse = T)) +
  coord_flip()+
  ggthemes::theme_tufte(base_size=12) + 
  labs( tag = "C)", y= "Oil content (relative %)")+ #title= "UFA vs SFA",
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        legend.position = "bottom", 
        legend.margin = margin(t = 0, unit='cm'),
        legend.box.margin = unit(c(0, 0,0,0), "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0, 0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"))-> fig3C;fig3C

# PCA for Panel D and E : PCA variables + species #####
# JOIN FAME composition (percentages relative to the total content of oil) merge with total oil for PCA
str(oil_data_full)
oil_data_full[, 4:33]%>%
  FactoMineR::PCA() -> pca_oil

pca_oil$var$contrib
pca_oil$eig
oil_data_full[, 4:33]%>%  cor() # remove C12:0, C20:1n9, C20:2n6 (hihgly correlated with other variables with higher contributions to axes)

# PCa with only the FA with more than 3% relative proportion 
oil_data_full%>%
  dplyr::select(C16.0, C18.0, C18.1n7c, C18.1n9c, C18.2n6c, C18.3n3, 
                C18.3n6, C20.1n9,C20.3n6, C20.4n6, C22.1n9, C24.1n9,
                oil.content, ratio)%>% #FAME >3% relative proportion
  FactoMineR::PCA() -> pca_oil
x11()
pca_oil$var$contrib
pca_oil$eig
oil_data_full%>%
  dplyr::select(C16.0, C18.0, C18.1n7c, C18.1n9c, C18.2n6c, C18.3n3, 
                C18.3n6, C20.1n9,C20.3n6, C20.4n6, C22.1n9, C24.1n9,
                oil.content, ratio) %>% cor()

# looks like 2 groups of variables highly correlated
# 1. C20.4n6, C20.3n6  >0.8
# 2. C22.1n9, C24.1n9  >0.8
 

cbind(oil_data_full, data.frame(pca_oil$ind$coord[, 1:2])) %>%
  mutate(species = factor(Taxon))%>%
  mutate(family = factor(family))-> pcaInds

pca_oil$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") -> pcaVars


### PANEL D) Plot PCA species#########
#order like family fct relevel
poales <- c("#78E208", "#45A747", "#396F3E")
rosids <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220","#fcb622", "darkgoldenrod3",  "#f68b08","#ff7125" ,"#c87107",  "#8E5005") #, 
asterids <- c("#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f")
col <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220","#fcb622", "darkgoldenrod3", "#f68b08", "#ff7125" ,,"#c87107",  "#8E5005",
         "#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f",
         "#78E208", "#45A747", "#396F3E")
x11()
pcaInds%>%
  as.data.frame()%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))->pcaInds 
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  #coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(fill = family), color = "black", show.legend = T, size = 4, shape = 21) + # family
  #geom_text_repel (data = pcaInds, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4, max.overlaps = 5) +
  ggthemes::theme_tufte(base_size=12) + 
  scale_fill_manual(values=col)+
  guides(fill=guide_legend(ncol=1, keywidth=0.1,keyheight=0.1,default.unit="cm")) +
  labs( tag = "D)")+#title= "PCA species",
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
                                  "% variance explained)", sep = "")) -> fig3D;fig3D#, limits = c(-4, 4)
### PANEL E) Plot PCA variables#########
x11()
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  #coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2), arrow = arrow()) +
  #geom_label(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 4) +
  geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 3) +
  #geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 4, segment.size= 1,
   #point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  ggthemes::theme_tufte(base_size=12) + 
  labs( tag = "E)")+#title= "PCA variables",
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
                                  "% variance explained)", sep = ""))-> fig3E;fig3E #, limits = c(-4, 4)

# combine D + E panels ###
library(patchwork)
fig3D + fig3E + 
  plot_layout(guides = 'auto')-> fig3DE;fig3DE


# combine panels ####

ggpubr::ggarrange(fig3A, fig3B, fig3C,ncol =3, nrow= 1,common.legend = FALSE, widths = c(1.5,1,1),align = "h")->fig3ABC;fig3ABC

ggpubr::ggarrange(fig3ABC, fig3DE,ncol =1, nrow= 2,common.legend = FALSE,heights = c(1.5,1), align = "h")->fig3;fig3


ggsave(filename = "oil exploration.png", plot =fig3 , path = "results/figures", 
       device = "png", dpi = 600)
