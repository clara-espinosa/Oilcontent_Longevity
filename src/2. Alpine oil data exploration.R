library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan);library(ggpubr)
library(patchwork):library(paletteer)

# descriptive statistics
View (oil_alpine_data) # from script 1 format_oil_data

read.csv("data/oil_alpinedata.csv", sep= ",")%>%  
  filter(units == "percentage")%>%
  gather(oil_type, oil_PER, 5:30)%>%
  group_by (oil_type)%>%
  get_summary_stats(oil_PER)%>%
  dplyr::select(oil_type, min, max, mean, se)%>%
  rename(Fa_types = oil_type)%>%
  write.csv("results/supplementary/Table S2. FA types summary.csv", row.names = F)

read.csv("data/oil_alpinedata.csv", sep= ",")%>%  
  filter(units == "percentage")%>%
  dplyr:: select (Taxon,C18.2n6c, C18.1n9c, C18.3n3, C16.0)%>%
  gather(oil_type, oil_PER, 2:5)%>%
  group_by(Taxon)%>%
  summarise(oil_PER= sum(oil_PER))%>%
  get_summary_stats()

oil_alpine_data%>%
  dplyr::select(Taxon, family, ratio, UFA, SFA)%>%
  group_by(family)%>%
  get_summary_stats(ratio)


#### PANEL A) Total oil content (%) per species ####
# order like order fct relevel 
col_order <- c( "#8E5005","darkgoldenrod3", "gold1", "orange","darkorange2", "orangered1" ,
         "skyblue1","deepskyblue",  "#2d7faf", "#275381", "#42346f",
         "#78E208")

x11()
read.csv("data/species_header.csv")%>%
  merge (oil_alpine_data )%>%
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
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales",
                              "Gentianales" ,"Ericales","Caryophyllales",  
                              "Fabales","Malpighiales","Brassicales", "Malvales", 
                              "Saxifragales", "Poales"))%>% 
  group_by(Taxon, family, order)%>%
  summarise(oil.content = mean(oil.content))%>%
  ggplot(aes(x=Taxon, y= oil.content, fill=order) )+ #family
  geom_bar(stat = "identity", show.legend = F, color="black")+ 
  coord_flip()+
  scale_fill_manual (values=col_order)+
  theme_classic(base_size=12) + 
  labs( tag = "(a)", y= "Oil content (%)")+ 
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        plot.tag = element_text(face="bold"),
        legend.position = "", 
        plot.background = element_blank(),
        plot.margin = unit(c(0, 0.2,0,0), "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black", face= "italic"))->fig2A;fig2A

#### PANEL B) oil types percentatges per sp (% specific oil /total oil ID) ####
paletteer_d("colorBlindness::Blue2DarkOrange12Steps")
read.csv("data/species_header.csv")%>%
  merge(oil_alpine_data)%>%
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
  ggplot(aes(x=Taxon, y= oil_PER, fill= oil_type))+
  geom_bar(position = position_stack(reverse=T), stat = "identity",  color = "black", show.legend = T)+ #position = "stack", 
  coord_flip()+
  scale_fill_manual(name= "", values = c("#1E8E99FF", "#51C3CCFF", "#99F9FFFF", "#CCFEFFFF", "#FFE5CCFF", "#FFCA99FF", "#FFAD65FF"), 
                    guide = guide_legend (title.position = "top",direction = "horizontal")) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "(b)", y= "Oil content (relative %)")+# title= "FA types  (>3%)", 
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        legend.position = "bottom", 
        legend.title = element_blank(),
        plot.tag = element_text(face="bold"),
        legend.key.size = unit(1,"line"),
        legend.spacing.x = unit(0.1, "lines"),
        legend.text = element_text(size = 9, color = "black"),
        legend.margin = margin(t = 0, unit='cm'),
        legend.box.margin = unit(c(0, 0,0,0), "cm"),
        plot.margin = unit(c(0, 0.2,0,0), "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black")) -> fig2B;fig2B

### PANEL C) unsaturated vs saturated fatty acids composition per species ####
read.csv("data/species_header.csv")%>%
  merge(oil_alpine_data )%>%
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
                        guide = guide_legend (title.position = "top",direction = "horizontal", nrow =2, reverse = T)) +
  coord_flip()+
  ggthemes::theme_tufte(base_size=12) + 
  labs( tag = "(c)", y= "Oil content (relative %)")+ #title= "UFA vs SFA",
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        plot.tag = element_text(face="bold"),
        legend.position = "bottom", 
        legend.key.size = unit(1,"line"),
        legend.margin = margin(t = 0, unit='cm'),
        legend.box.margin = unit(c(0, 0,0,0), "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 9, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0, 0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"))-> fig2C;fig2C

# PCA for Panel D and E : PCA variables + species #####
str(oil_alpine_data)
oil_alpine_data [, 7:36]%>%  cor() 
# remove C12:0, C20:1n9, C20:2n6 (hihgly correlated with other variables with higher contributions to axes)
oil_alpine_data %>%
  dplyr::select(C16.0, C18.0, C18.1n7c, C18.1n9c, C18.2n6c, C18.3n3, 
                C18.3n6, C20.1n9, C20.4n6, C20.3n6, C22.1n9, C24.1n9,
                oil.content, ratio) %>% cor()
# looks like 2 groups of variables highly correlated
# 1. C20.4n6, C20.3n6  >0.8
# 2. C22.1n9, C24.1n9  >0.8

# PCa with only the FA with more than 3% relative proportion and less than 0.8 correlation
oil_alpine_data %>%
  dplyr::select(C16.0, C18.0, C18.1n7c, C18.1n9c, C18.2n6c, C18.3n3, 
                C18.3n6, C20.1n9, C20.4n6, C22.1n9, 
                oil.content, ratio)%>% 
  FactoMineR::PCA() -> pca_oil
x11()
pca_oil$var$contrib
pca_oil$eig

cbind(oil_alpine_data, data.frame(pca_oil$ind$coord[, 1:2])) %>%
  mutate(species = factor(Taxon))%>%
  mutate(order= factor(order))%>%
  mutate(family = factor(family))-> pcaInds

pca_oil$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") -> pcaVars


### PANEL D) Plot PCA species#########

col_order <- c( "#8E5005","darkgoldenrod3", "gold1", "orange","darkorange2", "orangered1" ,
                "skyblue1","deepskyblue",  "#2d7faf", "#275381", "#42346f",
                "#78E208")
x11()
pcaInds%>%
  as.data.frame()%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales",
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", 
                             "Saxifragales", "Poales")) ->pcaInds 
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(fill = order), color = "black", show.legend = T, size = 4, shape = 21) + # family
  ggthemes::theme_tufte(base_size=12) + 
  scale_fill_manual(values=col_order)+
  guides(fill=guide_legend(ncol=1, keywidth=0.1,keyheight=0.1,default.unit="cm")) +
  labs( tag = "(d)")+#title= "PCA species",
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        plot.tag = element_text(face="bold"),
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
                                  "% variance explained)", sep = "")) -> fig2D;fig2D#, limits = c(-4, 4)

### PANEL E) Plot PCA variables#########
x11()
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2), arrow = arrow()) +
  geom_label_repel(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 3) +
  ggthemes::theme_tufte(base_size=12) + 
  labs( tag = "(e)")+#title= "PCA variables",
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        plot.tag = element_text(face="bold"),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        plot.margin = unit(c(0, 0,0,0), "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10, color = "black"))+
  scale_x_continuous(name = paste("Axis 1 (", round(pca_oil$eig[1, 2], 0),
                                  "% variance explained)", sep = "") ) + #,limits = c(-5, 5)
  scale_y_continuous(name = paste("Axis 2 (", round(pca_oil$eig[2, 2], 0), 
                                  "% variance explained)", sep = ""))-> fig2E;fig2E #, limits = c(-4, 4)


# combine panels ####

# combine A+B+C panels ###
fig2A + fig2B + fig2C+
  plot_layout(widths = c(1.5,1,1))->fig2ABC;fig2ABC

# combine D + E panels ###

fig2D + fig2E + 
  plot_layout()-> fig2DE;fig2DE

ggpubr::ggarrange(fig2ABC, fig2DE,ncol =1, nrow= 2,common.legend = FALSE,heights = c(1.8,1))->fig2;fig2

ggsave(filename = "fig 2. oil exploration.png", plot =fig2 , path = "results/figures", 
       device = "png", dpi = 600)
