library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan)

#### OWN oil data handling #####
# oil data mg oil/g sample transformation to long format
read.csv("data/oil_data.csv", sep= ";")%>%
  #dplyr:: select(Taxon:C24.0)%>%
  #mutate(species=make.cepnames(species))%>%# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
  convert_as_factor(Taxon,  family, community)%>%
  filter(units == "mg")%>%
  gather(oil_type, value, 7:32) %>%
  group_by(Taxon, community) %>%
  summarise(oil.content= (sum(value)/10))%>% # transformation mg of oil/ gr of sample to % of oil content 
  merge(read.csv("data/oil_data.csv", sep= ";"), by = c("Taxon", "community"))%>%
  filter(units=="percentage")%>%
  dplyr::select(Taxon,  family, community, oil.content, C12.0:C24.1n9) -> oil_data_per


read.csv("data/oil_data.csv", sep= ";")%>%  
  filter(units == "percentage")%>%
  gather(oil_type, oil_PER, 7:32)%>%
  merge(read.csv("data/FA_types.csv"))%>%
  dplyr::select(Taxon, community,oil_PER, saturation_type)%>%
  group_by(Taxon, community, saturation_type)%>%
  summarise(oil_proportion = sum(oil_PER))%>%
  spread(saturation_type, oil_proportion)%>%
  mutate(ratio=UFA/SFA)%>%
  merge(oil_data_per, by = c("Taxon", "community"))%>%
  dplyr::select(Taxon,  family, community, oil.content, ratio, UFA, SFA, C12.0:C24.1n9)%>%
  merge(read.csv("data/species_oil.csv"), by = c("Taxon", "community", "family"))%>%
  dplyr::select(Taxon,  family, community, year_analisis, ecology, oil.content, ratio, UFA, SFA, C12.0:C24.1n9)-> oil_data_full


#### PCA EXPLORATION ####
# check correlation of FAMEs with <3% relative proportion
oil_data_full%>%
  dplyr::select(C16.0, C18.0, C18.1n7c, C18.1n9c, C18.2n6c, C18.3n3, 
                C18.3n6, C20.1n9,C20.3n6, C20.4n6, C22.1n9, C24.1n9,
                oil.content, ratio) %>% cor()#FAME >3% relative proportion

# looks like 2 groups of variables highly correlated
# 1. C20.4n6, C20.3n6  >0.8 keep C20.4n6 higher contribution
# 2. C22.1n9, C24.1n9  >0.8 keep C22.1n9 higher contribution

oil_data_full%>%
  dplyr::select(C16.0, C18.0, C18.1n7c, C18.1n9c, C18.2n6c, C18.3n3, 
                C18.3n6, C20.1n9, C20.4n6, C22.1n9,
                oil.content, ratio)%>% #FAME >3% relative proportion
  FactoMineR::PCA() -> pca_oil
x11()
pca_oil$var$contrib
pca_oil$eig

oil_data_full%>%
  dplyr::select(C16.0, C18.0, C18.1n7c, C18.1n9c, C18.2n6c, C18.3n3, 
                C18.3n6, C20.1n9, C20.4n6, C22.1n9, 
                oil.content, ratio) %>% cor()

cbind(oil_data_full[, 1:5], data.frame(pca_oil$ind$coord[, 1:2])) %>%
  mutate(species = factor(Taxon))%>%
  mutate(family = factor(family))-> pcaInds

pca_oil$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") -> pcaVars

#order like family fct relevel
poales <- c("#78E208", "#45A747", "#396F3E")
rosids <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220","#fcb622", "darkgoldenrod3",  "#f68b08","#ff7125" ,"#c87107",  "#8E5005") #, 
asterids <- c("#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f")
col <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220","#fcb622", "darkgoldenrod3", "#f68b08", "#ff7125" ,"#c87107",  "#8E5005",
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
  geom_text_repel (data = pcaInds, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4, max.overlaps = 5) +
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
                                  "% variance explained)", sep = ""))
