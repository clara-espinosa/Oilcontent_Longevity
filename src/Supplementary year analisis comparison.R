library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan);library(ggplot2)

# check potential differences between years analysis 
## COMPARISON between years of oil analysis NOT SIGNIFICANT DIFFERENCES ####
oil_alpine_data%>% # from script 1 format oil data
  dplyr::select(Taxon, community, family, oil.content, ratio, year_analisis)%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia, year_analisis) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content),
         Lratio=log(ratio))%>%
  as.data.frame()-> oil_year

str(oil_year)

MCMCglmm::MCMCglmm(Lratio ~ year_analisis , #  
                   random = ~ animal,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = oil_year,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> ratio_y
plot(oil_y)   
summary(oil_y) # No differences between years

plot(ratio_y)   
summary(ratio_y) # no differences between years

oil_alpine_data%>%
  dplyr::select(Taxon, community, family, oil.content, ratio, year_analisis)%>%
  convert_as_factor(Taxon, year_analisis) %>%
  gather(trait, value, oil.content:ratio)%>%
  ggplot()+
  geom_boxplot(aes(x= year_analisis, y = value, fill = year_analisis))+
  facet_grid (~trait)+
  scale_fill_manual(name = "Year of analisis", labels = c("2022 n = 5", "2023 n = 31", "2024 n = 13"), 
                    values = c("salmon", "lightgreen", "lightblue"))+
  theme_classic(base_size = 10)

oil_alpine_data%>%
  dplyr::select(Taxon, community, family, oil.content, ratio, year_analisis)%>%
  group_by(year_analisis)%>%
  tally()
#### PCA EXPLORATION 
# check correlation of FAMEs with <3% relative proportion
oil_alpine_data%>%
  dplyr::select(C16.0, C18.0, C18.1n7c, C18.1n9c, C18.2n6c, C18.3n3, 
                C18.3n6, C20.1n9,C20.3n6, C20.4n6, C22.1n9, C24.1n9,
                oil.content, ratio) %>% cor()#FAME >3% relative proportion

# looks like 2 groups of variables highly correlated
# 1. C20.4n6, C20.3n6  >0.8 keep C20.4n6 higher contribution
# 2. C22.1n9, C24.1n9  >0.8 keep C22.1n9 higher contribution

oil_alpine_data%>%
  dplyr::select(C16.0, C18.0, C18.1n7c, C18.1n9c, C18.2n6c, C18.3n3, 
                C18.3n6, C20.1n9, C20.4n6, C22.1n9,
                oil.content, ratio)%>% #FAME >3% relative proportion
  FactoMineR::PCA() -> pca_oil
x11()
pca_oil$var$contrib
pca_oil$eig

oil_alpine_data%>%
  dplyr::select(C16.0, C18.0, C18.1n7c, C18.1n9c, C18.2n6c, C18.3n3, 
                C18.3n6, C20.1n9, C20.4n6, C22.1n9, 
                oil.content, ratio) %>% cor()

cbind(oil_alpine_data[, 1:5], data.frame(pca_oil$ind$coord[, 1:2])) %>%
  mutate(species = factor(Taxon))%>%
  mutate(family = factor(family))-> pcaInds

pca_oil$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") -> pcaVars

ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  #coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(fill = as.factor(year_analisis)), color = "black", show.legend = T, size = 4, shape = 21) + # family
  geom_text_repel (data = pcaInds, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4, max.overlaps = 5) +
  ggthemes::theme_tufte(base_size=12) + 
  #scale_fill_manual(values=col)+
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

## COMPARISON between years of longevity analysis NOT SIGNIFICANT DIFFERENCES ####
# longevity (p50)####
# read data
read.csv("data/genstat all.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  dplyr::select(Taxon,  family, community, year, p50)%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia, year) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID) %>%
  mutate(Lp50 = log(p50), 
         SRp50 = sqrt(p50), )%>%
  na.omit()%>%
  as.data.frame()%>%
  droplevels()-> p50_year #n=41
str(p50_year)
unique(p50_year$familia) # 16 families
unique(p50_year$Taxon) #33 levels (silene ciliata and thymus in both communities

hist(p50_year$p50) # quite normally distributed
hist(p50_year$SRp50) 
# model
MCMCglmm::MCMCglmm(p50 ~ year , # 
                   random = ~ animal,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = p50_year,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> p50_y

x11()
plot(p50_y)
summary(p50_y)

p50_year%>%
  ggplot()+
  geom_boxplot(aes(x= year, y = p50, fill = year))+
  scale_fill_manual(name = "Year of analisis", labels = c("2022 n = 25", "2024 n = 11"), 
                    values = c("salmon",  "lightblue"))+
  theme_classic(base_size = 10)
