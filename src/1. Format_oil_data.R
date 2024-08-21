library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan)

#### oil data handling #####
# oil data mg oil/g sample transformation to long format
read.csv("data/oil_data.csv")%>%
  #dplyr:: select(Taxon:C24.0)%>%
  #mutate(species=make.cepnames(species))%>%# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
  convert_as_factor(Taxon,  family, community)%>%
  filter(units == "mg")%>%
  gather(oil_type, value, 7:26) %>%
  group_by(Taxon, community) %>%
  summarise(oil.content= (sum(value)/10))%>% # transformation mg of oil/ gr of sample to % of oil content 
  merge(read.csv("data/oil_data.csv"), by = c("Taxon", "community"))%>%
  filter(units=="percentage")%>%
  dplyr::select(Taxon,  family, community, oil.content, C12.0:C24.0) -> oil_data_pca


read.csv("data/oil_data.csv")%>%  
  filter(units == "percentage")%>%
  gather(oil_type, oil_PER, 7:26)%>%
  merge(read.csv("data/FA_types.csv"))%>%
  dplyr::select(Taxon, community,oil_PER, saturation_type)%>%
  group_by(Taxon, community, saturation_type)%>%
  summarise(oil_proportion = sum(oil_PER))%>%
  spread(saturation_type, oil_proportion)%>%
  mutate(ratio=UFA/SFA)%>%
  merge(oil_data_pca, by = c("Taxon", "community"))%>%
  dplyr::select(Taxon,  family, community, oil.content, ratio, UFA, SFA, C12.0:C24.0)-> oil_data_full


