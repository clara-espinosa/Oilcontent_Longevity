library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan)

#### OWN oil data handling #####
# oil data mg oil/g sample transformation to long format
read.csv("data/oil_alpinedata.csv", sep= ",")%>%
  #dplyr:: select(Taxon:C24.0)%>%
  #mutate(species=make.cepnames(species))%>%# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
  convert_as_factor(Taxon,community)%>%
  filter(units == "mg")%>%
  gather(oil_type, value, 7:30) %>%
  group_by(Taxon, community) %>%
  summarise(oil.content= (sum(value)/10))%>% # transformation mg of oil/ gr of sample to % of oil content 
  merge(read.csv("data/oil_alpinedata.csv", sep= ","), by = c("Taxon", "community"))%>%
  filter(units=="percentage")%>%
  dplyr::select(Taxon, community, oil.content, C12.0:C24.1n9) -> oil_data_per # df with oil content and FA types percentages

read.csv("data/oil_alpinedata.csv", sep= ",")%>%  
  filter(units == "percentage")%>%
  gather(oil_type, oil_PER, 5:30)%>%
  merge(read.csv("data/FA_types.csv"))%>%
  dplyr::select(Taxon, community,oil_PER, saturation_type)%>%
  group_by(Taxon, community, saturation_type)%>%
  summarise(oil_proportion = sum(oil_PER))%>%
  spread(saturation_type, oil_proportion)%>%
  mutate(ratio=UFA/SFA)%>% # calculation of saturated and unsaturated fatty acids, and the ratio between them.
  merge(oil_data_per, by = c("Taxon", "community"))%>%
  dplyr::select(Taxon, community, oil.content, ratio, UFA, SFA, C12.0:C24.1n9)%>%
  merge(read.csv("data/species_header.csv"), by = c("Taxon", "community"))%>%
  dplyr::select(Taxon,  family, order, community, year_analisis, ecology, oil.content, ratio, UFA, SFA, C12.0:C24.1n9)-> oil_alpine_data # data frame with all oil data together in wide format
