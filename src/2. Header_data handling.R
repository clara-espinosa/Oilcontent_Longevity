library(tidyverse);library (rstatix);library (stringr)
library(ggpattern); library (vegan) ;library (ggrepel)

### read species oil content data ####
oil_data_full%>% # from format_oil_data script
  select(community, Taxon,family, ecology, oil.content, ratio)-> oil_data # 49 species
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

### read germ data (t50/EHS) ####
read.csv("data/t50_EHS_data.csv", sep= ";")->t50  #45 species with t50 trait
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
  left_join(p50,by= c("Taxon", "community") )%>%
  left_join(seedmass, by= c("Taxon", "community"))%>%
  left_join(t50, by= c("Taxon", "community", "family"))%>%
  left_join(sp_pref, by = c("Taxon", "community"))%>%
  dplyr::select(Taxon, community, family, ecology, 
                oil.content, ratio, mass_50, p50, T50_F, T50_S,EHS_F, EHS_S , GDD, FDD, Snw)%>%
  write.csv("results/species_traits_summary.csv")
