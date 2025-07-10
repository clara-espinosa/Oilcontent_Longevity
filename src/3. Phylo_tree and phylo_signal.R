library(tidyverse);library(V.PhyloMaker); library(scales);library(V.PhyloMaker)
library(phylosignal);library(phylobase);library(ape);library(tidytree)

# Phylo tree regional oil data set#####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/oil_regionaldata.csv", sep =",") %>% 
  dplyr::select (Taxon, family) %>%
  unique %>%
  separate(Taxon, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  arrange(species) %>%
  na.omit %>%
  select(species, genus, family)-> 
  ranks1

#devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
x11()
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/tree_regionaldata.tree")

# Phylogenetic signal with regional data using phylosignal package
ape::read.tree("results/tree_regionaldata.tree")-> tree
x11()
plot(tree)
read.csv("data/oil_regionaldata.csv", sep =",") %>% 
  group_by(Taxon, ecology)%>%
  summarise(oil_content = mean(oil_content),
            seed_mass = mean (X50seed_mass_mg))%>% 
  select(Taxon, oil_content, seed_mass, ecology)%>%
  mutate(label= gsub(" ","_", Taxon))%>%
  remove_rownames %>% 
  column_to_rownames(var="label") %>% 
  select(oil_content, seed_mass)-> oil_phylo_regional

obj_oil_regional <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo_regional), match.data=TRUE)
# calculate and test lambda phylogenetic signal
phyloSignal(obj_oil_regional)

# Phylo tree oil alpine data  #####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/species_header.csv", sep =",") %>%
  select (Taxon, family) %>%
  unique %>%
  separate(Taxon, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  arrange(species) %>%
  na.omit %>%
  select(species, genus, family)-> 
  ranks1

#devtools::install_github("jinyizju/V.PhyloMaker")
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
x11()
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/tree_alpinedata.tree")

# Phylogenetic signal of traits with alpine data using phylosignal package ####
ape::read.tree("results/tree_alpinedata.tree")-> tree
x11()
plot(tree)

read.csv("data/species_traits.csv")%>%
  group_by(Taxon)%>%
  summarise(oil.content = mean(oil.content),
            ratio = mean(ratio),
            seed.mass = mean(mass_50))%>%
  select(Taxon,oil.content, ratio, seed.mass)%>%
  mutate(label= gsub(" ","_", Taxon))%>%
  remove_rownames %>% 
  column_to_rownames(var="label") %>% 
  select(oil.content, ratio, seed.mass)-> oil_phylo_alpine 

obj_oil_alpine <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo_alpine), match.data=TRUE)
# calculation of lambda 
phyloSignal(obj_oil_alpine)

# phylogenetic signal for p50 trait (less species available)####
read.csv("data/species_header.csv", sep =",") %>%
  merge(read.csv("data/species_traits.csv"))%>%
  select (Taxon, family, p50) %>%
  na.omit()%>%
  select(Taxon, family)%>%
  unique %>%
  separate(Taxon, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  arrange(species) %>%
  na.omit %>%
  select(species, genus, family)-> 
  ranks1

#devtools::install_github("jinyizju/V.PhyloMaker")
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
x11()
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/tree_p50.tree")
ape::read.tree("results/tree_p50.tree")-> tree
x11()
plot(tree)

read.csv("data/species_traits.csv")%>%
  group_by(Taxon)%>%
  summarise(p50 = mean(p50),
            oil.content = mean(oil.content),
            ratio = mean(ratio))%>%
  select(Taxon, p50, oil.content, ratio)%>%
  mutate(label= gsub(" ","_", Taxon))%>%
  remove_rownames %>% 
  column_to_rownames(var="label") %>% 
  na.omit()%>%
  select(p50, oil.content, ratio)-> oil_phylo_p50

obj_oil_p50 <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo_p50), match.data=TRUE)
# calculation of lambda 
phyloSignal(obj_oil_p50)

# phylogenetic signal for EHS trait (less species available)####
read.csv("data/species_header.csv", sep =",") %>%
  merge(read.csv("data/species_traits.csv"))%>%
  select (Taxon, family, EHS_mean) %>%
  na.omit()%>%
  select(Taxon, family)%>%
  unique %>%
  separate(Taxon, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  arrange(species) %>%
  na.omit %>%
  select(species, genus, family)-> 
  ranks1

#devtools::install_github("jinyizju/V.PhyloMaker")
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
x11()
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/tree_EHS.tree")
ape::read.tree("results/tree_EHS.tree")-> tree
x11()
plot(tree)

read.csv("data/species_traits.csv")%>%
  group_by(Taxon)%>%
  summarise(EHS = mean(EHS_mean),
            oil.content = mean(oil.content),
            ratio = mean(ratio))%>%
  select(Taxon, EHS, oil.content, ratio)%>%
  mutate(label= gsub(" ","_", Taxon))%>%
  remove_rownames %>% 
  column_to_rownames(var="label") %>% 
  na.omit()%>%
  select(EHS, oil.content, ratio)-> oil_phylo_EHS

obj_oil_EHS <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo_EHS), match.data=TRUE)
# calculation of lambda 
phyloSignal(obj_oil_EHS)
