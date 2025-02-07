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

#https://www.francoiskeck.fr/phylosignal/demo_plots.html

obj_oil_regional <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo_regional), match.data=TRUE)

#visualization oil content and seed mass with phylogenetic tree
barplot.phylo4d(obj_oil_regional, tree.ratio = 0.2,  center=F, bar.col = rainbow(80), #bar.col = mat.col,  error.bar.sup = mat.e, 
                trait.bg.col = "white", show.box = F, grid.vertical = F,
                grid.horizontal = F, tip.labels = NULL, tree.ladderize =T) #rainbow(37)
barplot(obj_oil_regional ,tree.type = "fan", tip.cex = 0.6, tree.open.angle = 160, trait.cex = 0.6)
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
            seed.mass = mean(mass_50),
            p50 = mean(p50),
            EHS= mean(EHS_mean))%>%
  select(Taxon,oil.content, ratio, seed.mass, p50, EHS)%>%
  mutate(label= gsub(" ","_", Taxon))%>%
  remove_rownames %>% 
  column_to_rownames(var="label") %>% 
  select(oil.content, ratio, seed.mass, p50, EHS)-> oil_phylo_alpine 

#https://www.francoiskeck.fr/phylosignal/demo_plots.html

obj_oil_alpine <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo_alpine), match.data=TRUE)
# visualization of tree with traits
barplot.phylo4d(obj_oil_alpine, tree.ratio = 0.2,  center=F, bar.col = rainbow(47),#bar.col = mat.col,  error.bar.sup = mat.e
                trait.bg.col = "white", show.box = T, grid.vertical = F,
                grid.horizontal = F, tip.labels = NULL, tree.ladderize =T) #rainbow(37)
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

#https://www.francoiskeck.fr/phylosignal/demo_plots.html

obj_oil_p50 <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo_p50), match.data=TRUE)
# visualization of tree with traits
barplot.phylo4d(obj_oil_p50, tree.ratio = 0.2,  center=F, bar.col = rainbow(47),#bar.col = mat.col,  error.bar.sup = mat.e
                trait.bg.col = "white", show.box = T, grid.vertical = F,
                grid.horizontal = F, tip.labels = NULL, tree.ladderize =T) #rainbow(37)
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

#https://www.francoiskeck.fr/phylosignal/demo_plots.html

obj_oil_EHS <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo_EHS), match.data=TRUE)
# visualization of tree with traits
barplot.phylo4d(obj_oil_EHS, tree.ratio = 0.2,  center=F, bar.col = rainbow(47),#bar.col = mat.col,  error.bar.sup = mat.e
                trait.bg.col = "white", show.box = T, grid.vertical = F,
                grid.horizontal = F, tip.labels = NULL, tree.ladderize =T) #rainbow(37)
# calculation of lambda 
phyloSignal(obj_oil_EHS)
