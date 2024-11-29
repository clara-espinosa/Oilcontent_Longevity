library(tidyverse);library(V.PhyloMaker); library(scales);library(V.PhyloMaker)
library(phylosignal);library(phylobase);library(ape);library(tidytree)

### Phylo tree full oil data set#####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/oil_fulldataset.csv", sep =";") %>% 
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

write.tree(tree$scenario.3, file = "results/tree_oil_fulldataset.tree")

# phylosignal tree graph 
ape::read.tree("results/tree_oil_fulldataset.tree")-> tree
x11()
plot(tree)
read.csv("data/oil_fulldataset.csv", sep =";") %>% 
  group_by(Taxon, ecology)%>%
  summarise(oil_content = mean(oil_content),
            seed_mass = mean (X50seed_mass_mg))%>% 
  select(Taxon, oil_content, seed_mass, ecology)%>%
  mutate(label= gsub(" ","_", Taxon))%>%
  remove_rownames %>% 
  column_to_rownames(var="label") %>% 
  select(oil_content, seed_mass)-> oil_phylo 

#https://www.francoiskeck.fr/phylosignal/demo_plots.html

obj_oil <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo), match.data=TRUE)
#mat.col <- ifelse(tdata(obj_oil , "tip") < 1,  "darkred","darkgrey")
#mat.e <- matrix(abs(rnorm(19 * 3, 0, 0.5)), ncol = 3,
# dimnames = list(tipLabels(p4d), names(tdata(p4d))))
barplot.phylo4d(obj_oil, tree.ratio = 0.2,  center=F, bar.col = rainbow(80), #bar.col = mat.col,  error.bar.sup = mat.e, 
                trait.bg.col = "white", show.box = F, grid.vertical = F,
                grid.horizontal = F, tip.labels = NULL, tree.ladderize =T) #rainbow(37)

barplot(obj_oil ,tree.type = "fan", tip.cex = 0.6, tree.open.angle = 160, trait.cex = 0.6)
dotplot.phylo4d (obj_oil , tree.type= "cladogram")
gridplot.phylo4d (obj_oil ,  tip.cex = 1.5, show.trait = T) #tree.type = "fan",
phyloSignal(obj_oil )
phyloCorrelogram (obj_oil )


### Phylo tree oil own data  #####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/species_oil.csv", sep =",") %>%
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

write.tree(tree$scenario.3, file = "results/tree_oil.tree")

# phylosignal tree graph ###
ape::read.tree("results/tree_oil.tree")-> tree
x11()
plot(tree)

read.csv("data/species_traits_summary.csv")%>%
  group_by(Taxon, ecology)%>%
  summarise(oil.content = mean(oil.content),
            ratio = mean(ratio),
            seed.mass = mean(mass_50))%>%
  select(Taxon,ecology, oil.content, ratio, seed.mass)%>%
  mutate(label= gsub(" ","_", Taxon))%>%
  remove_rownames %>% 
  column_to_rownames(var="label") %>% 
  select(oil.content, ratio, seed.mass)-> oil_phylo2 

#https://www.francoiskeck.fr/phylosignal/demo_plots.html

obj_oil2 <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo2), match.data=TRUE)
#mat.col <- ifelse(tdata(obj_oil , "tip") < 1,  "darkred","darkgrey")
#mat.e <- matrix(abs(rnorm(19 * 3, 0, 0.5)), ncol = 3,
               # dimnames = list(tipLabels(p4d), names(tdata(p4d))))
barplot.phylo4d(obj_oil2, tree.ratio = 0.2,  center=F, bar.col = rainbow(47),#bar.col = mat.col,  error.bar.sup = mat.e
                trait.bg.col = "white", show.box = T, grid.vertical = F,
                grid.horizontal = F, tip.labels = NULL, tree.ladderize =T) #rainbow(37)

barplot(obj_oil2 ,tree.type = "fan", tip.cex = 0.6, tree.open.angle = 160, trait.cex = 0.6)
dotplot.phylo4d (obj_oil2 , tree.type= "cladogram")
gridplot.phylo4d (obj_oil2 ,  tip.cex = 1.5, show.trait = T) #tree.type = "fan",
phyloSignal(obj_oil2 )
phyloCorrelogram (obj_oil2 )

### Phylo tree only species with longevity data #####

read.csv("data/longevity/species.csv", sep =",") %>%
  separate(Taxon, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  #filter (community == "Mediterranean")%>%
  select(species, genus, family) %>%
  unique %>%
  #mutate(family() = fct_recode(family, 
  #                            "Asteraceae" = "Compositae",
  #                           "Fabaceae" = "Leguminosae",
  #                          "Asphodelaceae" = "Xanthorrhoeaceae")) %>%
  arrange(species) %>%
  na.omit -> 
  ranks1
#devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)

phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/treelongevity.tree")
