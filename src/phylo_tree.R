library(tidyverse);library(V.PhyloMaker)
library(phylosignal);library(phylobase);library(ape);library(tidytree)
library (tidyverse); library(scales)

### Phylo tree both communitites #####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/species.csv", sep =",") %>%
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

write.tree(tree$scenario.3, file = "results/tree.tree")

# phylosignal tree graph ###
ape::read.tree("results/tree.tree")-> tree
x11()
plot(tree)
read.csv("data/species.csv", sep =",") %>%
  select(Taxon, PERoil, ratio, seedmass)%>%
  group_by(Taxon) %>%
  summarise(PERoil= mean(PERoil), 
            ratio= mean(ratio),
            seedmass = mean(seedmass))%>%
  select(Taxon, PERoil, ratio, seedmass)%>%
  mutate(label= gsub(" ","_", Taxon))%>%
  remove_rownames %>% 
  column_to_rownames(var="label") %>% 
  select(PERoil, ratio, seedmass)-> oil_phylo 


obj_oil <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo), match.data=TRUE)
#mat.col <- ifelse(tdata(obj_oil , "tip") < 1,  "darkred","darkgrey")
barplot.phylo4d(obj_oil, tree.ratio = 0.2,  center=T, #bar.col = mat.col,
                trait.bg.col = "white", show.box = T, grid.vertical = F,
                grid.horizontal = F, tip.labels = NULL, tree.ladderize =T) #rainbow(37)

barplot(obj_oil ,tree.type = "fan", tip.cex = 0.6, tree.open.angle = 160, trait.cex = 0.6)
dotplot.phylo4d (obj_oil , tree.type= "cladogram")
gridplot.phylo4d (obj_oil , tree.type = "fan", tip.cex = 0.6, show.trait = FALSE)
phyloSignal(obj_oil )
phyloCorrelogram (obj_oil )
