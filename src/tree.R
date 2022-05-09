library(tidyverse)

### Phylo tree

read.csv("data/species.csv") %>%
  separate(Specie, into = c("Genus", "Species"), sep = " ") %>%
  mutate(species = paste(Genus, Species),
         genus = Genus,
         family = Family) %>%
  select(species, genus, family) %>%
  unique %>%
  mutate(family = fct_recode(family, 
                             "Asteraceae" = "Compositae",
                             "Fabaceae" = "Leguminosae",
                             "Asphodelaceae" = "Xanthorrhoeaceae")) %>%
  arrange(species) %>%
  na.omit -> 
  ranks1

library(V.PhyloMaker)

phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/tree.tree")
