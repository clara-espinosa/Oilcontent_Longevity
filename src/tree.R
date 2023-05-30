library(tidyverse)

### Phylo tree 2022 #####

read.csv("data/2022/species22.csv", sep =";") %>%
  separate(species, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = familia) %>%
  select(species, genus, family) %>%
  #unique %>%
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

write.tree(tree$scenario.3, file = "results/tree22.tree")

## Phylo tree 2023 #####

read.csv("data/2023/species23.csv", sep =";") %>%
  separate(species, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus) %>%
  select(species, genus, family) %>%
  #unique %>%
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

write.tree(tree$scenario.3, file = "results/tree23.tree")
