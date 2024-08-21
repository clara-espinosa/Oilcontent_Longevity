library(tidyverse)

### Phylo tree longevity #####

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

