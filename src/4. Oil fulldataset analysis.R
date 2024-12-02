library(tidyverse);library(readxl);library(rstatix)
library(vegan);library(glmmTMB);library(DHARMa)
library(phylosignal);library(phylobase);library(ape);library(tidytree)

# oil content data (own + literature)
read.csv("data/oil_fulldataset.csv", sep = ";")%>% rename(seed_mass= X50seed_mass_mg)-> oil_fulldataset # 80 species
unique(oil_fulldataset$family) # 19 families

# data exploration and visualization
hist(oil_fulldataset$oil_content) # clearly not normal
hist(oil_fulldataset$seed_mass)# clearly not normal

# Data format preparation for mcmc###

read.csv("data/oil_fulldataset.csv", sep = ";")%>%
  rename(seed_mass= X50seed_mass_mg)%>%
  group_by(Taxon, family, ecology)%>%
  summarise(oil_content = mean(oil_content), seed_mass = mean(seed_mass))%>%
  rename(familia = family)%>%
  convert_as_factor(Taxon, familia, ecology) %>%
  mutate(ecology= fct_relevel(ecology, "Strict lowland", "Generalist", "Strict alpine"))%>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil_content = log(oil_content),
         Lseed_mass=log(seed_mass))%>%
  na.omit()%>%
  as.data.frame()-> full_oil_analysis

# descriptive statistics
full_oil_analysis%>% 
  dplyr::select(Taxon, ecology, oil_content, seed_mass)%>%
  #group_by(ecology)%>%
  get_summary_stats(seed_mass)

hist(full_oil_analysis$Loil_content) # normally distributed
hist(full_oil_analysis$Lseed_mass)# normally distributed

##### USE MCMC to take into account phylogeny! 
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree_oil_fulldataset.tree")), 
                    ape::read.tree("results/tree_oil_fulldataset.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL
plot(ape::read.tree("results/tree_oil_fulldataset.tree"))
### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        #G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)

# 1. COMPARISON between ecology, NO DIFFERENCES RELATED TO ECOLOGY #####
# correct glm? ASK EDUARDO!!!
str(full_oil_analysis)
MCMCglmm::MCMCglmm(Loil_content ~ ecology, 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = full_oil_analysis,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g1

plot(g1)
summary(g1) 
x11()
read.csv("data/oil_fulldataset.csv", sep = ";")%>%
  convert_as_factor(family, ecology)%>%
  mutate(ecology= fct_relevel(ecology, "Strict lowland", "Generalist", "Strict alpine"))%>%
  ggplot()+
  geom_boxplot(aes(x=ecology, y=oil_content, fill=ecology), color= "black")+
  geom_point(aes(x=ecology, y=oil_content, fill=ecology), color= "black", alpha = 0.5, position = "jitter", shape = 21, size= 4, show.legend = F)+
  scale_fill_manual (name= "Ecology distribution", values = c("orange3", "darkcyan", "orchid4"), 
                     guide = guide_legend (title.position = "top",direction = "horizontal"))+
  ggthemes::theme_tufte(base_size=12) + 
  labs( title= "A) Species oil content (%)", y= "Oil content (%)")+ #",
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        legend.position = c(0.5,-0.06), 
        plot.margin = unit(c(0, 0,0.2,0.1), "cm"),
        legend.title = element_text(size = 12, color = "black", hjust = 0.5),
        legend.text = element_text(size = 10, color = "black"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.box.margin = unit(c(0, 0,0,0), "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"))-> fig2A;fig2A
  
# 2. Relationship between oil content and seed mass (trend but not significant, too much variation) ####
str(full_oil_analysis)
MCMCglmm::MCMCglmm(Loil_content ~ Lseed_mass, 
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = full_oil_analysis,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g2

plot(g2)
summary(g2) # no significant effect of seed mass in oil content (or viceversa)

col <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220", "#fcb622","darkgoldenrod3", "#f68b08", "#ff7125" ,"#c87107",  "#8E5005",
         "#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f",
         "#78E208", "#45A747", "#396F3E")
read.csv("data/oil_fulldataset.csv", sep = ";")%>%
  convert_as_factor(family, ecology)%>%
  mutate(ecology= fct_relevel(ecology, "Strict lowland", "Generalist", "Strict alpine"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>% 
  rename(seed_mass= X50seed_mass_mg)%>%
  #filter(source == "own_data")%>%
  ggplot()+
  geom_point(aes(x=seed_mass, y=oil_content, fill=family), color= "black", shape = 21, size= 4)+
  #geom_smooth(aes(x=seed_mass, y=oil_content),method = "glm", color = "black")+
  scale_fill_manual (name= "Families", values=col)+
  ggthemes::theme_tufte(base_size=12) + 
  labs( title = "B) Oil content vs seed mass relationship", y= "Oil content (%)", x= "Seed mass (mg)")+ 
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 14),
        legend.position = "right", 
        plot.margin = unit(c(0, 0,0.2,0.1), "cm"),
        legend.title = element_text(size = 12, color = "black", face = "bold"),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))-> fig2B;fig2B

# combine plots
library(patchwork)
fig2A + fig2B + plot_layout(widths= c(0.45, 0.55))-> fig2;fig2

# save plots
ggsave(filename = "regional patterns.png", plot =fig2 , path = "results/figures", 
        device = "png", dpi = 600) #scale = 1, height = 150, width = 180, units = "mm",
