############################
##### PHYLO ANALYSIS #######
############################

# Script 2/4

# This script runs phylogenetic analysis reported in the Schroeder et al.
# 'topHeaviness' project. INPUT 'Top_Heaviness' is produced by Script 1 (avalailable
# upon request) and is loaded with other needed data in first steps below.
# OUTPUT of this script is the results of phylo analysis.

# The script should run as-is if the data files and are present
# in your working directory - otherwise, 'read_csv' calls may need to be amended
# to locate data files in your local file structure. Please refer to manuscript 
# for details on data collection and interpretation.

rm(list = ls())
library(caper) # phylogenetic analysis
library(dplyr)
library(ape) # phylogenetic analysis
library(phytools) # phylogenetic analysis
library(ggplot2)
library(readr)
library(ggtree)

#### loading data ####
# reminder to set your working directory w/ setwd()

# load 
load('Top_Heaviness.RData') # focal trait data
ant.tree.moreau<-read.nexus("MoreauTree2016") # genus-level ant tree 

#### estimating phylogenetic signal ####
# get genus for tree
Top_Heaviness <- Top_Heaviness %>%
  mutate(genus = sub(" .*", "", Species))

#### Full phylo tree ####
ant.tree.moreau<-read.nexus("MoreauTree2016")

## Prune tree to the genera included in analysis
# ant.tree.moreau$tip.label
tips<-c("Acromyrmex_versicolor",
        "Aphaenogaster_occidentalis_NW",
        "Camponotus_maritimus",
        "Dorymyrmex_bicolor",
        "Ectatomma_opaciventre",
        "Formica_moki",
        "Harpegnathos_saltator",
        "Monomorium_pharaonis",
        "Mycetagroicus_triangularis",
        "Mycetarotes_acutus",
        "Mycetophylax_conformis",
        "Mycocepurus_goeldii",
        "Odontomachus_coquereli",
        "Pachycondyla_harpax",
        "Pheidole_longispinosa",
        "Prenolepis_imparis",
        "Sericomyrmex_Sp",
        "Trachymyrmex_arizonensis",
        "Veromessor_andrei"
)

# new tree with only represented genera
pruned.tree.all<-drop.tip(ant.tree.moreau,ant.tree.moreau$tip.label[-match(tips, ant.tree.moreau$tip.label)])

#Rename tip labels to create genus level tree
for (i in 1:length(pruned.tree.all$tip.label)) {
  split.tips<-strsplit(pruned.tree.all$tip.label[i],"_")
  genus.tip<-split.tips[[1]]
  pruned.tree.all$tip.label[i]<-genus.tip[1]
}

# species summary of top-heaviness
sp_Top_Heaviness <- Top_Heaviness %>%
  group_by(Species) %>%
  summarize(mean_topheavy=mean(topHeaviness),
            n = n(),
            genus = genus[1]) #find n which is the number of nests from that species 

### prep trait values for tree, phylo signal estimation
# filter to 'representative' dataset - the combination of species w/ most nests
# in it (ties just broken by alphabetical)
sp_Top_Heaviness_rep <- sp_Top_Heaviness %>%
  group_by(genus) %>%
  slice_max(n, with_ties = FALSE) %>%
  ungroup()
# make df to set rownames
sp_Top_Heaviness_rep <- as.data.frame(sp_Top_Heaviness_rep)

## make sp_Top_Heaviness_rep how phylosig wants it (matrix w/ genus rownames)
#set rownames as genera to match to tree
rownames(sp_Top_Heaviness_rep) <- sp_Top_Heaviness_rep$genus
# genus column is now redundant
sp_Top_Heaviness_rep <- sp_Top_Heaviness_rep %>%
  dplyr::select(-genus)

# make matrix
repTH_forMap <- as.matrix(sp_Top_Heaviness_rep)[,'mean_topheavy']


# w/ genus rownames
repTH_forMap <- setNames(as.numeric(repTH_forMap), names(repTH_forMap))
 
# put data in readable frame
repTH_df <- data.frame(
  genus = names(repTH_forMap),
  Top_Heaviness = repTH_forMap
)
# make sure top-heaviness is numeric
repTH_df$Top_Heaviness <- as.numeric(repTH_df$Top_Heaviness)

# Plot
p <- ggtree(pruned.tree.all) %<+% repTH_df +
  geom_tippoint(aes(fill = Top_Heaviness), size = 5, stroke = 0.5, shape = 21, color = 'black') +
  geom_tiplab(aes(label = label),                          # genus/species labels
              hjust = -0.1,                                  # offset a bit to right of tips
              size = 6,
              fontface = 'italic') +                                  # text size
    scale_fill_viridis_c(option = "plasma", name = "Top-heaviness") +
  theme(
    legend.position.inside = c(0.05, 0.95),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 50, 10, 10)                    # give room for labels
  ) +
  xlim(0, max(ggtree(pruned.tree.all)$data$x) + 50)  # increase this number if still cut off


p2 <- p +
  theme(
    legend.position = c(0.14, 0.94), 
    legend.justification = c(0, 1), 
    legend.background = element_rect(fill = "white", color = "white", size = 0.2),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(10, 10, 10, 10)      # reduce big right margin so legend can overlap
  ) +
  coord_cartesian(clip = "off")               # allow legend to draw outside panel clipping

p2

# p value
phylo.res.rep <- data.frame(phylo.lambda=numeric(),
                            phylo.lambda.p=numeric(),
                            phylo.K=numeric(),
                            phylo.K.p=numeric())

phylo.lambda.rep<-phylosig(pruned.tree.all, repTH_forMap, method = "lambda", test=T)
# w/ K
phylo.K.rep<-phylosig(pruned.tree.all, repTH_forMap, method = "K", test=T)
## write to output object phylo.res
phylo.res.rep[1,"phylo.lambda"]<-phylo.lambda.rep$lambda
phylo.res.rep[1,"phylo.lambda.p"]<-phylo.lambda.rep$P
phylo.res.rep[1,"phylo.K"]<-phylo.K.rep$K
phylo.res.rep[1,"phylo.K.p"]<-phylo.K.rep$P


### Get phylogenetic signal for every combination of species (1 / genus)
### Make object with 1 row for each combination of species
# make object with just each species and its genus
sp_gn <- sp_Top_Heaviness %>% 
  dplyr::select(Species,genus)

# make df w/ species+genus for repeated genera
# empty df to hold app rows
rep_sp_gn<-data.frame(Species=character(),
                      genus=character())

# fill w/ repeats
rep_sp_gn <- sp_gn %>%
  group_by(genus) %>%
  filter(n() > 1) %>%
  ungroup()

# make df w/ species+genus for species in every combo (only one species per genus)
# empty df to hold app rows
single_sp_gn<-data.frame(Species=character(),
                         genus=character())

# fill w/ singles
single_sp_gn <- sp_gn %>%
  group_by(genus) %>%
  filter(n() == 1) %>%
  ungroup()
# single_sp_gn <- as.vector(single_sp_gn)

tbl=table(rep_sp_gn$genus)
all_sp_poss = matrix(NA, ncol=length(unique(rep_sp_gn$genus)))
unq_gn=unique(rep_sp_gn$genus)

for (i in 1:tbl[1]){
  for (j in 1:tbl[2]){
    for (k in 1:tbl[3]){
      for (m in 1:tbl[4]){
        for (n in 1:tbl[5]) {
          for (o in 1:tbl[6]){
            gni=unq_gn[1]
            spi=rep_sp_gn$Species[rep_sp_gn$genus==gni][i]
            gnj=unq_gn[2]
            spj=rep_sp_gn$Species[rep_sp_gn$genus==gnj][j]
            gnk=unq_gn[3]
            spk=rep_sp_gn$Species[rep_sp_gn$genus==gnk][k]
            gnm=unq_gn[4]
            spm=rep_sp_gn$Species[rep_sp_gn$genus==gnm][m]
            gnn=unq_gn[5]
            spn=rep_sp_gn$Species[rep_sp_gn$genus==gnn][n]
            gno=unq_gn[6]
            spo=rep_sp_gn$Species[rep_sp_gn$genus==gno][o]
            all_sp_poss=rbind(all_sp_poss,c(spi,spj,spk,spm,spn,spo))
          }
        }
      }
    }
  }
}

#remove NA row at top
all_sp_poss<-all_sp_poss[-1,]

## perform phylo signal estimation on all combos made in all_sp_poss
#make empty matrix to be filled w/ all combos
all_subs<-matrix(data=NA,nrow=288,ncol=19)
# combine single rep species (single_sp_gn) w/ all_sp_poss to have matrix w/ each combo (all_aa_subs)
for (i in 1:nrow(all_sp_poss)) {
  all_subs[i,]<-c(all_sp_poss[i,], single_sp_gn$Species)
}


## Use each row in a phylosig call to get results for every combination of rep species
# Start df
phylo.res<-data.frame(topheavy_sum=numeric(),
                      phylo.lambda=numeric(),
                      phylo.lambda.p=numeric(),
                      phylo.K=numeric(),
                      phylo.K.p=numeric())

#turn into df
sp_Top_Heaviness<- as.data.frame(sp_Top_Heaviness)

# Write results for each subset
for (i in 1:nrow(all_subs)) {
  ## save ith subset to df
  topheavy_subi<-sp_Top_Heaviness %>% filter(Species %in% all_subs[i,])
  rownames(topheavy_subi)<-topheavy_subi$genus 
  ## make object for phylosig call
  topheavy.trait<-setNames(topheavy_subi[,2],rownames(topheavy_subi)) 
  
  # add sum of all topheavy values as a column to check that loop is running
  phylo.res[i,"topheavy_sum"] <- sum(topheavy.trait)
  
  ## test for phylogenetic signal
  # w/ lambda
  phylo.lambda<-phylosig(pruned.tree.all, topheavy.trait, method = "lambda", test=T)
  # w/ K
  phylo.K<-phylosig(pruned.tree.all, topheavy.trait, method = "K", test=T)
  ## write to output object phylo.res
  phylo.res[i,"phylo.lambda"]<-phylo.lambda$lambda
  phylo.res[i,"phylo.lambda.p"]<-phylo.lambda$P
  phylo.res[i,"phylo.K"]<-phylo.K$K
  phylo.res[i,"phylo.K.p"]<-phylo.K$P
}

##plot results of trees
par(mfrow = c(1,2),
    mar = c(4.5,4,1,3),
    mgp=c(2.3,1,0))
boxplot(phylo.res$phylo.lambda.p,
        ylab = "Likelihood ratio test p-value",
        xlab = "Nest Feature",
        col = "white",
        ylim = c(0,1),
        cex.axis = 0.7,
        cex.lab = 1.1,
        las = 1)
boxplot(phylo.res$phylo.lambda,
        ylab = "Lambda (Phylogenetic signal)",
        xlab = "Nest Feature",
        col = "white",
        ylim = c(0,2),
        cex.axis = 0.7,
        cex.lab = 1.1,
        las = 1)



