############################
#### FORAGING ANALYSIS #####
############################

# Script 4/4

# This script runs the testing of the relationship between foraging strategies and
# top-heaviness reported in the Schroeder et al. 'topHeaviness' project.
# INPUT 'Top_Heaviness' should be loaded by Script 1 prior
# to running this script, and other needed data are loaded in first section.
# OUTPUT of this script is the results of foraging analysis.

# The script should run as-is if the data files and are present
# in your working directory - otherwise, 'read_csv' calls may need to be amended
# to locate data files in your local file structure. Please refer to manuscript 
# for details on data collection and interpretation.

rm(list = ls())
library(dplyr)
library(ggplot2)
library(readr)
library(car)
library(FSA)
library(multcompView)
library(lme4)
library(lmerTest) 
library(emmeans)   

#### load and clean data ####

# reminder to set working directory w/ setwd()

# load data
load('Top_Heaviness.RData') # focal trait data
Foraging_Strategy_Raw <- read_csv("Foraging-Strategy_topHeaviness.csv") # foraging strategy of each species

## clean foraging strategy sheet
# take only needed columns, rename
Foraging_Strategy <- Foraging_Strategy_Raw %>%
  dplyr::select(`Species Name`,`Foraging Strategy`) %>%
  rename(Species = `Species Name`,
         ForagingStrategy = `Foraging Strategy`)


# fix spelling so that they are all the same formatting
Foraging_Strategy[Foraging_Strategy == "Mass Recruitment by Pheromone Trail"] <- "Mass Recruitment"
Foraging_Strategy[Foraging_Strategy == "Mass recruitment by pheromone trail"] <- "Mass Recruitment"
Foraging_Strategy[Foraging_Strategy == "Solitary foraging"] <- "Solitary Foraging"
# Group like categories
Foraging_Strategy[Foraging_Strategy == "Mass recruitment by foraging trail"] <- "Mass Recruitment"
Foraging_Strategy[Foraging_Strategy == "Tandem Running"] <- "Group Recruitment"
Foraging_Strategy[Foraging_Strategy == "Stable Trunk Trail"] <- "Stable Trail"
Foraging_Strategy[Foraging_Strategy == "Stable Trunk Trails"] <- "Stable Trail"
Foraging_Strategy[Foraging_Strategy == "Long-term Trail Network"] <- "Stable Trail"

# # how many species have each of the 4 strategies?
# table(Foraging_Strategy$ForagingStrategy)

# attach foraging strat to focal trait data to get object of analysis
Top_Heaviness_foraging <- left_join(Top_Heaviness, Foraging_Strategy, by = 'Species')
Top_Heaviness_foraging <- na.omit(Top_Heaviness_foraging)
# table(Top_Heaviness_foraging$ForagingStrategy)

# order foraging categories
Top_Heaviness_foraging <- Top_Heaviness_foraging %>%
  mutate(
       ForagingStrategy = factor(
       ForagingStrategy,
       levels = c("Solitary Foraging", "Group Recruitment", "Stable Trail", "Mass Recruitment")
    )
  )

### modeling
# fit the model
mod_foraging <- lmer(topHeaviness ~ ForagingStrategy + (1|Species), 
                     data = Top_Heaviness_foraging)

# ANOVA-style test for overall effect of ForagingStrategy
anova(mod_foraging)

# Nakagawa for conditional R2 - marginal R2
r2_nakagawa(mod_foraging)


# Plot
ggplot(Top_Heaviness_foraging, aes(x = ForagingStrategy, y = topHeaviness)) +
  geom_boxplot(fill = "gray", color = "black", outlier.shape = 21, 
               outlier.size = 2, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Foraging Strategy",
       y = "Top-heaviness") +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(size = 19),
    axis.title = element_text(size = 22, face = "bold"),
    legend.position = "none"
  )


