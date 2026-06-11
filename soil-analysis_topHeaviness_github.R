############################
###### SOIL ANALYSIS #######
############################

# Script 3/4

# This script runs the testing of the relationship between soil properties and
# top-heaviness reported in the Schroeder et al. 'topHeaviness' project.
# INPUT 'Top_Heaviness' should be loaded by Script 1 (available upon request) prior
# to running this script, and other needed data are loaded in first section.
# OUTPUT of this script is the results of soil analysis.

# The script should run as-is if the data files and are present
# in your working directory - otherwise, 'read_csv' calls may need to be amended
# to locate data files in your local file structure. Please refer to manuscript 
# for details on data collection and interpretation.

# clear environment
rm(list = ls())
# load packages
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(lmerTest)
library(sjPlot)
library(purrr)
library(GGally)
library(ggeffects)
library(performance)

# reminder to set working directory w/ setwd()

# load data
load('Top_Heaviness.RData') # focal trait data
allSoil_byNest <- read.csv("allSoil_byNest.csv")

#### Calculate per-nest mean and log-ratios ####

# Define the band breaks and their ranges
bands <- list(
  "0.5" = c(0, 5),
  "5.15" = c(5, 15),
  "15.30" = c(15, 30),
  "30.60" = c(30, 60),
  "60.100" = c(60, 100),
  "100.200" = c(100, 200)
)

band_labels <- names(bands)
soil_properties <- c("sand", "silt", "clay", "cfvo", "ocd", "bdod")

compute_soil_metrics_row <- function(...) {
  row <- list(...)
  depth <- row$NestDepthcm
  midpoint <- row$midNestDepthcm
  
  result <- list()
  
  for (prop in soil_properties) {
    total_weighted_sum <- 0
    total_thickness <- 0
    above_sum <- 0
    above_thickness <- 0
    below_sum <- 0
    below_thickness <- 0
    
    for (i in seq_along(bands)) {
      band <- bands[[i]]
      band_name <- band_labels[i]
      band_col <- paste0(prop, "_", band_name, "cm_mean")
      
      if (!(band_col %in% names(row))) next
      value <- as.numeric(row[[band_col]])
      if (is.na(value)) next
      
      band_start <- band[1]
      band_end <- band[2]
      
      ## -------- Whole depth --------
      if (i == length(bands) && depth > band_end) {
        eff_end <- depth
      } else {
        eff_end <- min(band_end, depth)
      }
      eff_start <- max(band_start, 0)
      
      if (eff_end > eff_start) {
        thickness <- eff_end - eff_start
        total_weighted_sum <- total_weighted_sum + value * thickness
        total_thickness <- total_thickness + thickness
      }
      
      ## -------- Above midpoint --------
      if (i == length(bands) && midpoint > band_start) {
        eff_end_top <- min(depth, midpoint)   # cap at midpoint not depth
      } else {
        eff_end_top <- min(band_end, midpoint)
      }
      eff_start_top <- max(band_start, 0)
      
      if (eff_end_top > eff_start_top) {
        thickness_top <- eff_end_top - eff_start_top
        above_sum <- above_sum + value * thickness_top
        above_thickness <- above_thickness + thickness_top
      }
      
      ## -------- Below midpoint --------
      if (i == length(bands) && depth > band_end) {
        eff_end_bot <- depth
      } else {
        eff_end_bot <- min(band_end, depth)
      }
      eff_start_bot <- max(band_start, midpoint)
      
      if (eff_end_bot > eff_start_bot) {
        thickness_bot <- eff_end_bot - eff_start_bot
        below_sum <- below_sum + value * thickness_bot
        below_thickness <- below_thickness + thickness_bot
      }
    }
    
    # Compute weighted means
    prop_mean <- if (total_thickness > 0) total_weighted_sum / total_thickness else NA
    prop_above <- if (above_thickness > 0) above_sum / above_thickness else NA
    prop_below <- if (below_thickness > 0) below_sum / below_thickness else NA
    prop_diff <- if (!is.na(prop_above) && !is.na(prop_below) && prop_below > 0) {
      log(prop_above / prop_below)
    } else {
      NA
    }
    
    result[[paste0(prop, "_MeanToDepth")]] <- prop_mean
    result[[paste0(prop, "_Above")]] <- prop_above
    result[[paste0(prop, "_Below")]] <- prop_below
    result[[paste0(prop, "_Diff")]] <- prop_diff
  }
  
  return(result)
}


# Apply to dataframe using purrr::pmap
soil_metrics_df <- allSoil_byNest %>%
  mutate(row_id = row_number()) %>%  # To maintain row order
  pmap(compute_soil_metrics_row) %>%
  bind_rows()

# Join back to original
allSoil_byNest <- bind_cols(allSoil_byNest, soil_metrics_df)


### join soil metrics to topHeaviness
Top_Heaviness_soil <- Top_Heaviness %>%
  left_join(
    allSoil_byNest %>%
      dplyr::select(NestID, ends_with("Above"), ends_with("Below"), ends_with("Diff"), ends_with("MeanToDepth")),
    by = "NestID"
  )

# Create sand/clay ratio
Top_Heaviness_soil <- Top_Heaviness_soil %>%
  mutate(claySand_ratio = log(clay_MeanToDepth / sand_MeanToDepth))

#### Top_Heaviness_soil - object of analysis ####
Top_Heaviness_soil <- na.omit(Top_Heaviness_soil)


#### Modeling soil ####
## scaling / centering variables to standardize for comparison
Top_Heaviness_soil <- Top_Heaviness_soil %>%
  mutate(
    cfvo_z = as.numeric(scale(cfvo_MeanToDepth)),
    ocd_z = as.numeric(scale(ocd_MeanToDepth)),
    bdod_z = as.numeric(scale(bdod_MeanToDepth)),
    claySand_z = as.numeric(scale(claySand_ratio))
  )

# Pairplot
ggpairs(Top_Heaviness_soil[,c('cfvo_z','ocd_z','bdod_z','claySand_z')])

ggpairs(Top_Heaviness_soil[,c('cfvo_MeanToDepth','ocd_MeanToDepth','bdod_MeanToDepth','claySand_ratio')])

# bdod + claySand_ratio (cor is 0.704, but they get at different things)
model_vars <- lmer(topHeaviness ~ bdod_z + claySand_z + (1 | Species), data = Top_Heaviness_soil)
check_model(model_vars)
summary(model_vars)


# bdod + claySand_ratio unscaled (cor is 0.704, but they get at different things)
model_vars_unscaled <- lmer(topHeaviness ~ bdod_MeanToDepth + claySand_ratio + (1 | Species), data = Top_Heaviness_soil)
check_model(model_vars_unscaled)
summary(model_vars_unscaled)

r2_nakagawa(model_vars_unscaled)


#### Visualizing bulk density  ####
# Get model-predicted effects for bdod_MeanToDepth
pred_bdod <- ggpredict(model_vars_unscaled, terms = "bdod_MeanToDepth [all]")

# Plot
ggplot() +
  # raw data points, colored by Species for transparency
  geom_point(data = Top_Heaviness_soil,
             aes(x = bdod_MeanToDepth, y = topHeaviness),
             alpha = 0.6, size = 2) +
  
  # model prediction line + 95% CI ribbon
  geom_ribbon(data = pred_bdod,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "grey70", alpha = 0.3) +
  geom_line(data = pred_bdod,
            aes(x = x, y = predicted),
            color = "black", linewidth = 1.2) +
  
  # labels
  labs(x = "Bulk Density (to depth, g/cm³)",
       y = "Top-heaviness") +
  
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(size = 19),
    axis.title = element_text(size = 22, face = "bold"),
    legend.position = "none"
  )

####