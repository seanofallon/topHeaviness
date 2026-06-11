############################
###### SOIL ANALYSIS #######
############################

# Script 3/4

# This script runs the testing of the relationship between soil properties and
# top-heaviness reported in the Schroeder et al. 'topHeaviness' project.
# INPUT 'Top_Heaviness' should be loaded by Script 1 prior
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

# set working directory 
setwd('C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/TopHeaviness/topHeaviness')
# load data
load('Top_Heaviness.RData') # focal trait data
load('Nest_Summaries.RData') # ref for nests
# attach soil data - HWSD attached w/ lat/lon via SiteID (matches NestSummaries)
HWSDSoil_Sheet_Raw <- read.csv("C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/Soil/HSWDsoil-sheet_w-siteid.csv")
# SoilGrids to be attached by matching lat/lon
SoilGrids_AllSites <- read.csv("C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/Soil/SoilGrids_AllSites.csv")
SoilGrids_BufferedExtraction1500m <- read.csv("C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/Soil/SoilGrids_BufferedExtraction1500m.csv")

#### attaching HWSD Data ####
## HWSDSoil_Sheet should have one row for each SiteID - raw datasheet has multiple 
## rows w/ same SiteID for two cases: 1) multi-species sites - more than one species was 
## excavated at the site, 2) HWSD border sites - when coordinates appear to be on
## the border between two HWSD soil mapping units, we record both 

# resolve HWSD border sites by averaging numeric columns, replacing character columns
# w/ 'border site' when the two entries differ
HWSDSoil_Sheet <- HWSDSoil_Sheet_Raw %>%
  group_by(SiteID, Species) %>%
  summarise(across(everything(), ~ {
    if (is.numeric(.x)) {
      mean(.x, na.rm = TRUE)
    } else {
      vals <- unique(.x)
      if (length(vals) == 1) vals else "border site"
    }
  }), .groups = "drop")

# resolve multi-species sites by keeping just the first row w/ each SiteID - 
# Species column should be dropped as it is not used for merges, just for reference in raw sheet
HWSDSoil_Sheet <- HWSDSoil_Sheet %>%
  arrange(SiteID) %>%
  distinct(SiteID, .keep_all = TRUE) %>%
  dplyr::select(-Species)


Soil_byNest <- Nest_Summaries[Nest_Summaries$NestID %in% Top_Heaviness$NestID,]


Soil_byNest <- Soil_byNest %>%
  left_join(
    HWSDSoil_Sheet %>% dplyr::select(SiteID, Field.Site.Lat., Field.Site.Lon.,
                                     AWC.for.rootable.soil.depth, Coarse.fragments....0.20.,
                                     Sand....0.20., Silt....0.20., Clay....0.20.,
                                     Bulk.Density..0.20., Reference.Bulk.Density..0.20.,
                                     Organic.Carbon.content....0.20.),
    by = "SiteID"
  )


Soil_byNest <- Soil_byNest[!is.na(Soil_byNest$SiteID), ]


# Step 4.1: Identify ambiguous SiteIDs and extract publication codes
ambiguousPubs <- Soil_byNest %>%
  filter(str_ends(SiteID, "_ambiguous")) %>%
  distinct(pubCode = str_remove(SiteID, "_ambiguous"))

# If there are any ambiguous SiteIDs, ensure lat2 and lon2 columns exist
if (nrow(ambiguousPubs) > 0) {
  if (!"lat2" %in% names(Soil_byNest)) Soil_byNest$lat2 <- NA
  if (!"lon2" %in% names(Soil_byNest)) Soil_byNest$lon2 <- NA
}

# Step 4.2: For each ambiguous pubCode...
for (i in seq_len(nrow(ambiguousPubs))) {
  pubCode <- ambiguousPubs$pubCode[i]
  ambiguous_id <- paste0(pubCode, "_ambiguous")
  
  # Get HWSD rows matching this pubCode
  matches <- HWSDSoil_Sheet %>%
    filter(str_starts(SiteID, paste0(pubCode, "_")))
  
  if (nrow(matches) != 2) {
    warning(paste("Expected 2 HWSD rows for", pubCode, "but found", nrow(matches)))
    next
  }
  
  # Get coordinate sets
  lat1 <- matches$Field.Site.Lat.[1]
  lon1 <- matches$Field.Site.Lon.[1]
  lat2 <- matches$Field.Site.Lat.[2]
  lon2 <- matches$Field.Site.Lon.[2]
  
  # put coords into Soil_byNest rows w/ matching SiteID
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"Field.Site.Lat."] <- lat1
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"Field.Site.Lon."] <- lon1
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"lat2"] <- lat2
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"lon2"] <- lon2
  
  # put means of HWSD vars of interest into Soil_byNest rows w/ matching SiteID
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"AWC.for.rootable.soil.depth"] <-
    mean(matches$AWC.for.rootable.soil.depth, na.rm = TRUE)
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"Coarse.fragments....0.20."] <-
    mean(matches$Coarse.fragments....0.20., na.rm = TRUE)
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"Sand....0.20."] <-
    mean(matches$Sand....0.20., na.rm = TRUE)
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"Silt....0.20."] <-
    mean(matches$Silt....0.20., na.rm = TRUE)
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"Clay....0.20."] <-
    mean(matches$Clay....0.20., na.rm = TRUE)
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"Bulk.Density..0.20."] <-
    mean(matches$Bulk.Density..0.20., na.rm = TRUE)
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"Reference.Bulk.Density..0.20."] <-
    mean(matches$Reference.Bulk.Density..0.20., na.rm = TRUE)
  Soil_byNest[Soil_byNest$SiteID==ambiguous_id,"Organic.Carbon.content....0.20."] <-
    mean(matches$Organic.Carbon.content....0.20., na.rm = TRUE)
}


# rename columns to not be annoying
Soil_byNest <- Soil_byNest %>%
  rename(lat = Field.Site.Lat.,
         lon = Field.Site.Lon.,
         SandHWSD = Sand....0.20.,
         SiltHWSD = Silt....0.20.,
         ClayHWSD = Clay....0.20.,
         BulkDensityHWSD = Bulk.Density..0.20.,
         CoarseFragmentsHWSD = Coarse.fragments....0.20.,
         OrganicCarbonContentHWSD = Organic.Carbon.content....0.20.,
         AWCforRootableSoilDepthHWSD = AWC.for.rootable.soil.depth,
         ReferenceBulkDensityHWSD = Reference.Bulk.Density..0.20.)

# fix formatting of lat/lon to numeric
Soil_byNest$lat <- as.numeric(Soil_byNest$lat)
Soil_byNest$lon <- as.numeric(Soil_byNest$lon)
Soil_byNest$lat2 <- as.numeric(Soil_byNest$lat2)
Soil_byNest$lon2 <- as.numeric(Soil_byNest$lon2)


#### attach SoilGrids data ####
SoilGrids_BufferedExtraction1500m <- SoilGrids_BufferedExtraction1500m %>%
  rename_with(
    ~ str_remove(., "_mean$"),
    .cols = ends_with("_mean")
  )

# Step 1: Join on lat and lon
SoilGrids_joined <- SoilGrids_AllSites %>%
  left_join(SoilGrids_BufferedExtraction1500m, by = c("lat", "lon"), suffix = c("", "_buffer"))

# Step 2: Fill in NAs from buffer version for columns that exist in both
shared_cols <- intersect(names(SoilGrids_AllSites), names(SoilGrids_BufferedExtraction1500m))
shared_cols <- setdiff(shared_cols, c("lat", "lon", "system.index",
                                      "site_id" , ".geo"))  # exclude join keys

for (col in shared_cols) {
  buffer_col <- paste0(col, "_buffer")
  SoilGrids_joined[[col]] <- ifelse(
    is.na(SoilGrids_joined[[col]]),
    SoilGrids_joined[[buffer_col]],
    SoilGrids_joined[[col]]
  )
}

# drop buffer columns 
SoilGrids <- SoilGrids_joined %>%
  dplyr::select(-ends_with("_buffer"))

# put lat/lon columns at the front
SoilGrids <- SoilGrids %>%
  dplyr::select(lat, lon, everything())

# remove columns not needed for attachment to Soils_byNest
SoilGrids <- SoilGrids %>%
  dplyr::select(-ends_with("_count"))
SoilGrids <- SoilGrids %>%
  dplyr:: select(-c('system.index','site_id','.geo'))

# split Soil_byNest into two groups for attaching SoilGrids
# nests w/ a single set of coords
singleSite_nests <- Soil_byNest %>%
  filter(is.na(lat2) | is.na(lon2))  # no second coordinate
# nests w/ ambiguous field sites - 2 sets of coords
ambiguousNests <- Soil_byNest %>%
  filter(!is.na(lat2) & !is.na(lon2))  # has both lat/lon and lat2/lon2
# fix formatting so coord columns are numeric
# single site
singleSite_nests$lat <- as.numeric(singleSite_nests$lat)
singleSite_nests$lon <- as.numeric(singleSite_nests$lon)
singleSite_nests$lat2<- as.numeric(singleSite_nests$lat2)
singleSite_nests$lon2 <- as.numeric(singleSite_nests$lon2)
# ambiguous site
ambiguousNests$lat <- as.numeric(ambiguousNests$lat)
ambiguousNests$lon <- as.numeric(ambiguousNests$lon)
ambiguousNests$lat2 <- as.numeric(ambiguousNests$lat2)
ambiguousNests$lon2 <- as.numeric(ambiguousNests$lon2)


# First, identify the SoilGrids columns to summarize (e.g., those ending in "_mean")
soil_cols <- names(SoilGrids)[grepl("_mean$", names(SoilGrids))]

# Prepare an empty list to collect results
results_list <- vector("list", nrow(ambiguousNests))

# Loop over each row of ambiguousNests
for (i in seq_len(nrow(ambiguousNests))) {
  # Extract the 4 coordinate pairs
  row <- ambiguousNests[i, ]
  coord_pairs <- list(
    c(row$lat, row$lon),
    c(row$lat2, row$lon2)
  )
  
  # Find all matching rows in SoilGrids
  matching <- SoilGrids %>%
    dplyr::filter(paste(lat, lon) %in% sapply(coord_pairs, function(x) paste(x[1], x[2])))
  
  # Compute column means for the selected _mean columns
  if (nrow(matching) > 0) {
    means <- matching %>%
      dplyr::select(all_of(soil_cols)) %>%
      summarise(across(everything(), mean, na.rm = TRUE))
  } else {
    # If no match, fill with NA
    means <- as.data.frame(matrix(NA, nrow = 1, ncol = length(soil_cols)))
    names(means) <- soil_cols
  }
  
  # Combine the original row with the computed means
  results_list[[i]] <- cbind(row, means)
}

# Bind all rows into a single data frame
ambiguousNests_withSoil <- do.call(rbind, results_list)



# join soil data to single site nests
singleSite_nests <- singleSite_nests %>%
  left_join(SoilGrids, by = c("lat", "lon"))

#
allSoil_byNest <- bind_rows(singleSite_nests, ambiguousNests_withSoil)


#### Calculate per-nest mean and log-ratios ####
### (should work w/ just allSoil_byNest and topHeaviness)
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
  mutate(sandClay_ratio = log(sand_MeanToDepth / clay_MeanToDepth),
         claySand_ratio = log(clay_MeanToDepth / sand_MeanToDepth))

#### Top_Heaviness_soil - object of analysis ####
Top_Heaviness_soil <- na.omit(Top_Heaviness_soil)


#### Modeling soil ####
## scaling / centering variables to standardize for comparison
Top_Heaviness_soil <- Top_Heaviness_soil %>%
  mutate(
    cfvo_z = as.numeric(scale(cfvo_MeanToDepth)),
    ocd_z = as.numeric(scale(ocd_MeanToDepth)),
    bdod_z = as.numeric(scale(bdod_MeanToDepth)),
    sandClay_z = as.numeric(scale(sandClay_ratio)),
    claySand_z = as.numeric(scale(claySand_ratio))
  )

# Pairplot
ggpairs(Top_Heaviness_soil[,c('cfvo_z','ocd_z','bdod_z','sandClay_z','claySand_z')])

ggpairs(Top_Heaviness_soil[,c('cfvo_MeanToDepth','ocd_MeanToDepth','bdod_MeanToDepth','sandClay_ratio','claySand_ratio')])

# bdod + claySand_ratio (cor is 0.704, but they get at different things)
model_vars <- lmer(topHeaviness ~ bdod_z + claySand_z + (1 | Species), data = Top_Heaviness_soil)
check_model(model_vars)
summary(model_vars)

# bdod + sandClay_ratio unscaled (cor is 0.704, but they get at different things)
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