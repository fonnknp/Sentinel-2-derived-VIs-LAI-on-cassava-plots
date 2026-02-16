# Load libraries
library(tidyverse)    # For data manipulation & visualization

# --- 0. Data preparation ----

# Ground-LAI data
ALL_LAI <- read_csv("DATA/Ground-LAI-Clean.csv")
#unique(ALL_LAI$MAP) # MAP at 2,3,4,5,6,7,8
#unique(ALL_LAI$Plot) # 47 Plots

# All index data
ALL_INDEX <- read_csv("DATA/all_index_zonal_median.csv")
#unique(ALL_INDEX$MAP) # MAP at 2,3,4,5,7
#unique(ALL_INDEX$Plot) # 47 Plots

# Ground-LAI was averaged by plot and MAP (mean LAI)
PLOT_MAP_LAI <- ALL_LAI %>%
  filter(!MAP %in% c(6, 8)) %>% 
  group_by(Plot, MAP) %>%
  summarise(mean_Index = mean(LAI)) %>%
  mutate(Index = "Ground-LAI") %>%
  select(Plot, MAP, Index, mean_Index) %>%
  ungroup()

# Each index was averaged by plot & MAP (mean index)
PLOT_MAP_ALL_INDEX <- ALL_INDEX %>%
  mutate(Plot = recode(Plot, E7 = "E07", E8 = "E08", E9 = "E09")) %>% 
  group_by(Plot, MAP, Index) %>% 
  summarise(mean_Index = mean(Value)) %>%
  filter(Index %in% c("BNDVI","CIG","DVI","EVI","GNDVI","GRVI","NDVI","NDWI","RVI","SAVI","SeLI","VIG","TCARI")) %>%
  ungroup()

# Combine mean VIs and mean ground-LAI
ALL_INDEX_LAI <- rbind(PLOT_MAP_LAI, PLOT_MAP_ALL_INDEX)

# Define plots without trees and assign plot type
without_tree_plots <- c("S32","S15","S14","S23","S17","S09","S20","S18","S19","S22","S26","S27","S28","S33","S29","S30","K13","K14","K15")

ALL_INDEX_LAI <- ALL_INDEX_LAI %>%
  mutate(Plot_type = ifelse(Plot %in% without_tree_plots, "Without tree", "With tree"))

length(unique(ALL_INDEX_LAI[ALL_INDEX_LAI$Plot_type=="Without tree",]$Plot)) # 19 plots without trees
length(unique(ALL_INDEX_LAI[ALL_INDEX_LAI$Plot_type=="With tree",]$Plot)) # 28 plots with trees

# ===================================================================================================================
# Figure 3: Time-series of vegetation indices (VIs) and ground-LAI
# on cassava plots by Month After Planting (MAP)
# ===================================================================================================================

# Define facet ordering based on RMSE ranking from mixed-effects models (Section 3.2)
rmse.rank <- c("GNDVI", "NDWI", "SeLI", "SAVI", "NDVI", "DVI", "EVI",
               "BNDVI", "GRVI", "CIG", "TCARI", "RVI", "VIG")
rmse.rank <- append(rmse.rank, "Ground-LAI")

# Calculate mean of each index across MAP
Mean_INDEX <- ALL_INDEX_LAI %>%
  group_by(MAP, Index) %>%
  summarise(MAP_INDEX = mean(mean_Index)) %>%
  ungroup()

# Plot Time-series of VIs and ground-LAI by MAP with facet ordering based on RMSE ranking
ggplot() +
  geom_jitter(data = ALL_INDEX_LAI,
              aes(x = factor(MAP), y = mean_Index),
              width = 0.2, alpha = 0.2, size = 0.9) +
  geom_line(data = Mean_INDEX,
            aes(x = factor(MAP), y = MAP_INDEX, group = Index),
            linewidth = 0.5) +
  facet_wrap(~factor(Index, rmse.rank), scale = "free_y") +
  theme_bw() +
  theme(
    axis.title  = element_text(face = "bold", size = 15),
    strip.text  = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 15),
    legend.text  = element_text(size = 15),
    legend.position = "none"
  ) +
  labs(x = "Month After Planting (MAP)", y = "Value") +
  guides(alpha = "none") +
  NULL

dev.off()

# =====================================================================================================================
# Figure 4: Comparison of vegetation indices (VIs) with ground–LAI for cassava plots by Month After Planting (MAP)
# =====================================================================================================================
ALL_INDEX_LAI_wider <- ALL_INDEX_LAI %>% 
  group_by(MAP, Index) %>%
  summarise(Value = mean(mean_Index)) %>% 
  pivot_wider(names_from = Index, values_from = Value) %>% 
  ungroup()

# Function for min-max normalization
min_max_normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

all_index_norm <- ALL_INDEX_LAI_wider %>% 
  mutate(lai_norm=min_max_normalize(`Ground-LAI`),
         bndvi_norm=min_max_normalize(BNDVI),
         cig_norm=min_max_normalize(CIG),
         dvi_norm=min_max_normalize(DVI),
         evi_norm=min_max_normalize(EVI),
         gndvi_norm=min_max_normalize(GNDVI),
         grvi_norm=min_max_normalize(GRVI),
         ndvi_norm=min_max_normalize(NDVI),
         ndwi_norm=min_max_normalize(NDWI),
         rvi_norm=min_max_normalize(RVI),
         savi_norm=min_max_normalize(SAVI),
         seli_norm=min_max_normalize(SeLI),
         tcari_norm=min_max_normalize(TCARI),
         vig_norm=min_max_normalize(VIG))

# Create a function to prepare the tidy data
prepare_comparison_data <- function(data) {
  # Extract ground LAI for comparison
  ground_lai <- data %>%
    dplyr::select(MAP, lai_norm) %>%
    rename(value = lai_norm) %>%
    mutate(series = "ground_lai_norm") %>% 
    mutate(index = "ground_lai")
  
  # Process each vegetation index to include ground LAI comparison
  bndvi_comparison <- data %>%
    dplyr::select(MAP, bndvi_norm) %>%
    rename(value = bndvi_norm) %>%
    mutate(series = "bndvi_norm") %>%
    mutate(index = "bndvi") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "BNDVI with Ground-LAI")
  cig_comparison <- data %>%
    dplyr::select(MAP, cig_norm) %>%
    rename(value = cig_norm) %>%
    mutate(index = "cig") %>%
    mutate(series = "cig_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "CIG with Ground-LAI")
  dvi_comparison <- data %>%
    dplyr::select(MAP, dvi_norm) %>%
    mutate(index = "dvi") %>%
    rename(value = dvi_norm) %>%
    mutate(series = "dvi_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "DVI with Ground-LAI")
  evi_comparison <- data %>%
    dplyr::select(MAP, evi_norm) %>%
    mutate(index = "evi") %>%
    rename(value = evi_norm) %>%
    mutate(series = "evi_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "EVI with Ground-LAI")
  gndvi_comparison <- data %>%
    dplyr::select(MAP, gndvi_norm) %>%
    mutate(index = "gndvi") %>%
    rename(value = gndvi_norm) %>%
    mutate(series = "gndvi_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "GNDVI with Ground-LAI")
  grvi_comparison <- data %>%
    dplyr::select(MAP, grvi_norm) %>%
    mutate(index = "grvi") %>%
    rename(value = grvi_norm) %>%
    mutate(series = "grvi_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "GRVI with Ground-LAI")
  ndvi_comparison <- data %>%
    dplyr::select(MAP, ndvi_norm) %>%
    mutate(index = "ndvi") %>%
    rename(value = ndvi_norm) %>%
    mutate(series = "ndvi_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "NDVI with Ground-LAI")
  ndwi_comparison <- data %>%
    dplyr::select(MAP, ndwi_norm) %>%
    mutate(index = "ndwi") %>%
    rename(value = ndwi_norm) %>%
    mutate(series = "ndwi_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "NDWI with Ground-LAI")
  rvi_comparison <- data %>%
    dplyr::select(MAP, rvi_norm) %>%
    mutate(index = "rvi") %>%
    rename(value = rvi_norm) %>%
    mutate(series = "rvi_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "RVI with Ground-LAI")
  savi_comparison <- data %>%
    dplyr::select(MAP, savi_norm) %>%
    mutate(index = "savi") %>%
    rename(value = savi_norm) %>%
    mutate(series = "savi_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "SAVI with Ground-LAI")
  seli_comparison <- data %>%
    dplyr::select(MAP, seli_norm) %>%
    mutate(index = "seli") %>%
    rename(value = seli_norm) %>%
    mutate(series = "seli_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "SeLI with Ground-LAI")
  tcari_comparison <- data %>%
    dplyr::select(MAP, tcari_norm) %>%
    mutate(index = "tcari") %>%
    rename(value = tcari_norm) %>%
    mutate(series = "tcari_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "TCARI with Ground-LAI")
  vig_comparison <- data %>%
    dplyr::select(MAP, vig_norm) %>%
    mutate(index = "vig") %>%
    rename(value = vig_norm) %>%
    mutate(series = "vig_norm") %>%
    bind_rows(ground_lai) %>%
    mutate(panel = "VIG with Ground-LAI")
  
  # Combine all comparisons
  combined_data <- bind_rows(bndvi_comparison,cig_comparison,dvi_comparison,evi_comparison,
                             gndvi_comparison,grvi_comparison,ndvi_comparison,ndwi_comparison,
                             rvi_comparison,savi_comparison,seli_comparison,tcari_comparison,vig_comparison)
  
  return(combined_data)
}

# Prepare the data for comparison plots
comparison_data <- prepare_comparison_data(all_index_norm)

# Create the rank from 'm4' model performance (mixed-effects models) in 'Regression.R'
rmse.rank2 <- c("GNDVI", "NDWI", "SeLI", "SAVI", "NDVI", "DVI", "EVI",
               "BNDVI", "GRVI", "CIG", "TCARI", "RVI", "VIG")
rmse.rank2 <- paste(rmse.rank2, " with Ground-LAI", sep = "")

# Plot comparison of VIs with ground-LAI by MAP with facet ordering based on RMSE ranking from mixed-effects models
ggplot(comparison_data, aes(x = factor(MAP), y = value, group = series, linetype = index)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~factor(panel, rmse.rank2), scale = "free_y") +
  scale_linetype_manual(values = c(
    "evi" = "solid", "seli" = "solid", "dvi" = "solid", "grvi" = "solid",
    "cig" = "solid", "savi" = "solid", "ndvi" = "solid", "rvi" = "solid",
    "vig" = "solid", "tcari" = "solid", "gndvi" = "solid",
    "ndwi" = "solid", "bndvi" = "solid", "ground_lai" = "dashed")) +
  labs(x = "Month After Planting (MAP)", y = "Normalized Value (0-1 scale)") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 12)) +
  guides(linetype = "none") +
  theme(legend.position = "bottom")

dev.off()
