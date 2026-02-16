# ==============================================================================
# Section 3.4. Model Performance Comparison by Plot Types 
# ==============================================================================

# To evaluate the influence of tree presence on model performance, we addressed the imbalanced distribution of plot types (with/without trees) 
# through a rigorous bootstrap validation approach. Given that plots with and without large interspersed trees num-bered 28 and 19 respectively,
# we implemented bootstrapping with 100 iterations to prevent bias toward the more abundant plot type. In each iteration,
# we randomly resampled data from plots with trees to match the sample size of plots without trees (n = 19) prior to model fitting.

# Load libraries
library(tidyverse)    # For data manipulation & visualization
library(broom)        # For linear models (lm)
library(broom.mixed)  # For linear mixed-effects models (lmer)
library(lmerTest)     # For linear mixed-effects models with p-values
library(MuMIn)        # For marginal & conditional R²
library(performance)  # For model evaluation metric - RMSE
library(tidytext)     # For reorder in faceted plots
library(ggridges)     # For ridge density plots

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

# Prepare ground-LAI data: remove the Index column and rename mean_Index to PLOT_MAP_LAI
LAI_DATA <- PLOT_MAP_LAI %>% select(-Index) %>% rename(PLOT_MAP_LAI = mean_Index)

# Prepare vegetation index data: rename mean_Index to PLOT_MAP_INDEX
VI_DATA <- PLOT_MAP_ALL_INDEX %>% rename(PLOT_MAP_INDEX = mean_Index)

# Merge LAI and VI data by Plot and MAP (each row = one plot × one MAP × one VI)
ALL_DATA <- left_join(LAI_DATA, VI_DATA, by = c("Plot", "MAP"))

# Define plots without trees and assign plot type
without_tree_plots <- c("S32","S15","S14","S23","S17","S09","S20","S18","S19","S22","S26","S27","S28","S33","S29","S30","K13","K14","K15")

# Classify plots into "With tree" and "Without tree" categories
PLOT_TYPE <- ALL_DATA %>%
  mutate(Plot_type = ifelse(Plot %in% without_tree_plots, "Without tree", "With tree"))

# Function 1: Bootstrap Sampling ----
create_bootstrap_datasets <- function(data, n_bootstrap = 1000) {
  
  # Get unique plot IDs that belong to the "With tree" category
  with_tree_plots <- unique(data$Plot[data$Plot_type == "With tree"])
  
  # Initialize an empty list to store bootstrap datasets
  bootstrap_datasets <- list()
  
  # Loop through each bootstrap iteration
  for(i in 1:n_bootstrap) {
    
    cat("Creating bootstrap dataset:", i, "/", n_bootstrap, "\r")
    
    # Randomly sample 19 plots (without replacement) from the "With tree" group
    sampled_with_tree_plots <- sample(with_tree_plots, 19, replace = FALSE)
    
    # Filter data to keep only the sampled "With tree" plots
    balanced_with_tree <- data %>% 
      filter(Plot %in% sampled_with_tree_plots, Plot_type == "With tree")
    
    # Keep all "Without tree" plots (no subsampling needed)
    balanced_without_tree <- data %>%
      filter(Plot_type == "Without tree")
    
    # Combine both subsets into a balanced dataset
    BALANCED_INDEX_LAI <- rbind(balanced_with_tree, balanced_without_tree)
    
    # Store the balanced dataset with its bootstrap iteration ID
    bootstrap_datasets[[i]] <- BALANCED_INDEX_LAI %>%
      mutate(bootstrap_id = i)
  }
  
  cat("\n")
  return(bootstrap_datasets)
}

# Function 2: Mixed Effects Analysis
run_mixed_effects_analysis <- function(balanced_data) {
  
  # lme1: All plots
  lme1 <- balanced_data %>% 
    group_by(Index) %>% 
    nest() %>%
    mutate(
      model = map(data, ~ lmerTest::lmer(PLOT_MAP_LAI ~ PLOT_MAP_INDEX + (1|Plot), data = .)),    
      model_stats = map(model, ~glance(.x)),
      fixed_effects = map(model, ~tidy(.x, effects = "fixed")),
      random_effects = map(model, ~tidy(.x, effects = "ran_pars")),
      r_squared = map(model, ~MuMIn::r.squaredGLMM(.x))
    ) %>%
    mutate(
      # fixed effect - coefficient, std. error, degree of freedom, p-value
      n_obs = map_int(data, nrow),
      fixed_df = map_dbl(fixed_effects, ~.x$df[.x$term == "PLOT_MAP_INDEX"]), 
      fixed_coef = map_dbl(fixed_effects, ~.x$estimate[.x$term == "PLOT_MAP_INDEX"]),
      fixed_std.error = map_dbl(fixed_effects, ~.x$std.error[.x$term == "PLOT_MAP_INDEX"]),
      fixed_p_value = map_dbl(fixed_effects, ~.x$p.value[.x$term == "PLOT_MAP_INDEX"]),
      sig_p = case_when(
        fixed_p_value < 0.001 ~ "< 0.001",
        fixed_p_value < 0.01 ~ "< 0.01",
        fixed_p_value < 0.05 ~ "< 0.05",
        fixed_p_value < 0.1 ~ "< 0.1",
        TRUE ~ "ns"),
      #r squared - marginal & conditional
      R2_marginal = map_dbl(r_squared, ~.x[1, "R2m"]), # R2 for fixed effects only
      R2_conditional = map_dbl(r_squared, ~.x[1, "R2c"]), # R2 for fixed + random effects
      #random effect - standard deviation & variance
      random_sd = map_dbl(random_effects, ~.x$estimate[.x$term == "sd__(Intercept)"]),
      random_var = map_dbl(random_effects, ~(.x$estimate[.x$term == "sd__(Intercept)"])^2),
      AIC = map_dbl(model, AIC),
      BIC = map_dbl(model, BIC),
      RMSE = map_dbl(model, ~ performance::performance_rmse(.x, normalized = FALSE))) %>%
    mutate(MAP = "All MAP", Plot_type = "All plots", Model = "Mixed effect (Plot as random)") %>%
    select(Index, MAP, Plot_type, Model, n_obs, R2_marginal, fixed_p_value, sig_p, RMSE) %>% arrange(RMSE)
  
  # lme2: Without tree
  lme2 <- balanced_data %>% 
    filter(Plot_type == "Without tree") %>% 
    group_by(Index) %>% 
    nest() %>% 
    mutate(
      model = map(data, ~ lmerTest::lmer(PLOT_MAP_LAI ~ PLOT_MAP_INDEX + (1|Plot), data = .)),    
      model_stats = map(model, ~glance(.x)),
      fixed_effects = map(model, ~tidy(.x, effects = "fixed")),
      random_effects = map(model, ~tidy(.x, effects = "ran_pars")),
      r_squared = map(model, ~MuMIn::r.squaredGLMM(.x))
    ) %>%
    mutate(
      # fixed effect - coefficient, std. error, degree of freedom, p-value
      n_obs = map_int(data, nrow),
      fixed_df = map_dbl(fixed_effects, ~.x$df[.x$term == "PLOT_MAP_INDEX"]), 
      fixed_coef = map_dbl(fixed_effects, ~.x$estimate[.x$term == "PLOT_MAP_INDEX"]),
      fixed_std.error = map_dbl(fixed_effects, ~.x$std.error[.x$term == "PLOT_MAP_INDEX"]),
      fixed_p_value = map_dbl(fixed_effects, ~.x$p.value[.x$term == "PLOT_MAP_INDEX"]),
      sig_p = case_when(
        fixed_p_value < 0.001 ~ "< 0.001",
        fixed_p_value < 0.01 ~ "< 0.01",
        fixed_p_value < 0.05 ~ "< 0.05",
        fixed_p_value < 0.1 ~ "< 0.1",
        TRUE ~ "ns"),
      #r squared - marginal & conditional
      R2_marginal = map_dbl(r_squared, ~.x[1, "R2m"]), # R2 for fixed effects only
      R2_conditional = map_dbl(r_squared, ~.x[1, "R2c"]), # R2 for fixed + random effects
      #random effect - standard deviation & variance
      random_sd = map_dbl(random_effects, ~.x$estimate[.x$term == "sd__(Intercept)"]),
      random_var = map_dbl(random_effects, ~(.x$estimate[.x$term == "sd__(Intercept)"])^2),
      AIC = map_dbl(model, AIC),
      BIC = map_dbl(model, BIC),
      RMSE = map_dbl(model, ~ performance::performance_rmse(.x, normalized = FALSE))) %>%
    mutate(MAP = "All MAP", Plot_type = "Without tree", Model = "Mixed effect (Plot as random)") %>%
    select(Index, MAP, Plot_type, Model, n_obs, R2_marginal, fixed_p_value, sig_p, RMSE) %>% arrange(RMSE)
  
  # lme3: With tree
  lme3 <- balanced_data %>% 
    filter(Plot_type == "With tree") %>% 
    group_by(Index) %>% 
    nest() %>% 
    mutate(
      model = map(data, ~ lmerTest::lmer(PLOT_MAP_LAI ~ PLOT_MAP_INDEX + (1|Plot), data = .)),    
      model_stats = map(model, ~glance(.x)),
      fixed_effects = map(model, ~tidy(.x, effects = "fixed")),
      random_effects = map(model, ~tidy(.x, effects = "ran_pars")),
      r_squared = map(model, ~MuMIn::r.squaredGLMM(.x))
    ) %>%
    mutate(
      # fixed effect - coefficient, std. error, degree of freedom, p-value
      n_obs = map_int(data, nrow),
      fixed_df = map_dbl(fixed_effects, ~.x$df[.x$term == "PLOT_MAP_INDEX"]), 
      fixed_coef = map_dbl(fixed_effects, ~.x$estimate[.x$term == "PLOT_MAP_INDEX"]),
      fixed_std.error = map_dbl(fixed_effects, ~.x$std.error[.x$term == "PLOT_MAP_INDEX"]),
      fixed_p_value = map_dbl(fixed_effects, ~.x$p.value[.x$term == "PLOT_MAP_INDEX"]),
      sig_p = case_when(
        fixed_p_value < 0.001 ~ "< 0.001",
        fixed_p_value < 0.01 ~ "< 0.01",
        fixed_p_value < 0.05 ~ "< 0.05",
        fixed_p_value < 0.1 ~ "< 0.1",
        TRUE ~ "ns"),
      #r squared - marginal & conditional
      R2_marginal = map_dbl(r_squared, ~.x[1, "R2m"]), # R2 for fixed effects only
      R2_conditional = map_dbl(r_squared, ~.x[1, "R2c"]), # R2 for fixed + random effects
      #random effect - standard deviation & variance
      random_sd = map_dbl(random_effects, ~.x$estimate[.x$term == "sd__(Intercept)"]),
      random_var = map_dbl(random_effects, ~(.x$estimate[.x$term == "sd__(Intercept)"])^2),
      AIC = map_dbl(model, AIC),
      BIC = map_dbl(model, BIC),
      RMSE = map_dbl(model, ~ performance::performance_rmse(.x, normalized = FALSE))) %>%
    mutate(MAP = "All MAP", Plot_type = "With tree", Model = "Mixed effect (Plot as random)") %>%
    select(Index, MAP, Plot_type, Model, n_obs, R2_marginal, fixed_p_value, sig_p, RMSE) %>% arrange(RMSE)
  
  # Combine results from all three model subsets into one data frame
  combined_results <- bind_rows(lme1, lme2, lme3)
  
  return(list(
    lme1 = lme1,
    lme2 = lme2, 
    lme3 = lme3,
    combined = combined_results
  ))
}

# Function 3: Run Bootstrap Analysis ----
run_bootstrap_analysis <- function(data, n_bootstrap = 1000) {
  
  # Step 1: Create bootstrap datasets
  cat("Step 1: Creating bootstrap datasets...\n")
  bootstrap_datasets <- create_bootstrap_datasets(data, n_bootstrap)
  
  # Step 2: Run analysis on each bootstrap dataset
  cat("Step 2: Running mixed effects analysis on bootstrap datasets...\n")
  all_results <- list()
  
  for(i in 1:length(bootstrap_datasets)) {
    cat("Analyzing bootstrap dataset:", i, "/", length(bootstrap_datasets), "\r")
    
    # Run mixed effects analysis on the current bootstrap dataset
    analysis_results <- run_mixed_effects_analysis(bootstrap_datasets[[i]])
    
    # Append the bootstrap iteration ID to the combined results
    analysis_results$combined <- analysis_results$combined %>%
      mutate(bootstrap_id = i)
    
    all_results[[i]] <- analysis_results$combined
  }
  
  cat("\n")
  
  # Combine all bootstrap iteration results into a single data frame
  bootstrap_results <- bind_rows(all_results)
  
  return(bootstrap_results)
}

# Function 4: Summarize Bootstrap Results ----
summarize_bootstrap_results <- function(bootstrap_results) {
  
  summary <- bootstrap_results %>%
    group_by(Index, Plot_type) %>%
    summarise(
      n_bootstrap = n(),
      mean_RMSE = mean(RMSE, na.rm = TRUE),
      median_RMSE = median(RMSE, na.rm = TRUE),
      sd_RMSE = sd(RMSE, na.rm = TRUE),
      RMSE_ci_lower = quantile(RMSE, 0.025, na.rm = TRUE),
      RMSE_ci_upper = quantile(RMSE, 0.975, na.rm = TRUE),
      RMSE_cv = sd_RMSE / mean_RMSE,  # Coefficient of variation
      .groups = "drop"
    ) %>%
    arrange(Plot_type, mean_RMSE)
  
  return(summary)
}

# Main Execution ----
set.seed(1234)
bootstrap_results <- run_bootstrap_analysis(PLOT_TYPE, n_bootstrap = 100)
summary_results <- summarize_bootstrap_results(bootstrap_results)

# Filter summary results for "Without tree" plots only
summary_without_tree <- summary_results %>% filter(Plot_type == "Without tree") 

# Prepare ridge plot data: exclude "All plots", remove bootstrap ID, and reorder indices
ridge_data <- bootstrap_results %>% 
  filter(Plot_type != "All plots") %>% 
  select(-bootstrap_id) %>%
  mutate(Index = factor(Index, levels = index_order)) %>% 
  distinct()

# Define index ordering based on descending mean RMSE
index_order <- ridge_data %>%
  group_by(Index) %>%
  summarise(mean_rmse = mean(RMSE, na.rm = TRUE)) %>%
  arrange(desc(mean_rmse)) %>%
  pull(Index)

# Prepare vertical line (arrow) data for "Without tree" mean RMSE markers
vline_data <- summary_without_tree %>%
  mutate(Index = factor(Index, levels = index_order)) %>% 
  mutate(y_pos = as.numeric(as.factor(Index)))

# Create the ridge density plot comparing RMSE distributions by vegetation index
ggplot() +
  geom_density_ridges(
    data = ridge_data,
    aes(x = RMSE, y = Index, fill = "With tree"),
    color = "white",
    alpha = 0.7, 
    scale = 1,
    rel_min_height = 0.01,
    bandwidth = 0.01,
    jittered_points = F,
    point_alpha = 0.3,
    point_size = 0.4
  ) +
  geom_segment(
    data = vline_data,
    aes(x = mean_RMSE, xend = mean_RMSE, 
        y = y_pos + 0.1, yend = y_pos, color = "Without tree"),
    arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
    size = 0.1
  ) +
  scale_fill_manual(name = NULL, values = c("With tree" = "#05595B")) +
  scale_color_manual(values = c("Without tree" = "black")) +
  labs(
    x = "RMSE (Root Mean Square Error)",
    y = "Vegetation Indices (VIs)",
    fill = "Plot type", color = "Plot type") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 16))

dev.off()
