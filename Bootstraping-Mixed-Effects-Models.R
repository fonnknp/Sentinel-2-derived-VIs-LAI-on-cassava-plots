# Bootstrap Data First, Then Analysis
library(tidyverse)
library(lmerTest)
library(broom)
library(broom.mixed)
library(performance)

# Import data----
ALL_LAI <- read_csv("Data/all-index-pran/Ground-LAI-Clean.csv")
PLOT_MAP_LAI <- ALL_LAI %>%
  group_by(Plot,MAP) %>%
  summarise(PLOT_MAP_LAI = mean(LAI)) %>%
  ungroup() %>%
  dplyr::select(Plot, MAP, PLOT_MAP_LAI)

ALL_INDEX <- read_csv("Data/all_index_zonal_median.csv")
PLOT_MAP_ALL_INDEX <- ALL_INDEX %>%
  mutate(Plot = recode(Plot,E7="E07",E8="E08",E9="E09")) %>% 
  group_by(Plot,MAP,Index) %>% 
  summarise(PLOT_MAP_INDEX = mean(Value)) %>%
  filter(Index %in% c("BNDVI","CIG","DVI","EVI","GNDVI","GRVI","NDVI","NDWI","RVI","SAVI","SeLI","VIG","TCARI")) %>%
  ungroup()

ALL_INDEX_LAI <- left_join(PLOT_MAP_LAI,PLOT_MAP_ALL_INDEX, by=c("Plot", "MAP"))

without <- ALL_INDEX_LAI %>% 
  filter(Plot %in% c("S32","S15","S14","S23","S17","S09","S20","S18","S19","S22","S26","S27","S28","S33","S29","S30","K13","K14","K15")) %>% 
  mutate(Plot_type = "Without tree")

with <- ALL_INDEX_LAI %>% 
  filter(!Plot %in% c("S32","S15","S14","S23","S17","S09","S20","S18","S19","S22","S26","S27","S28","S33","S29","S30","K13","K14","K15")) %>% 
  mutate(Plot_type = "With tree")

ALL_INDEX_LAI <- rbind(without,with)

# Function 1: Bootstrap Sampling ----
create_bootstrap_datasets <- function(data, n_bootstrap = 1000) {
  
  with_tree_plots <- unique(data$Plot[data$Plot_type == "With tree"])
  
  bootstrap_datasets <- list()
  
  for(i in 1:n_bootstrap) {
    
    cat("Creating bootstrap dataset:", i, "/", n_bootstrap, "\r")
    
    sampled_with_tree_plots <- sample(with_tree_plots, 19, replace = FALSE)
    
    balanced_with_tree <- data %>% 
      filter(Plot %in% sampled_with_tree_plots, Plot_type == "With tree")
    
    balanced_without_tree <- data %>%
      filter(Plot_type == "Without tree")
    
    BALANCED_INDEX_LAI <- rbind(balanced_with_tree, balanced_without_tree)
    
    bootstrap_datasets[[i]] <- BALANCED_INDEX_LAI %>%
      mutate(bootstrap_id = i)
  }
  
  cat("\n")
  return(bootstrap_datasets)
}

# Function 2: Mixed Effects Analysis----
run_mixed_effects_analysis <- function(balanced_data) {
  
  # lme1: All plots
  lme1 <- balanced_data %>% 
    drop_na() %>%
    filter(MAP < 8) %>% 
    group_by(Index) %>% 
    nest() %>% 
    mutate(
      model = map(data, ~ safely(lmerTest::lmer)(PLOT_MAP_LAI ~ PLOT_MAP_INDEX + (1|Plot), data = .)),
      model_success = map_lgl(model, ~ is.null(.x$error)),
      model = map(model, "result")
    ) %>%
    filter(model_success) %>%
    mutate(
      model_stats = map(model, ~glance(.x)),
      fixed_effects = map(model, ~tidy(.x, effects = "fixed")),
      random_effects = map(model, ~tidy(.x, effects = "ran_pars"))
    ) %>%
    mutate(
      n_obs = map_int(data, nrow),
      fixed_df = map_dbl(fixed_effects, ~.x$df[.x$term == "PLOT_MAP_INDEX"]), 
      fixed_coef = map_dbl(fixed_effects, ~.x$estimate[.x$term == "PLOT_MAP_INDEX"]),
      fixed_std.error = map_dbl(fixed_effects, ~.x$std.error[.x$term == "PLOT_MAP_INDEX"]),
      fixed_p_value = map_dbl(fixed_effects, ~.x$p.value[.x$term == "PLOT_MAP_INDEX"]),
      random_sd = map_dbl(random_effects, ~.x$estimate[.x$term == "sd__(Intercept)"]),
      random_var = map_dbl(random_effects, ~(.x$estimate[.x$term == "sd__(Intercept)"])^2),
      AIC = map_dbl(model, AIC),
      BIC = map_dbl(model, BIC),
      RMSE = map_dbl(model, ~ performance::performance_rmse(.x, normalized = FALSE))
    ) %>%
    mutate(MAP = "All MAP", Plot_type = "All plots", Model = "Mixed effect (Plot as random)") %>%
    select(Index, MAP, Plot_type, Model, RMSE) %>% 
    arrange(RMSE)
  
  # lme2: Without tree
  lme2 <- balanced_data %>% 
    drop_na() %>% 
    filter(Plot_type == "Without tree") %>% 
    filter(MAP < 8) %>% 
    group_by(Index) %>% 
    nest() %>% 
    mutate(
      model = map(data, ~ safely(lmerTest::lmer)(PLOT_MAP_LAI ~ PLOT_MAP_INDEX + (1|Plot), data = .)),
      model_success = map_lgl(model, ~ is.null(.x$error)),
      model = map(model, "result")
    ) %>%
    filter(model_success) %>%
    mutate(
      model_stats = map(model, ~glance(.x)),
      fixed_effects = map(model, ~tidy(.x, effects = "fixed")),
      random_effects = map(model, ~tidy(.x, effects = "ran_pars"))
    ) %>%
    mutate(
      n_obs = map_int(data, nrow),
      fixed_df = map_dbl(fixed_effects, ~.x$df[.x$term == "PLOT_MAP_INDEX"]), 
      fixed_coef = map_dbl(fixed_effects, ~.x$estimate[.x$term == "PLOT_MAP_INDEX"]),
      fixed_se = map_dbl(fixed_effects, ~.x$std.error[.x$term == "PLOT_MAP_INDEX"]),
      fixed_p_value = map_dbl(fixed_effects, ~.x$p.value[.x$term == "PLOT_MAP_INDEX"]),
      random_sd = map_dbl(random_effects, ~.x$estimate[.x$term == "sd__(Intercept)"]),
      random_var = map_dbl(random_effects, ~(.x$estimate[.x$term == "sd__(Intercept)"])^2),
      AIC = map_dbl(model, AIC),
      BIC = map_dbl(model, BIC),
      RMSE = map_dbl(model, ~ performance::performance_rmse(.x, normalized = FALSE))
    ) %>%
    mutate(MAP = "All MAP", Plot_type = "Without tree", Model = "Mixed effect (Plot as random)") %>% 
    select(Index, MAP, Plot_type, Model, RMSE) %>% 
    arrange(RMSE)
  
  # lme3: With tree
  lme3 <- balanced_data %>% 
    drop_na() %>% 
    filter(Plot_type == "With tree") %>% 
    filter(MAP < 8) %>% 
    group_by(Index) %>% 
    nest() %>% 
    mutate(
      model = map(data, ~ safely(lmerTest::lmer)(PLOT_MAP_LAI ~ PLOT_MAP_INDEX + (1|Plot), data = .)),
      model_success = map_lgl(model, ~ is.null(.x$error)),
      model = map(model, "result")
    ) %>%
    filter(model_success) %>%
    mutate(
      model_stats = map(model, ~glance(.x)),
      fixed_effects = map(model, ~tidy(.x, effects = "fixed")),
      random_effects = map(model, ~tidy(.x, effects = "ran_pars"))
    ) %>%
    mutate(
      n_obs = map_int(data, nrow),
      fixed_df = map_dbl(fixed_effects, ~.x$df[.x$term == "PLOT_MAP_INDEX"]), 
      fixed_coef = map_dbl(fixed_effects, ~.x$estimate[.x$term == "PLOT_MAP_INDEX"]),
      fixed_se = map_dbl(fixed_effects, ~.x$std.error[.x$term == "PLOT_MAP_INDEX"]),
      fixed_p_value = map_dbl(fixed_effects, ~.x$p.value[.x$term == "PLOT_MAP_INDEX"]),
      random_sd = map_dbl(random_effects, ~.x$estimate[.x$term == "sd__(Intercept)"]),
      random_var = map_dbl(random_effects, ~(.x$estimate[.x$term == "sd__(Intercept)"])^2),
      AIC = map_dbl(model, AIC),
      BIC = map_dbl(model, BIC),
      RMSE = map_dbl(model, ~ performance::performance_rmse(.x, normalized = FALSE))
    ) %>%
    mutate(MAP = "All MAP", Plot_type = "With tree", Model = "Mixed effect (Plot as random)") %>% 
    select(Index, MAP, Plot_type, Model, RMSE) %>% 
    arrange(RMSE)
  
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

    analysis_results <- run_mixed_effects_analysis(bootstrap_datasets[[i]])
    

    analysis_results$combined <- analysis_results$combined %>%
      mutate(bootstrap_id = i)
    
    all_results[[i]] <- analysis_results$combined
  }
  
  cat("\n")
  
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

cat("=== BOOTSTRAP ANALYSIS WORKFLOW ===\n")
cat("This will take some time...\n\n")

# bootstrap analysis
bootstrap_results <- run_bootstrap_analysis(ALL_INDEX_LAI, n_bootstrap = 100) ทดสอบ

cat("Step 3: Summarizing results...\n")
summary_results <- summarize_bootstrap_results(bootstrap_results)

# Show results----
cat("\n=== BOOTSTRAP MIXED EFFECTS MODEL RESULTS ===\n")
cat("Summary across all bootstrap iterations\n\n")

cat("=== ALL PLOTS ===\n")
all_plots_summary <- summary_results %>% filter(Plot_type == "All plots")
print(all_plots_summary)

cat("\n=== WITHOUT TREE ===\n")
without_tree_summary <- summary_results %>% filter(Plot_type == "Without tree")
print(without_tree_summary)

cat("\n=== WITH TREE ===\n")
with_tree_summary <- summary_results %>% filter(Plot_type == "With tree")
print(with_tree_summary)

# Top performers by category
cat("\n=== TOP 3 PERFORMERS BY CATEGORY ===\n")
top_performers <- summary_results %>%
  group_by(Plot_type) %>%
  slice_min(mean_RMSE, n = 3) %>%
  select(Plot_type, Index, mean_RMSE, RMSE_ci_lower, RMSE_ci_upper, RMSE_cv)

top_performers

# Visualization ----
library(ggplot2)

comparison_plot <- summary_results %>%
  ggplot(aes(x = reorder(Index, mean_RMSE), y = mean_RMSE, color = Plot_type)) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = RMSE_ci_lower, ymax = RMSE_ci_upper),
                width = 0.2, position = position_dodge(width = 0.3)) +
  #facet_wrap(~Plot_type, scales = "free_y") +
  coord_flip() +
  labs(title = "Bootstrap Mixed Effects Model Performance",
       subtitle = "RMSE with 95% Bootstrap Confidence Intervals",
       x = "Vegetation Index", 
       y = "RMSE",
       color = "Plot Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

comparison_plot

# Stability plot (CV)
stability_plot <- summary_results %>%
  ggplot(aes(x = reorder(Index, RMSE_cv), y = RMSE_cv, fill = Plot_type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Model Stability (Coefficient of Variation)",
       subtitle = "Lower CV = More stable across bootstrap iterations",
       x = "Vegetation Index", 
       y = "Coefficient of Variation (CV)",
       fill = "Plot Type") +
  theme_minimal()

print(stability_plot)

cat("\n=== ANALYSIS COMPLETED ===\n")
cat("Bootstrap iterations completed:", max(bootstrap_results$bootstrap_id, na.rm = TRUE), "\n")
cat("Total models fitted:", nrow(bootstrap_results), "\n")

# Save results (optional)
# write_csv(bootstrap_results, "bootstrap_mixed_effects_results.csv")
# write_csv(summary_results, "bootstrap_summary_results.csv")

# Visualization - Distribution----
library(ggplot2)

# Function: Create Bootstrap Distribution Plots ----
create_bootstrap_distribution_plots <- function(bootstrap_results) {
  
  distribution_data <- bootstrap_results %>%
    filter(Plot_type %in% c("With tree", "Without tree")) %>%
    select(Index, Plot_type, RMSE, bootstrap_id)
  
  # 1. Density plots for all indices combined
  overall_density_plot <- distribution_data %>%
    ggplot(aes(x = RMSE, fill = Plot_type, color = Plot_type)) +
    geom_density(alpha = 0.6, size = 1) +
    facet_wrap(~Plot_type, scales = "free") +
    labs(title = "Bootstrap RMSE Distributions: With Tree vs Without Tree",
         subtitle = "Density curves showing distribution of RMSE across all vegetation indices",
         x = "RMSE",
         y = "Density",
         fill = "Plot Type",
         color = "Plot Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 2. Individual index distributions (top 6 indices)
  top_indices <- distribution_data %>%
    group_by(Index) %>%
    summarise(mean_rmse = mean(RMSE, na.rm = TRUE)) %>%
    slice_min(mean_rmse, n = 6) %>%
    pull(Index)
  
  individual_density_plot <- distribution_data %>%
    filter(Index %in% top_indices) %>%
    ggplot(aes(x = RMSE, fill = Plot_type, color = Plot_type)) +
    geom_density(alpha = 0.6, size = 0.8) +
    facet_grid(Index ~ Plot_type, scales = "free") +
    labs(title = "Bootstrap RMSE Distributions by Vegetation Index",
         subtitle = "Top 6 performing indices: With Tree vs Without Tree",
         x = "RMSE",
         y = "Density",
         fill = "Plot Type",
         color = "Plot Type") +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 9))
  
  # 3. Violin plots for comparison
  violin_plot <- distribution_data %>%
    filter(Index %in% top_indices) %>%
    ggplot(aes(x = Index, y = RMSE, fill = Plot_type)) +
    geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), 
                 alpha = 0.8, outlier.size = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2,
                 position = position_dodge(width = 0.8), color = "red") +
    coord_flip() +
    labs(title = "Bootstrap RMSE Distribution Comparison",
         subtitle = "Violin plots with boxplots (red diamonds = means)",
         x = "Vegetation Index",
         y = "RMSE",
         fill = "Plot Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 4. Ridge plots (if ggridges available)
  ridge_plot <- tryCatch({
    library(ggridges)
    distribution_data %>%
      filter(Index %in% top_indices) %>%
      mutate(Index_PlotType = paste(Index, "-", Plot_type)) %>%
      ggplot(aes(x = RMSE, y = Index_PlotType, fill = Plot_type)) +
      geom_density_ridges(alpha = 0.7, scale = 2) +
      labs(title = "Ridge Plot: Bootstrap RMSE Distributions",
           subtitle = "Comparing With Tree vs Without Tree across indices",
           x = "RMSE",
           y = "Index - Plot Type",
           fill = "Plot Type") +
      theme_minimal() +
      theme(legend.position = "bottom")
  }, error = function(e) {
    cat("ggridges package not available, skipping ridge plot\n")
    return(NULL)
  })
  
  return(list(
    overall = overall_density_plot,
    individual = individual_density_plot,
    violin = violin_plot,
    ridge = ridge_plot
  ))
}

# Function: Statistical Comparison ----
compare_distributions <- function(bootstrap_results) {
  
  cat("\n=== STATISTICAL COMPARISON OF DISTRIBUTIONS ===\n")
  

  comparison_data <- bootstrap_results %>%
    filter(Plot_type %in% c("With tree", "Without tree")) %>%
    select(Index, Plot_type, RMSE, bootstrap_id) %>%
    pivot_wider(names_from = Plot_type, values_from = RMSE, 
                names_prefix = "RMSE_", id_cols = c(Index, bootstrap_id)) %>%
    drop_na()
  

  statistical_tests <- comparison_data %>%
    group_by(Index) %>%
    summarise(
      # Descriptive statistics
      mean_with_tree = mean(`RMSE_With tree`, na.rm = TRUE),
      mean_without_tree = mean(`RMSE_Without tree`, na.rm = TRUE),
      mean_difference = mean_with_tree - mean_without_tree,
      
      # Standard deviations
      sd_with_tree = sd(`RMSE_With tree`, na.rm = TRUE),
      sd_without_tree = sd(`RMSE_Without tree`, na.rm = TRUE),
      
      # Effect size (Cohen's d)
      cohens_d = mean_difference / sqrt((sd_with_tree^2 + sd_without_tree^2) / 2),
      
      # Paired t-test
      t_test_p_value = tryCatch({
        t.test(`RMSE_With tree`, `RMSE_Without tree`, paired = TRUE)$p.value
      }, error = function(e) NA),
      
      # Wilcoxon signed-rank test (non-parametric alternative)
      wilcox_p_value = tryCatch({
        wilcox.test(`RMSE_With tree`, `RMSE_Without tree`, paired = TRUE)$p.value
      }, error = function(e) NA),
      
      n_observations = n(),
      .groups = "drop"
    ) %>%
    mutate(
      significant_t_test = t_test_p_value < 0.05,
      significant_wilcox = wilcox_p_value < 0.05,
      effect_size_interpretation = case_when(
        abs(cohens_d) < 0.2 ~ "Small",
        abs(cohens_d) < 0.5 ~ "Small",
        abs(cohens_d) < 0.8 ~ "Medium", 
        TRUE ~ "Large"
      )
    ) %>%
    arrange(abs(mean_difference)) %>%
    arrange(desc(abs(cohens_d)))
  
  print(statistical_tests %>%
          select(Index, mean_difference, cohens_d, effect_size_interpretation, 
                 significant_t_test, significant_wilcox))
  
  return(statistical_tests)
}

# Create distribution plots
cat("\nCreating bootstrap distribution visualizations...\n")
distribution_plots <- create_bootstrap_distribution_plots(bootstrap_results)

# Show plots
print(distribution_plots$overall)
print(distribution_plots$individual)
print(distribution_plots$violin)
if(!is.null(distribution_plots$ridge)) {
  print(distribution_plots$ridge)
}

# Statistical comparison
statistical_comparison <- compare_distributions(bootstrap_results)

# Original plots
# กราฟเปรียบเทียบ RMSE ระหว่าง Plot_type
comparison_plot <- summary_results %>%
  ggplot(aes(x = reorder(Index, mean_RMSE), y = mean_RMSE, color = Plot_type)) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = RMSE_ci_lower, ymax = RMSE_ci_upper), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  facet_wrap(~Plot_type, scales = "free_y") +
  coord_flip() +
  labs(title = "Bootstrap Mixed Effects Model Performance",
       subtitle = "RMSE with 95% Bootstrap Confidence Intervals",
       x = "Vegetation Index", 
       y = "RMSE",
       color = "Plot Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(comparison_plot)

# Stability plot (CV)
stability_plot <- summary_results %>%
  ggplot(aes(x = reorder(Index, RMSE_cv), y = RMSE_cv, fill = Plot_type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Model Stability (Coefficient of Variation)",
       subtitle = "Lower CV = More stable across bootstrap iterations",
       x = "Vegetation Index", 
       y = "Coefficient of Variation (CV)",
       fill = "Plot Type") +
  theme_minimal()

print(stability_plot)

cat("\n=== ANALYSIS COMPLETED ===\n")
cat("Bootstrap iterations completed:", max(bootstrap_results$bootstrap_id, na.rm = TRUE), "\n")
cat("Total models fitted:", nrow(bootstrap_results), "\n")

# Save results (optional)
# write_csv(bootstrap_results, "bootstrap_mixed_effects_results.csv")
# write_csv(summary_results, "bootstrap_summary_results.csv")
