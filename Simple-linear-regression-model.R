# ==============================================================================
# Section 3.3. Stage-Specific Performance of Vegetation Indices ----
# ==============================================================================

# To evaluate the stage-specific performance of each VI during cassava growth, we fitted separate ordinary linear regression models, 
# at five distinct stages: 2, 3, 4, 5, and 7 months after planting (MAP). At each growth stage, thirteen models (one per VI) were fitted to determine 
# which indices most accurately predicted ground–LAI at that particular stage of development.

# Load libraries
library(tidyverse)    # For data manipulation & visualization
library(broom)        # For linear models (lm)
library(broom.mixed)  # For linear mixed-effects models (lmer)
library(lmerTest)     # For linear mixed-effects models with p-values
library(MuMIn)        # For marginal & conditional R²
library(performance)  # For model evaluation metric - RMSE
library(tidytext)     # For reorder in faceted plots

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

# Fit linear mixed-effects models (LAI ~ VI with Plot as random intercept) for each vegetation index
m1 <- ALL_DATA %>% 
  group_by(Index) %>% 
  nest() %>%                # Nest data by Index for model fitting
  mutate(
    # Fit mixed-effects model: ground-LAI predicted by VI value, with Plot as random intercept
    model = map(data, ~ lmerTest::lmer(PLOT_MAP_LAI ~ PLOT_MAP_INDEX + (1|Plot), data = .)),    
    # Extract overall model statistics
    model_stats = map(model, ~glance(.x)),
    # Extract fixed effect estimates (intercept and slope)
    fixed_effects = map(model, ~tidy(.x, effects = "fixed")),
    # Extract random effect parameters (variance components)
    random_effects = map(model, ~tidy(.x, effects = "ran_pars")),
    # Compute marginal and conditional R-squared
    r_squared = map(model, ~MuMIn::r.squaredGLMM(.x))
  ) %>%
  mutate(
    # Fixed effect - coefficient, std. error, degree of freedom, p-value
    n_obs = map_int(data, nrow),
    fixed_df = map_dbl(fixed_effects, ~.x$df[.x$term == "PLOT_MAP_INDEX"]), 
    fixed_coef = map_dbl(fixed_effects, ~.x$estimate[.x$term == "PLOT_MAP_INDEX"]),
    fixed_std.error = map_dbl(fixed_effects, ~.x$std.error[.x$term == "PLOT_MAP_INDEX"]),
    fixed_p_value = map_dbl(fixed_effects, ~.x$p.value[.x$term == "PLOT_MAP_INDEX"]),
    # Categorize p-value into significance levels
    sig_p = case_when(
      fixed_p_value < 0.001 ~ "< 0.001",
      fixed_p_value < 0.01 ~ "< 0.01",
      fixed_p_value < 0.05 ~ "< 0.05",
      fixed_p_value < 0.1 ~ "< 0.1",
      TRUE ~ "ns"),
    # R-squared - marginal & conditional
    R2_marginal = map_dbl(r_squared, ~.x[1, "R2m"]),     # R² for fixed effects only
    R2_conditional = map_dbl(r_squared, ~.x[1, "R2c"]),   # R² for fixed + random effects
    # Random effect - standard deviation & variance of Plot intercept
    random_sd = map_dbl(random_effects, ~.x$estimate[.x$term == "sd__(Intercept)"]),
    random_var = map_dbl(random_effects, ~(.x$estimate[.x$term == "sd__(Intercept)"])^2),
    # Model fit criteria
    AIC = map_dbl(model, AIC),
    BIC = map_dbl(model, BIC),
    RMSE = map_dbl(model, ~ performance::performance_rmse(.x, normalized = FALSE))) %>%
  # Label metadata for this model set
  mutate(MAP = "All MAP", Plot_type = "All plots", Model = "Mixed effect (Plot as random)") %>%
  # Select key summary columns for reporting
  select(Index, MAP, Plot_type, Model, n_obs, R2_marginal, fixed_p_value, sig_p, RMSE) %>% 
  # Rank vegetation indices by RMSE (best performing first)
  arrange(RMSE)

# Fit simple linear regression models (LAI ~ VI) for each vegetation index at each MAP separately
m2 <- ALL_DATA %>%
  group_by(Index, MAP) %>%
  nest() %>%                # Nest data by Index and MAP for per-stage model fitting
  mutate(model = map(data, ~ lm(PLOT_MAP_LAI ~ PLOT_MAP_INDEX, data = .))) %>%
  # Extract overall model statistics (R², F-statistic, etc.)
  mutate(results = map(model, glance)) %>% 
  # Extract coefficient-level statistics (estimate, std.error, p-value)
  mutate(results_tidy = map(model, tidy)) %>%
  # Calculate RMSE from model residuals
  mutate(RMSE = map_dbl(model, ~ sqrt(mean(.x$residuals^2)))) %>%
  # Unnest model-level and coefficient-level results into columns
  unnest(results, names_sep = "_") %>%
  unnest(results_tidy, names_sep = "_") %>%
  # Keep only the slope term (VI predictor), exclude the intercept
  filter(results_tidy_term == "PLOT_MAP_INDEX") %>%
  # Convert MAP to character for consistency with m1
  mutate(MAP = as.character(MAP)) %>% 
  # Label metadata for this model set
  mutate(Plot_type = "All plots", Model = "Non-mixed effect") %>% 
  # Select and rename key summary columns for reporting
  select(Index, MAP, Plot_type, Model, Adjusted_R2 = results_adj.r.squared,
    p_value = results_tidy_p.value, RMSE) %>%
  ungroup() %>%
  # Categorize p-value into significance levels
  mutate(sig_p = case_when(
    p_value < 0.001 ~ "< 0.001",
    p_value < 0.01 ~ "< 0.01",
    p_value < 0.05 ~ "< 0.05",
    p_value < 0.1 ~ "< 0.1",
    TRUE ~ "ns")) %>%
  # Rank by MAP and RMSE (best performing first within each MAP)
  arrange(MAP, RMSE)
m2

# ==============================================================================
# Figure 6: Predictive model performance of vegetation indices (VIs) in predicting 
# ground LAI in all cassava plots (N = 47 plots) by Month After Planting (MAP) 
# and across all time points from 2 MAP to 7 MAP (ALL MAP) ----
# ==============================================================================

# Combine mixed-effects results (All MAP) with per-stage linear regression results
lmer_all_plots <- rbind(m1, m2)

# Define custom facet labels for MAP values
MAP_names <- as_labeller(c('2' = 'MAP 2', '3' = 'MAP 3', '4' = 'MAP 4', '5' = 'MAP 5', '7' = 'MAP 7','8' = 'MAP 8', 'All MAP' = "ALL MAP"))

# Plot Lollipop chart of RMSE by vegetation index, faceted by MAP
lmer_all_plots %>%
  group_by(MAP) %>%
  arrange(RMSE) %>%
  ungroup() %>% 
  # Reorder indices by descending RMSE within each MAP facet
  mutate(MAP = as.factor(MAP),
         Index = reorder_within(Index, desc(RMSE), MAP)) %>%
  ggplot(aes(x = Index, y = RMSE)) +
  # Draw lollipop stems from 0 to RMSE value
  geom_segment(aes(xend = Index, yend = 0)) +
  # Draw lollipop points at RMSE value
  geom_point(size = 4, color = "orange") +
  # Clean up reorder_within labels (remove "__MAP" suffix)
  scale_x_reordered() +
  facet_wrap(~MAP, scale = "free_y", labeller = MAP_names) +
  coord_flip() +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 16),
    axis.title = element_text(face = "bold", size = 18)) +
  labs(y = "Root Mean Squared Error (RMSE)", x = "Vegetation indices (VIs)")

dev.off()
