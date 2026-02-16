# ==============================================================================
# Section 3.2. Predictive Performance of Vegetation Indices Across All Growth Stages 
# ==============================================================================

# To evaluate the ability of Sentinel-2-derived-VIs in predicting ground–LAI, we constructed linear mixed-effects models. For each of the thirteen VIs, 
# an individual model was fitted using monthly VI values as fixed-effect predictors and monthly ground–LAI as the response variable.

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

# --- 1. Linear mixed effects model with 'plot' random effect ----

# Fit linear mixed-effects models (LAI ~ VI with Plot as random intercept) for each vegetation index
m1 <- ALL_DATA %>% 
  group_by(Index) %>% 
  nest() %>%             # Nest data by Index for model fitting
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
m1

# =============================================================================================================================================
# Figure 5: Vegetation indices (VIs) predictive capability for ground–LAI estimation in all cassava plots (N = 47 plots) across all time points
# =============================================================================================================================================

# Create the rank from 'm1' model performance (mixed-effects models)
rmse.rank <- m1$Index

# Plot graph
ggplot(ALL_DATA) +
  geom_point(aes(x = PLOT_MAP_INDEX, y = PLOT_MAP_LAI, color = factor(MAP)), alpha = 0.6)+
  geom_smooth(aes(x = PLOT_MAP_INDEX, y = PLOT_MAP_LAI), method = "lm", se = F, color = "black", linewidth = 0.5) +
  geom_text(data = m1, aes(x = -Inf, y = Inf,
            label = paste0("R² = ", round(R2_marginal, 3), ", p ", sig_p)),
            hjust = -0.02, vjust = 1.2, size = 4, fontface = "bold", 
            color = "black") +
  geom_text(data = m1, aes(x = -Inf, y = Inf, 
            label = paste0("RMSE = ", round(RMSE, 3))),
            hjust = -0.02, vjust = 3.0, size = 4, fontface = "bold",
            color = "black") +
  facet_wrap(~factor(Index, rmse.rank),scale="free_x") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "Vegetation indices (VIs) value", y = "Ground-LAI value", color = "MAP")

dev.off()
