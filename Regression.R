#Import library----
library(tidyverse)
library(tidytext)
library(lubridate)
library(ggpmisc)
library(broom)
library(lmerTest) #For lmer() - fitted mixed effect model
library(broom.mixed) #For tidy() และ glance() - 
library(MuMIn) #For r.squaredGLMM() - marginal & conditional r2
library(performance) #For evaluation metric - AIC,BIC,RMSE

#Ground-LAI data----
ALL_LAI <- read_csv("Data/all-index-pran/Ground-LAI-Clean.csv")
#unique(ALL_LAI$MAP) # MAP at 2,3,4,5,6,7,8
#unique(ALL_LAI$Plot) # 47 Plots

#Ground-LAI grouped by plot and MAP (mean LAI)
PLOT_MAP_LAI <- ALL_LAI %>%
  group_by(Plot,MAP) %>%
  summarise(PLOT_MAP_LAI = mean(LAI)) %>%
  ungroup() %>%
  dplyr::select(Plot, MAP, as.double(), PLOT_MAP_LAI)

#Summarize monthly ground-LAI
meanLAI <- PLOT_MAP_LAI %>% group_by(MAP) %>% summarize(meanLAI=mean(PLOT_MAP_LAI),sdLAI=sd(PLOT_MAP_LAI))

#Plot ground-LAI times-series
ggplot() +
  geom_jitter(data=PLOT_MAP_LAI,aes(x=MAP,y=PLOT_MAP_LAI),width = 0.25) +
  geom_line(data=meanLAI,aes(x=MAP,y=meanLAI)) + labs(x='MAP',y='Ground-LAI value')

#All index data----
ALL_INDEX <- read_csv("Data/all_index_zonal_median.csv")

#All index by averaging plot & MAP (mean index)
#unique(ALL_INDEX$Index) #Index names
PLOT_MAP_ALL_INDEX <- ALL_INDEX %>%
  mutate(Plot = recode(Plot,E7="E07",E8="E08",E9="E09")) %>% 
  group_by(Plot,MAP,Index) %>% 
  summarise(PLOT_MAP_INDEX = mean(Value)) %>%
  filter(Index %in% c("BNDVI","CIG","DVI","EVI","GNDVI","GRVI","NDVI","NDWI","RVI","SAVI","SeLI","VIG","TCARI")) %>%
  ungroup()

PLOT_MAP_ALL_INDEX %>% filter(Index %in% c('GNDVI','NDWI')) %>% select(Plot, MAP, Index,PLOT_MAP_INDEX) %>% pivot_wider(names_from = Index, values_from = PLOT_MAP_INDEX)

PLOT_MAP_ALL_INDEX %>% 
  ggplot(aes(x=Index,y=PLOT_MAP_INDEX))+
  geom_boxplot()

#Combine VIs and ground-LAI----
unique(PLOT_MAP_ALL_INDEX$Plot)
unique(PLOT_MAP_LAI$Plot)
unique(PLOT_MAP_ALL_INDEX$Index)

ALL_INDEX_LAI <- left_join(PLOT_MAP_LAI,PLOT_MAP_ALL_INDEX, by=c("Plot", "MAP"))

#Create column 'Plot_type' for cassava plots without tree (19 plots)
without <- ALL_INDEX_LAI %>% filter(Plot %in% c("S32","S15","S14","S23","S17","S09","S20","S18","S19","S22","S26","S27","S28","S33","S29","S30","K13","K14","K15")) %>% mutate(Plot_type = "Without tree")

#Create column 'Plot_type' for cassava plots with tree (28 plots)
with <- ALL_INDEX_LAI %>% filter(!Plot %in% c("S32","S15","S14","S23","S17","S09","S20","S18","S19","S22","S26","S27","S28","S33","S29","S30","K13","K14","K15")) %>% mutate(Plot_type = "With tree")

#Combine data with LAI and VIs
ALL_INDEX_LAI <- rbind(without,with)
write_csv(ALL_INDEX_LAI,"Data/all_index_lai_plot_map.csv")

ALL_INDEX_LAI %>% group_by(Index) %>% summarise(mean=mean(PLOT_MAP_INDEX),sd=sd(PLOT_MAP_INDEX))

mean(ALL_INDEX_LAI$PLOT_MAP_LAI)
sd(ALL_INDEX_LAI$PLOT_MAP_LAI)

#19 plots without trees
length(unique(ALL_INDEX_LAI[ALL_INDEX_LAI$Plot_type=="Without tree",]$Plot))
#28 plots with trees
length(unique(ALL_INDEX_LAI[ALL_INDEX_LAI$Plot_type=="With tree",]$Plot))

#Linear regression----
##Results (m1): all cassava plots----
m1 <- ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(MAP < 8) %>% 
  group_by(Index) %>%
  nest() %>% 
  mutate(model = map(data, ~ lm(PLOT_MAP_LAI ~ PLOT_MAP_INDEX, data = .))) %>%
  mutate(results = map(model, glance)) %>% 
  mutate(RMSE = map(model,~sqrt(mean(.x$residuals^2)))) %>% 
  unnest(c("results","RMSE")) %>% 
  arrange(desc(adj.r.squared)) %>%
  mutate(MAP = "All MAP", Plot_type = "All plots", Model = "Non-mixed effect") %>% 
  select(Index,MAP,Plot_type,Model,adj.r.squared,RMSE,p.value) %>% 
  ungroup()
m1
#write_csv(m1,"Results/results-non-mixed-effect-all-plots.csv")

#Create the rank from 'm4' model performance (mixed-effects models)
r2.rank <- m1$Index

ALL_INDEX_LAI %>% 
  drop_na() %>%
  filter(MAP < 8) %>% 
  ggplot() +
  geom_point(aes(x = PLOT_MAP_INDEX, y = PLOT_MAP_LAI, color = factor(MAP)), alpha = 0.6) +
  geom_smooth(aes(x = PLOT_MAP_INDEX, y = PLOT_MAP_LAI), method = "lm", se = F, color = "black", linewidth = 0.5) +
  geom_text(data = m1, 
            aes(x = -Inf, y = Inf, 
                label = paste("R² =", round(adj.r.squared, 2))), 
            hjust = -0.1, vjust = 1.5, size = 4.5, fontface = "bold", color = "black") +
  geom_text(data = m1, 
            aes(x = -Inf, y = Inf, 
                label = paste("RMSE =", round(RMSE, 2))), 
            hjust = -0.07, vjust = 3.0, size = 4.5, fontface = "bold", color = "black") +
  geom_text(data = m1, 
            aes(x = -Inf, y = Inf, 
                label = paste("p =", format(p.value, digits = 2, scientific = TRUE))), 
            hjust = -0.1, vjust = 4.5, size = 4.5, fontface = "bold", color = "black") +
  facet_wrap(~factor(Index, levels = r2.rank), scales = "free_x") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 15),
        strip.text = element_text(face = "bold", size = 18),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 15)) +
  labs(x = "Vegetation indices (VIs) value", 
       y = "Ground-LAI value", 
       color = "MAP")

#ggsave("vegetation_indices_performance.pdf", width = 12, height = 8, dpi = 300)

m12 <- ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(MAP < 8) %>% 
  group_by(Index, MAP) %>%
  nest() %>% 
  mutate(model = map(data, ~ lm(PLOT_MAP_LAI ~ PLOT_MAP_INDEX, data = .))) %>%
  mutate(
    adj.r.squared = map_dbl(model, ~ summary(.x)$adj.r.squared),
    RMSE = map_dbl(model, ~ sqrt(mean(.x$residuals^2))),
    p.value = map_dbl(model, ~ summary(.x)$coefficients[2, 4])  # p-value for slope coefficient
  ) %>%
  arrange(desc(adj.r.squared)) %>%
  mutate(Plot_type = "All plots", Model = "Non-mixed effect") %>% 
  select(Index,MAP,Plot_type,Model,adj.r.squared,RMSE,p.value) %>% ungroup()
m12

lm_all_plots <- rbind(m1,m12)

MAP_names <- as_labeller(c('2'='MAP 2','3'='MAP 3','4'='MAP 4','5'='MAP 5','7'='MAP 7','8'='MAP 8','All MAP'="ALL MAP"))

lm_all_plots %>%
  ungroup() %>% 
  mutate(MAP = as.factor(MAP),
         Index = reorder_within(Index,adj.r.squared,MAP)) %>% 
  ggplot(aes(x=Index,y=adj.r.squared))+
  geom_col(show.legend = F,fill="#619cff") +
  scale_x_reordered() +
  facet_wrap(~MAP,scale="free_y",labeller=MAP_names) +
  coord_flip() +
  theme_bw() +
  theme(axis.title.x = element_text(size=15,face = "bold"),
        axis.title.y = element_text(size=15,face = "bold"),
        axis.text = element_text(size=15),
        strip.text = element_text(face = "bold",size=18),
        legend.title = element_text(face = "bold",size=15),
        legend.text = element_text(size=15))+
  labs(y="Adjusted R-squared",x="Vegetation indices (VIs)")

#ggsave("vegetation_indices_MAP_performance.pdf", width = 12, height = 8, dpi = 300)

lm_all_plots %>% filter(MAP=="7") %>% arrange(desc(adj.r.squared))

##Results (m2): plots without trees----
m2 <- ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(Plot_type=="Without tree") %>% 
  filter(MAP<8) %>%
  group_by(Index) %>%
  nest() %>% 
  mutate(model=map(data,~lm(PLOT_MAP_LAI~PLOT_MAP_INDEX,data=.))) %>%
  mutate(results=map(model,glance)) %>% 
  mutate(RMSE = map(model,~sqrt(mean(.x$residuals^2)))) %>% 
  unnest(c("results","RMSE")) %>% 
  arrange(RMSE) %>%
  mutate(MAP = "All MAP", Plot_type = "Without tree",Model = "Non-mixed effect") %>% 
  select(Index,MAP,Plot_type,Model,RMSE)
m2
#write_csv(m2,"Results/results-non-mixed-effect-without-tree.csv")

##Results (m3): plots with trees----
m3 <- ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(Plot_type == "With tree") %>% 
  filter(MAP < 8) %>%
  group_by(Index) %>%
  nest() %>% 
  mutate(model = map(data, ~ lm(PLOT_MAP_LAI ~ PLOT_MAP_INDEX, data = .))) %>%
  mutate(results = map(model, glance)) %>% 
  mutate(RMSE = map(model,~sqrt(mean(.x$residuals^2)))) %>% 
  unnest(c("results","RMSE")) %>% 
  arrange(RMSE) %>%
  mutate(MAP = "All MAP", Plot_type = "With tree",Model = "Non-mixed effect") %>% 
  select(Index,MAP,Plot_type,Model,RMSE)
m3
#write_csv(m3,"Results/results-non-mixed-effect-with-tree.csv")

#Mixed effects model----
#Mixed-effects model with 'plot' random effect
##Results (m4): all cassava plots----
m4 <- ALL_INDEX_LAI %>% 
  drop_na() %>%
  filter(MAP < 8) %>% 
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
        p_value_sig = ,
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
  #select(Index,MAP,Plot_type,Model,n_obs,fixed_df,fixed_coef,fixed_std.error,RMSE) %>%
  select(Index,MAP,Plot_type,Model,n_obs,R2_marginal,fixed_p_value,RMSE) %>% 
  arrange(RMSE)
m4

m4 <- ALL_INDEX_LAI %>% 
  drop_na() %>%
  filter(MAP < 8) %>% 
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
  #select(Index,MAP,Plot_type,Model,n_obs,fixed_df,fixed_coef,fixed_std.error,RMSE) %>%
  select(Index,MAP,Plot_type,Model,n_obs,R2_marginal,fixed_p_value,sig_p,RMSE) %>% 
  arrange(RMSE)
#write_csv(m4,"Results/results-mixed-effect-all-plots.csv")

###Figure 4: Plot regression analysis graphs----
#png("Results/Figure4-lmer-regression-all-plot.png",width = 1000,height = 600,res = 90)
pdf("Results/Figure4-lmer-regression-all-plot.pdf",width = 11,height = 7)

#Create the rank from 'm4' model performance (mixed-effects models)
rmse.rank <- m4$Index
r2 <- m4$R2_marginal

ALL_INDEX_LAI %>% 
  drop_na() %>%
  filter(MAP < 8) %>% 
  ggplot()+
  geom_point(aes(x=PLOT_MAP_INDEX,y=PLOT_MAP_LAI,color=factor(MAP)),alpha=0.6)+
  geom_smooth(aes(x=PLOT_MAP_INDEX,y=PLOT_MAP_LAI), method = "lm",se=F, color="black", linewidth=0.5) +
  geom_text(data = m4, 
            aes(x = -Inf, y = Inf, 
                label = paste0("R² = ", round(R2_marginal, 3), ", p ", sig_p)), 
            hjust = -0.02, vjust = 1.2, size = 4, fontface = "bold", color = "black") +
  geom_text(data = m4, 
            aes(x = -Inf, y = Inf, 
                label = paste0("RMSE = ", round(RMSE, 3))),
            hjust = -0.02, vjust = 3.0, size = 4, fontface = "bold", color = "black") +
  facet_wrap(~factor(Index, rmse.rank),scale="free_x")+
  theme_bw()+
  theme(axis.title = element_text(face = "bold",size=12),
        strip.text = element_text(face = "bold",size=12),
        legend.title = element_text(face = "bold",size=12),
        legend.text = element_text(size=12)) +
  labs(x="Vegetation indices (VIs) value",y="Ground-LAI value",color="MAP")
dev.off()

##Results (m5): plots without trees----
m5 <- ALL_INDEX_LAI %>% 
  drop_na() %>% 
  filter(Plot_type == "Without tree") %>% 
  filter(MAP < 8) %>% 
  group_by(Index) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~ lmerTest::lmer(PLOT_MAP_LAI ~ PLOT_MAP_INDEX + (1|Plot), data = .)),    
    model_stats = map(model, ~glance(.x)),
    fixed_effects = map(model, ~tidy(.x, effects = "fixed")),
    random_effects = map(model, ~tidy(.x, effects = "ran_pars"))
  ) %>%
  mutate(
    # fixed effect - coefficient, std. error, degree of freedom, p-value
    n_obs = map_int(data, nrow),
    fixed_df = map_dbl(fixed_effects, ~.x$df[.x$term == "PLOT_MAP_INDEX"]), 
    fixed_coef = map_dbl(fixed_effects, ~.x$estimate[.x$term == "PLOT_MAP_INDEX"]),
    fixed_se = map_dbl(fixed_effects, ~.x$std.error[.x$term == "PLOT_MAP_INDEX"]),
    fixed_p_value = map_dbl(fixed_effects, ~.x$p.value[.x$term == "PLOT_MAP_INDEX"]),
    #random effect - standard deviation & variance
    random_sd = map_dbl(random_effects, ~.x$estimate[.x$term == "sd__(Intercept)"]),
    random_var = map_dbl(random_effects, ~(.x$estimate[.x$term == "sd__(Intercept)"])^2),
    AIC = map_dbl(model, AIC),
    BIC = map_dbl(model, BIC),
    RMSE = map_dbl(model, ~ performance::performance_rmse(.x, normalized = FALSE))) %>%
  mutate(MAP = "All MAP", Plot_type = "Without tree", Model = "Mixed effect (Plot as random)") %>% 
  select(Index,MAP,Plot_type,Model,RMSE,AIC) %>% arrange(RMSE)
m5
#write_csv(m5,"Results/results-mixed-effect-without-tree.csv")

##Results (m6): plots with trees----
m6 <- ALL_INDEX_LAI %>% 
  drop_na() %>% 
  filter(Plot_type == "With tree") %>% 
  filter(MAP < 8) %>% 
  group_by(Index) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~ lmerTest::lmer(PLOT_MAP_LAI ~ PLOT_MAP_INDEX + (1|Plot), data = .)),    
    model_stats = map(model, ~glance(.x)),
    fixed_effects = map(model, ~tidy(.x, effects = "fixed")),
    random_effects = map(model, ~tidy(.x, effects = "ran_pars"))
  ) %>%
  mutate(
    # fixed effect - coefficient, std. error, degree of freedom, p-value
    n_obs = map_int(data, nrow),
    fixed_df = map_dbl(fixed_effects, ~.x$df[.x$term == "PLOT_MAP_INDEX"]), 
    fixed_coef = map_dbl(fixed_effects, ~.x$estimate[.x$term == "PLOT_MAP_INDEX"]),
    fixed_se = map_dbl(fixed_effects, ~.x$std.error[.x$term == "PLOT_MAP_INDEX"]),
    fixed_p_value = map_dbl(fixed_effects, ~.x$p.value[.x$term == "PLOT_MAP_INDEX"]),
    #random effect - standard deviation & variance
    random_sd = map_dbl(random_effects, ~.x$estimate[.x$term == "sd__(Intercept)"]),
    random_var = map_dbl(random_effects, ~(.x$estimate[.x$term == "sd__(Intercept)"])^2),
    AIC = map_dbl(model, AIC),
    BIC = map_dbl(model, BIC),
    RMSE = map_dbl(model, ~ performance::performance_rmse(.x, normalized = FALSE))) %>%
  mutate(MAP = "All MAP", Plot_type = "With tree", Model = "Mixed effect (Plot as random)") %>% 
  select(Index,MAP,Plot_type,Model,RMSE) %>% arrange(RMSE)
m6
#write_csv(m6,"Results/results-mixed-effect-with-tree.csv")

#Grouped by Index & MAP-Linear----
##Results (m7): all cassava plots----
ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(MAP < 8) %>% 
  group_by(Index,MAP) %>% count(Index,MAP)

m7 <- ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(MAP < 8) %>% 
  group_by(Index, MAP) %>%
  nest() %>% 
  mutate(model = map(data, ~ lm(PLOT_MAP_LAI ~ PLOT_MAP_INDEX, data = .))) %>%
  mutate(results = map(model, glance)) %>% 
  mutate(results_tidy = map(model, tidy)) %>%
  mutate(RMSE = map_dbl(model, ~ sqrt(mean(.x$residuals^2)))) %>%
  unnest(results, names_sep = "_") %>%
  unnest(results_tidy, names_sep = "_") %>%
  filter(results_tidy_term == "PLOT_MAP_INDEX") %>%
  mutate(MAP = as.character(MAP)) %>% 
  mutate(Plot_type = "All plots", Model = "Non-mixed effect") %>% 
  select(
    Index,
    MAP,
    Plot_type,
    Model,
    Adjusted_R2 = results_adj.r.squared,
    p_value = results_tidy_p.value,
    RMSE
  ) %>%
  ungroup() %>%
  mutate(sig_p = case_when(
           p_value < 0.001 ~ "< 0.001",
           p_value < 0.01 ~ "< 0.01",
           p_value < 0.05 ~ "< 0.05",
           p_value < 0.1 ~ "< 0.1",
           TRUE ~ "ns")) %>%
  arrange(MAP, RMSE)
m7

###Plot regression analysis graphs by MAP----
####2 MAP----
MAP2 <- m7 %>% filter(MAP=="2") %>% arrange(RMSE)
rmse.rank.map2 <- MAP2$Index
ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(MAP == 2) %>% 
  group_by(Index,MAP) %>% 
  ggplot(aes(x=PLOT_MAP_INDEX,y=PLOT_MAP_LAI,color=Plot)) +
  geom_point() +
  geom_smooth(method = "lm",se=F, color="black", linewidth=0.5) +
  geom_text(data = MAP2, aes(x = -Inf, y = Inf, label = paste("RMSE =", round(RMSE, 3))), 
            hjust = -0.1, vjust = 1.5, size = 3.5, fontface = "bold", color = "black") +
  facet_wrap(~factor(Index, rmse.rank.map2),scale="free_x")+
  theme_bw()+
  theme(axis.title = element_text(face = "bold",size=12),
        strip.text = element_text(face = "bold",size=12),
        legend.title = element_text(face = "bold",size=12),
        legend.position = "none") +
  labs(x="Vegetation indices (VIs) value",y="Ground-LAI value",title = "At MAP 2 - Effectiveness of vegetation indices (VIs) for ground-LAI estimation in cassava plots")
dev.off()

####3 MAP----
MAP3 <- m7 %>% filter(MAP=="3") %>% arrange(RMSE)
rmse.rank.map3 <- MAP3$Index
ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(MAP == 3) %>% 
  group_by(Index,MAP) %>% 
  ggplot(aes(x=PLOT_MAP_INDEX,y=PLOT_MAP_LAI,color=Plot)) +
  geom_point() +
  geom_smooth(method = "lm",se=F, color="black", linewidth=0.5) +
  geom_text(data = MAP3, aes(x = -Inf, y = Inf, label = paste("RMSE =", round(RMSE, 3))), 
            hjust = -0.1, vjust = 1.5, size = 3.5, fontface = "bold", color = "black") +
  facet_wrap(~factor(Index, rmse.rank.map3),scale="free_x")+
  theme_bw()+
  theme(axis.title = element_text(face = "bold",size=12),
        strip.text = element_text(face = "bold",size=12),
        legend.title = element_text(face = "bold",size=12),
        legend.position = "none") +
  labs(x="Vegetation indices (VIs) value",y="Ground-LAI value",title = "At MAP 3 - Effectiveness of vegetation indices (VIs) for ground-LAI estimation in cassava plots")
dev.off()

####Example RVI----
RVI <- m7 %>% filter(Index == "RVI") %>% arrange(RMSE)
rmse.rank.rvi <- RVI$Index
ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(Index == "RVI") %>% 
  group_by(Index,MAP) %>% 
  ggplot(aes(x=PLOT_MAP_INDEX,y=PLOT_MAP_LAI)) +
  geom_point() +
  geom_smooth(method = "lm",se=F, color="black", linewidth=0.5) +
  facet_wrap(~factor(MAP),scale="free_x")+
  theme_bw()+
  theme(axis.title = element_text(face = "bold",size=12),
        strip.text = element_text(face = "bold",size=12),
        legend.title = element_text(face = "bold",size=12),
        legend.position = "none") +
  labs(x="Vegetation indices (VIs) value",y="Ground-LAI value", title = "Effectiveness of RVI for ground-LAI estimation in cassava plots")
dev.off()



###Export linear results by MAP file.csv----
MAP2 <- m7 %>% filter(MAP=="2") %>% arrange(RMSE)
MAP3 <- m7 %>% filter(MAP=="3") %>% arrange(RMSE)
MAP4 <- m7 %>% filter(MAP=="4") %>% arrange(RMSE)
MAP5 <- m7 %>% filter(MAP=="5") %>% arrange(RMSE)
MAP7 <- m7 %>% filter(MAP=="7") %>% arrange(RMSE)

#write_csv(MAP2,"Results/results-linear-MAP2-all-plots.csv")
#write_csv(MAP3,"Results/results-linear-MAP3-all-plots.csv")
#write_csv(MAP4,"Results/results-linear-MAP4-all-plots.csv")
#write_csv(MAP5,"Results/results-linear-MAP5-all-plots.csv")
#write_csv(MAP7,"Results/results-linear-MAP7-all-plots.csv")

lm_all_plots <- rbind(m7,m4) #combine data of each MAP and All MAP

#lm_all_plots %>% filter(MAP=="7") %>% arrange(RMSE)

MAP_names <- as_labeller(c('2'='MAP 2','3'='MAP 3','4'='MAP 4','5'='MAP 5','7'='MAP 7','8'='MAP 8','All MAP'="ALL MAP"))

lm_all_plots %>% filter(MAP=="All MAP") %>% arrange(RMSE)

##Results (m8): plot without trees----
ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(MAP < 8) %>%
  filter(Plot_type == "Without tree") %>% distinct(Plot) # 19 plots without trees

m8 <- ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(MAP < 8) %>%
  filter(Plot_type == "Without tree") %>% 
  group_by(Index,MAP) %>%
  nest() %>% 
  mutate(model = map(data, ~ lm(PLOT_MAP_LAI ~ PLOT_MAP_INDEX, data = .))) %>%
  mutate(results = map(model, glance)) %>% 
  mutate(RMSE = map(model,~sqrt(mean(.x$residuals^2)))) %>% 
  unnest(c("results","RMSE")) %>% 
  arrange(RMSE) %>%
  mutate(MAP = as.character(MAP)) %>% 
  mutate(Plot_type = "Without tree", Model = "Non-mixed effect") %>%
  select(Index,MAP,Plot_type,Model,RMSE)
m8

lm_without_tree <- rbind(m8,m2) #combine data of each MAP and All MAP

lm_without_tree %>% filter(MAP=='7') %>% arrange(RMSE)

##Results (m9): plot with trees----
ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(MAP < 8) %>%
  filter(Plot_type == "With tree") %>% distinct(Plot) #28 plots with tree

m9 <- ALL_INDEX_LAI %>%
  drop_na() %>%
  filter(MAP < 8) %>%
  filter(Plot_type == "With tree") %>% 
  group_by(Index,MAP) %>%
  nest() %>% 
  mutate(model = map(data, ~ lm(PLOT_MAP_LAI ~ PLOT_MAP_INDEX, data = .))) %>%
  mutate(results = map(model, glance)) %>% 
  mutate(RMSE = map(model,~sqrt(mean(.x$residuals^2)))) %>% 
  unnest(c("results","RMSE")) %>% 
  arrange(RMSE) %>%
  mutate(MAP = as.character(MAP)) %>% 
  mutate(Plot_type = "With tree", Model = "Non-mixed effect") %>%
  select(Index,MAP,Plot_type,Model,RMSE)
m9

lm_with_tree <- rbind(m9,m6) #combine data of each MAP and All MAP

MAP_names <- as_labeller(c('2'='MAP 2','3'='MAP 3','4'='MAP 4','5'='MAP 5','7'='MAP 7','8'='MAP 8','All MAP'="ALL MAP"))

lm_with_tree %>%
  group_by(MAP) %>%
  arrange(RMSE) %>%
  ungroup() %>% 
  mutate(MAP = as.factor(MAP),
         Index = reorder_within(Index,desc(RMSE),MAP)) %>%
  ggplot(aes(x=Index,y=RMSE))+
  geom_col(show.legend = F,fill="#9ecae1") +
  scale_x_reordered() +
  facet_wrap(~MAP,scale="free_y",labeller=MAP_names) +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(face = "bold",size=12),
        strip.text = element_text(face = "bold",size=12),
        legend.title = element_text(face = "bold",size=12),
        legend.text = element_text(size=12)) +
  labs(y="Root Mean Squared (RMSE)",x="Vehetation indices (VIs)") +
  #ggtitle("Predictive capability of vegetation indices (VIs) for ground-LAI estimation in 'cassava plots with tree' in each Month After Planting (MAP)") +
  NULL
dev.off()

#Comparison of the Best VI Performance across Plot Types in Each Month----
lmer_all_plots <- rbind(m7,m4) #combine data of each MAP and All MAP
lmer_without_tree <- rbind(m8,m5) #combine data of each MAP and All MAP
lmer_with_tree <- rbind(m9,m6) #combine data of each MAP and All MAP

lmer_all_plots %>% filter(MAP=='All MAP') %>% arrange(RMSE)
lmer_without_tree %>% filter(MAP=='All MAP') %>% arrange(RMSE)
lmer_with_tree %>% filter(MAP=='All MAP') %>% arrange(RMSE)

ALL_PLOTS <- lmer_all_plots %>%
  group_by(MAP) %>%
  filter(RMSE == min(RMSE)) 
  #filter(!(Index=="NDWI"&MAP=="ALL MAP")) %>%
  #mutate(Index=if_else(Index=="GNDVI"&MAP=="ALL MAP", "GNDVI & NDWI", Index)) %>% 
  #filter(!(Index=='NDWI'&MAP=='3')) %>% 
  #mutate(Index=if_else(Index=='GNDVI'&MAP=='3','GNDVI & NDWI',Index))

WITHOUT_TREES <- lmer_without_tree %>% 
  group_by(MAP) %>%
  filter(RMSE == min(RMSE))  
  #filter(!(Index=='NDWI'&MAP=='7')) %>% 
  #mutate(Index=if_else(Index=='GNDVI'&MAP=='7','GNDVI & NDWI',Index)) %>%
  #filter(!(Index=='NDWI'&MAP=='4')) %>% 
  #mutate(Index=if_else(Index=='GNDVI'&MAP=='4','GNDVI & NDWI',Index))

WITH_TREES <- lmer_with_tree %>% 
  group_by(MAP) %>%
  filter(RMSE == min(RMSE))
#filter(!(Index=='NDWI'&MAP=='7')) %>% 
#mutate(Index=if_else(Index=='GNDVI'&MAP=='7','GNDVI & NDWI',Index)) %>%
#filter(!(Index=='NDWI'&MAP=='4')) %>% 
#mutate(Index=if_else(Index=='GNDVI'&MAP=='4','GNDVI & NDWI',Index))

BEST_INDEX <- rbind(ALL_PLOTS,WITHOUT_TREES,WITH_TREES) %>% print(n="all")

unique(BEST_INDEX$MAP)

BEST_INDEX %>% 
  select(Plot_type,MAP,Index,RMSE) %>% 
  arrange(Plot_type,MAP) %>% print(n='all')

#Figure 5: Model performance of all plots----
png("Results/Figure/Figure5-Model-performance-rmse-by-MAP-and-All-MAP.png",width = 1000,height = 600,res = 90)
#pdf("Results/Figure/Figure5-Model-performance-rmse-by-MAP-and-All-MAP.pdf",width = 11,height = 7)

lmer_all_plots <- rbind(m7,m4) #combine data of each MAP and All MAP

MAP_names <- as_labeller(c('2'='MAP 2','3'='MAP 3','4'='MAP 4','5'='MAP 5','7'='MAP 7','8'='MAP 8','All MAP'='ALL MAP'))

lmer_all_plots %>%
  group_by(MAP) %>%
  arrange(RMSE) %>%
  ungroup() %>% 
  mutate(MAP = as.factor(MAP),
         Index = reorder_within(Index,desc(RMSE),MAP)) %>%
  ggplot(aes(x=Index,y=RMSE))+
  geom_col(show.legend = F,fill="#F8766D") +
  scale_x_reordered() +
  facet_wrap(~MAP,scale="free_y",labeller=MAP_names) +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(face = "bold",size=12),
        strip.text = element_text(face = "bold",size=12),
        legend.title = element_text(face = "bold",size=12),
        legend.text = element_text(size=12)) +
  labs(y="Root Mean Squared Error (RMSE)",x="Vegetation indices (VIs)")+
  #ggtitle("Predictive capability of vegetation indices (VIs) for ground-LAI estimation in 'all cassava plots' in each Month After Planting (MAP)") +
  NULL
dev.off()

#png("Results/Figure/Figure5-1-Model-performance-rmse-by-MAP-and-All-MAP.png",width = 1000,height = 600,res = 90)
pdf("Results/Figure/Figure5-1-Model-performance-rmse-by-MAP-and-All-MAP.pdf",width = 11,height = 7)
lmer_all_plots %>%
  group_by(MAP) %>%
  arrange(RMSE) %>%
  ungroup() %>% 
  mutate(MAP = as.factor(MAP),
         Index = reorder_within(Index,desc(RMSE),MAP)) %>%
  ggplot(aes(x=Index,y=RMSE))+
  geom_segment(aes(xend=Index, yend=0)) +
  geom_point(size=4, color="orange") +
  scale_x_reordered() +
  facet_wrap(~MAP,scale="free_y",labeller=MAP_names) +
  coord_flip() +
  theme_bw() +
  theme(
        strip.text = element_text(face = "bold",size=16),
        legend.text = element_text(size=12),
        axis.text = element_text(size = 16),
        axis.title = element_text(face = "bold",size=18)) +
  labs(y="Root Mean Squared Error (RMSE)",x="Vegetation indices (VIs)")
dev.off()

#Figure 6: Model performance of plots without trees----
#png("Results/Figure/Figure6-Model-performance-rmse-by-MAP-and-withoutTree.png",width = 1000,height = 600,res = 90)
#pdf("Results/Figure/Figure6-Model-performance-rmse-by-MAP-and-withoutTree.pdf",width = 11,height = 7)

lmer_without_tree <- rbind(m8,m5) #combine data of each MAP and All MAP

MAP_names <- as_labeller(c('2'='MAP 2','3'='MAP 3','4'='MAP 4','5'='MAP 5','7'='MAP 7','8'='MAP 8','All MAP'="ALL MAP"))

lmer_without_tree %>%
  ungroup() %>% 
  mutate(MAP = as.factor(MAP),
         Index = reorder_within(Index,desc(RMSE),MAP)) %>% 
  ggplot(aes(x=Index,y=RMSE))+
  geom_col(show.legend = F,fill="#619cff") +
  scale_x_reordered() +
  facet_wrap(~MAP,scale="free_y",labeller=MAP_names) +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(face = "bold",size=12),
        strip.text = element_text(face = "bold",size=12),
        legend.title = element_text(face = "bold",size=12),
        legend.text = element_text(size=12)) +
  labs(y="Root Mean Squared Error (RMSE)",x="Vegetation indices (VIs)")+
  #ggtitle("Predictive capability of vegetation indices (VIs) for ground-LAI estimation in 'cassava plots without tree' in each Month After Planting (MAP)") +
  NULL
dev.off()

#Figure 7: Model performance of plots with trees----
png("Results/Figure/Figure7-Model-performance-rmse-by-MAP-and-withTree.png",width = 1000,height = 600,res = 90)
pdf("Results/Figure/Figure7-Model-performance-rmse-by-MAP-and-withTree.pdf",width = 11,height = 7)

lmer_with_tree <- rbind(m9,m6) #combine data of each MAP and All MAP

MAP_names <- as_labeller(c('2'='MAP 2','3'='MAP 3','4'='MAP 4','5'='MAP 5','7'='MAP 7','8'='MAP 8','All MAP'="ALL MAP"))

lmer_with_tree %>%
  ungroup() %>% 
  mutate(MAP = as.factor(MAP),
         Index = reorder_within(Index,desc(RMSE),MAP)) %>% 
  ggplot(aes(x=Index,y=RMSE))+
  geom_col(show.legend = F,fill="#00BA38") +
  scale_x_reordered() +
  facet_wrap(~MAP,scale="free_y",labeller=MAP_names) +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(face = "bold",size=12),
        strip.text = element_text(face = "bold",size=12),
        legend.title = element_text(face = "bold",size=12),
        legend.text = element_text(size=12)) +
  labs(y="Root Mean Squared Error (RMSE)",x="Vegetation indices (VIs)")
dev.off()

#Figure 8: The best VIs----
png("Results/Figure/Figure8-the-best-VIs.png",width = 1000,height = 600,res = 90)
pdf("Results/Figure/Figure8-the-best-VIs.pdf",width = 11,height = 7)

BEST_INDEX %>% 
  select(MAP,Index,Plot_type,RMSE) %>%
  filter(Index != 'NDWI') %>% 
  ggplot(aes(x = RMSE, y = MAP, fill = factor(Plot_type, levels = c("All plots","With tree","Without tree")))) +
  geom_col(position = position_dodge(width = 0.7),width = 0.7) +
  geom_text(aes(label = Index),
            position = position_dodge(width = 0.7), hjust = -0.15, size = 4, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values=c("All plots" = "#F8766D", "With tree" = "#00BA38", "Without tree" = "#619cff")) +
  labs(x = "Root Mean Squared Error (RMSE)", 
       y = "MAP (Month After Planting)",
       fill = "Plot Types") + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold",size=12),
        axis.text = element_text(size = 13),
        strip.text = element_text(face = "bold",size=12),
        legend.title = element_text(face = "bold",size=12),
        legend.text = element_text(face = "bold",size=12),
        legend.position = "bottom") +
  #ggtitle('Comparison of The Best VI Performance across Plot Types in Each Month') +
  NULL
dev.off()

lmer_all_plots <- rbind(m7,m4) #combine data of each MAP and All MAP
lmer_without_tree <- rbind(m8,m5) #combine data of each MAP and All MAP
lmer_with_tree <- rbind(m9,m6) #combine data of each MAP and All MAP

lmer_all_plots %>% filter(MAP=='All MAP') %>% arrange(RMSE)
lmer_without_tree %>% filter(MAP=='All MAP') %>% arrange(RMSE)
lmer_with_tree %>% filter(MAP=='All MAP') %>% arrange(RMSE)

lmer_all_index <- rbind(lmer_all_plots,lmer_without_tree,lmer_with_tree)

lmer_gndvi <- lmer_all_index %>% filter(Index=="GNDVI")

lmer_all_index %>% 
  ggplot(aes(x = RMSE, y = MAP, fill = factor(Plot_type, levels = c("All plots","With tree","Without tree")))) +
  geom_col(position = position_dodge(width = 0.7),width = 0.7) +
  facet_wrap(~Index)+
  #geom_text(aes(label = Index), position = position_dodge(width = 0.7), hjust = -0.15, size = 4, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values=c("All plots" = "#F8766D", "With tree" = "#00BA38", "Without tree" = "#619cff")) +
  labs(x = "Root Mean Squared Error (RMSE)", 
       y = "MAP (Month After Planting)",
       fill = "Plot Types") + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold",size=12),
        axis.text = element_text(size = 13),
        strip.text = element_text(face = "bold",size=12),
        legend.title = element_text(face = "bold",size=12),
        legend.text = element_text(face = "bold",size=12),
        legend.position = "bottom") +
  ggtitle('Comparison of Performance Based on GNDVI across Plot Types in Each Month') +
  NULL

#Mixed effects models (2)----
#Mixed-effects model with 'plot' as random effect, 
#When 'plot type' & 'growth stages' as interaction 
##Results (m9): all cassava plots----
ALL_INDEX_LAI %>% 
  drop_na() %>%
  filter(MAP < 8) %>%
  mutate(MAP = factor(MAP)) %>% 
  group_by(Index) %>% 
  filter(Index == "GNDVI") %>% 
  lmerTest::lmer(PLOT_MAP_LAI ~ PLOT_MAP_INDEX * Plot_type * MAP + (1|Plot),data=.) %>% 
  summary()
  
m9 <- ALL_INDEX_LAI %>% 
  drop_na() %>%
  filter(MAP < 8) %>%
  mutate(MAP = factor(MAP)) %>% 
  group_by(Index) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~ lmerTest::lmer(PLOT_MAP_LAI ~ PLOT_MAP_INDEX * Plot_type * MAP + (1|Plot), data = .)),
    model_stats = map(model, ~glance(.x)),
    fixed_effects = map(model, ~tidy(.x, effects = "fixed")),
    random_effects = map(model, ~tidy(.x, effects = "ran_pars"))
  ) %>%
  mutate(
    # fixed effect - coefficient, std. error, degree of freedom, p-value
    n_obs = map_int(data, nrow),
    fixed_df = map_dbl(fixed_effects, ~.x$df[.x$term == "PLOT_MAP_INDEX"]), 
    fixed_coef = map_dbl(fixed_effects, ~.x$estimate[.x$term == "PLOT_MAP_INDEX"]),
    fixed_std.error = map_dbl(fixed_effects, ~.x$std.error[.x$term == "PLOT_MAP_INDEX"]),
    fixed_p_value = map_dbl(fixed_effects, ~.x$p.value[.x$term == "PLOT_MAP_INDEX"]),
    #random effect - standard deviation & variance
    random_sd = map_dbl(random_effects, ~.x$estimate[.x$term == "sd__(Intercept)"]),
    random_var = map_dbl(random_effects, ~(.x$estimate[.x$term == "sd__(Intercept)"])^2),
    AIC = map_dbl(model, AIC),
    BIC = map_dbl(model, BIC),
    RMSE = map_dbl(model, ~ performance::performance_rmse(.x, normalized = FALSE))) %>%
  mutate(MAP = "All MAP", Plot_type = "All plots", Model = "Mixed effect (Plot as random)") %>%
  #select(Index,MAP,Plot_type,Model,n_obs,fixed_df,fixed_coef,fixed_std.error,RMSE) %>%
  select(Index,MAP,Plot_type,Model,RMSE) %>% 
  arrange(RMSE)
m9

#Examine the effect of plot types based on 'GNDVI' index----
lmer_all_plots <- rbind(m7,m4) #combine data of each MAP and All MAP
lmer_without_tree <- rbind(m8,m5) #combine data of each MAP and All MAP
lmer_with_tree <- rbind(m9,m6) #combine data of each MAP and All MAP
