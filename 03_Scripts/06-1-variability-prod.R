##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create total result figure of variability and cumulative PP

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-data.R")

#### Load/wrangle simulated data ####

results_phase_df <- import_data(path = "02_Data/result-phase.rds")

results_noise_df <- import_data(path = "02_Data/result-noise.rds")

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))

#### Split data #### 

results_combined_list <- dplyr::filter(results_combined_df, measure == "gamma",
                                       part %in% c("ag_production", "bg_production", 
                                                   "ttl_production")) |>
  dplyr::group_by(scenario, part) |>
  dplyr::group_split()

#### Fit regression model ####

full_lm_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
  
  df_temp_stand <- dplyr::select(df_temp, value.prod, value.cv, biotic, abiotic,
                                 value.biomass, nutrient_input) |> 
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) |>
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
  
  purrr::map_dfr(c(cv = "value.cv", separated = "biotic * abiotic * value.biomass * nutrient_input"), 
                 y = df_temp_stand, function(x, y) {
    
    lm_temp <- lm(as.formula(paste0("value.prod ~ ", paste0(x))), data = y) 
    
    broom::tidy(lm_temp) |> 
      dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                    measure = unique(df_temp$measure), .before = term) |> 
      dplyr::mutate(p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                                     p.value < 0.05 ~ "*", p.value >= 0.05 ~ ""), 
                    r2 = summary(lm_temp)$adj.r.squared, 
                    direction = dplyr::case_when(estimate < 0.0 ~ "decrease", 
                                                 estimate > 0.0 ~ "increase"))}, 
    .id = "explanatory")}) |> 
  dplyr::mutate(explanatory = factor(explanatory, levels = c("cv", "separated")))

#### Save linear regression results ####

purrr::walk(c("phase", "noise"), function(i) {
  dplyr::filter(full_lm_df, scenario == i, explanatory == "separated") |> 
    suppoRt::save_rds(filename = paste0("lm_variability_prod_", i, ".rds"), 
                      path = "05_Results/", overwrite = FALSE)
})

#### Relative importance R2 ####

rel_importance_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
  
  message("\r> Progress: scenario=", unique(df_temp$scenario), "; part=", unique(df_temp$part), "; measure=",
          unique(df_temp$measure), "\t\t", appendLF = FALSE)
  
  df_temp_stand <- dplyr::select(df_temp, value.prod, biotic, abiotic, 
                                 value.biomass, nutrient_input) |> 
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) |>
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
  
  lm_temp <- lm(value.prod ~ biotic + abiotic + value.biomass + nutrient_input, 
                data = df_temp_stand) 
  
  rel_r2 <- relaimpo::boot.relimp(lm_temp, type = "lmg", b = 500, level = 0.95, fixed = FALSE) |>
    relaimpo::booteval.relimp(bty = "basic")
  
  
  tibble::tibble(
    scenario = unique(df_temp$scenario), part = unique(df_temp$part), measure = unique(df_temp$measure),
    beta = c(rel_r2@namen[-1], "residual"), mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)), 
    lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA), r2 = summary(lm_temp)$adj.r.squared) |> 
    dplyr::mutate(dplyr::across(mean:higher, ~ .x * 100))}) |> 
  dplyr::mutate(beta = factor(beta))

#### Save relative importance results ####

purrr::walk(c("phase", "noise"), function(i) {
  dplyr::filter(rel_importance_df, scenario == i) |> 
    suppoRt::save_rds(filename = paste0("relimp_variability_prod_", i, ".rds"), 
                      path = "05_Results/", overwrite = FALSE)
})

#### Setup ggplot ####

# color_term <- c("biotic" = "#00af73", "abiotic" = "#006d9a", "biotic:abiotic" = "#ffaa3a", 
#                 "cv" = "#ff3b18")
# 
# color_term <- c("biotic" = "#00af73", "abiotic" = "#006d9a",
#                 "biotic:abiotic" = "#ffaa3a", "residuals" = "grey")
# 
# size_point <- 5
# size_line <- 0.75
# size_text <- 2.5
# size_base <- 10.0
# 
# alpha <- 0.5
# 
# width_pos <- 0.65

#### Create ggplot model parameters ####

# gg_coef_scenario <- purrr::map(c("phase", "noise"), function(scenario_i) {
#   
#   ggplot(data = dplyr::filter(regression_df, scenario == scenario_i, 
#                               term %in% c("biotic", "abiotic"))) + 
#     
#     # adding geoms
#     geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
#     geom_linerange(aes(x = pop_n, ymin = 0.0, ymax = estimate, color = term),
#                    position = position_dodge(width = width_pos), size = size_line) +
#     geom_point(aes(x = pop_n, y = estimate, fill = term, color = term, shape = direction),
#                position = position_dodge(width = width_pos), size = size_point * 0.75) +
#     geom_text(aes(x = pop_n, y = estimate, label = p.value, group = term), position = position_dodge(width = width_pos), 
#               vjust = 0.75, size = size_text, color = "white") +
#     
#     # facet gridding
#     facet_grid(rows = dplyr::vars(part), # cols = dplyr::vars(scenario),
#                labeller = labeller(part = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
#                                             "ttl_production" = "Total"), 
#                                    scenario = c("phase" = "Phase scenario","noise" = "Noise scenario"))) +
#     
#     # set scales and labs
#     scale_color_manual(name = "", values = color_term[c(1,2)], 
#                        labels = c("biotic" = "Consumer behavior", "abiotic" =  "Abiotic subsidies", 
#                                   "biotic:abiotic" =  "Interaction Behavior:Subsidies", "residuals" = "Residuals")) +
#     scale_fill_manual(name = "", values =  color_term[c(1,2)],
#                       labels = c("biotic" = "Consumer behavior", "abiotic" =  "Abiotic subsidies", 
#                                  "biotic:abiotic" =  "Interaction Behavior:Subsidies", "residuals" = "Residuals")) +
#     scale_shape_manual(name = "", values = c("decrease" = 25, "increase" = 24)) +
#     coord_flip() +
#     scale_x_discrete(limits = rev(levels(regression_df$pop_n))) +
#     scale_y_continuous(limits = function(x) range(x), labels = function(x) round(x, 2),
#                        breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4)) +
#     
#     # setup theme
#     labs(x = "Population size", y = "Parameter estimate") +
#     guides(color = guide_legend(order = 1, override.aes = list(shape = 17)), 
#            fill = "none", shape = guide_legend(order = 2)) +
#     theme_classic(base_size = size_base) + 
#     theme(axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA),
#           strip.background = element_blank(), strip.text = element_text(hjust = 0.5),
#           legend.position = "bottom",
#           plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt"))
#   
# })
# 
# names(gg_coef_scenario) <- c("phase", "noise")

#### Create ggplot relative importance ####

# gg_relimp_scenario <- purrr::map(c("phase", "noise"), function(scenario_i) {
#   
#   ggplot(data = dplyr::filter(importance_df, scenario == scenario_i)) + 
#     
#   # relative importance bars
#   geom_col(aes(x = pop_n, y = mean * 100, fill = beta)) + 
#   
#   # set scales
#   scale_fill_manual(name = "", values = color_term) +
#   scale_y_continuous(labels = function(x) paste0(x, "%")) + 
#   
#   # facet gridding
#   facet_grid(rows = dplyr::vars(part), # cols = dplyr::vars(scenario),
#              labeller = labeller(part = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
#                                           "ttl_production" = "Total"), 
#                                  scenario = c("phase" = "Phase scenario","noise" = "Noise scenario"))) +
#     
#   # setup theme
#   labs(x = "Population size", y = "Parameter estimate") +
#   theme_classic(base_size = size_base) + 
#   theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.5), 
#         axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA),
#         legend.position = "bottom",
#         plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt"))
# 
# })
# 
# names(gg_relimp_scenario) <- c("phase", "noise")

#### Save plot ####

# overwrite <- FALSE
# 
# suppoRt::save_ggplot(plot = gg_coef_scenario$noise, filename = paste0("Figure-4", extension),
#                      path = "04_Figures/", width = width * 0.6, height = height * 0.65,
#                      units = units, dpi = dpi, overwrite = overwrite)
# 
# Rel importance
# 
# suppoRt::save_ggplot(plot = gg_coef_scenario$phase, filename = paste0("Figure-A5", extension),
#                      path = "04_Figures/Appendix/", width = width * 0.65, height = height * 0.75,
#                      units = units, dpi = dpi, overwrite = overwrite)
# 
# suppoRt::save_ggplot(plot = gg_relimp_scenario$noise, filename = paste0("Figure-A2", extension),
#                      path = "04_Figures/Appendix/", width = width * 0.65, height = height * 0.65,
#                      units = units, dpi = dpi, overwrite = overwrite)
# 
# suppoRt::save_ggplot(plot = gg_relimp_scenario$phase, filename = paste0("Figure-A7", extension),
#                      path = "04_Figures/Appendix/", width = width * 0.65, height = height * 0.75,
#                      units = units, dpi = dpi, overwrite = overwrite)
