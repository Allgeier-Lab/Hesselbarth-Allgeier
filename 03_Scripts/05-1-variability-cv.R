##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create total result figure of variability and cv

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

results_combined_list <- dplyr::filter(results_combined_df, measure %in% c("alpha", "gamma"),
                                       part %in% c("ag_production", "bg_production", 
                                                   "ttl_production")) |>
  dplyr::select(-c(pop_n, value.prod, value.move)) |> 
  dplyr::group_by(scenario, part, measure) |>
  dplyr::group_split()

#### Fit regression model ####

full_lm_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
  
  df_temp_stand <- dplyr::select(df_temp, value.cv, value.cv.sd, value.cv.mn, 
                                 biotic, abiotic, value.biomass, nutrient_input) |> 
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) |>
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
  
  purrr::map_dfr(c(cv = "value.cv", sd = "value.cv.sd", mn = "value.cv.mn"), y = df_temp_stand, function(x, y) {
    
    lm_temp <- lm(as.formula(paste0(x, " ~ biotic * abiotic * value.biomass * nutrient_input")), data = y) 

    broom::tidy(lm_temp) |> 
      dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                    measure = unique(df_temp$measure), .before = term) |> 
      dplyr::mutate(p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                                     p.value < 0.05 ~ "*", p.value >= 0.05 ~ ""), 
                    r2 = summary(lm_temp)$adj.r.squared)
    }, .id = "response")}) |> 
  dplyr::mutate(response = factor(response, levels = c("cv", "mn", "sd")), 
                term = factor(term), p.value.class = factor(p.value.class, levels = c("*", "**", "***", "")))

#### Save linear regression results ####

purrr::walk(c("phase", "noise"), function(i) {
  dplyr::filter(full_lm_df, scenario == i, response == "cv") |> 
    suppoRt::save_rds(filename = paste0("lm_variability_cv_", i, ".rds"), 
                      path = "05_Results/", overwrite = FALSE)
})

purrr::walk(c("phase", "noise"), function(i) {
  dplyr::filter(full_lm_df, scenario == i, response != "cv") |> 
    suppoRt::save_rds(filename = paste0("lm_variability_cv_seperated-", i, ".rds"), 
                      path = "05_Results/", overwrite = FALSE)
})

#### Relative importance R2 ####

rel_importance_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
  
  message("\r> Progress: scenario=", unique(df_temp$scenario), "; part=", unique(df_temp$part), "; measure=",
          unique(df_temp$measure), "\t\t", appendLF = FALSE)
  
  df_temp_stand <- dplyr::select(df_temp, value.cv, value.cv.sd, value.cv.mn, 
                                 biotic, abiotic, value.biomass, nutrient_input) |> 
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) |>
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
  
  purrr::map_dfr(c(cv = "value.cv", sd = "value.cv.sd", mn = "value.cv.mn"), y = df_temp_stand, function(x, y) {
    
    lm_temp <- lm(as.formula(paste0(x, " ~ biotic + abiotic + value.biomass + nutrient_input")), 
                  data = y) 
    
    rel_r2 <- relaimpo::boot.relimp(lm_temp, type = "lmg", b = 500, level = 0.95, fixed = FALSE) |>
      relaimpo::booteval.relimp(bty = "basic")
    
    tibble::tibble(
      scenario = unique(df_temp$scenario), part = unique(df_temp$part), measure = unique(df_temp$measure),
      beta = c(rel_r2@namen[-1], "residual"), mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)), 
      lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA), r2 = summary(lm_temp)$adj.r.squared) |> 
      dplyr::mutate(dplyr::across(mean:higher, ~ .x * 100))
    
    }, .id = "response")}) |> 
  dplyr::mutate(response = factor(response, levels = c("cv", "mn", "sd")), beta = factor(beta))

#### Save relative importance results ####

purrr::walk(c("phase", "noise"), function(i) {
  dplyr::filter(rel_importance_df, scenario == i, response == "cv") |> 
    suppoRt::save_rds(filename = paste0("relimp_variability_cv_", i, ".rds"), 
                      path = "05_Results/", overwrite = FALSE)
})

purrr::walk(c("phase", "noise"), function(i) {
  dplyr::filter(rel_importance_df, scenario == i, response != "cv") |> 
    suppoRt::save_rds(filename = paste0("relimp_variability_cv_seperated-", i, ".rds"), 
                      path = "05_Results/", overwrite = FALSE)
})

#### Setup ggplot ####

# size_base <- 10.0
# size_text <- 2.0
# size_point <- 3.5
# 
# width_doge <- 0.5
# 
# color_scale <- c("biotic" = "#41b282", "abiotic" = "#007aa1", "biotic:abiotic" = "#fcb252", "residual" = "grey")
# 
# gg_dummy <- data.frame(beta = c("biotic", "abiotic", "biotic:abiotic", "residual"),
#                        mean = c(1, 1, 1, 1)) |> 
#   dplyr::mutate(beta = factor(beta, levels = c("biotic", "abiotic", "biotic:abiotic", "residual"))) |> 
#   dplyr::filter(beta != "biotic:abiotic") |> 
#   ggplot() + 
#   geom_col(aes(x = beta, y = mean, fill = beta)) + 
#   scale_fill_manual(name = "", values = color_scale[-3], 
#                     labels = c("biotic" = "Consumer behavior", "abiotic" = "Abiotic subsidies", 
#                                "biotic:abiotic" = "Interaction Behavior:Subsidies", "residual" = "Residuals")) +
#   guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
#   theme_classic(base_size = size_base) + 
#   theme(legend.position = "bottom")

#### Create ggplot cv ####

# y_range_cv <- dplyr::filter(regression_df, response == "cv") |> 
#   dplyr::pull(estimate) |> 
#   range()
# 
# gg_scenario_cv <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
#   
#   gg_cv_list <- purrr::map(c("alpha", "gamma"), function(measure_i) {
#     
#     purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
#       
#       regression_df_temp <- dplyr::filter(regression_df, scenario == scenario_i,
#                                           measure == measure_i, part == part_i, response == "cv", 
#                                           term %in% c("biotic", "abiotic"))
#       
#       importance_df_temp <- dplyr::filter(importance_df,  scenario == scenario_i,
#                                           measure == measure_i, part == part_i, response == "cv", 
#                                           beta %in% c("biotic", "abiotic", "residual"))
#       
#       
#       if (measure_i == "gamma") { 
#         if (part_i == "ag_production") {
#           label_part <-  "Aboveground"
#         } else if (part_i == "bg_production") {
#           label_part <-  "Belowground"
#         }  else {
#           label_part <- "Total      "
#         }
#       } else {
#         label_part <- ""
#       }
# 
#       gg_regression <- ggplot(data = regression_df_temp) + 
#         
#         # zero line
#         geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
#         
#         # Lines
#         geom_line(aes(x = pop_n, y = estimate, group = term, color = term),
#                   alpha = 0.5, position = position_dodge(width = width_doge)) +
#         
#         # Points
#         geom_point(aes(x = pop_n, y = estimate, color = term),
#                    size = size_point, shape = 19, position = position_dodge(width = width_doge)) +
#         
#         # Text
#         geom_text(aes(x = pop_n, y = estimate, label = p.value, group = term), color = "white",
#                   size = size_text, position = position_dodge(width = width_doge), vjust = 0.75) +
#         
#         # set scales
#         scale_color_manual(name = "Scale", values = color_scale) +
#         scale_y_continuous(limits = y_range_cv, breaks = seq(y_range_cv[[1]], y_range_cv[[2]], length.out = 4), 
#                            labels = function(x) round(x, digits = 2)) + 
#         coord_cartesian(clip = "off") +
#         
#         # labels and themes
#         labs(x = "", y = "") +
#         theme_classic(base_size = size_base) + 
#         theme(strip.background = element_blank(), strip.text = element_blank(), 
#               axis.line = element_blank(), panel.border = element_rect(size = 0.25, fill = NA),
#               legend.position = "none",
#               plot.margin = margin(t = 5.5, r = 0.5, b = 5.5, l = 5.5, unit = "pt"))
#       
#       gg_relimp <- ggplot(data = importance_df_temp) + 
#         
#         # relative importance bars
#         geom_col(aes(x = pop_n, y = mean * 100, fill = beta)) + 
#         
#         # set scales
#         scale_fill_manual(name = "", values = color_scale) +
#         scale_y_continuous(labels = function(x) paste0(x, "%")) + 
#         
#         # labels and themes
#         labs(x = "", y = "") +
#         theme_classic(base_size = size_base) + 
#         theme(strip.background = element_blank(), strip.text = element_blank(), 
#               axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA),
#               legend.position = "none",
#               plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.5, unit = "pt"))
#       
#       cowplot::plot_grid(gg_regression, gg_relimp, ncol = 2, rel_widths = c(0.5, 0.5)) |> 
#         cowplot::ggdraw(xlim = c(0, 1.05)) +
#         cowplot::draw_label(label = label_part, x = 1.0, y = 0.75, vjust = -0.5, angle = 270, 
#                             size = size_base * 0.85)
#       
#     })
#   }) |> purrr::flatten()
#   
#   gg_cv_combined <- cowplot::plot_grid(plotlist = gg_cv_list, nrow = 3, ncol = 2, byrow = FALSE)
#   
#   gg_cv_combined <- cowplot::ggdraw(gg_cv_combined, xlim = c(-0.05, 1.0), ylim = c(-0.05, 1.05)) + 
#     cowplot::draw_label("Population size", x = 0.5, y = 0, angle = 0, size = size_base) + 
#     cowplot::draw_label("Parameter estimate / Relative importance [%]", x = 0.0, y = 0.5,
#                         angle = 90, size = size_base) +
#     cowplot::draw_label("Local", x = 0.25, y = 1.0, vjust = -0.35, angle = 0, size = size_base * 1.0) + 
#     cowplot::draw_label("Meta-ecosystem", x = 0.75, y = 1.0, vjust = -0.35, size = size_base * 1.0)
#   
#   cowplot::plot_grid(gg_cv_combined, cowplot::get_legend(gg_dummy),
#                      nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))
#   
# })

#### Create ggplot sd and mean ####

# y_range_frac <- dplyr::filter(regression_df, response != "cv") |> 
#   dplyr::pull(estimate) |> 
#   range()
# 
# gg_scenario_frac <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
#   
#   gg_frac_list <- purrr::map(c("alpha", "gamma"), function(measure_i) {
#     
#     purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
#       
#       regression_df_temp <- dplyr::filter(regression_df, scenario == scenario_i,
#                                           measure == measure_i, part == part_i, response != "cv", 
#                                           term %in% c("biotic", "abiotic"))
#     
#       importance_df_temp <- dplyr::filter(importance_df, scenario == scenario_i,
#                                           measure == measure_i, part == part_i, response != "cv", 
#                                           beta %in% c("biotic", "abiotic", "residual"))
#       
#       if (measure_i == "gamma") { 
#         if (part_i == "ag_production") {
#           label_part <-  "Aboveground"
#         } else if (part_i == "bg_production") {
#           label_part <-  "Belowground"
#         }  else {
#           label_part <- "Total      "
#         }
#       } else {
#         label_part <- ""
#       }
#      
#       gg_regression <- ggplot(data = regression_df_temp) + 
#         
#         # zero line
#         geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
#         
#         # Lines
#         geom_line(aes(x = pop_n, y = estimate, group = term, color = term),
#                   alpha = 0.5, position = position_dodge(width = width_doge)) +
#         
#         # Points
#         geom_point(aes(x = pop_n, y = estimate, color = term),
#                    size = size_point * 0.75, shape = 19, position = position_dodge(width = width_doge)) +
#         
#         # Text
#         geom_text(aes(x = pop_n, y = estimate, label = p.value, group = term), color = "white",
#                   size = size_text, position = position_dodge(width = width_doge), vjust = 0.75) +
#         
#         # wrap lm models
#         facet_wrap(. ~ response, nrow = 2, labeller = labeller(response = c("mn" = "Mean", "sd" = "Standard deviation"))) +
#         
#         # set scales
#         scale_color_manual(name = "Scale", values = color_scale) +
#         scale_y_continuous(limits = y_range_frac, breaks = seq(y_range_frac[[1]], y_range_frac[[2]], length.out = 4), 
#                            labels = function(x) round(x, digits = 2)) + 
#         coord_cartesian(clip = "off") +
#         
#         # labels and themes
#         labs(x = "", y = "") +
#         theme_classic(base_size = size_base) + 
#         theme(strip.background = element_blank(), strip.text = element_text(hjust = 0), 
#               axis.line = element_blank(), panel.border = element_rect(size = 0.25, fill = NA),
#               legend.position = "none",
#               plot.margin = margin(t = 5.5, r = 0.5, b = 5.5, l = 5.5, unit = "pt"))
#       
#       gg_relimp <- ggplot(data = importance_df_temp) + 
#         
#         # relative importance bars
#         geom_col(aes(x = pop_n, y = mean * 100, fill = beta)) + 
#         
#         # set scales
#         scale_fill_manual(name = "", values = color_scale) +
#         scale_y_continuous(labels = function(x) paste0(x, "%")) + 
#         
#         # wrap lm models
#         facet_wrap(. ~ response, nrow = 2, labeller = labeller(response = c("mn" = "", "sd" = ""))) + 
#         
#         # labels and themes
#         labs(x = "", y = "") +
#         theme_classic(base_size = size_base) + 
#         theme(strip.background = element_blank(), strip.text = element_text(hjust = 0), 
#               axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA),
#               legend.position = "none",
#               plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.5, unit = "pt"))
#       
#       cowplot::plot_grid(gg_regression, gg_relimp, ncol = 2, rel_widths = c(0.5, 0.5)) |> 
#         cowplot::ggdraw(xlim = c(0, 1.05)) +
#         cowplot::draw_label(label = label_part, x = 1.0, y = 0.75, vjust = -0.5, angle = 270, 
#                             size = size_base * 0.85)
#       
#     })
#   }) |> purrr::flatten()
#   
#   gg_frac_combined <- cowplot::plot_grid(plotlist = gg_frac_list, nrow = 3, ncol = 2, byrow = FALSE)
#   
#   gg_frac_combined <- cowplot::ggdraw(gg_frac_combined, xlim = c(-0.05, 1.0), ylim = c(-0.05, 1.05)) + 
#     cowplot::draw_label("Population size", x = 0.5, y = 0, angle = 0, size = size_base) + 
#     cowplot::draw_label("Parameter estimate / Relative importance [%]", x = 0.0, y = 0.5,
#                         angle = 90, size = size_base) +
#     cowplot::draw_label("Local", x = 0.25, y = 1.0, vjust = -0.35, angle = 0, size = size_base * 1.0) + 
#     cowplot::draw_label("Meta-ecosystem", x = 0.75, y = 1.0, vjust = -0.35, size = size_base * 1.0)
#   
#   gg_frac_combined <- cowplot::plot_grid(gg_frac_combined, cowplot::get_legend(gg_dummy),
#                                          nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))
#   
# })

#### Save plot ####

# overwrite <- FALSE
# 
# # CV
# 
# suppoRt::save_ggplot(plot = gg_scenario_cv$noise, filename = paste0("Figure-3", extension),
#                      path = "04_Figures/", width = height, height = width * 0.75,
#                      units = units, dpi = dpi, overwrite = overwrite)
# 
# suppoRt::save_ggplot(plot = gg_scenario_cv$phase, filename = paste0("Figure-A2", extension),
#                      path = "04_Figures/Appendix", width = height, height = width * 0.7,
#                      units = units, dpi = dpi, overwrite = overwrite)

# separated
# 
# suppoRt::save_ggplot(plot = gg_scenario_frac$noise, filename = paste0("Figure-A1", extension),
#                      path = "04_Figures/Appendix/", width = width, height = height * 0.85,
#                      units = units, dpi = dpi, overwrite = overwrite)
# 
# 
# suppoRt::save_ggplot(plot = gg_scenario_frac$phase, filename = paste0("Figure-A4", extension),
#                      path = "04_Figures/Appendix/", width = width, height = height * 0.85,
#                      units = units, dpi = dpi, overwrite = overwrite)
