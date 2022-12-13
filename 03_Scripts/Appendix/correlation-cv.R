##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Correlation between alpha amd gamma scale vc

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-cv.R")

#### Load/wrangle simulated data ####

n <- 5

results_phase_df <- import_cv(path = "02_Data/result-phase.rds", near = FALSE)

results_noise_df <- import_cv(path = "02_Data/result-noise.rds", near = FALSE)

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") %>% 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))

#### Reshape data ####

result_wide_df <- dplyr::select(results_combined_df, scenario, part, measure, pop_n, value.cv) |> 
  dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), measure %in% c("alpha", "gamma")) |> 
  dplyr::mutate(row_id = rep(x = seq(from = 1, to = dplyr::n() / 2), each = 2)) |> 
  tidyr::pivot_wider(names_from = measure, values_from = value.cv) |> 
  dplyr::select(-row_id)

corr_df <- dplyr::group_by(result_wide_df, scenario, part, pop_n) |> 
  dplyr::group_split() |> 
  purrr::map_dfr(function(df_temp) {
    
    cor_test <- cor.test(df_temp$alpha, df_temp$gamma, alternative = "two.sided", 
                         method = "kendall")
    
    tibble::tibble(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                   pop_n = unique(df_temp$pop_n),
                   cor = cor_test$estimate, p_value = cor_test$p.value)}) |> 
  dplyr::mutate(p_value = dplyr::case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", 
                                           p_value < 0.05 ~ "*", TRUE ~ "n.s."))

#### Create ggplot ####

colors_pop <- c("8" = "#92c6de", "16" = "#9cde81", "32" = "#c9a6cf", 
                "64" = "#ffba68", "128" = "#ff908f")

gg_cv_corr <- purrr::map(c("phase", "noise"), function(i) {
  
  data_temp <- dplyr::filter(result_wide_df, scenario == i)
  
  ggplot(data = data_temp, aes(x = alpha, y = gamma)) + 
    
    # adding geoms
    geom_point(shape = 1, alpha = 0.75) + 
    geom_smooth(method = "lm", formula = "y ~ x", se = FALSE, color = "#32b2da", size = 0.5) +
    # ggpubr::stat_regline_equation(aes(label = ..rr.label..), color = "#32b2da", 
    #                               geom = "label", size = 2.5) +
    
    ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = ..adj.rr.label..),  color = "#32b2da",
                          parse = TRUE, rr.digits = 2, size = 2.5, geom = "label_npc") +
    
    # facet wraping
    facet_wrap(. ~ part + pop_n, nrow = 3, ncol = 5, scales = "free",
               labeller = labeller(part = c("ag_production" = "Aboveground", "bg_production" = "Belowground",
                                            "ttl_production" = "Total"),
                                   pop_n = function(x) paste("Pop. size", x))) +
    
    # scales
    labs(x = "Local scale", y = "Meta-ecosystem scale") +
    scale_x_continuous(limits = function(x) range(x), 
                       breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4), 
                       labels = function(x) round(x, 3)) +
    scale_y_continuous(limits = function(x) range(x), 
                       breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4), 
                       labels = function(x) round(x, 3)) +
    
    # themes
    theme_classic(base_size = 10) + 
    theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.0), 
          legend.position = "bottom")
  
})

names(gg_cv_corr) <- c("phase", "noise")

#### Save ggplot ####

overwrite <- FALSE

# suppoRt::save_ggplot(plot = gg_cv_corr$phase, filename = paste0("Figure-A6", extension),
#                      path = "04_Figures/Appendix/", width = height, height = width,
#                      units = units, dpi = dpi, overwrite = overwrite)
# 
# suppoRt::save_ggplot(plot = gg_cv_corr$noise, filename = paste0("Figure-A7", extension),
#                      path = "04_Figures/Appendix/", width = height, height = width,
#                      units = units, dpi = dpi, overwrite = overwrite)
