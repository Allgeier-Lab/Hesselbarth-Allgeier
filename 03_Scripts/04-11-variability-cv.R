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
source("01_Functions/import-cv.R")
source("01_Functions/import-mortality.R")

#### Load/wrangle simulated data ####

results_phase_df <- import_cv(path = "02_Data/result-phase.rds")

results_noise_df <- import_cv(path = "02_Data/result-noise.rds")

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))

#### Load mortality data ####

mortality_phase_df <- import_mortality(path = "02_Data/result-phase.rds") |> 
  dplyr::filter(pop_n > 0)

mortality_noise_df <- import_mortality(path = "02_Data/result-noise.rds") |> 
  dplyr::filter(pop_n > 0)

mortality_combined_df <- dplyr::bind_rows(phase = mortality_phase_df, noise = mortality_noise_df,
                                          .id = "scenario") |>
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise"))) |>
  dplyr::group_by(scenario, row_id, pop_n, nutrient_input) |>
  dplyr::summarise(died_total = mean(died_total), .groups = "drop")

#### Filter abundances/mortality #### 

results_combined_df <- dplyr::left_join(x = results_combined_df, y = mortality_combined_df, 
                                        by = c("scenario", "row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(pop_n == 128 & nutrient_input == "low" &
                                             died_total > threshold_mort ~ "no", 
                                           TRUE ~ "yes"))

#### Split data #### 

results_combined_list <- dplyr::filter(results_combined_df, measure %in% c("alpha", "gamma"),
                                       part %in% c("ag_production", "bg_production", "ttl_production"), 
                                       pop_n != 0, biotic != 0.0, abiotic != 0.0, include == "yes") |>
  dplyr::group_by(scenario, part, measure, pop_n, nutrient_input) |>
  dplyr::group_split()

#### Fit regression model ####

full_lm_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
  
  message("\r> Progress: scenario=", unique(df_temp$scenario), "; part=", unique(df_temp$part), 
          "; measure=", unique(df_temp$measure), "; pop_n=", unique(df_temp$pop_n), 
          "; nutrient_input=", unique(df_temp$nutrient_input), appendLF = FALSE)
  
  df_temp_stand <- dplyr::select(df_temp, value.cv, abiotic, biotic) # |>
  # dplyr::mutate(value.cv = log(value.cv))
  # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x)))
  # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
  
  lm_temp <- lm(value.cv ~ abiotic * biotic, data = df_temp_stand) 
  
  broom::tidy(lm_temp) |> 
    dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                  measure = unique(df_temp$measure), pop_n = unique(df_temp$pop_n), 
                  nutrient_input = unique(df_temp$nutrient_input), .before = term) |> 
    dplyr::mutate(p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                                   p.value < 0.05 ~ "*", p.value >= 0.05 ~ " "), 
                  r2 = summary(lm_temp)$adj.r.squared)}) |> 
  dplyr::mutate(term = factor(term), p.value.class = factor(p.value.class, levels = c("*", "**", "***", " ")))

#### Model assumptions ####
# 
# assumption_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
# 
#   message("\r> Progress: scenario=", unique(df_temp$scenario), "; part=", unique(df_temp$part),
#           "; measure=", unique(df_temp$measure), "; pop_n=", unique(df_temp$pop_n),
#           "; nutrient_input=", unique(df_temp$nutrient_input), appendLF = FALSE)
# 
#   df_temp_stand <- dplyr::select(df_temp, value.cv, abiotic, biotic) # |>
#   # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) |>
#   # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
# 
#   lm_temp <- lm(value.cv ~ abiotic * biotic, data = df_temp_stand)
# 
#   residuals_temp <- resid(lm_temp)
# 
#   tibble::tibble(scenario = unique(df_temp$scenario), part = unique(df_temp$part),
#                  measure = unique(df_temp$measure), pop_n = unique(df_temp$pop_n),
#                  nutrient_input = unique(df_temp$nutrient_input),
#                  shapiro =  shapiro.test(residuals_temp)[[2]], lillie = nortest::lillie.test(residuals_temp)[[2]])})
# 
# dplyr::filter(assumption_df, shapiro < 0.05 & lillie < 0.05)

#### Relative importance R2 ####

# rel_importance_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
#   
#   message("\r> Progress: scenario=", unique(df_temp$scenario), "; part=", unique(df_temp$part), 
#           "; measure=", unique(df_temp$measure), "; pop_n=", unique(df_temp$pop_n), 
#           "; nutrient_input=", unique(df_temp$nutrient_input), appendLF = FALSE)
#   
#   df_temp_stand <- dplyr::select(df_temp, value.cv, abiotic, biotic) # |>
#   # dplyr::mutate(value.cv = log(value.cv))
#   # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x)))
#   # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
#   
#   lm_temp <- lm(value.cv ~ abiotic * biotic, data = df_temp_stand) 
#   
#   rel_r2 <- relaimpo::boot.relimp(lm_temp, type = "lmg", b = 100, level = 0.95, fixed = FALSE) |>
#     relaimpo::booteval.relimp(bty = "basic")
#   
#   tibble::tibble(
#     scenario = unique(df_temp$scenario), part = unique(df_temp$part), measure = unique(df_temp$measure),
#     pop_n = unique(df_temp$pop_n), nutrient_input = unique(df_temp$nutrient_input),
#     beta = c(rel_r2@namen[-1], "residual"), mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)), 
#     lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA), r2 = summary(lm_temp)$adj.r.squared) |> 
#     dplyr::mutate(dplyr::across(mean:higher, ~ .x * 100))}) |> 
#   dplyr::mutate(beta = factor(beta))

#### Setup ggplot ####

size_base <- 10.0
size_text <- 2.0
size_point <- 2.5

width_doge <- 0.5

color_beta <- c("0" = "grey", "8" = "#0D0887FF", "16" = "#7E03A8FF", "32" = "#CC4678FF", "64" = "#F89441FF", "128" = "#F0F921FF")
color_regression <- c("abiotic" = "#007aa1", "biotic" = "#41b282", "abiotic:biotic" = "#fcb252", "(Intercept)" = "#e76254")

gg_dummy_beta <- data.frame(pop_n = factor(c(0, 8, 16, 32, 64, 128), ordered = TRUE), 
                            mean = c(1, 1, 1, 1, 1, 1)) |>
  ggplot() +
  geom_col(aes(x = pop_n, y = mean, fill = pop_n)) +
  scale_fill_manual(name = "Population size", values = color_beta) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_classic(base_size = size_base) +
  theme(legend.position = "bottom")

gg_dummy_regression <- data.frame(beta = c("abiotic", "biotic", "abiotic:biotic", "(Intercept)"), 
                                  mean = c(1, 1, 1, 1)) |>
  dplyr::mutate(beta = factor(beta, levels = c("abiotic", "biotic", "abiotic:biotic", "(Intercept)"))) |>
  dplyr::filter(beta != "(Intercept)") |> 
  ggplot() +
  geom_point(aes(x = beta, y = mean, color = beta), size = 3.5) +
  scale_color_manual(name = "", values = color_regression,
                     labels = c("abiotic" = "Abiotic subsidies", "biotic" = "Consumer connectivity",
                                "abiotic:biotic" = "Subsidies:Connectivity", "residual" = "Residuals")) +
  guides(colour = guide_legend(nrow = 2)) +
  theme_classic(base_size = size_base) +
  theme(legend.position = "bottom")

#### Create ggplot cv ####

gg_scenario <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  purrr::map(c(ag = "ag_production", bg =  "bg_production", ttl = "ttl_production"), function(part_i) {
    
    regression_temp <- dplyr::filter(full_lm_df, scenario == scenario_i, part == part_i, 
                                     term != "(Intercept)")
      
    beta_temp <- dplyr::filter(results_combined_df, scenario == scenario_i, part == part_i, 
                               measure == "beta", include == "yes") |> 
      dplyr::mutate(type = dplyr::case_when(pop_n == 0 ~ "a", abiotic == 0.0 ~ "b", TRUE ~ "ab"), 
                    type = factor(type, levels = c("a", "ab", "b"), 
                                  labels = c("Abiotic subsidies only", "Abiotic and connectivity", 
                                             "Consumer connectivity only")))
    
    y_upper <- dplyr::group_by(beta_temp, nutrient_input, type, pop_n) |> 
      dplyr::summarise(upper = quantile(value.cv, probs = 0.75), .groups = "drop") |> 
      dplyr::pull(upper) |> 
      max()

    gg_regression <- ggplot(data = regression_temp) +
      
      # zero line
      geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
      
      # lines
      geom_line(aes(x = pop_n, y = estimate, group = term, color = term),
                alpha = 0.25, position = position_dodge(width = width_doge)) +
      
      # points
      geom_point(aes(x = pop_n, y = estimate, color = term, shape = p.value.class),
                 shape = 19, position = position_dodge(width = width_doge), size = size_point) +
      
      # text
      geom_text(aes(x = pop_n, y = estimate, label = p.value.class, group = term), color = "grey",
                size = size_text, position = position_dodge(width = width_doge), vjust = 0.75) +
      
      # facet
      # facet_wrap(. ~ measure + nutrient_input, ncol = 2, nrow = 3) +
      facet_grid(rows = vars(nutrient_input), cols = vars(measure), 
                 labeller = labeller(nutrient_input = function(x) paste("Nutr. input:", x))) +
      
      # set scales
      scale_color_manual(name = "Scale", values = color_regression) +

      # labels and themes
      labs(x = "", y = "Parameter estimate") +
      theme_classic(base_size = size_base) +
      theme(strip.background = element_blank(), strip.text.x  = element_blank(),
            axis.line = element_blank(), panel.border = element_rect(linewidth = 0.25, fill = NA),
            legend.position = "none")
    
    gg_regression <- cowplot::ggdraw(gg_regression, ylim = c(0.0, 1.025)) +
      cowplot::draw_label("Local scale", x = 0.3, y = 1, angle = 0, size =  size_base * 2/3) + 
      cowplot::draw_label("Meta-ecosystem scale", x = 0.755, y = 1, angle = 0, size = size_base * 2/3)
    
    gg_regression <- cowplot::plot_grid(gg_regression, cowplot::get_legend(gg_dummy_regression),
                                        nrow = 2, ncol = 1, rel_heights = c(0.9, 0.1))
    
    # ggplot
    gg_beta <- ggplot(beta_temp, aes(x = type, y = value.cv, fill = pop_n)) + 
      
      # geoms
      geom_boxplot(outlier.shape = NA) +
      geom_hline(yintercept = 1.0, linetype = 2, color = "grey") +
      
      scale_fill_manual(name = "Population size", values = c("0" = "grey", "8" = "#0D0887FF", "16" = "#7E03A8FF", 
                                                             "32" = "#CC4678FF", "64" = "#F89441FF", "128" = "#F0F921FF")) +
      coord_cartesian(ylim = c(1.0, y_upper * 1.2)) +
      
      # facet
      facet_wrap(. ~ nutrient_input, ncol = 1, nrow = 3) +
      
      # themes
      labs(x = "", y = "Portfolio effect") +
      theme_classic(base_size = 10) +
      theme(strip.background = element_blank(), strip.text = element_blank(),
            axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA), 
            legend.position = "none")
    
    gg_beta <- cowplot::ggdraw(gg_beta, ylim = c(0.0, 1.025))
    
    gg_beta <- cowplot::plot_grid(gg_beta, cowplot::get_legend(gg_dummy_beta),
                                  nrow = 2, ncol = 1, rel_heights = c(0.9, 0.10))
    
    cowplot::plot_grid(gg_beta, gg_regression, ncol = 2)
  
  })
})

# gg_scenario$noise$ag

#### Save ggplot #### 

suppoRt::save_ggplot(plot = gg_scenario$noise$ag, filename = "Figure-2.pdf",
                     path = "04_Figures/", width = height * 0.8, height = width * 0.75,
                     units = units, dpi = dpi, overwrite = T)

