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
  dplyr::summarise(died_total = sum(died_total), .groups = "drop")

#### Filter quantiles mortality #### 

results_combined_df <- dplyr::group_by(mortality_combined_df, scenario, pop_n, nutrient_input) |> 
  dplyr::summarise(thres = quantile(died_total, probs = 0.75), .groups = "drop") |> 
  dplyr::right_join(y = mortality_combined_df, by = c("scenario", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(pop_n == 128 & nutrient_input == "low" & died_total > thres ~ "no", 
                                           TRUE ~ "yes")) |> 
  dplyr::select(scenario, row_id, include) |> 
  dplyr::right_join(y = results_combined_df, by = c("scenario", "row_id"))

# dplyr::filter(results_combined_df, include == "yes") |> 
#   dplyr::group_by(scenario, pop_n, nutrient_input) |> 
#   dplyr::summarise(n = dplyr::n() / (4 * 3)) |>  # divide by four because of measures and three because of parts
#   dplyr::arrange(n)

#### Split data #### 

results_combined_list <- dplyr::filter(results_combined_df, measure %in% c("alpha", "gamma"),
                                       part %in% c("ag_production", "bg_production", "ttl_production"), 
                                       pop_n != 0, include == "yes") |>
  dplyr::group_by(scenario, part, measure, pop_n, nutrient_input) |>
  dplyr::group_split()

#### Fit regression model ####

full_lm_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
  
  message("\r> Progress: scenario=", unique(df_temp$scenario), "; part=", unique(df_temp$part), 
          "; measure=", unique(df_temp$measure), "; pop_n=", unique(df_temp$pop_n), 
          "; nutrient_input=", unique(df_temp$nutrient_input), appendLF = FALSE)
  
  df_temp_stand <- dplyr::select(df_temp, value.cv, value.sd, value.mn, biotic, abiotic) |> 
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) |>
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
  
  purrr::map_dfr(c(cv = "value.cv", sd = "value.sd", mn = "value.mn"), y = df_temp_stand, function(x, y) {
    
    lm_temp <- lm(as.formula(paste0(x, " ~ biotic * abiotic")), data = y) 

    broom::tidy(lm_temp) |> 
      dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                    measure = unique(df_temp$measure), pop_n = unique(df_temp$pop_n), 
                    nutrient_input = unique(df_temp$nutrient_input), .before = term) |> 
      dplyr::mutate(p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                                     p.value < 0.05 ~ "*", p.value >= 0.05 ~ " "), 
                    r2 = summary(lm_temp)$adj.r.squared)
    }, .id = "response")}) |> 
  dplyr::mutate(response = factor(response, levels = c("cv", "mn", "sd")), 
                term = factor(term), p.value.class = factor(p.value.class, levels = c("*", "**", "***", " ")))

#### Relative importance R2 ####

rel_importance_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
  
  message("\r> Progress: scenario=", unique(df_temp$scenario), "; part=", unique(df_temp$part), 
          "; measure=", unique(df_temp$measure), "; pop_n=", unique(df_temp$pop_n), 
          "; nutrient_input=", unique(df_temp$nutrient_input), appendLF = FALSE)
  
  df_temp_stand <- dplyr::select(df_temp, value.cv, value.sd, value.mn, biotic, abiotic) |> 
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) |>
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
  
  purrr::map_dfr(c(cv = "value.cv", sd = "value.sd", mn = "value.mn"), y = df_temp_stand, function(x, y) {
    
    lm_temp <- lm(as.formula(paste0(x, " ~ biotic * abiotic")), data = y) 
    
    rel_r2 <- relaimpo::boot.relimp(lm_temp, type = "lmg", b = 100, level = 0.95, fixed = FALSE) |>
      relaimpo::booteval.relimp(bty = "basic")
    
    tibble::tibble(
      scenario = unique(df_temp$scenario), part = unique(df_temp$part), measure = unique(df_temp$measure),
      pop_n = unique(df_temp$pop_n), nutrient_input = unique(df_temp$nutrient_input),
      beta = c(rel_r2@namen[-1], "residual"), mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)), 
      lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA), r2 = summary(lm_temp)$adj.r.squared) |> 
      dplyr::mutate(dplyr::across(mean:higher, ~ .x * 100))}, .id = "response")}) |> 
  dplyr::mutate(response = factor(response, levels = c("cv", "mn", "sd")), beta = factor(beta))

#### Setup ggplot ####

size_base <- 10.0
size_text <- 2.0
size_point <- 3.5

width_doge <- 0.5

color_scale <- c("biotic" = "#41b282", "abiotic" = "#007aa1", "biotic:abiotic" = "#fcb252", "residual" = "grey")

gg_dummy <- data.frame(beta = c("biotic", "abiotic", "biotic:abiotic", "residual"),
                       mean = c(1, 1, 1, 1)) |>
  dplyr::mutate(beta = factor(beta, levels = c("abiotic", "biotic", "biotic:abiotic", "residual"))) |>
  # dplyr::filter(beta != "biotic:abiotic") |>
  ggplot() +
  geom_col(aes(x = beta, y = mean, fill = beta)) +
  scale_fill_manual(name = "", values = color_scale,
                    labels = c("biotic" = "Consumer behavior", "abiotic" = "Abiotic subsidies",
                               "biotic:abiotic" = "Subsidies:Behavior", "residual" = "Residuals")) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  theme_classic(base_size = size_base) +
  theme(legend.position = "bottom")

#### Create ggplot cv ####

y_range_cv <- dplyr::filter(full_lm_df, response == "cv") |>
  dplyr::pull(estimate) |>
  range(na.rm = TRUE)

gg_scenario_cv <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  gg_part <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {

    gg_cv_list <- purrr::map(c("alpha", "gamma"), function(measure_i) {
  
      purrr::map(levels(results_combined_df$nutrient_input), function(nutrient_i) {
  
        regression_df_temp <- dplyr::filter(full_lm_df, scenario == scenario_i,
                                            measure == measure_i, part == part_i, 
                                            nutrient_input == nutrient_i, response == "cv",
                                            term %in% c("biotic", "abiotic", "biotic:abiotic"))
  
        importance_df_temp <- dplyr::filter(rel_importance_df, scenario == scenario_i,
                                            measure == measure_i, part == part_i, 
                                            nutrient_input == nutrient_i, response == "cv",
                                            beta %in% c("biotic", "abiotic", "biotic:abiotic", "residual"))
        
        label_part <- paste0("Nutr. input: ", nutrient_i)
        
        gg_regression <- ggplot(data = regression_df_temp) +
  
          # zero line
          geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
  
          # Lines
          geom_line(aes(x = pop_n, y = estimate, group = term, color = term),
                    alpha = 0.25, position = position_dodge(width = width_doge)) +
  
          # Points
          geom_point(aes(x = pop_n, y = estimate, color = term, size = term),
                     shape = 19, position = position_dodge(width = width_doge)) +
  
          # Text
          geom_text(aes(x = pop_n, y = estimate, label = p.value.class, group = term), color = "white",
                    size = size_text, position = position_dodge(width = width_doge), vjust = 0.75) +
  
          # set scales
          scale_color_manual(name = "Scale", values = color_scale) +
          scale_y_continuous(limits = y_range_cv, breaks = seq(y_range_cv[[1]], y_range_cv[[2]], length.out = 4),
                             labels = function(x) round(x, digits = 2)) +
          scale_size_manual(values = c("biotic" = size_point, "abiotic" = size_point, "biotic:abiotic" = size_point / 2)) +
          # scale_shape_manual(values = c("***" = 19, "**" = 19, "*" = 19, " " = 1)) +
          coord_cartesian(clip = "off") +
  
          # labels and themes
          labs(x = "", y = "") +
          theme_classic(base_size = size_base) +
          theme(strip.background = element_blank(), strip.text = element_blank(), plot.title = element_text(size = size_base),
                axis.line = element_blank(), panel.border = element_rect(linewidth = 0.25, fill = NA),
                legend.position = "none",
                plot.margin = margin(t = 5.5, r = 0.5, b = 5.5, l = 5.5, unit = "pt"))
  
        gg_relimp <- ggplot(data = importance_df_temp) +
  
          # relative importance bars
          geom_col(aes(x = pop_n, y = mean, fill = factor(beta, levels = rev(levels(beta))))) +
  
          # set scales
          scale_fill_manual(name = "", values = color_scale) +
          scale_y_continuous(labels = function(x) paste0(x, "%")) +
  
          # labels and themes
          labs(x = "", y = "") +
          theme_classic(base_size = size_base) +
          theme(strip.background = element_blank(), strip.text = element_blank(),
                axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA),
                legend.position = "none",
                plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.5, unit = "pt"))
  
        cowplot::plot_grid(gg_regression, gg_relimp, ncol = 2, rel_widths = c(0.5, 0.5)) |>
          cowplot::ggdraw(xlim = c(0, 1.05)) +
          cowplot::draw_label(label = label_part, x = 1.0, y = 0.65, angle = 270, size = size_base * 0.85)
  
      })
    }) |> purrr::flatten()
  
    gg_cv_combined <- cowplot::plot_grid(plotlist = gg_cv_list, nrow = 3, 
                                         ncol = 2, byrow = FALSE)
  
    gg_cv_combined <- cowplot::ggdraw(gg_cv_combined, xlim = c(-0.05, 1.0), ylim = c(-0.05, 1.05)) +
      cowplot::draw_label("Population size", x = 0.5, y = 0, angle = 0, size = size_base) +
      cowplot::draw_label("Parameter estimate / Relative importance [%]", x = 0.0, y = 0.5,
                          angle = 90, size = size_base) +
      cowplot::draw_label("Local", x = 0.25, y = 1.0, vjust = -0.35, angle = 0, size = size_base * 1.0) +
      cowplot::draw_label("Meta-ecosystem", x = 0.75, y = 1.0, vjust = -0.35, size = size_base * 1.0)
  
    cowplot::plot_grid(gg_cv_combined, cowplot::get_legend(gg_dummy),
                       nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))
  
  })
  
  names(gg_part) <- c("ag_production", "bg_production", "ttl_production")
  gg_part
  
})

#### Save ggplot #### 

suppoRt::save_ggplot(plot = gg_scenario_cv$noise$ag_production, filename = "Figure-2.pdf",
                     path = "04_Figures/", width = height * 0.8, height = width * 0.75,
                     units = units, dpi = dpi, overwrite = T)
