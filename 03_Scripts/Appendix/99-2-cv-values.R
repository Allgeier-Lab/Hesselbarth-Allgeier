##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: 

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-cv.R")
source("01_Functions/import-mortality.R")

#### Load/wrangle simulated data ####

near <- FALSE

results_phase_df <- import_cv(path = "02_Data/result-phase.rds", near = near)

results_noise_df <- import_cv(path = "02_Data/result-noise.rds", near = near)

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario")
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

results_final_df <- dplyr::left_join(x = results_combined_df, y = mortality_combined_df, 
                                     by = c("scenario", "row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(died_total > threshold_mort ~ "no", TRUE ~ "yes")) |> 
  dplyr::mutate(treatment = dplyr::case_when(abiotic != 0.0 & biotic == 0.0 ~ "subsidies", 
                                             abiotic == 0.0 & biotic != 0.0 ~ "connectivity", 
                                             abiotic != 0.0 & biotic != 0.0 ~ "combined")) |>
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")),
                part = factor(part, labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                               "ttl_productio" = "Total")),
                treatment = factor(treatment, levels = c("subsidies", "connectivity", "combined")))

#### Setup ggplot ####

base_size <- 10.0

w <- 0.5

color_means <- c(low = "#007895", medium = "#f2d445", high = "#e93c26")

gg_dummy <- data.frame(x = c(1, 2, 3), y = c(0.5, 0.5, 0.5), z = c("low", "medium", "high")) |>
  dplyr::mutate(z = factor(z, levels = c("low", "medium", "high"))) |>
  ggplot(aes(x = x, y = y, color = z, fill = z)) +
  geom_boxplot(alpha = 0.25) +
  scale_color_manual(name = "Enrichment treatment", values = color_means) +
  scale_fill_manual(name = "Enrichment treatment", values = color_means) +
  theme_classic(base_size = base_size) +
  theme(legend.position = "bottom")

#### Create ggplot: Main ####

gg_cv <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  ranges_temp <- dplyr::filter(results_final_df, scenario == scenario_i, measure %in% c("alpha", "gamma"), 
                               include == "yes", treatment == "combined", part != "Total") |> 
    dplyr::group_by(part) |> 
    dplyr::summarise(y_low = min(value.cv), y_max = max(value.cv))
  
  gg_scale <- purrr::map(c(alpha = "alpha", gamma = "gamma"), function(measure_i) {
    
    if (measure_i == "alpha") y_lab <- expression(italic(Coefficient~of~variation[alpha])) else 
      y_lab <- expression(italic(Coefficient~of~variation[gamma]))
    
    if (measure_i == "alpha") labels_temp <- c("a)", "c)") else labels_temp <- c("b)", "d)")
    
    gg_part <- purrr::map(c(Aboveground = "Aboveground", Belowground = "Belowground"), function(part_i) {
      
      if (part_i == "Aboveground") x_axis <- element_blank() else x_axis <- NULL
      
      cv_temp <- dplyr::filter(results_final_df, scenario == scenario_i, part == part_i, measure == measure_i,
                               include == "yes", treatment == "combined")
      
      ggplot(data = dplyr::filter(cv_temp, part == part_i), aes(x = pop_n, y = value.cv)) + 
        
        # adding points
        geom_boxplot(aes(color = nutrient_input, fill = nutrient_input), alpha = 0.25) +
        
        # facet_wrap(. ~ measure, scales = "fixed") +
        
        # change scales
        scale_color_manual(name = "", values = color_means) +
        scale_fill_manual(name = "", values = color_means) +
        scale_y_continuous(limits = c(dplyr::filter(ranges_temp, part == part_i)$y_low, 
                                      dplyr::filter(ranges_temp, part == part_i)$y_max)) +
        
        # adding box
        annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
        annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
        annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
        annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
        
        # change theme
        labs(y = y_lab) +
        theme_classic(base_size = base_size) + 
        theme(legend.position = "none", axis.line = element_blank(), 
              axis.title = element_blank(), axis.text.x = x_axis)
      
    })
    
    gg_part <- cowplot::plot_grid(plotlist = gg_part, nrow = 2, ncol = 1, 
                                  labels = labels_temp)
    
    cowplot::ggdraw(gg_part, xlim = c(-0.025, 1.0)) +
      cowplot::draw_label(y_lab, x = -0.01, y = 0.5, angle = 90, size = base_size)
    
  })
  
  gg_scale <- cowplot::plot_grid(plotlist = gg_scale, nrow = 1, ncol = 2)
  
  gg_total <- cowplot::ggdraw(gg_scale, ylim = c(-0.025, 1.0)) +
    cowplot::draw_label("Population size", y = -0.01, x = 0.5, angle = 0, size = base_size)
  
  cowplot::plot_grid(gg_total, cowplot::get_legend(gg_dummy), nrow = 2, ncol = 1, 
                     rel_heights = c(0.9, 0.1))
  
})

#### Save ggplot #### 

suppoRt::save_ggplot(plot = gg_cv$noise, filename = "Figure-A2.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.45,
                     units = units, dpi = dpi, overwrite = FALSE)
