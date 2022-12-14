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
source("01_Functions/import-cv.R")
source("01_Functions/import-mortality.R")

#### Load/wrangle simulated data ####

results_phase_df <- import_cv(path = "02_Data/result-phase.rds", near = TRUE)

results_noise_df <- import_cv(path = "02_Data/result-noise.rds", near = TRUE)

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
                                           TRUE ~ "yes")) |> 
  dplyr::filter(abiotic != 0.0, measure == "gamma")

obs_n_df <- dplyr::filter(results_combined_df, include == "yes") |>
  dplyr::group_by(scenario, pop_n, nutrient_input) |>
  dplyr::summarise(n = dplyr::n() / (3 * 4), .groups = "drop") |>  # divide by four because of measures and three because of parts
  dplyr::arrange(n, pop_n, scenario, nutrient_input)

#### Density PP ####

size_base <- 10.0

gg_pp <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  gg_part <- purrr::map(c(ag = "ag_production", bg =  "bg_production"), function(part_i) {
    
    df_temp <- dplyr::filter(results_combined_df, scenario == scenario_i, 
                             part %in% part_i, include == "yes")
    
    obs_n_temp <- dplyr::filter(obs_n_df, scenario == scenario_i)
    
    # ggplot
    ggplot(df_temp, aes(x = pop_n, y = value.prod)) + 
      
      # geoms
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha = 0.1) +
      # geom_text(data = obs_n_temp, aes(x = pop_n, y = 0.5, label = paste0("n=", n)), size = 2.5) +
      
      # scales
      scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4), 
                         labels = function(x) round(x, digits = 2)) +
      
      # facet
      facet_wrap(. ~ nutrient_input, scales = "fixed", nrow = 3, ncol = 1) +
      
      # themes
      labs(x = "", y = "") +
      theme_classic(base_size = 10) +
      theme(strip.background = element_blank(), strip.text = element_blank(),
            axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA))
    
  })
  
  gg_combined <- cowplot::plot_grid(plotlist = gg_part, ncol = 2)
  
  cowplot::ggdraw(gg_combined, xlim = c(-0.025, 1.025), ylim = c(-0.025,  1.025)) +
    cowplot::draw_label("Population size", x = 0.5, y = 0, angle = 0, size = size_base) +
    cowplot::draw_label(expression(paste("Primary production [", gDW~d^-1~m^-2, "]")), x = 0.0, y = 0.5, angle = 90, size = size_base) +
    cowplot::draw_label("Aboveground", x = 0.265, y = 1.0, angle = 0, size = size_base * 0.75) +
    cowplot::draw_label("Belowground", x = 0.755, y = 1.0, angle = 0, size = size_base * 0.75) +
    cowplot::draw_label("Nutr. input: low", x = 1.0, y = 0.85, angle = 270, size = size_base * 0.65) +
    cowplot::draw_label("Nutr. input: medium", x = 1.0, y = 0.525, angle = 270, size = size_base * 0.65) +
    cowplot::draw_label("Nutr. input: high", x = 1.0, y = 0.225, angle = 270, size = size_base * 0.65)
  
})

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_pp$noise, filename = "Figure-A9.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = FALSE)
