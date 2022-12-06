##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Figure of CV values

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

#### Filter data ####

results_combined_df <- dplyr::left_join(x = results_combined_df, y = mortality_combined_df, 
                                        by = c("scenario", "row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(pop_n == 128 & nutrient_input == "low" &
                                             died_total > threshold_mort ~ "no", 
                                           TRUE ~ "yes")) |> 
  dplyr::mutate(type = dplyr::case_when(pop_n == 0 ~ "a", abiotic == 0.0 ~ "b", 
                                        TRUE ~ "ab"), 
                type = factor(type, levels = c("a", "ab", "b"), labels = c("Abiotic subsidies only", 
                                                                           "Abiotic and connectivity", 
                                                                           "Connectivity only")))

obs_n_df <- dplyr::filter(results_combined_df, include == "yes") |>
  dplyr::group_by(scenario, pop_n, nutrient_input) |>
  dplyr::summarise(n = dplyr::n() / (3 * 4), .groups = "drop") |>  # divide by four because of measures and three because of parts
  dplyr::arrange(n, pop_n, scenario, nutrient_input)

#### CV densities ####

gg_cv <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  purrr::map(c(ag = "ag_production", bg =  "bg_production", ttl = "ttl_production"), function(part_i) {
    
    df_temp <- dplyr::filter(results_combined_df, scenario == scenario_i, part == part_i, 
                             measure == "beta", include == "yes")

    # y_max <- quantile(df_temp$value.cv, probs = 0.99)
    
    # obs_n_temp <- dplyr::filter(obs_n_df, scenario == scenario_i)
    
    # ggplot
    ggplot(df_temp, aes(x = type, y = value.cv, fill = pop_n)) + 
      
      # geoms
      geom_boxplot(outlier.shape = NA) +
      geom_hline(yintercept = 1.0, linetype = 2, color = "grey") +

      facet_grid(rows = vars(nutrient_input), scales = "fixed",
                 labeller = labeller(nutrient_input = function(x) paste("Nutr. input:", x))) +
    
      scale_fill_manual(name = "Population size", values = c("0" = "grey", "8" = "#0D0887FF", "16" = "#7E03A8FF", 
                                                             "32" = "#CC4678FF", "64" = "#F89441FF", "128" = "#F0F921FF")) +
      coord_cartesian(ylim = c(1.0, 15.0)) +
      # coord_flip() +
      
      # themes
      labs(x = "", y = "Portfolio effect") +
      theme_classic(base_size = 10) +
      theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.5),
            axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA), 
            legend.position = "bottom")
    
  })
})

gg_cv$noise$ag

#### Save plots ####

overwrite <- FALSE

suppoRt::save_ggplot(plot = gg_cv$noise$ag, filename = "Figure-A2.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)

# suppoRt::save_ggplot(plot = gg_cv$noise$bg, filename = "Figure-A3.pdf",
#                      path = "04_Figures/Appendix/", width = width, height = height * 0.5,
#                      units = units, dpi = dpi, overwrite = overwrite)
