##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create density of CV values

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-cv.R")
source("01_Functions/import-mortality.R")
source("01_Functions/import-abundance.R")

#### Load/wrangle simulated data ####

results_phase_df <- import_cv(path = "02_Data/result-phase.rds")

results_noise_df <- import_cv(path = "02_Data/result-noise.rds")

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))

#### Load abundance data ####

abundance_phase_df <- import_abundance(path = "02_Data/result-phase.rds")

abundance_noise_df <- import_abundance(path = "02_Data/result-noise.rds")

abundance_combined_df <- dplyr::bind_rows(phase = abundance_phase_df, noise = abundance_noise_df,
                                          .id = "scenario") |>
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise"))) |> 
  dplyr::group_by(scenario, row_id, pop_n, nutrient_input) |>
  dplyr::summarise(abundance_max = max(mean), .groups = "drop")

results_combined_df <- dplyr::left_join(x = results_combined_df, y = abundance_combined_df, 
                                        by = c("scenario", "row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(pop_n == 128 & nutrient_input == "low" &
                                             abundance_max > threshold_abundance ~ "no", 
                                           TRUE ~ "yes"))

obs_n_df <- dplyr::filter(results_combined_df, include == "yes") |>
  dplyr::group_by(scenario, pop_n, nutrient_input) |>
  dplyr::summarise(n = dplyr::n() / (4 * 3), .groups = "drop") |>  # divide by four because of measures and three because of parts
  dplyr::arrange(n)

#### CV densities ####

gg_cv <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  purrr::map(c(ag = "ag_production", bg =  "bg_production", ttl = "ttl_production"), function(part_i) {
    
    df_temp <- dplyr::filter(results_combined_df, scenario == scenario_i, part == part_i, 
                             measure %in% c("alpha", "gamma"), include == "yes")
    
    df_extra <- dplyr::filter(results_combined_df, scenario == scenario_i, part == part_i, 
                              measure %in% c("alpha", "gamma"), include == "no")
    
    obs_n_temp <- dplyr::filter(obs_n_df, scenario == scenario_i)
    
    # ggplot
    ggplot(df_temp, aes(x = pop_n, y = value.cv)) + 
      
      # geoms
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha = 0.1) +
      geom_jitter(data = df_extra, aes(x = pop_n, y = value.cv), col = "#ef8a47", alpha = 0.25) +
      geom_text(data = obs_n_temp, aes(x = pop_n, y = 0.0, label = paste0("n=", n)), size = 2.5) +
      
      # facet
      facet_grid(rows = vars(nutrient_input), cols = vars(measure), 
                 labeller = labeller(measure = c("alpha" = "Local scale",
                                                 "gamma" = "Meta-ecosystem scale"), 
                                     nutrient_input = function(x) paste("Nutr. input:", x))) +
      
      # themes
      labs(x = "Population size", y = "Coefficient of variation") +
      theme_classic(base_size = 10) +
      theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.5),
            axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA))
    
  })
})

#### Save plots ####

overwrite <- FALSE

suppoRt::save_ggplot(plot = gg_cv$noise$ag, filename = "Figure-A2.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_cv$noise$bg, filename = "Figure-A3.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)
