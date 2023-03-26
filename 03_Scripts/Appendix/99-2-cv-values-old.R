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

near <- FALSE

results_phase_df <- import_cv(path = "02_Data/result-phase.rds", near = near)

results_noise_df <- import_cv(path = "02_Data/result-noise.rds", near = near)

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))

#### Load mortality data ####

mortality_phase_df <- import_mortality(path = "02_Data/result-phase.rds")

mortality_noise_df <- import_mortality(path = "02_Data/result-noise.rds")

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
  dplyr::filter(abiotic != 0.0)

obs_n_df <- dplyr::filter(results_combined_df, include == "yes") |>
  dplyr::group_by(scenario, pop_n, nutrient_input) |>
  dplyr::summarise(n = dplyr::n() / (3 * 4), .groups = "drop") |>  # divide by four because of measures and three because of parts
  dplyr::arrange(n, pop_n, scenario, nutrient_input)

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
      
      # zero line
      geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
      
      # boxplot
      geom_boxplot(outlier.shape = NA) +
      
      # jitter points included
      geom_jitter(alpha = 0.1) +
      
      # jitter points excluded
      geom_jitter(data = df_extra, aes(x = pop_n, y = value.cv), col = "#ef8a47", alpha = 0.25) +
      
      # number of obs
      geom_text(data = obs_n_temp, aes(x = pop_n, y = 0.0, label = paste0("n=", n)), size = 2.5) +
      
      # facet
      facet_grid(rows = vars(nutrient_input), cols = vars(measure), scales = "free_y",
                 labeller = labeller(measure = c("alpha" = "Local scale",
                                                 "gamma" = "Meta-ecosystem scale",
                                                 "beta" = "Portfolio effect"),
                                     nutrient_input = function(x) paste("Nutr. input:", x))) +
      
      # themes
      labs(x = "Population size", y = "Coefficient of variation") +
      theme_classic(base_size = 12) +
      theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.5),
            axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA))
    
  })
})

#### Save plots ####

overwrite <- FALSE

suppoRt::save_ggplot(plot = gg_cv$noise$ag, 
                     filename = ifelse(near, yes = "Figure-A2-near.pdf", no = "Figure-A2.pdf"),
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_cv$noise$bg, 
                     filename = ifelse(near, yes = "Figure-A3-near.pdf", no = "Figure-A3.pdf"),
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)
