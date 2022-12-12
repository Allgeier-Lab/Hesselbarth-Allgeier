##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Model runs with with variation (abiotic or biotic) only

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
                                             died_total> threshold_mort ~ "no", 
                                           TRUE ~ "yes")) |> 
  dplyr::filter(include == "yes")

#### Abiotic variation only ####

gg_abiotic <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  purrr::map(c(ag_production = "ag_production", bg_production = "bg_production", ttl_production = "ttl_production"), function(part_i) {
    
      df_temp <- dplyr::filter(results_combined_df, scenario == scenario_i, 
                               part == part_i, pop_n == 0, include == "yes") |> 
        dplyr::select(row_id, pop_n, nutrient_input, measure, value.cv, abiotic, biotic)
      
      alpha_gamma_temp <- dplyr::filter(df_temp, measure %in% c("alpha", "gamma")) |> 
        tidyr::pivot_wider(names_from = measure, values_from = value.cv)
        
      beta_temp <- dplyr::filter(df_temp, measure == "beta")
      
      gg_alpha_gamma <- ggplot(data = alpha_gamma_temp, aes(x = alpha, y = gamma)) + 
        
        geom_point(pch = 1, alpha = 0.5) + 
        geom_abline(color = "grey", linetype = 2) +
        geom_smooth(method = 'loess', formula = "y ~ x", se = FALSE, color = "#6ad5e8") + 
        
        facet_grid(rows = dplyr::vars(nutrient_input), scales = "fixed") + 
        
        labs(x = bquote(CV[Local~scale]), y = bquote(CV[Meta-ecosystem~scale])) +
        theme_classic() + 
        theme(strip.background = element_blank(), strip.text = element_blank(),
              axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA),
              legend.position = "none")
      
      gg_beta <- ggplot(data = beta_temp, aes(x = abiotic, y = value.cv)) + 
        
        geom_point(pch = 1, alpha = 0.5) + 
        geom_smooth(method = 'loess', formula = "y ~ x", se = FALSE, color = "#6ad5e8") + 
        
        facet_grid(rows = dplyr::vars(nutrient_input), scales = "fixed", 
                   labeller = labeller(nutrient_input = function(x) paste0("Nutr. input: ", x))) + 
        
        labs(x = bquote(Variation[Nutrient~subsidies]), y = bquote(CV[beta])) +
        theme_classic() + 
        theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.5),
              axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA),
              legend.position = "none")
      
      cowplot::plot_grid(gg_alpha_gamma, gg_beta, ncol = 2)
    
  })
})

#### Biotic connectivity only ####

gg_dummy <- ggplot(data = data.frame(pop_n = c(8, 16, 32, 64, 128),  x = c(1, 2, 3, 4, 5),
                                     y = c(5, 4, 3, 2, 1)), 
                   aes(x = x, y = y, color = factor(pop_n))) +
  geom_point(size = 3.5) + 
  scale_color_viridis_d(name = "Population size", option = "C") + 
  theme_classic() +
  theme(legend.position = "bottom")

gg_biotic <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  purrr::map(c(ag_production = "ag_production", bg_production = "bg_production", ttl_production = "ttl_production"), function(part_i) {
    
    df_temp <- dplyr::filter(results_combined_df, scenario == scenario_i, 
                             part == part_i, abiotic == 0.0, include == "yes") |> 
      dplyr::select(row_id, pop_n, nutrient_input, measure, value.cv, abiotic, biotic)
    
    alpha_gamma_temp <- dplyr::filter(df_temp, measure %in% c("alpha", "gamma")) |> 
      tidyr::pivot_wider(names_from = measure, values_from = value.cv)
    
    beta_temp <- dplyr::filter(df_temp, measure == "beta")
    
    gg_alpha_gamma <- ggplot(data = alpha_gamma_temp, aes(x = alpha, y = gamma, color = pop_n)) + 
      
      geom_point(pch = 1, alpha = 0.5) + 
      geom_abline(color = "grey", linetype = 2) +
      geom_smooth(method = 'loess', formula = "y ~ x", se = FALSE) + 
      
      facet_grid(rows = dplyr::vars(nutrient_input), scales = "fixed") + 
      
      scale_color_viridis_d(option = "C") +
      
      labs(x = bquote(CV[Local~scale]), y = bquote(CV[Meta-ecosystem~scale])) +
      theme_classic() + 
      theme(strip.background = element_blank(), strip.text = element_blank(),
            axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA),
            legend.position = "none")
    
    gg_beta <- ggplot(data = beta_temp, aes(x = biotic, y = value.cv, color = pop_n)) +
      
      geom_point(pch = 1, alpha = 0.5) + 
      geom_smooth(method = 'loess', formula = "y ~ x", se = FALSE) + 
      
      facet_grid(rows = dplyr::vars(nutrient_input), scales = "fixed", 
                 labeller = labeller(nutrient_input = function(x) paste0("Nutr. input: ", x))) + 
      
      scale_color_viridis_d(option = "C") +
      
      labs(x = "Connectivity", y = bquote(CV[beta])) +
      theme_classic() + 
      theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.5),
            axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA),
            legend.position = "none")
    
    gg_combined <- cowplot::plot_grid(gg_alpha_gamma, gg_beta, ncol = 2)
    
    cowplot::plot_grid(gg_combined, cowplot::get_legend(gg_dummy),
                       nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))
    
  })
})

#### Save ggplot ####

overwrite <- FALSE

suppoRt::save_ggplot(plot = gg_abiotic$noise$ag_production, filename = "Figure-A5.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_abiotic$noise$bg_production, filename = "Figure-A6.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_biotic$noise$ag_production, filename = "Figure-A7.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_biotic$noise$bg_production, filename = "Figure-A8.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)
