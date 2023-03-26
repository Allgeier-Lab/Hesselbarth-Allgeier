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


#### Plot Total vs AG/BG ####

base_size <- 10

df_pp <- dplyr::filter(results_final_df, measure == "gamma", include == "yes") |> 
  dplyr::select(row_id, scenario, part, value.prod, pop_n, nutrient_input) |> 
  tidyr::pivot_wider(names_from = part, values_from = value.prod) |> 
  tidyr::pivot_longer(c(Aboveground, Belowground), names_to = "Part", values_to = "Value")

gg_prod <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  df_temp <- dplyr::filter(df_pp, scenario == scenario_i)
  
  ggplot(data = df_temp, aes(x = Value, y = Total, color = Part)) + 
    
    # point geom
    geom_point(shape = 19, alpha = 0.65) + 

    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    
    # set scales
    scale_x_continuous(limits = c(0, max(c(df_temp$Total, df_temp$Value)))) +
    scale_y_continuous(limits = c(0, max(c(df_temp$Total, df_temp$Value)))) +
    scale_color_manual(name = "", values = c(Aboveground = "#00496f", Belowground = "#dd4124")) + 
    coord_equal() +
    
    # set theme
    labs(x = "Aboveground/Belowground", y = "Total", title = expression(paste("Primary production [gDW ",m^-2, d^-1, "]"))) +
    theme_classic(base_size = base_size) + 
    theme(axis.line = element_blank(), legend.position = "bottom")
  
})

#### Plot all PE values ####

df_pe <- dplyr::left_join(x = dplyr::filter(results_final_df, measure == "beta", 
                                            include == "yes") |> 
                            dplyr::select(scenario, row_id, part, pop_n, nutrient_input, biotic, abiotic, value.cv),
                          y = dplyr::filter(results_final_df, measure == "gamma",
                                            include == "yes") |> 
                            dplyr::select(scenario, row_id, part, pop_n, nutrient_input, biotic, abiotic, value.prod), 
                          by = c("scenario", "row_id", "part", "pop_n", "nutrient_input", "biotic", "abiotic"))

gg_pe <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  df_temp <- dplyr::filter(df_pe, scenario == scenario_i)
  
  ggplot(data = df_temp, aes(x = value.cv, y = value.prod, color = part)) + 
    
    # point geom
    geom_point(shape = 19, alpha = 0.65) + 
    
    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    
    # set scales
    scale_color_manual(name = "", values = c(Aboveground = "#00496f", Belowground = "#dd4124", 
                                             Total = "#edd746")) + 

    # set theme
    labs(x = expression(paste("Portfolio effect ", italic(CV), italic(beta))), 
         y = expression(paste("Primary production [gDW ",m^-2, d^-1, "]"))) +
    theme_classic(base_size = base_size) + 
    theme(axis.line = element_blank(), legend.position = "bottom")
  
})

#### Save figures #### 

overwrite <- TRUE

suppoRt::save_ggplot(plot = gg_prod$noise, filename = "Figure-A2.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.45,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_prod$noise, filename = "Figure-A2.jpg",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.45,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_pe$noise, filename = "Figure-A3.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.45,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_pe$noise, filename = "Figure-A3.jpg",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.45,
                     units = units, dpi = dpi, overwrite = overwrite)
