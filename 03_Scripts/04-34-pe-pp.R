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
source("01_Functions/plot-lm.R")

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

#### Join PE and PP values ####

results_pe_pp_df <- dplyr::left_join(x = dplyr::filter(results_final_df, measure == "beta", treatment == "combined", 
                                                       include == "yes", part != "Belowground") |> 
                                       dplyr::select(scenario, row_id, part, pop_n, nutrient_input, biotic, abiotic, value.cv),
                                     y = dplyr::filter(results_final_df, measure == "gamma", treatment == "combined", 
                                                       include == "yes", part != "Belowground") |> 
                                       dplyr::select(scenario, row_id, part, pop_n, nutrient_input, biotic, abiotic, value.prod), 
                                     by = c("scenario", "row_id", "part", "pop_n", "nutrient_input", "biotic", "abiotic"))

#### Setup ggplot ####

base_size <- 10.0

#### Create ggplot ####

gg_pe_pp <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  gg_part <- purrr::map(c(Aboveground = "Aboveground", Total = "Total"), function(part_i) {
    
    if (part_i == "Aboveground") {
      
      legend_position <- "none"
      x_lab <- " "
      lab_temp <- labeller(nutrient_input = c(low = "low", medium = "medium", high = "high"))
      
      
    } else if (part_i == "Total") {
      
      legend_position <- "bottom"
      x_lab <- expression(paste("Portfolio effect ", italic(CV), italic(beta)))
      lab_temp <- labeller(nutrient_input = c(low = " ", medium = " ", high = " "))
    }
    
    df_temp <- dplyr::filter(results_pe_pp_df, scenario == scenario_i, part == part_i) |> 
      dplyr::mutate(pe_class = cut(value.cv, breaks = 50)) |> 
      dplyr::group_by(pe_class) |> 
      dplyr::summarise(value.prod = mean(value.prod), 
                       n = dplyr::n())
    
    ggplot(data = df_temp) +
      
      # adding tiles with PP
      geom_col(aes(x = pe_class, y = n, fill = value.prod)) +
      
      # adding box
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set scales 
      scale_fill_viridis_c(name = "Primary production") + 
      # scale_color_manual(name = "Effect", values = c("black", "grey", "white"), 
      #                    limits = c("Connectivity", "Variation", "Portfolio")) +
      # guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
      
      # theming
      theme_classic(base_size = base_size) + 
      theme(legend.position = "right", axis.line = element_blank(), axis.title = element_blank(),
            axis.text.x = element_blank())
    
  })
  
  cowplot::plot_grid(gg_part$Aboveground, gg_part$Total, 
                     nrow = 2, labels = c("a)", "b)"), hjust = -0.25,
                     rel_heights = c(0.5, 0.5)) |> 
    cowplot::ggdraw(xlim = c(-0.025, 1.0), ylim = c(-0.025, 1.0)) + 
    cowplot::draw_label("Count", x = -0.015, y = 0.5, angle = 90, size = base_size) +
    cowplot::draw_label("Portfolio effect", x = 0.5, y = 0.015, angle = 0, size = base_size)
  
})

suppoRt::save_ggplot(plot = gg_pe_pp$noise, filename = "Figure-4-alt.pdf",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_pe_pp$noise, filename = "Figure-4-alt.jpg",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = FALSE)
