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
source("01_Functions/import-abundance.R")

#### Load/wrangle simulated data ####

results_phase_df <- import_cv(path = "02_Data/result-phase.rds")

results_noise_df <- import_cv(path = "02_Data/result-noise.rds")

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))

#### Load abundance data ####

abundance_phase_df <- import_abundance(path = "02_Data/result-phase.rds") |> 
  dplyr::filter(pop_n > 0)

abundance_noise_df <- import_abundance(path = "02_Data/result-noise.rds") |> 
  dplyr::filter(pop_n > 0)

abundance_combined_df <- dplyr::bind_rows(phase = abundance_phase_df, noise = abundance_noise_df,
                                          .id = "scenario") |>
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise"))) |> 
  dplyr::group_by(scenario, row_id, pop_n, nutrient_input) |>
  dplyr::summarise(abundance_max = max(mean), .groups = "drop")

#### Filter abundances/mortality #### 

results_combined_df <- dplyr::left_join(x = results_combined_df, y = abundance_combined_df, 
                                        by = c("scenario", "row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(pop_n == 128 & nutrient_input == "low" &
                                             abundance_max > threshold_abundance ~ "no", 
                                           TRUE ~ "yes"))

#### Split data #### 

results_combined_list <- dplyr::filter(results_combined_df, measure == "gamma",
                                       part %in% c("ag_production", "bg_production", "ttl_production"), 
                                       pop_n != 0, include == "yes") |>
  dplyr::group_by(scenario, part, pop_n, nutrient_input) |>
  dplyr::group_split()

#### Fit regression model ####

full_lm_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
  
  df_temp_stand <- dplyr::select(df_temp, value.prod, biotic, abiotic) |> 
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) |>
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
  
 lm_temp <- lm(value.prod ~ biotic * abiotic, data = df_temp_stand)
  
  broom::tidy(lm_temp) |> 
    dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                  measure = unique(df_temp$measure), pop_n = unique(df_temp$pop_n),
                  nutrient_input = unique(df_temp$nutrient_input), .before = term) |> 
    dplyr::mutate(p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                                   p.value < 0.05 ~ "*", p.value >= 0.05 ~ ""), 
                  r2 = summary(lm_temp)$adj.r.squared, 
                  direction = dplyr::case_when(estimate < 0.0 ~ "decrease", 
                                               estimate > 0.0 ~ "increase"))
})

#### Save linear regression results ####

# purrr::walk(c("phase", "noise"), function(i) {
#   dplyr::filter(full_lm_df, scenario == i) |> 
#     suppoRt::save_rds(filename = paste0("lm_variability_prod_", i, ".rds"), 
#                       path = "05_Results/", overwrite = FALSE)
# })

#### Setup ggplot ####

color_term <- c("biotic" = "#00af73", "abiotic" = "#006d9a", "biotic:abiotic" = "#ffaa3a")

color_term <- c("biotic" = "#00af73", "abiotic" = "#006d9a",
                "biotic:abiotic" = "#ffaa3a", "residuals" = "grey")

size_point <- 5 #* 0.75
size_line <- 0.75
size_text <- 2.5
size_base <- 10.0

alpha <- 0.5

width_pos <- 0.65

# dplyr::filter(full_lm_df, term %in% c("biotic", "abiotic"))

#### Create ggplot model parameters ####

gg_coef_scenario <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  range_estimate <- dplyr::filter(full_lm_df, scenario == scenario_i, 
                                  term %in% c("biotic", "abiotic", "biotic:abiotic")) |> 
    dplyr::pull(estimate) |> 
    range() |> 
    abs() |> 
    max()
  
  gg_part <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_j) {
  
    dplyr::filter(full_lm_df, scenario == scenario_i, 
                  part == part_j, term %in% c("biotic", "abiotic", "biotic:abiotic")) |> 
      
      ggplot() +
      
      # adding geoms
      geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
      geom_linerange(aes(x = pop_n, ymin = 0.0, ymax = estimate, color = term),
                     alpha = 0.25, position = position_dodge(width = width_pos), linewidth = size_line) +
      geom_point(aes(x = pop_n, y = estimate, fill = term, color = term, size = term, shape = direction),
                 position = position_dodge(width = width_pos)) +
      geom_text(aes(x = pop_n, y = estimate, label = p.value.class, group = term), 
                position = position_dodge(width = width_pos),vjust = 0.75, size = size_text, color = "white") +

      # facet gridding
      facet_grid(rows = dplyr::vars(nutrient_input), labeller = labeller(nutrient_input = function(x) paste0("Nutr. input: ", x))) +

      # set scales and labs
      scale_color_manual(name = "", values = color_term,
                         labels = c("biotic" = "Consumer behavior", "abiotic" =  "Abiotic subsidies",
                                    "biotic:abiotic" = "Interaction Behavior:Subsidies", "residuals" = "Residuals")) +
      scale_fill_manual(name = "", values =  color_term,
                        labels = c("biotic" = "Consumer behavior", "abiotic" =  "Abiotic subsidies",
                                   "biotic:abiotic" = "Subsidies:Behavior", "residuals" = "Residuals")) +
      scale_shape_manual(name = "", values = c("decrease" = 25, "increase" = 24)) +
      scale_size_manual(values = c("biotic" = size_point, "abiotic" = size_point, 
                                   "biotic:abiotic" = size_point / 2)) +
      coord_flip() +
      scale_x_discrete(limits = rev(levels(full_lm_df$pop_n)[-1])) +
      scale_y_continuous(limits = c(-range_estimate, range_estimate), labels = function(x) round(x, 2),
                         breaks = seq(-range_estimate, range_estimate, length.out = 5)) +

      # setup theme
      labs(x = "Population size", y = "Parameter estimate") +
      guides(color = guide_legend(order = 1, override.aes = list(shape = 17, size = size_point)),
             shape = guide_legend(order = 2, override.aes = list(size = size_point)), 
             fill = "none", size = "none") +
      theme_classic(base_size = size_base) +
      theme(axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA),
            strip.background = element_blank(), strip.text = element_text(hjust = 0.1), plot.title = element_text(size = size_base),
            legend.position = "bottom",
            plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt"))
    
  })
  
  names(gg_part) <- c("ag_production", "bg_production", "ttl_production")
  
  gg_part
  
})

#### Save ggplot #### 

suppoRt::save_ggplot(plot = gg_coef_scenario$noise$ag_production, filename = "Figure-3.pdf",
                     path = "04_Figures/", width = width, height = height * 0.75,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_coef_scenario$noise$bg_production, filename = "Figure-A5.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.75,
                     units = units, dpi = dpi, overwrite = FALSE)
