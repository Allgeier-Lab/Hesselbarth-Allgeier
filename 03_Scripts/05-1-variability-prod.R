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
                                             died_total > threshold_mort ~ "no", 
                                           TRUE ~ "yes"))

#### Split data #### 

results_combined_list <- dplyr::left_join(x = dplyr::filter(results_combined_df, measure == "beta") |> 
                                            dplyr::select(scenario, row_id, part, pop_n, nutrient_input, include, value.cv),
                                          y = dplyr::filter(results_combined_df, measure == "gamma") |> 
                                            dplyr::select(scenario, row_id, part, pop_n, nutrient_input, include, value.prod), 
                                          by = c("scenario", "row_id", "part", "pop_n", "nutrient_input", "include")) |> 
  dplyr::filter(include == "yes") |> 
  dplyr::group_by(scenario, part, pop_n, nutrient_input) |>
  dplyr::group_split()
  
#### Fit regression model ####

full_lm_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
  
  df_temp_stand <- dplyr::select(df_temp, value.cv, value.prod) |> 
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) |> 
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
  
  lm_temp <- lm(value.prod ~ value.cv, data = df_temp_stand)
  
  broom::tidy(lm_temp) |> 
    dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                  pop_n = unique(df_temp$pop_n), nutrient_input = unique(df_temp$nutrient_input), .before = term) |> 
    dplyr::mutate(p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                                   p.value < 0.05 ~ "*", p.value >= 0.05 ~ ""), 
                  r2 = summary(lm_temp)$adj.r.squared, 
                  direction = dplyr::case_when(estimate < 0.0 ~ "decrease", 
                                               estimate > 0.0 ~ "increase"))
})

#### Setup ggplot ####

color_part <- c("ag_production" = "#c97644", "bg_production" = "#541f1b")

width_pos <- 0.35

#### Create ggplot model parameters ####

gg_coef_scenario <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  dplyr::filter(full_lm_df, scenario == scenario_i, part %in% c("ag_production", "bg_production"), term == "value.cv") |> 
    
    ggplot(aes(x = pop_n, y = estimate, color = part, fill = part, group = part)) +
    
    # adding geoms
    geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
    geom_linerange(aes(ymin = 0.0, ymax = estimate), alpha = 0.75, linewidth = 0.75, 
                   position = position_dodge(width = width_pos)) +
    geom_point(aes(shape = direction), size = 5.0, 
               position = position_dodge(width = width_pos)) +
    geom_text(aes(label = p.value.class), vjust = 0.75, size = 2.5, color = "white",
              position = position_dodge(width = width_pos)) +
      
    # facet grid
    facet_grid(rows = dplyr::vars(nutrient_input), scales = "fixed",
               labeller = labeller(nutrient_input = function(x) paste0("Nutr. input: ", x))) +
  
    # set scales and labs
    scale_color_manual(name = "", values = color_part) +
    scale_fill_manual(name = "", values = color_part, 
                      labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground")) +
    # scale_shape_manual(name = "", values = c("decrease" = 25, "increase" = 24)) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(full_lm_df$pop_n))) +
    scale_y_continuous(breaks = function(x) seq(min(x), max(x), length.out = 5), 
                       labels = function(x) round(x, 2)) +
      
    # setup theme
    labs(x = "Population size", y = "Parameter estimate") +
    guides(shape = "none", color = "none", 
           fill = guide_legend(override.aes = list(shape = 24, color = "white"))) +
    theme_classic(base_size = 10.0) +
    theme(axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA),
          strip.background = element_blank(), strip.text = element_text(hjust = 0.5),
          legend.position = "bottom")
  
})

#### Save ggplot ####

overwrite <- TRUE

suppoRt::save_ggplot(plot = gg_coef_scenario$noise, filename = "Figure-3.pdf",
                     path = "04_Figures/", width = width, height = height * 0.65,
                     units = units, dpi = dpi, overwrite = overwrite)
