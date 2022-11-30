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

#### Load/wrangle simulated data ####

results_phase_df <- import_cv(path = "02_Data/result-phase.rds")

results_noise_df <- import_cv(path = "02_Data/result-noise.rds")

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise"))) |> 
  dplyr::filter(pop_n > 0)

#### Load mortality data ####

mortality_phase_df <- import_mortality(path = "02_Data/result-phase.rds") |> 
  dplyr::filter(pop_n > 0)

mortality_noise_df <- import_mortality(path = "02_Data/result-noise.rds") |> 
  dplyr::filter(pop_n > 0)

mortality_combined_df <- dplyr::bind_rows(phase = mortality_phase_df, noise = mortality_noise_df, 
                                          .id = "scenario") |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise"))) |> 
  dplyr::group_by(scenario, row_id, pop_n, nutrient_input) |> 
  dplyr::summarise(died_total = sum(died_total), .groups = "drop")

#### Filter quantiles mortality #### 

results_combined_df <- dplyr::group_by(mortality_combined_df, scenario, pop_n, nutrient_input) |> 
  dplyr::summarise(thres = quantile(died_total, probs = 0.75), .groups = "drop") |> 
  dplyr::right_join(y = mortality_combined_df, by = c("scenario", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(pop_n == 128 & nutrient_input == "low" & died_total > thres ~ "no", 
                                           TRUE ~ "yes")) |> 
  dplyr::select(scenario, row_id, include) |> 
  dplyr::right_join(y = results_combined_df, by = c("scenario", "row_id"))
  
#### CV densities ####

size_base <- 10

gg_cv <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  purrr::map(c(ag = "ag_production", bg =  "bg_production", ttl = "ttl_production"), function(part_i) {
    
    df_temp <- dplyr::filter(results_combined_df, scenario == scenario_i, part == part_i, 
                             measure %in% c("alpha", "gamma"), include == "yes") 
    
    df_extra <- dplyr::filter(results_combined_df, scenario == scenario_i, part == part_i, 
                              measure %in% c("alpha", "gamma"), include == "no") 
    
    # ggplot
    ggplot(df_temp, aes(x = pop_n, y = value.cv)) + 
      
      # geoms
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha = 0.25) +
      geom_jitter(data = df_extra, aes(x = pop_n, y = value.cv), col = "#ef8a47", alpha = 0.25) +
      
      # facet
      facet_grid(rows = vars(nutrient_input), cols = vars(measure), 
                 labeller = labeller(measure = c("alpha" = "Local scale",
                                                 "gamma" = "Meta-ecosystem scale"), 
                                     nutrient_input = function(x) paste("Nutr. input:", x))) +
      
      # themes
      labs(x = "Population size", y = "Coefficient of variation") +
      theme_classic(base_size = size_base) +
      theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.0),
            axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA))
    
  })
})

suppoRt::save_ggplot(plot = gg_cv$noise$ag, filename = "Figure-cv.pdf",
                     path = "04_Figures/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = FALSE)
