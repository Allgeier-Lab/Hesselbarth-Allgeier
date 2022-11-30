##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create resulting figure of variability and synchrony

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-cv.R")
source("01_Functions/import-abundance.R")

#### Load simulated data ####

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

#### Run total regression model ####

full_lm_df <- dplyr::filter(results_combined_df, measure == "synchrony", 
                            part %in% c("ag_production", "bg_production", "ttl_production"), 
                            pop_n != 0, include == "yes") |> 
  dplyr::select(-c(measure, fish_biomass, value.sd, value.mn, value.prod, value.connectivity)) |> 
  dplyr::group_by(scenario, part, pop_n, nutrient_input) |> 
  dplyr::group_split() |> 
  purrr::map_dfr(function(df_temp) {
    
    df_temp_stand <- dplyr::select(df_temp, value.cv, biotic, abiotic) |> 
      dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) |> 
      dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
    
    lm_temp <- lm(value.cv ~ biotic * abiotic, data = df_temp_stand)
    
    broom::tidy(lm_temp) |> 
      dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part),
                    pop_n = unique(df_temp$pop_n), nutrient_input = unique(df_temp$nutrient_input), 
                    .before = "term") |> 
      dplyr::mutate(p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                                     p.value < 0.05 ~ "*", p.value >= 0.05 ~ ""), 
                    r2 = summary(lm_temp)$adj.r.squared)}) |> 
  dplyr::mutate(term = factor(term), p.value.class = factor(p.value.class, levels = c("*", "**", "***", ""))) |> 
  dplyr::filter(nutrient_input %in% c("low", "medium", "high"))

dplyr::filter(full_lm_df, scenario == "noise", part == "ag_production", term != "(Intercept)") |> 
  View()

#### Save results ####

purrr::walk(c("phase", "noise"), function(i) {
  dplyr::filter(full_lm_df, scenario == i) |> 
    suppoRt::save_rds(filename = paste0("lm_variability_synchrony_", i, ".rds"), 
                      path = "05_Results/", overwrite = FALSE)
})
