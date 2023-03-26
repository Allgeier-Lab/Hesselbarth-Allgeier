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

#### Split data #### 

results_final_list <- dplyr::group_by(results_pe_pp_df, scenario, part) |>
  dplyr::group_split()

names_list <- purrr::map(results_final_list, function(i) c(unique(i$scenario), unique(i$part)))

#### Fit regression model and dredge ####

best_lm_list <- vector(mode = "list", length = length(results_final_list))

lm_summary_list <- vector(mode = "list", length = length(results_final_list))

for(i in 1:length(results_final_list)) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  df_temp <- results_final_list[[i]] |>
    dplyr::mutate(value.cv = (value.cv - mean(value.cv)) / sd(value.cv),
                  biotic = (biotic - mean(biotic)) / sd(biotic),
                  abiotic = (abiotic - mean(abiotic)) / sd(abiotic))
  
  if (unique(df_temp$part) == "Aboveground") {
    
    df_temp$value.prod <- sqrt(df_temp$value.prod)
    
  } else if (unique(df_temp$part) == "Total") {
    
    df_temp$value.prod <- log(df_temp$value.prod)
    
  }
  
  lm_temp <- lm(value.prod ~ abiotic + biotic + nutrient_input + pop_n,
                data = df_temp, na.action = "na.fail")
  
  lm_dredge <- MuMIn::dredge(lm_temp, extra = c("R^2"))
  
  lm_summary_list[[i]] <- subset(lm_dredge, subset = 1:3) |>
    tibble::as_tibble() |>
    dplyr::mutate(scenario = names_list[[i]][[1]], part = names_list[[i]][[2]], 
                  .before = "(Intercept)")
  
  best_lm_list[[i]] <- get.models(lm_dredge, subset = 1)
  
}

# convert to data.frame
lm_summary_df <- dplyr::bind_rows(lm_summary_list) |> 
  dplyr::filter(scenario == "noise") |>
  dplyr::select(-scenario, -logLik, -delta, -weight, -df) |>
  dplyr::select(part, `(Intercept)`, abiotic, biotic, nutrient_input, pop_n, 
                `R^2`, AICc) |> 
  dplyr::rename("Part" = "part", "Intercept" = "(Intercept)", "Spatial variation" = "abiotic",
                "Connectivity" = "biotic", "Enrichment" = "nutrient_input", 
                "Population size" = "pop_n") |>
  dplyr::mutate_if(is.numeric, round, digits = 3)

#### Check models ####

# gg_assumptions <- purrr::map(best_lm_list, function(i) {
#   
#   plot_lm(i[[1]])
#   
# })

#### Save results ####

overwrite <- FALSE

if (overwrite) readr::write_csv2(x = lm_summary_df, file = "04_Figures/Table-2-alt.csv")
