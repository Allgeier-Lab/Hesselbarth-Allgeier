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

### Alternative model ####

results_final_list <- dplyr::filter(results_final_df, measure %in% c("alpha", "gamma"), include == "yes",
                                    part %in% c("Aboveground", "Total"), treatment == "combined") |>
  dplyr::group_by(scenario, part, measure) |>
  dplyr::group_split()

names_list <- purrr::map(results_final_list, function(i) c(unique(i$scenario), unique(i$part)))

best_lm_list <- vector(mode = "list", length = length(results_final_list))

lm_summary_list <- vector(mode = "list", length = length(results_final_list))

for(i in 1:length(results_final_list)) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  df_temp <- results_final_list[[i]] |> 
    dplyr::mutate(value.cv = log(value.cv), abiotic = log(abiotic), biotic = log(biotic)) |>
    dplyr::mutate(biotic = (biotic - mean(biotic)) / sd(biotic),
                  abiotic = (abiotic - mean(abiotic)) / sd(abiotic))
  
  lm_temp <- lm(value.cv ~ nutrient_input + abiotic + pop_n + biotic + 
                  nutrient_input*pop_n + nutrient_input*biotic + pop_n*abiotic + 
                  abiotic*biotic, data = df_temp, na.action = "na.fail")
  
  lm_dredge <- MuMIn::dredge(lm_temp, extra = c("R^2"))
  
  lm_summary_list[[i]] <- subset(lm_dredge, subset = 1:3) |>
    tibble::as_tibble() |>
    dplyr::mutate(scenario = names_list[[i]][[1]], part = names_list[[i]][[2]], 
                  measure = unique(df_temp$measure), .before = "(Intercept)")
  
  best_lm_list[[i]] <- get.models(lm_dredge, subset = 1)
  
}

# convert to data.frame
lm_summary_df <- dplyr::bind_rows(lm_summary_list) |> 
  dplyr::filter(scenario == "noise") |>
  dplyr::select(-scenario, -logLik, -delta, -weight, -df) |>
  dplyr::rename("Part" = "part", "Intercept" = "(Intercept)", "Spatial variation" = "abiotic",
                "Connectivity" = "biotic", "Enrichment" = "nutrient_input", "Population size" = "pop_n",
                "Variation:Connectivity" = "abiotic:biotic", "Variation:Population" = "abiotic:pop_n",
                "Enrichment:Connectivity" = "biotic:nutrient_input", "Enrichment:Population" = "nutrient_input:pop_n") |>
  dplyr::mutate_if(is.numeric, round, digits = 3)

#### Save table and ggplot #### 

overwrite <- FALSE

if (overwrite) readr::write_csv2(x = lm_summary_df, file = "04_Figures/Table-1-alt.csv")


