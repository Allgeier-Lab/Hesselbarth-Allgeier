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

#### Join PE and PP values ####

results_pe_pp_df <- dplyr::left_join(x = dplyr::filter(results_final_df, measure == "beta", treatment == "combined", 
                                                       include == "yes", part != "Total") |> 
                                       dplyr::select(scenario, row_id, part, pop_n, nutrient_input, biotic, abiotic, value.cv),
                                     y = dplyr::filter(results_final_df, measure == "gamma", treatment == "combined", 
                                                       include == "yes", part != "Total") |> 
                                       dplyr::select(scenario, row_id, part, pop_n, nutrient_input, biotic, abiotic, value.prod), 
                                     by = c("scenario", "row_id", "part", "pop_n", "nutrient_input", "biotic", "abiotic"))

#### Split data #### 

results_final_list <- dplyr::group_by(results_pe_pp_df, scenario, part) |>
  dplyr::group_split()

names_list <- purrr::map(results_final_list, function(i) c(unique(i$scenario), unique(i$part)))

#### Fit regression model and dredge ####

best_lm_list <- purrr::map(seq_along(results_final_list), function(i) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  df_temp <- results_final_list[[i]]
  
  lm_temp <- lm(value.prod ~ value.cv + abiotic + biotic + nutrient_input + pop_n,
                data = df_temp, na.action = "na.fail")
  
  lm_dredge <- MuMIn::dredge(lm_temp)
  
  get.models(lm_dredge, delta == 0.0)[[1]]
  
})

#### Summarize model and rel importance

lm_summary_df <- purrr::map_dfr(seq_along(results_final_list), function(i) {
  
  broom::tidy(best_lm_list[[i]]) |> 
    dplyr::mutate(r2 = summary(best_lm_list[[i]])$adj.r.squared) |> 
    dplyr::mutate(scenario = names_list[[i]][[1]], part = names_list[[i]][[2]], 
                  .before = term)}) |> 
  dplyr::mutate(p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ ""))

#### Setup ggplot ####

base_size <- 10.0

color_pop_n <- MetBrewer::met.brewer(name = "Isfahan2", n = 5, type = "discrete")

#### Create ggplot: Main ####

gg_treatments <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  df_temp <- dplyr::filter(results_pe_pp_df, scenario == scenario_i)
  
  ggplot(data = df_temp, aes(x = log(value.cv), y = log(value.prod), 
                                         color = pop_n, shape = nutrient_input)) +
    
    # adding geoms
    geom_point(alpha = 1, size = 1) +
    geom_smooth(aes(linetype = nutrient_input), 
                method = "lm", se = FALSE, formula = "y ~ x", linewidth = 0.5) +
    
    # scales
    scale_color_manual(name = "Population size", values = color_pop_n) +
    scale_shape_manual(name = "Amount subsidies", values = c("low" = 0, "medium" = 1, "high" = 2)) + 
    scale_linetype_manual(name = "Amount subsidies", values = c("low" = "solid", "medium" = "dotted", "high" = "dashed")) + 
    
    # facet wrap by part
    facet_wrap(. ~ part, ncol = 2, scales = "free", 
               labeller = labeller(part = c("Aboveground" = "a)", "Belowground" = "b)"))) +
    
    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    
    # theming
    labs(x = "log(Portfolio effect)", y = "log(Total primary production)") +
    theme_classic(base_size = base_size) + 
    theme(legend.position = "bottom", axis.line = element_blank(), strip.background = element_blank(), 
          strip.text = element_text(hjust = 0))

  })

#### Save ggplots ####

dplyr::filter(lm_summary_df, scenario == "noise") |> 
  dplyr::select(-scenario, -r2) |> 
  dplyr::mutate_if(is.numeric, round, digits = 3) |> 
  readr::write_csv(file = "04_Figures/Table-2.csv")

suppoRt::save_ggplot(plot = gg_treatments$noise, filename = "Figure-4.pdf",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = FALSE)
