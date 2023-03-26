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

results_final_list <- dplyr::group_by(results_pe_pp_df, scenario, part, nutrient_input, 
                                      pop_n) |>
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
  
  lm_temp <- lm(value.prod ~ value.cv, data = df_temp, na.action = "na.fail")
  
  lm_dredge <- MuMIn::dredge(lm_temp, extra = c("R^2"))
  
  lm_summary_list[[i]] <- subset(lm_dredge, subset = 1:3) |>
    tibble::as_tibble() |>
    dplyr::mutate(scenario = names_list[[i]][[1]], part = names_list[[i]][[2]], 
                  nutrient_input = unique(df_temp$nutrient_input), pop_n = unique(df_temp$pop_n),
                  .before = "(Intercept)")
  
  best_lm_list[[i]] <- get.models(lm_dredge, subset = 1)
  
}

# convert to data.frame
lm_summary_df <- dplyr::bind_rows(lm_summary_list)

#### Check models ####

# gg_assumptions <- purrr::map(best_lm_list, function(i) {
#   
#   plot_lm(i[[1]])
#   
# })

#### Setup ggplot ####

base_size <- 10.0

#### Create ggplot ####

gg_pe_pp <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  gg_part <- purrr::map(c(Aboveground = "Aboveground", Total = "Total"), function(part_i) {
    
    lm_temp <- dplyr::filter(lm_summary_df, scenario == scenario_i, part == part_i) |> 
      dplyr::group_by(scenario, part, nutrient_input, pop_n) |> 
      dplyr::top_n(n = 1, wt = -AICc) |> 
      dplyr::ungroup()
    
    df_temp <- dplyr::filter(results_pe_pp_df, scenario == scenario_i, part == part_i) |> 
      dplyr::group_by(pop_n, nutrient_input) |> 
      dplyr::summarise(value.prod = mean(value.prod), .groups = "drop")
    
    if (part_i == "Aboveground") {
    
      legend_position <- "none"
      x_lab <- " "
      lab_temp <- labeller(nutrient_input = c(low = "low", medium = "medium", high = "high"))
      
      
    } else if (part_i == "Total") {
      
      legend_position <- "bottom"
      x_lab <- expression(paste("Portfolio effect ", italic(CV), italic(beta)))
      lab_temp <- labeller(nutrient_input = c(low = " ", medium = " ", high = " "))
    }
    
    ggplot(data = df_temp) +
      
      # adding tiles with PP
      geom_raster(aes(x = pop_n, y = nutrient_input, fill = value.prod)) + 
      
      geom_text(data = lm_temp, aes(x = pop_n, y = nutrient_input, label = round(value.cv, 3))) +
      
      # adding box
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set scales 
      scale_fill_viridis_c(name = "Primary production") + 
      
      # theming
      theme_classic(base_size = base_size) + 
      theme(legend.position = "right", axis.line = element_blank(), axis.title = element_blank(),
            strip.background = element_blank(), strip.text = element_text(hjust = 0.0))
      
  })
  
  cowplot::plot_grid(gg_part$Aboveground, gg_part$Total, 
                     nrow = 2, labels = c("a)", "b)"), hjust = -0.25,
                     rel_heights = c(0.5, 0.5)) |> 
    cowplot::ggdraw(xlim = c(-0.025, 1.0), ylim = c(-0.025, 1.0)) + 
    cowplot::draw_label("Enrichment treatment", x = -0.015, y = 0.5, angle = 90, size = base_size) +
    cowplot::draw_label("Population size", x = 0.5, y = 0.015, angle = 0, size = base_size)
  
})

suppoRt::save_ggplot(plot = gg_pe_pp$noise, filename = "Figure-4.pdf",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_pe_pp$noise, filename = "Figure-4.jpg",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = FALSE)
