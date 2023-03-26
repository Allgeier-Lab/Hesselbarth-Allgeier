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

#### Split data #### 

results_final_list <- dplyr::filter(results_final_df, measure == "beta", include == "yes",
                                    part %in% c("Aboveground", "Total"), treatment == "combined") |>
  dplyr::group_by(scenario, part) |>
  dplyr::group_split()

names_list <- purrr::map(results_final_list, function(i) c(unique(i$scenario), unique(i$part)))

#### Fit regression model and dredge ####

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
                  .before = "(Intercept)")
  
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

#### Relative importance ####

rel_importance_df <- purrr::map_dfr(seq_along(results_final_list), function(i) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  rel_r2 <- relaimpo::boot.relimp(best_lm_list[[i]][[1]], type = "lmg", b = 500, level = 0.95, 
                                  fixed = FALSE) |>
    relaimpo::booteval.relimp(bty = "basic")
  
  tibble::tibble(
    scenario = names_list[[i]][[1]], part = names_list[[i]][[2]],
    beta = c(rel_r2@namen[-1], "residual"), mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)),
    lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA)) |>
    dplyr::mutate(dplyr::across(mean:higher, ~ .x * 100), 
                  lower = dplyr::case_when(lower < 0 ~ 0.0, TRUE ~ lower))}) |> 
  dplyr::mutate(beta = factor(beta, levels = c("abiotic", "biotic", "nutrient_input", "pop_n", 
                                               "nutrient_input:pop_n", "biotic:nutrient_input", 
                                               "abiotic:pop_n", "abiotic:biotic", "residual"), 
                              labels = c("abiotic" = "Spatial variation", "biotic" = "Connectivity", 
                                         "nutrient_input" = "Enrichment", "pop_n" = "Population size",
                                         "nutrient_input:pop_n" = "Enrichment:Population", "biotic:nutrient_input" = "Enrichment:Connectivity", 
                                         "abiotic:pop_n" = "Variation:Population", "abiotic:biotic" = "Variation:Connectivity", 
                                         "residual" = "Residuals")))

#### Marginal means ####

marginal_means <- purrr::map_dfr(seq_along(results_final_list), function(i) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  modelbased::estimate_means(best_lm_list[[i]][[1]], at = c("pop_n", "nutrient_input"), 
                             ci = 0.95) |> 
    tibble::as_tibble() |> 
    dplyr::mutate(scenario = names_list[[i]][[1]], part = names_list[[i]][[2]], 
                  .before = "pop_n")

})

#### Create dummy plot #### 

base_size <- 12.0

w <- 0.5

color_imp <- c("Spatial variation" = "#377EB8", "Connectivity" = "#4DAF4A",
               "Enrichment" = "#984EA3", "Population size" = "#A65628",
               "Enrichment:Population" = "#FF7F00", "Enrichment:Connectivity" = "#FFFF33",
               "Variation:Population" = "#E41A1C", "Variation:Connectivity" = "#F781BF",
               "Residuals" = "grey85")

color_means <- c(low = "#007895", medium = "#f2d445", high = "#e93c26")

gg_dummy <- data.frame(nutrient_input = c("low", "low", "low", "medium", "medium", "medium", "high", "high"), 
                       beta = unique(rel_importance_df$beta)[-9],
                       x = 1:8, y = 1:8) |> 
  dplyr::mutate(nutrient_input = factor(nutrient_input, levels = c("low", "medium", "high"))) |> 
  
  # create ggplot
  ggplot(aes(x = x, y = y)) +
  
  # adding geoims
  geom_point(aes(color = nutrient_input)) + 
  geom_errorbar(aes(color = nutrient_input, ymin = y - 1, ymax = y + 1), width = w) +
  geom_col(aes(fill = beta)) + 
  
  # change scales
  scale_color_manual(name = "Nutr.\nenrichment", values = color_means) +
  scale_fill_manual(name = "Rel. imp", values = color_imp) +
  
  # themeing
  guides(color = guide_legend(order = 1, nrow = 2, byrow = TRUE), 
         fill = guide_legend(order = 2, nrow = 3, byrow = FALSE)) +
  theme_classic(base_size = base_size * 0.9) +
  theme(legend.position = "bottom")

#### Create ggplot ####

gg_relimp <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  gg_part <- purrr::map(c(Aboveground = "Aboveground", Total = "Total"), function(part_i) {
    
    mean_temp <- dplyr::filter(marginal_means, scenario == scenario_i, part == part_i) 
    
    importance_temp <- dplyr::filter(rel_importance_df, scenario == scenario_i, part == part_i) |>
      dplyr::mutate(part = " ")
    
    pe_max <- dplyr::filter(marginal_means, scenario == scenario_i) |> 
      dplyr::pull(CI_high) |> 
      max() |> 
      exp()
    
    gg_means <- ggplot(data = mean_temp) +
      
      # adding mean and error bars
      geom_point(aes(x = pop_n, y = exp(Mean), color = nutrient_input)) + 
      geom_errorbar(aes(x = pop_n, ymin = exp(CI_low), ymax = exp(CI_high), color = nutrient_input), 
                    width = 0.0) +
      
      geom_line(aes(x = pop_n, y = exp(Mean), color = nutrient_input, group = nutrient_input), 
                alpha = 0.25) + 
      
      # adding box
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set colors
      scale_color_manual(name = "Amount nutrient subsidies", values = color_means) +
      scale_y_continuous(limits = c(1.0, pe_max)) +
      
      # change theme
      labs(y = "", x = "") +
      theme_classic(base_size = base_size) +
      theme(legend.position = "none", axis.line = element_blank(),
            strip.background = element_blank(), strip.text = element_text(hjust = 0.0))
    
    gg_importance <- ggplot(data = importance_temp) +
      
      # zero line
      geom_hline(yintercept = 0.0, linetype = 2, color = "grey95") +
      
      # relative importance bars
      geom_col(aes(x = part, y = mean, fill = forcats::fct_rev(beta)), position = "stack") +
      
      # adding box
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set scales
      scale_fill_manual(name = "", values = color_imp) +
      scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100),
                         breaks = seq(0, 100, 25)) +
      coord_flip() +

      # labels and themes
      labs(y = "", x = "") +
      theme_classic(base_size = base_size) +
      theme(legend.position = "none", axis.line = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_text(angle = 90))
    
    # cowplot::plot_grid(gg_importance, gg_means, nrow = 2, rel_heights = c(0.3, 0.7))
    
    list(imp = gg_importance, means = gg_means)
    
  })
  
  gg_imp <- cowplot::plot_grid(gg_part$Aboveground$imp, gg_part$Total$imp,
                               ncol = 2, labels = c("a)", "b)")) |> 
    cowplot::ggdraw(xlim = c(-0.025, 1.0), ylim = c(-0.025, 1.0)) + 
    cowplot::draw_label(label = "Relative importance [%]", x = 0.5, y = 0.15, size = base_size)
  
  gg_means <- cowplot::plot_grid(gg_part$Aboveground$means, gg_part$Total$means, ncol = 2) |> 
    cowplot::ggdraw(xlim = c(-0.025, 1.0), ylim = c(-0.025, 1.0)) + 
    cowplot::draw_label(label = "Population size", x = 0.5, y = 0.05, angle = 0, size = base_size) + 
    cowplot::draw_label(expression(paste("Portfolio effect ", italic(cv), italic(beta))),
                        x = -0.01, y = 0.5, angle = 90, size = base_size)
  
  gg_combined <- cowplot::plot_grid(gg_imp, gg_means, nrow = 2, rel_heights = c(0.25, 0.75))
  
  cowplot::plot_grid(gg_combined, cowplot::get_legend(gg_dummy), nrow = 2, ncol = 1, 
                     rel_heights = c(0.85, 0.15))
    
})

#### Save table and ggplot #### 

overwrite <- FALSE

if (overwrite) readr::write_csv2(x = lm_summary_df, file = "04_Figures/Table-1.csv")

suppoRt::save_ggplot(plot = gg_relimp$noise, filename = "Figure-3.pdf",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_relimp$noise, filename = "Figure-3.jpg",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = overwrite)

