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
                part = factor(part, levels = c("ag_production", "bg_production", "ttl_production"),
                              labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                         "ttl_production" = "Total")),
                treatment = factor(treatment, levels = c("subsidies", "connectivity", "combined")))

#### Split data #### 

results_final_list <- dplyr::filter(results_final_df, scenario == "noise", measure != "synchrony", 
                                    include == "yes", part %in% c("Aboveground", "Total"), 
                                    treatment == "combined") |>  
  dplyr::group_by(part, measure) |>
  dplyr::group_split()

names_list <- purrr::map(results_final_list, function(i) as.character(c(unique(i$scenario), unique(i$part), unique(i$measure))))

#### Fit regression model and dredge ####

best_lm_list <- vector(mode = "list", length = length(results_final_list))

lm_summary_list <- vector(mode = "list", length = length(results_final_list))

for (i in 1:length(results_final_list)) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  df_temp <- results_final_list[[i]] |> 
    dplyr::mutate(value.cv = log(value.cv), abiotic = log(abiotic), biotic = log(biotic)) |>
    dplyr::mutate(value.cv = (value.cv - mean(value.cv)) / sd(value.cv),
                  biotic = (biotic - mean(biotic)) / sd(biotic),
                  abiotic = (abiotic - mean(abiotic)) / sd(abiotic))
  
  lm_temp <- lm(value.cv ~ nutrient_input + abiotic + pop_n + biotic + 
                  nutrient_input*pop_n + nutrient_input*biotic + pop_n*abiotic + 
                  abiotic*biotic, data = df_temp, na.action = "na.fail")
  
  lm_dredge <- MuMIn::dredge(lm_temp, extra = c("R^2"))
  
  lm_summary_list[[i]] <- subset(lm_dredge, subset = 1:3) |>
    tibble::as_tibble() |>
    dplyr::mutate(scenario = names_list[[i]][[1]], part = names_list[[i]][[2]], measure = names_list[[i]][[3]],
                  .before = "(Intercept)")
  
  best_lm_list[[i]] <- get.models(lm_dredge, subset = 1)
  
}

# convert to data.frame
lm_summary_df <- dplyr::bind_rows(lm_summary_list) |> 
  dplyr::rename("intercept" = "(Intercept)", "spatialVariation" = "abiotic",
                "connectivity" = "biotic", "enrichment" = "nutrient_input", "popSize" = "pop_n",
                "var:connect" = "abiotic:biotic", "var:pop" = "abiotic:pop_n",
                "enrich:connect" = "biotic:nutrient_input", "enrich:pop" = "nutrient_input:pop_n") |>
  dplyr::mutate(measure = factor(measure, levels = c("alpha", "gamma" , "beta"), 
                                 labels = c("alpha" = "Local scale", "gamma" = "Meta-ecosystem scale", 
                                            "beta" = "Portfolio effect"))) |> 
  dplyr::mutate_if(is.numeric, round, digits = 3)

#### Relative importance ####

rel_importance_df <- purrr::map_dfr(seq_along(results_final_list), function(i) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  rel_r2 <- relaimpo::boot.relimp(best_lm_list[[i]][[1]], type = "lmg", b = 500, level = 0.95, 
                                  fixed = FALSE) |>
    relaimpo::booteval.relimp(bty = "basic")
  
  tibble::tibble(
    scenario = names_list[[i]][[1]], part = names_list[[i]][[2]], measure = names_list[[i]][[3]],
    beta = c(rel_r2@namen[-1], "residual"), mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)),
    lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA)) |>
    dplyr::mutate(dplyr::across(mean:higher, ~ .x * 100), 
                  lower = dplyr::case_when(lower < 0 ~ 0.0, TRUE ~ lower))}) |> 
  dplyr::mutate(measure = factor(measure, levels = c("alpha", "gamma" , "beta"), 
                                 labels = c("alpha" = "Local scale", "gamma" = "Meta-ecosystem scale", 
                                            "beta" = "Portfolio effect")),
                beta = factor(beta, levels = c("abiotic", "biotic", "nutrient_input", "pop_n", 
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
    dplyr::mutate(scenario = names_list[[i]][[1]], part = names_list[[i]][[2]], measure = names_list[[i]][[3]],
                  .before = "pop_n")}) |> 
  dplyr::mutate(measure = factor(measure, levels = c("alpha", "gamma" , "beta"), 
                                 labels = c("alpha" = "Local scale", "gamma" = "Meta-ecosystem scale", 
                                            "beta" = "Portfolio effect")))

#### Create dummy plot #### 

base_size <- 12.0

w <- 0.5

# color_scale <- c("Spatial variation" = "#377EB8", "Connectivity" = "#4DAF4A",
#                "Enrichment" = "#984EA3", "Population size" = "#A65628",
#                "Enrichment:Population" = "#FF7F00", "Enrichment:Connectivity" = "#FFFF33",
#                "Variation:Population" = "#E41A1C", "Variation:Connectivity" = "#F781BF",
#                "Residuals" = "grey85")

# color_scale <- c("Local scale" = "#4DAF4A", "Meta-ecosystem scale" = "#F781BF", "Portfolio effect" = "#A65628")
color_scale <- c("Local scale" = "#0f7ba2", "Meta-ecosystem scale" = "#43b284", "Portfolio effect" = "#fab255")

color_enrich <- c(low = "#007895", medium = "#f2d445", high = "#e93c26")

gg_dummy <- data.frame(measure = unique(rel_importance_df$measure), x = 1:3, y = 1:3) |>
  
  # create ggplot
  ggplot(aes(x = x, y = y)) +
  
  # adding geom
  geom_col(aes(fill = measure)) + 
  
  # change scales
  scale_fill_manual(name = "", values = color_scale) +
  
  # themeing
  guides(fill = guide_legend(nrow = 1, byrow = FALSE)) +
  theme_classic() +
  theme(legend.position = "bottom", legend.text = element_text(size = base_size * 0.75))

#### Create ggplot ####

rel_imp_max <- max(rel_importance_df$mean)
mar_means_max <- max(marginal_means$Mean)

gg_part <- gg_part <- purrr::map(c(Aboveground = "Aboveground", Total = "Total"), function(part_i) {
  
  mean_temp <- dplyr::filter(marginal_means, part == part_i, measure == "Portfolio effect") 
    
  importance_temp <- dplyr::filter(rel_importance_df, part == part_i) |>
    dplyr::filter(beta != "Residuals")
  
  lm_temp <- dplyr::filter(lm_summary_df, part == part_i) |> 
    dplyr::group_by(measure) |> 
    dplyr::slice(1) |> 
    dplyr::select(-c(scenario, intercept, `R^2`, df, logLik, AICc, delta, weight)) |> 
    dplyr::mutate(spatialVariation = dplyr::case_when(spatialVariation < 0 ~ "-", spatialVariation > 0 ~ "+"), 
                  connectivity  = dplyr::case_when(connectivity  < 0 ~ "-", connectivity  > 0 ~ "+"),
                  `var:connect` = dplyr::case_when(`var:connect` < 0 ~ "-", `var:connect` > 0 ~ "+")) |> 
    tidyr::pivot_longer(-c(part, measure), names_to = "beta", values_to = "direction") |> 
    tidyr::replace_na(list(direction = "NA")) |> 
    dplyr::mutate(beta = factor(beta, levels = c("spatialVariation", "connectivity", "enrichment", "popSize", 
                                                 "enrich:pop", "enrich:connect", "var:pop", "var:connect"), 
                                labels = c("spatialVariation" = "Spatial variation", "connectivity" = "Connectivity", 
                                           "enrichment" = "Enrichment", "popSize" = "Population size",
                                           "enrich:pop" = "Enrichment:Population", "enrich:connect" = "Enrichment:Connectivity", 
                                           "var:pop" = "Variation:Population", "var:connect" = "Variation:Connectivity")))
  
  if (part_i == "Aboveground") {
    
    lgd_pos <- c(0.3, 0.8)
    label_y <- NULL
    
  } else {
    
    lgd_pos <- "none"
    label_y <- element_blank()
    
  }
    
  gg_importance <- ggplot(data = importance_temp) +
      
    # relative importance bars
    geom_col(aes(x = forcats::fct_rev(beta), y = mean, fill = forcats::fct_rev(measure)), 
             color = "black", position = position_dodge(width = 0.6), width = 0.45) +
    geom_text(data = lm_temp, position = position_dodge(width = 0.75),
              aes(x = forcats::fct_rev(beta), y = rel_imp_max * 1.05, label = direction, 
                  color = forcats::fct_rev(measure)), size = 3.0) +
    
    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    
    # set scales
    scale_fill_manual(name = "", values = color_scale) +
    scale_color_manual(name = "", values = color_scale) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, rel_imp_max * 1.1)) +
    
    coord_flip() +
    
    # labels and themes
    labs(y = "", x = "") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "none", axis.line = element_blank(), 
          axis.text.y = label_y)
    
  gg_means <- ggplot(data = mean_temp) +
      
    # adding mean and error bars
    geom_point(aes(x = pop_n, y = exp(Mean), color = nutrient_input), size = 3.5) + 
    geom_errorbar(aes(x = pop_n, ymin = exp(CI_low), ymax = exp(CI_high), color = nutrient_input), 
                  width = 0, linewidth = 1.5) +
    
    geom_line(aes(x = pop_n, y = exp(Mean), color = nutrient_input, group = nutrient_input), 
              alpha = 0.25, linewidth = 1.5) + 
    
    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    
    # set colors
    scale_color_manual(name = "", values = color_enrich) +
    scale_y_continuous(limits = c(0.0, exp(mar_means_max))) +
    guides(color = guide_legend(nrow = 1)) +
      
    # change theme
    labs(y = "", x = "") +
    theme_classic(base_size = base_size) +
    theme(legend.position = lgd_pos, legend.text = element_text(size = base_size * 0.5),  
          axis.line = element_blank())
    
  list(imp = gg_importance, means = gg_means)
    
})
  
gg_imp <- cowplot::plot_grid(gg_part$Aboveground$imp, gg_part$Total$imp, 
                             ncol = 2, labels = c("a)", "b)"), rel_widths = c(0.575, 0.425))

gg_imp <- cowplot::plot_grid(gg_imp, cowplot::get_legend(gg_dummy), nrow = 2, ncol = 1,
                             rel_heights = c(0.9, 0.1))

gg_means <- cowplot::plot_grid(gg_part$Aboveground$means, gg_part$Total$means, ncol = 2) |>
  cowplot::ggdraw(xlim = c(-0.025, 1.0), ylim = c(-0.025, 1.0)) +
  cowplot::draw_label(label = "Population size", x = 0.5, y = 0.05, angle = 0, size = base_size) +
  cowplot::draw_label(expression(paste("Portfolio effect ", italic(cv), italic(beta))),
                      x = -0.01, y = 0.5, angle = 90, size = base_size)

gg_final <- cowplot::plot_grid(gg_imp, gg_means, nrow = 2, rel_heights = c(0.55, 0.45))

#### Save table and ggplot #### 

overwrite <- FALSE

write.table(x = lm_summary_df, file = "04_Figures/Table-Q2.csv", sep  = ";", dec = ".", row.names = FALSE)

suppoRt::save_ggplot(plot = gg_final, filename = "Figure-3.pdf",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_final, filename = "Figure-3.jpg",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = overwrite)

