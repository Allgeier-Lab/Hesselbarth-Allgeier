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

results_pe_pp_df <- dplyr::left_join(x = dplyr::filter(results_final_df, measure != "synchrony", treatment == "combined", 
                                                       include == "yes", part != "Belowground") |> 
                                       dplyr::select(scenario, measure, row_id, part, pop_n, nutrient_input, biotic, abiotic, value.cv),
                                     y = dplyr::filter(results_final_df, measure == "gamma", treatment == "combined", 
                                                       include == "yes", part != "Belowground") |> 
                                       dplyr::select(scenario, row_id, part, pop_n, nutrient_input, biotic, abiotic, value.prod), 
                                     by = c("scenario", "row_id", "part", "pop_n", "nutrient_input", "biotic", "abiotic"))

#### Split data #### 

results_final_list <-  dplyr::filter(results_pe_pp_df, scenario  == "noise") |>
  dplyr::group_by(part, measure) |>
  dplyr::group_split()

names_list <- purrr::map(results_final_list, function(i) as.character(c(unique(i$scenario), unique(i$part), unique(i$measure))))

#### Fit regression model and dredge ####

best_lm_list <- vector(mode = "list", length = length(results_final_list))

lm_summary_list <- vector(mode = "list", length = length(results_final_list))

for(i in 1:length(results_final_list)) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  df_temp <- results_final_list[[i]]
  
  if (unique(df_temp$part) == "Aboveground") {
    
    df_temp$value.prod <- sqrt(df_temp$value.prod)
    
  } else if (unique(df_temp$part) == "Total") {
    
    df_temp$value.prod <- log(df_temp$value.prod)
    
  }
  
  df_temp <- dplyr::mutate(df_temp, value.cv = (value.cv - mean(value.cv)) / sd(value.cv),
                           biotic = (biotic - mean(biotic)) / sd(biotic),
                           abiotic = (abiotic - mean(abiotic)) / sd(abiotic))
  
  lm_temp <- lm(value.prod ~ value.cv + abiotic + biotic + nutrient_input + pop_n,
                data = df_temp, na.action = "na.fail")
  
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
                "connectivity" = "biotic", "enrichment" = "nutrient_input", 
                "popSize" = "pop_n", "portfolioEffect" = "value.cv") |> 
  dplyr::mutate_if(is.numeric, round, digits = 5)

#### Rel imp ####

rel_importance_df <- purrr::map_dfr(seq_along(best_lm_list), function(i) {
  
  message("> Progress: ", i, "/", length(best_lm_list))
  
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
                beta = factor(beta, levels = c("abiotic", "biotic", "nutrient_input", "pop_n", "value.cv",
                                               "residual"),
                              labels = c("abiotic" = "Spatial variation", "biotic" = "Connectivity",
                                         "nutrient_input" = "Enrichment", "pop_n" = "Population size",
                                         "value.cv" = "Stability", "residual" = "Residuals")))

#### Setup ggplot ####

base_size <- 12.0

#### Create ggplot ####

gg_part <- purrr::map(c(Aboveground = "Aboveground", Total = "Total"), function(part_i) {
  
  df_temp <- dplyr::filter(results_pe_pp_df, part == part_i, measure == "gamma") |> 
    dplyr::mutate(pe_class = cut(value.cv, breaks = 50)) |> 
    dplyr::group_by(pe_class) |> 
    dplyr::summarise(value.prod = mean(value.prod), n = dplyr::n())
  
  rel_imp_temp <- dplyr::filter(rel_importance_df, part == part_i, 
                                measure == "Meta-ecosystem scale", beta != "Residuals")
  
  lm_temp <- dplyr::filter(lm_summary_df, part == part_i, measure == "gamma") |> 
    dplyr::slice(1) |> 
    dplyr::select(part, spatialVariation, connectivity, enrichment, popSize, portfolioEffect) |> 
    dplyr::mutate(spatialVariation = dplyr::case_when(spatialVariation < 0 ~ "-",
                                                      spatialVariation > 0 ~ "+"), 
                  connectivity  = dplyr::case_when(connectivity  < 0 ~ "-", 
                                                   connectivity  > 0 ~ "+"),
                  portfolioEffect = dplyr::case_when(portfolioEffect < 0 ~ "-", 
                                                     portfolioEffect > 0 ~ "+")) |> 
    tidyr::pivot_longer(-part, names_to = "beta", values_to = "direction") |> 
    dplyr::filter(!is.na(direction)) |> 
    dplyr::mutate(beta = factor(beta, levels = c("spatialVariation", "connectivity", 
                                                 "enrichment", "popSize", "portfolioEffect"),
                                labels = c("spatialVariation" = "Spatial variation", 
                                           "connectivity" = "Connectivity",
                                           "enrichment" = "Enrichment", "popSize" = "Population size",
                                           "portfolioEffect" = "Stability")))
  
  if (part_i == "Aboveground") {
    
    label_temp <- "a)"
    axis_text <- element_blank()
    
  } else if (part_i == "Total") {
    
    label_temp <- "b)"
    axis_text <- NULL
    
  }
  
  y_pos <- max(df_temp$n) * 1
  
  gg_pp <- ggplot(data = df_temp) +
    
    # adding tiles with PP
    geom_col(aes(x = pe_class, y = n, fill = value.prod)) +
    
    # adding a) b) labels
    annotate("text", x = 1.5, y = y_pos, label = label_temp) +
    
    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    
    # set scales 
    # scale_fill_gradient2(name = "Primary production", low =  "white", mid = "#A5C68C", high = "gare") +
    scale_fill_viridis_c(name = "Primary production", direction = 1) +
    scale_x_discrete(breaks = df_temp$pe_class[seq(1, length(df_temp$pe_class), by = 8)]) +
    
    # theming
    theme_classic(base_size = base_size) + 
    theme(legend.position = "right", axis.line = element_blank(), axis.title = element_blank(), 
          axis.text.x = axis_text)
  
  gg_inset <- ggplot(data = rel_imp_temp, aes(x = beta, y = mean)) + 
    geom_col(fill = NA, color = "black", width = 0.5) +
    geom_text(data = lm_temp, aes(x = beta, y = max(rel_imp_temp$mean) * 1.1, label = direction),
              size = 3.5) +
    scale_fill_discrete(name = "") +
    scale_y_continuous(limits = c(0, 100)) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_classic(base_size = base_size * 0.5) + 
    labs(y = "Rel. importance [%]") +
    theme(legend.position = "none", axis.title.x = element_blank())
  
  ggdraw(gg_pp) + 
    draw_plot(gg_inset, x = 0.4, y = 0.35, width = 0.35, height = 0.5)
  
})

gg_pe_pp <- cowplot::plot_grid(gg_part$Aboveground, gg_part$Total, nrow = 2) |> 
  cowplot::ggdraw(xlim = c(-0.025, 1.0), ylim = c(-0.025, 1.0)) +
  cowplot::draw_label("Count", x = -0.015, y = 0.5, angle = 90, size = base_size) +
  cowplot::draw_label(expression(paste("Coeffiecent of variation ", italic(gamma))), 
                      x = 0.425, y = 0.00, angle = 0, size = base_size)

#### Save ggplots ####

overwrite <- FALSE

suppoRt::save_ggplot(plot = gg_pe_pp, filename = "Figure-4-gamma.pdf",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_pe_pp, filename = "Figure-4-gamma.jpg",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = overwrite)
