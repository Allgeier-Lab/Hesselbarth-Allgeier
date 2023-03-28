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
  
  lm_temp <- lm(value.prod ~ value.cv + abiotic + biotic + nutrient_input + pop_n,
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
  dplyr::select(-logLik, -delta, -weight, -df) |>
  dplyr::rename("Variation" = "abiotic", "Connectivity" = "biotic", 
                "Enrichment" = "nutrient_input", "Population" = "pop_n", "PE" = "value.cv")

#### Setup ggplot ####

base_size <- 12.0

#### Create ggplot ####

gg_pe_pp <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  gg_part <- purrr::map(c(Aboveground = "Aboveground", Total = "Total"), function(part_i) {
    
    df_temp <- dplyr::filter(results_pe_pp_df, scenario == scenario_i, part == part_i) |> 
      dplyr::mutate(pe_class = cut(value.cv, breaks = 50)) |> 
      dplyr::group_by(pe_class) |> 
      dplyr::summarise(value.prod = mean(value.prod), 
                       n = dplyr::n())
    
    lm_temp <- dplyr::filter(lm_summary_df, scenario == scenario_i, part == part_i) |>
      dplyr::slice(1) |> 
      dplyr::select(-c(`(Intercept)`, `R^2`, AICc)) |> 
      # dplyr::mutate(Enrichment = dplyr::case)
      tidyr::pivot_longer(-c(scenario, part, Enrichment, Population)) |> 
      dplyr::mutate(value = tidyr::replace_na(value, 0)) |> 
      dplyr::mutate(name = factor(name, levels = c("PE", "Connectivity", "Variation")))
    
    df_temp$pe_label <- " "
    
    df_temp$pe_label[2] <- "low"
    df_temp$pe_label[floor(nrow(df_temp)/2)] <- "mid"
    df_temp$pe_label[nrow(df_temp) - 1] <- "high"
    
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
      scale_fill_gradient2(name = "Primary production", low =  "white", mid = "#A5C68C", high = "#1E5320") +
      scale_x_discrete(labels =  df_temp$pe_label) +
      
      # theming
      theme_classic(base_size = base_size) + 
      theme(legend.position = "right", axis.line = element_blank(), axis.title = element_blank(), 
            axis.text.x = axis_text)
    
    gg_inset <- ggplot(data = lm_temp, aes(x = name, y = value)) + 
      geom_col(aes(fill = name), color = "black", width = 0.5) + 
      scale_fill_discrete(name = "") +
      theme_classic(base_size = base_size * 0.75) + 
      labs(y = "Effect size") +
      theme(legend.position = "top", axis.text.x = element_blank(), axis.title.x = element_blank())
    
    ggdraw(gg_pp) + 
      draw_plot(gg_inset, x = 0.2, y = 0.65, width = 0.55, height = 0.25)
    
    
  })
  
  cowplot::plot_grid(gg_part$Aboveground, gg_part$Total, nrow = 2) |> 
    cowplot::ggdraw(xlim = c(-0.025, 1.0), ylim = c(-0.025, 1.0)) +
    cowplot::draw_label("Count", x = -0.015, y = 0.5, angle = 90, size = base_size) +
    cowplot::draw_label("Portfolio effect", x = 0.425, y = 0.00, angle = 0, size = base_size)
  
})

#### Save ggplots ####

overwrite <- TRUE

suppoRt::save_ggplot(plot = gg_pe_pp$noise, filename = "Figure-4.pdf",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_pe_pp$noise, filename = "Figure-4.jpg",
                     path = "04_Figures/", width = width, height = height * 1/2,
                     units = units, dpi = dpi, overwrite = overwrite)
