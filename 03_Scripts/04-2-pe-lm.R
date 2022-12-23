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

#### Split data #### 

results_final_list <- dplyr::filter(results_final_df, measure == "beta", include == "yes",
                                    part %in% c("Aboveground", "Belowground"), treatment == "combined") |>
  dplyr::group_by(scenario, part) |>
  dplyr::group_split()

names_list <- purrr::map(results_final_list, function(i) c(unique(i$scenario), unique(i$part)))

#### Fit regression model and dredge ####

best_lm_list <- purrr::map(seq_along(results_final_list), function(i) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  df_temp <- results_final_list[[i]]
  
  lm_temp <- lm(value.cv ~ nutrient_input + abiotic + pop_n + biotic + 
                  nutrient_input*pop_n + nutrient_input*biotic + pop_n*abiotic + 
                  abiotic*biotic, data = df_temp, na.action = "na.fail")
  
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

rel_importance_df <- purrr::map_dfr(seq_along(results_final_list), function(i) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  rel_r2 <- relaimpo::boot.relimp(best_lm_list[[i]], type = "lmg", b = 1000, level = 0.95, 
                                  fixed = FALSE) |>
    relaimpo::booteval.relimp(bty = "basic")
  
  tibble::tibble(
    scenario = names_list[[i]][[1]], part = names_list[[i]][[2]],
    beta = c(rel_r2@namen[-1], "residual"), mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)),
    lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA)) |>
    dplyr::mutate(dplyr::across(mean:higher, ~ .x * 100), 
                  lower = dplyr::case_when(lower < 0 ~ 0.0, TRUE ~ lower))}) |> 
  dplyr::mutate(beta = factor(beta, levels = c("abiotic", "biotic", "nutrient_input", "pop_n", 
                                               "nutrient_input:pop_n", "abiotic:biotic", "abiotic:pop_n",
                                               "biotic:nutrient_input", "residual"), 
                              labels = c("abiotic" = "Variation subsidies", "biotic" = "Connectivity", 
                                         "nutrient_input" = "Amount subsidies", "pop_n" = "Population size", 
                                         "nutrient_input:pop_n" = "Amount:Population", "abiotic:biotic" = "Variation:Connectivity", 
                                         "abiotic:pop_n" = "Variation:Population", "biotic:nutrient_input" = "Connectivity:Amount", 
                                         "residual" = "Residuals")))

#### Create ggplot #### 

base_size <- 10.0

w <- 0.5

color_beta <- c("Variation subsidies" = "#e76254", "Connectivity" = "#ef8a47", 
                "Amount subsidies" = "#f7aa58", "Population size" = "#ffd06f", 
                "Amount:Population" = "#ffe6b7", "Variation:Connectivity" = "#aadce0", 
               "Variation:Population" = "#72bcd5", "Connectivity:Amount" = "#376795", 
                "Residuals" = "grey")

gg_relimp <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
    importance_temp <- dplyr::filter(rel_importance_df, scenario == scenario_i) |>
      dplyr::mutate(label = dplyr::case_when(mean > 10.0 & beta != "Rexsiduals" ~ as.character(beta), TRUE ~ ""))
      # dplyr::arrange(part, desc(beta)) |> 
      # dplyr::group_by(part) |>
      # dplyr::mutate(y_label = cumsum(mean))
    
    max_imp <- ceiling(max(importance_temp$higher) / 10) * 10
    
    summary_temp <- dplyr::filter(lm_summary_df, scenario == scenario_i)
    
    r2_temp <- dplyr::group_by(summary_temp, part) |> 
      dplyr::summarise(r2 = unique(r2)) |> 
      dplyr::mutate(r2 = round(r2, 3), x = length(unique(importance_temp$beta)) - 1, 
                    y = c(max_imp, max_imp * 0.95), label = paste0(part, " R^2=", r2))
    
    # ggplot(data = importance_temp) +
    #   
    #   # zero line
    #   geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
    #   
    #   # relative importance bars
    #   geom_col(aes(x = beta, y = mean, fill = part), position = position_dodge(w), 
    #            width = w * 0.95, alpha = 0.5) +
    # 
    #   geom_errorbar(aes(x = beta, ymin = lower, ymax = higher, color = part),
    #                 position = position_dodge(w), width = w * 0.5) +
    #   
    #   geom_text(data = r2_temp, aes(x = x, y = y, label = label, color = part), 
    #             size = 2.5) +
    # 
    #   # facets
    #   annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    #   annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    #   annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    #   annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    #   
    #   # set scales
    #   scale_color_manual(name = "", values = c("Aboveground" = "#32b2da", "Belowground" = "#a62f00")) +
    #   scale_fill_manual(name = "", values = c("Aboveground" = "#32b2da", "Belowground" = "#a62f00")) +
    #   scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, max_imp), 
    #                      breaks = seq(0, max_imp, 5)) +
    # 
    #   # labels and themes
    #   labs(y = "Relative importance [%]") +
    #   theme_classic(base_size = base_size) +
    #   theme(axis.title.x = element_blank(), axis.line = element_blank(), 
    #         axis.text.x = element_text(angle = 45, hjust = 1),
    #         legend.position = "none")
    
    ggplot(data = importance_temp) +
      
      # zero line
      geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
      
      # relative importance bars
      geom_col(aes(x = part, y = mean, fill = forcats::fct_rev(beta)), position = "stack", 
               width = w * 0.95) +
      geom_text(aes(x = part, y = mean, label = label, group = forcats::fct_rev(beta)),
                position = position_stack(vjust = 0.5), size = 2) +

      # facets
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set scales
      scale_fill_manual(name = "", values = color_beta) +


      scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100),
                         breaks = seq(0, 100, 25)) +
      
      # labels and themes
      labs(y = "Relative importance [%]") +
      theme_classic(base_size = base_size) +
      theme(axis.title.x = element_blank(), axis.line = element_blank(), 
            # axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right")
    
  })

#### Save table and ggplot #### 

dplyr::filter(lm_summary_df, scenario == "noise") |> 
  dplyr::select(-scenario, -r2) |> 
  dplyr::mutate_if(is.numeric, round, digits = 3) |> 
  readr::write_csv(file = "04_Figures/Table-1.csv")

suppoRt::save_ggplot(plot = gg_relimp$noise, filename = "Figure-3.pdf",
                     path = "04_Figures/", width = width * 0.75, height = height * 1/3,
                     units = units, dpi = dpi, overwrite = FALSE)
