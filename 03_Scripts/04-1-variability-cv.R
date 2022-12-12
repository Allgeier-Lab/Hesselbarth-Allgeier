##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create total result figure of variability and cv

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-cv.R")
source("01_Functions/import-mortality.R")

#### Load/wrangle simulated data ####

results_phase_df <- import_cv(path = "02_Data/result-phase.rds")

results_noise_df <- import_cv(path = "02_Data/result-noise.rds")

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")), 
                type = dplyr::case_when(pop_n == 0 ~ "a", abiotic == 0.0 ~ "b", TRUE ~ "ab"), 
                type = factor(type, levels = c("a", "ab", "b"), 
                              labels = c("Abiotic subsidies only", "Abiotic and connectivity", 
                                         "Consumer connectivity only")))

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

results_combined_df <- dplyr::left_join(x = results_combined_df, y = mortality_combined_df, 
                                        by = c("scenario", "row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(pop_n == 128 & nutrient_input == "low" &
                                             died_total > threshold_mort ~ "no", 
                                           TRUE ~ "yes"))

#### Split data #### 

results_combined_list <- dplyr::filter(results_combined_df, measure %in% c("alpha", "gamma"),
                                       part %in% c("ag_production", "bg_production", "ttl_production"), 
                                       pop_n != 0, biotic != 0.0, abiotic != 0.0, include == "yes") |>
  dplyr::group_by(scenario, part, measure, pop_n, nutrient_input) |>
  dplyr::group_split()

#### Fit regression model ####

full_lm_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
  
  message("\r> Progress: scenario=", unique(df_temp$scenario), "; part=", unique(df_temp$part), 
          "; measure=", unique(df_temp$measure), "; pop_n=", unique(df_temp$pop_n), 
          "; nutrient_input=", unique(df_temp$nutrient_input), appendLF = FALSE)
  
  df_temp_stand <- dplyr::select(df_temp, value.cv, abiotic, biotic) # |>
  # dplyr::mutate(value.cv = log(value.cv))
  # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x)))
  # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
  
  lm_temp <- lm(value.cv ~ abiotic * biotic, data = df_temp_stand) 
  
  broom::tidy(lm_temp) |> 
    dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                  measure = unique(df_temp$measure), pop_n = unique(df_temp$pop_n), 
                  nutrient_input = unique(df_temp$nutrient_input), .before = term) |> 
    dplyr::mutate(p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                                   p.value < 0.05 ~ "*", p.value >= 0.05 ~ " "), 
                  r2 = summary(lm_temp)$adj.r.squared)}) |> 
  dplyr::mutate(term = factor(term), p.value.class = factor(p.value.class, levels = c("*", "**", "***", " ")))

#### Model assumptions ####
# 
# assumption_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
# 
#   message("\r> Progress: scenario=", unique(df_temp$scenario), "; part=", unique(df_temp$part),
#           "; measure=", unique(df_temp$measure), "; pop_n=", unique(df_temp$pop_n),
#           "; nutrient_input=", unique(df_temp$nutrient_input), appendLF = FALSE)
# 
#   df_temp_stand <- dplyr::select(df_temp, value.cv, abiotic, biotic) # |>
#   # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) |>
#   # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
# 
#   lm_temp <- lm(value.cv ~ abiotic * biotic, data = df_temp_stand)
# 
#   residuals_temp <- resid(lm_temp)
# 
#   tibble::tibble(scenario = unique(df_temp$scenario), part = unique(df_temp$part),
#                  measure = unique(df_temp$measure), pop_n = unique(df_temp$pop_n),
#                  nutrient_input = unique(df_temp$nutrient_input),
#                  shapiro =  shapiro.test(residuals_temp)[[2]], lillie = nortest::lillie.test(residuals_temp)[[2]])})
# 
# dplyr::filter(assumption_df, shapiro < 0.05 & lillie < 0.05)

#### Relative importance R2 ####

rel_importance_df <- purrr::map_dfr(results_combined_list, function(df_temp) {

  message("\r> Progress: scenario=", unique(df_temp$scenario), "; part=", unique(df_temp$part),
          "; measure=", unique(df_temp$measure), "; pop_n=", unique(df_temp$pop_n),
          "; nutrient_input=", unique(df_temp$nutrient_input), appendLF = FALSE)

  df_temp_stand <- dplyr::select(df_temp, value.cv, abiotic, biotic) # |>
  # dplyr::mutate(value.cv = log(value.cv))
  # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x)))
  # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))

  lm_temp <- lm(value.cv ~ abiotic * biotic, data = df_temp_stand)

  rel_r2 <- relaimpo::boot.relimp(lm_temp, type = "lmg", b = 100, level = 0.95, fixed = FALSE) |>
    relaimpo::booteval.relimp(bty = "basic")

  tibble::tibble(
    scenario = unique(df_temp$scenario), part = unique(df_temp$part), measure = unique(df_temp$measure),
    pop_n = unique(df_temp$pop_n), nutrient_input = unique(df_temp$nutrient_input),
    beta = c(rel_r2@namen[-1], "residual"), mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)),
    lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA), r2 = summary(lm_temp)$adj.r.squared) |>
    dplyr::mutate(dplyr::across(mean:higher, ~ .x * 100))}) |>
  dplyr::mutate(beta = factor(beta))

#### Setup ggplot ####

size_base <- 10.0
size_text <- 2.0
size_point <- 2.5

width_doge <- 0.5

color_regression <- c("abiotic" = "#007aa1", "biotic" = "#41b282", "abiotic:biotic" = "#fcb252")
color_relimp <- c("abiotic" = "#007aa1", "biotic" = "#41b282", "abiotic:biotic" = "#fcb252", "residual" = "grey")
color_beta <- c("0" = "grey", "8" = "#0D0887FF", "16" = "#7E03A8FF", "32" = "#CC4678FF", "64" = "#F89441FF", "128" = "#F0F921FF")

gg_dummy_upper <- data.frame(beta = c("abiotic", "biotic", "abiotic:biotic", "residual"),
                              mean = c(1, 1, 1, 1)) |>
  dplyr::mutate(beta = factor(beta, levels = c("abiotic", "biotic", "abiotic:biotic", "residual"))) |>
  ggplot() +
  geom_col(aes(x = beta, y = mean, fill = beta)) +
  scale_fill_manual(name = "", values = color_relimp,
                    labels = c("abiotic" = "Abiotic subsidies", "biotic" = "Consumer connectivity",
                               "abiotic:biotic" = "Subsidies:Connectivity", "residual" = "Residuals")) +
  theme_classic(base_size = size_base) +
  theme(legend.position = "bottom")

#### Create ggplot cv ####

gg_scenario <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  purrr::map(c(ag = "ag_production", bg =  "bg_production", ttl = "ttl_production"), function(part_i) {
    
    regression_temp <- dplyr::filter(full_lm_df, scenario == scenario_i, part == part_i, 
                                     term != "(Intercept)")
    
    importance_temp <- dplyr::filter(rel_importance_df, scenario == scenario_i, part == part_i)
      
    beta_temp <- dplyr::filter(results_combined_df, scenario == scenario_i, part == part_i, 
                               measure == "beta", include == "yes")
    
    label_temp <- dplyr::group_by(beta_temp, nutrient_input, type) |> 
      dplyr::group_split() |> 
      purrr::map_dfr(function(i) {
        
        if (any(i$pop_n != 0)) {
        
        anova <- aov(value.cv ~ pop_n, data = i)
        
        tukey <- TukeyHSD(anova)
        
        label_letters <- multcompView::multcompLetters4(anova, tukey)
        
        data.frame(nutrient_input = unique(i$nutrient_input), type = unique(i$type),
                   pop_n = names(label_letters$pop_n$Letters),
                   letter = unname(label_letters$pop_n$Letters))
        
        } else {
        
            data.frame(nutrient_input = unique(i$nutrient_input), type = unique(i$type),
                       pop_n = unique(i$pop_n), letter = "")
        
        }
      }) |> 
      dplyr::mutate(pop_n = factor(as.numeric(pop_n), ordered = TRUE)) |> 
      dplyr::arrange(nutrient_input, type, pop_n)
    
    y_lims <- range(regression_temp$estimate)

    gg_upper <- purrr::map(c(alpha = "alpha", beta = "gamma"), function(scale_i) {
      
      gg_regression <- ggplot(data = dplyr::filter(regression_temp, measure == scale_i)) +
        
        # zero line
        geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
        
        # lines
        geom_line(aes(x = pop_n, y = estimate, group = term, color = term),
                  alpha = 0.25, position = position_dodge(width = width_doge)) +
        
        # points
        geom_point(aes(x = pop_n, y = estimate, color = term, shape = p.value.class),
                   shape = 19, position = position_dodge(width = width_doge), size = size_point) +
        
        # text
        geom_text(aes(x = pop_n, y = estimate, label = p.value.class, group = term), color = "black",
                  size = size_text, position = position_dodge(width = width_doge), vjust = 0.75) +
        
        # facet
        facet_grid(rows = vars(nutrient_input)) +
        
        # set scales
        scale_color_manual(name = "Scale", values = color_regression) +
        scale_y_continuous(limits = y_lims) +
        
        # labels and themes
        theme_classic(base_size = size_base) +
        theme(strip.background = element_blank(), strip.text = element_blank(),
              axis.title = element_blank(), axis.line = element_blank(), 
              panel.border = element_rect(linewidth = 0.25, fill = NA),
              legend.position = "none")
      
      gg_relimp <- ggplot(data = dplyr::filter(importance_temp, measure == scale_i)) +
        
        # relative importance bars
        geom_col(aes(x = pop_n, y = mean, fill = factor(beta, levels = rev(levels(beta))))) +
        
        # facets
        facet_grid(rows = vars(nutrient_input)) +
        
        # set scales
        scale_fill_manual(name = "", values = color_relimp) +
        scale_y_continuous(labels = function(x) paste0(x, "%")) +
        
        # labels and themes
        theme_classic(base_size = size_base) +
        theme(strip.background = element_blank(), strip.text = element_blank(),
              axis.title = element_blank(), axis.line = element_blank(), 
              panel.border = element_rect(linewidth = 0.25, fill = NA),
              legend.position = "none")
      
      cowplot::plot_grid(gg_regression, gg_relimp, ncol = 2)
      
    })
    
    gg_upper <- cowplot::plot_grid(plotlist = gg_upper, ncol = 2, labels = c("a)", "b)"), 
                                   hjust = -0.0, vjust = 0.25)
  
    gg_upper <- cowplot::ggdraw(gg_upper, xlim = c(-0.025, 1.025), ylim = c(-0.025,  1.05)) +
      cowplot::draw_label("Population size", x = 0.5, y = 0, angle = 0, size = size_base) +
      cowplot::draw_label("Parameter estimate / Relative importance [%]", x = 0.0, y = 0.5, angle = 90, size = size_base) +
      cowplot::draw_label("Nutr. input: low", x = 1.0, y = 0.85, angle = 270, size = size_base * 0.65) +
      cowplot::draw_label("Nutr. input: medium", x = 1.0, y = 0.525, angle = 270, size = size_base * 0.65) +
      cowplot::draw_label("Nutr. input: high", x = 1.0, y = 0.225, angle = 270, size = size_base * 0.65)
    
    gg_upper <- cowplot::plot_grid(gg_upper, cowplot::get_legend(gg_dummy_upper), 
                                   nrow = 2, ncol = 1, rel_heights = c(0.9, 0.1))
    
    # ggplot
    gg_beta <- ggplot(beta_temp, aes(x = type, y = value.cv, fill = pop_n)) + 
      
      # geoms
      geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
      geom_hline(yintercept = 1.0, linetype = 2, color = "grey") +
      geom_text(data = label_temp, aes(x = type, y = max(beta_temp$value.cv), group = pop_n, label = letter), 
                 position = position_dodge(0.8), size = size_text) +
      
      # scales
      scale_fill_manual(name = "Population size", values = color_beta) +
      scale_color_manual(name = "", values = color_beta) +
      coord_cartesian(ylim = c(1.0, max(beta_temp$value.cv))) +
      
      # facet
      facet_grid(rows = vars(nutrient_input),
                 labeller = labeller(nutrient_input = function(x) paste("Nutr. input:", x))) +
      
      # themes
      labs(x = "", y = "Portfolio effect") +
      guides(fill = guide_legend(nrow = 1), color = "none") +
      theme_classic(base_size = 10) +
      theme(strip.background = element_blank(), strip.text = element_text(size = size_base * 0.65),
            axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA), 
            legend.position = "bottom")
    
   cowplot::plot_grid(gg_upper, gg_beta, nrow = 2, rel_heights = c(0.45, 0.55),
                      labels = c("", "c)"), hjust = -0.25, vjust = 0.5)
  
  })
})

#### Save ggplot #### 

overwrite <- FALSE

suppoRt::save_ggplot(plot = gg_scenario$noise$ag, filename = "Figure-2.pdf",
                     path = "04_Figures/", width = width, height = height * 0.75,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_scenario$noise$bg, filename = "Figure-A4.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.75,
                     units = units, dpi = dpi, overwrite = overwrite)
