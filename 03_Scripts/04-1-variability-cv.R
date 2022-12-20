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

near <- FALSE

results_phase_df <- import_cv(path = "02_Data/result-phase.rds", near = near)

results_noise_df <- import_cv(path = "02_Data/result-noise.rds", near = near)

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
  dplyr::mutate(include = dplyr::case_when(died_total > threshold_mort ~ "no", TRUE ~ "yes"))

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
  
  lm_temp <- lm(value.cv ~ abiotic * biotic, data = df_temp) 
  
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
#   lm_temp <- lm(value.cv ~ abiotic * biotic, data = df_temp)
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

  lm_temp <- lm(value.cv ~ abiotic * biotic, data = df_temp)

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
# color_beta <- c("0" = "grey", "8" = "#0D0887FF", "16" = "#7E03A8FF", "32" = "#CC4678FF", "64" = "#F89441FF", "128" = "#F0F921FF")
color_beta <- c("Abiotic subsidies only" = "#007aa1", "Consumer connectivity only" = "#41b282", 
                "Abiotic and connectivity" = "#fcb252")

gg_dummy <- data.frame(beta = c("abiotic", "biotic", "abiotic:biotic", "residual"),
                           mean = c(1, 1, 1, 1)) |>
  dplyr::mutate(beta = factor(beta, levels = c("abiotic", "biotic", "abiotic:biotic", "residual"))) |>
  ggplot() +
  geom_col(aes(x = beta, y = mean, fill = beta)) +
  scale_fill_manual(name = "", values = color_relimp,
                    labels = c("abiotic" = "Abiotic subsidies", "biotic" = "Consumer connectivity",
                               "abiotic:biotic" = "Subsidies:Connectivity", "residual" = "Residuals")) +
  theme_classic(base_size = size_base) +
  theme(legend.position = "bottom")

#### Create ggplot: Main ####

gg_alpha <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  purrr::map(c(ag = "ag_production", bg =  "bg_production", ttl = "ttl_production"), function(part_i) {
    
    regression_temp <- dplyr::filter(full_lm_df, scenario == scenario_i, part == part_i, 
                                     term != "(Intercept)", measure == "alpha")
    
    importance_temp <- dplyr::filter(rel_importance_df, scenario == scenario_i, part == part_i, 
                                     measure == "alpha")
      
    beta_temp <- dplyr::filter(results_combined_df, scenario == scenario_i, part == part_i, 
                               measure == "beta", include == "yes")
    
    gg_regression <- ggplot(data = regression_temp) +
      
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
      facet_grid(rows = vars(nutrient_input), scales = "free") +
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set scales
      scale_color_manual(name = "Scale", values = color_regression) +
      scale_y_continuous(limits = function(x) c(-max(abs(regression_temp$estimate)), 
                                                max(abs(regression_temp$estimate))), 
                         breaks = function(x) seq(-quantile(abs(regression_temp$estimate), 0.99), 
                                                  quantile(abs(regression_temp$estimate), 0.99), length.out = 4),
                         labels = function(x) round(x, 2)) +

      # labels and themes
      labs(y = "Parameter estimate") +
      theme_classic(base_size = size_base) +
      theme(strip.background = element_blank(), strip.text = element_blank(),
            axis.line = element_blank(), axis.title.x = element_blank(), legend.position = "none")
      
    gg_relimp <- ggplot(data = importance_temp) +
      
      # zero line
      geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
      
      # relative importance bars
      geom_col(aes(x = pop_n, y = mean, fill = factor(beta, levels = rev(levels(beta))))) +
      
      # facets
      facet_grid(rows = vars(nutrient_input), scales = "fixed") +
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set scales
      scale_fill_manual(name = "", values = color_relimp) +
      scale_y_continuous(labels = function(x) paste0(x, "%")) +
      
      # labels and themes
      labs(y = "Relative importance [%]") +
      theme_classic(base_size = size_base) +
      theme(strip.background = element_blank(), strip.text = element_blank(),
            axis.title.x = element_blank(), axis.line = element_blank(), legend.position = "none")
    
    # ggplot
    gg_beta <- ggplot(beta_temp, aes(x = pop_n, y = value.cv, color = type, fill = type)) +

      # geoms
      geom_boxplot(outlier.shape = NA, position = position_dodge2(width = 1, preserve = "single"),
                   alpha = 0.5) +
      geom_hline(yintercept = 1.0, linetype = 2, color = "grey") +

      # scales
      scale_fill_manual(name = "", values = color_beta) +
      scale_color_manual(name = "", values = color_beta) +
      scale_y_continuous(sec.axis = dup_axis()) +

      # facet
      facet_grid(rows = vars(nutrient_input), cols = vars(type), scales = "free_x", space = "free",
                 labeller = labeller(nutrient_input = function(x) paste("Nutr. input:", x))) +
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +

      # themes
      labs(y = "Portfolio effect") +
      guides(fill = guide_legend(nrow = 1), color = "none") +
      theme_classic(base_size = size_base) +
      theme(strip.background = element_blank(), strip.text = element_blank(),
            panel.spacing.x = unit(0.0, "mm"), axis.title.x = element_blank(),
            axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
            axis.title.y.right = element_blank(), axis.line.y = element_line(linewidth = 0.5, color = "grey"), 
            axis.line.x = element_blank(), legend.position = "none")
 
   gg_total <- cowplot::plot_grid(gg_regression, gg_relimp, gg_beta, 
                                  ncol = 3, rel_widths = c(1/3, 1/3, 1/3),
                                  labels = c("a)", "b)", "c)"))
   
   gg_total <- cowplot::ggdraw(gg_total, xlim = c(0, 1.025), ylim = c(-0.025,  1.0)) +
     cowplot::draw_label("Population size", x = 0.55, y = 0, angle = 0, size = size_base * 0.65) + 
     cowplot::draw_label("Nutr. input: low", x = 1.0, y = 0.8, angle = 270, size = size_base * 0.65, hjust = 1) +
     cowplot::draw_label("Nutr. input: medium", x = 1.0, y = 1/2, angle = 270, size = size_base * 0.65) +
     cowplot::draw_label("Nutr. input: high", x = 1.0, y = 0.15, angle = 270, size = size_base * 0.65, hjust = 1)
   
   gg_total <- cowplot::plot_grid(gg_total, cowplot::get_legend(gg_dummy),
                                  nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))
   
   gg_total
   
  })
})

#### Create ggplot: Appendix ####

gg_gamma <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  purrr::map(c(ag = "ag_production", bg =  "bg_production", ttl = "ttl_production"), function(part_i) {
    
    regression_temp <- dplyr::filter(full_lm_df, scenario == scenario_i, part == part_i, 
                                     term != "(Intercept)", measure == "gamma")
    
    importance_temp <- dplyr::filter(rel_importance_df, scenario == scenario_i, part == part_i, 
                                     measure == "gamma")
    
    beta_temp <- dplyr::filter(results_combined_df, scenario == scenario_i, part == part_i, 
                               measure == "beta", include == "yes", measure == "gamma")
    
    gg_regression <- ggplot(data = regression_temp) +
      
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
      facet_grid(rows = vars(nutrient_input), scales = "free") +
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set scales
      scale_color_manual(name = "Scale", values = color_regression) +
      scale_y_continuous(limits = function(x) c(-max(abs(regression_temp$estimate)), 
                                                max(abs(regression_temp$estimate))), 
                         breaks = function(x) seq(-quantile(abs(regression_temp$estimate), 0.99), 
                                                  quantile(abs(regression_temp$estimate), 0.99), length.out = 4),
                         labels = function(x) round(x, 2)) +
      
      # labels and themes
      labs(y = "Parameter estimate") +
      theme_classic(base_size = size_base) +
      theme(strip.background = element_blank(), strip.text = element_blank(),
            axis.title.x = element_blank(), axis.line = element_blank(), legend.position = "none")
    
    gg_relimp <- ggplot(data = importance_temp) +
      
      # zero line
      geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
      
      # relative importance bars
      geom_col(aes(x = pop_n, y = mean, fill = factor(beta, levels = rev(levels(beta))))) +
      
      # facets
      facet_grid(rows = vars(nutrient_input), scales = "fixed") +
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set scales
      scale_fill_manual(name = "", values = color_relimp) +
      scale_y_continuous(labels = function(x) paste0(x, "%")) +
      
      # labels and themes
      labs(y = "Relative importance [%]") +
      theme_classic(base_size = size_base) +
      theme(strip.background = element_blank(), strip.text = element_blank(),
            axis.line = element_blank(), legend.position = "none")
    
    gg_total <- cowplot::plot_grid(gg_regression, gg_relimp, 
                                   ncol = 2, rel_widths = c(0.5, 0.5),
                                   labels = c("a)", "b)"))
    
    gg_total <- cowplot::ggdraw(gg_total, xlim = c(0, 1.025), ylim = c(-0.025,  1.0)) +
      cowplot::draw_label("Population size", x = 0.55, y = 0, angle = 0, size = size_base) +
      cowplot::draw_label("Nutr. input: low", x = 1.0, y = 0.8, angle = 270, size = size_base * 0.65, hjust = 1) +
      cowplot::draw_label("Nutr. input: medium", x = 1.0, y = 1/2, angle = 270, size = size_base * 0.65) +
      cowplot::draw_label("Nutr. input: high", x = 1.0, y = 0.15, angle = 270, size = size_base * 0.65, hjust = 1)
    
    gg_total <- cowplot::plot_grid(gg_total, cowplot::get_legend(gg_dummy),
                                   nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))
    
    gg_total
    
  })
})

#### Save ggplot #### 

overwrite <- FALSE

suppoRt::save_ggplot(plot = gg_alpha$noise$ag, 
                     filename = ifelse(near, yes = "Figure-2-near.pdf", no = "Figure-2.pdf"),
                     path = "04_Figures/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_gamma$noise$ag, 
                     filename = ifelse(near, yes = "Figure-A4-near.pdf", no = "Figure-A4.pdf"),
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)

# Belowground

suppoRt::save_ggplot(plot = gg_alpha$noise$bg, 
                     filename = ifelse(near, yes = "Figure-2-near.pdf", no = "Figure-2.pdf"),
                     path = "04_Figures/Belowground/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_gamma$noise$bg, 
                     filename = ifelse(near, yes = "Figure-A4-near.pdf", no = "Figure-A4.pdf"),
                     path = "04_Figures/Belowground//", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)

# suppoRt::save_ggplot(plot = gg_scenario$noise$bg, filename = "Figure-A4.pdf",
#                      path = "04_Figures/Appendix/", width = width, height = height * 0.65,
#                      units = units, dpi = dpi, overwrite = overwrite)
