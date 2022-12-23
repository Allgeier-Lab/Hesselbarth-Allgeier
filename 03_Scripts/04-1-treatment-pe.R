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
                                    part %in% c("Aboveground", "Belowground")) |>
  dplyr::group_by(scenario, part) |>
  dplyr::group_split()

#### Fit regression model ####

anova_df <- purrr::map_dfr(results_final_list, function(df_temp) {
  
  anova <- aov(value.cv ~ treatment, data = df_temp)
  
  tukey <- TukeyHSD(anova)

  # compact letter display
  cld <- multcompView::multcompLetters4(anova, tukey)
  
  data.frame(treatment = names(cld$treatment$Letters), label = as.character(cld$treatment$Letters)) |> 
    dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                  .before = treatment)}) |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")), 
                treatment = factor(treatment, levels = c("subsidies", "connectivity", "combined")))
                
#### Setup ggplot ####

base_size <- 10.0

w <- 0.5

color_treatment <- c("subsidies" = "#007aa1", "connectivity" = "#41b282", "combined" = "#fcb252")

gg_dummy <- data.frame(x = c(1, 2, 3), y = c(0.5, 0.5, 0.5), z = c("subsidies", "connectivity", "combined")) |> 
  dplyr::mutate(z = factor(z, levels = c("subsidies", "connectivity", "combined"))) |> 
  ggplot(aes(x = x, y = y, color = z)) + 
  geom_point() + geom_errorbar(aes(ymin = 0, ymax = 1)) +
  scale_color_manual(name = "", values = color_treatment, 
                     labels = c("subsidies" = "Variation abiotic subsidies", 
                                "connectivity" = "Consumer connectivity", 
                                "combined" = "Variation and connectivity combined")) +
  theme_classic(base_size = base_size) +
  theme(legend.position = "bottom")

#### Create ggplot: Main ####

gg_tukey <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
    anova_temp <- dplyr::filter(anova_df, scenario == scenario_i)
    
    pe_temp <- dplyr::filter(results_final_df, scenario == scenario_i, measure == "beta", 
                             part != "Total", include == "yes")
    
    pe_sum <- dplyr::group_by(pe_temp, part, treatment) |> 
      dplyr::summarise(mean = mean(value.cv), lower = mean - sd(value.cv), upper = mean + sd(value.cv), 
                       .groups = "drop")
    
    cv_temp <- dplyr::filter(results_final_df, scenario == scenario_i, measure %in% c("alpha", "gamma"),
                             part != "Total", include == "yes") |> 
      dplyr::select(row_id, part, treatment, measure, value.cv) |> 
      tidyr::pivot_wider(names_from = measure, values_from = value.cv)
    
    gg_pe <- ggplot(data = pe_sum, aes(x = part, color = treatment, fill = treatment)) + 
      
      # adding hline at 1
      geom_hline(yintercept = 1, color = "grey", linetype = 2) +
      
      # adding transparent jitter
      geom_point(data = pe_temp, aes(x = part, y = value.cv, color = treatment, group = treatment), 
                 alpha = 0.1, size = 0.75,
                 position = position_jitterdodge(dodge.width = w, jitter.width = w * 0.75)) +
      
      # adding mean +- sd
      geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(w), 
                    width = w * 0.75, linewidth = 1.0) +
      geom_point(aes(y = mean), position = position_dodge(w), size = 1.5) +
      geom_text(data = anova_temp, position = position_dodge(w),
                aes(x = part, y = max(pe_sum$upper) * 1.1, 
                    label = label, group = treatment, color = treatment)) +

      # change scales
      scale_color_manual(name = "", values = color_treatment) +

      # adding box
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # change theme
      labs(x = "", y = expression(paste("Portfolio effect ", italic(cv[beta])))) +
      theme_classic(base_size = base_size) + 
      theme(legend.position = "none", axis.line = element_blank())
    
    gg_cv <- ggplot(data = cv_temp, aes(x = gamma, y = alpha, color = treatment)) + 
      
      geom_abline(slope = 1, intercept = 0, color = "grey", linetype = 2) +
      
      # adding points
      geom_point(alpha = 0.25) +
      
      # facet wrappung
      facet_wrap(. ~ part, scales = "free", nrow = 2) +
      
      # # change scales
      # coord_equal(ratio = 1) +
      scale_color_manual(name = "", values = color_treatment) + 
      
      # adding box
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # change theme
      labs(x = expression(paste("Meta-ecosystem scale ", italic(cv[gamma]))), 
           y = expression(paste("Local scale ", italic(cv[alpha])))) +
      theme_classic(base_size = base_size) + 
      theme(legend.position = "none", axis.line = element_blank(),
            strip.background = element_blank(), strip.text = element_text(hjust = 0))
    
    cowplot::plot_grid(cowplot::plot_grid(gg_pe, gg_cv, ncol = 2, labels = c("a)", "b)")), 
                       cowplot::get_legend(gg_dummy), nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))

})

#### Save ggplot #### 

suppoRt::save_ggplot(plot = gg_tukey$noise, filename = "Figure-2.pdf",
                     path = "04_Figures/", width = width, height = height * 1/3,
                     units = units, dpi = dpi, overwrite = FALSE)
