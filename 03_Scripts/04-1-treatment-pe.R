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
  dplyr::mutate(treatment = dplyr::case_when(abiotic != 0.0 & biotic == 0.0 ~ "enrichment", 
                                             abiotic == 0.0 & biotic != 0.0 ~ "connectivity", 
                                             abiotic != 0.0 & biotic != 0.0 ~ "combined")) |>
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")),
                part = factor(part, levels = c("ag_production", "bg_production", "ttl_production"),
                              labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                         "ttl_production" = "Total")),
                treatment = factor(treatment, levels = c("enrichment", "connectivity", "combined"))) |> 
  dplyr::filter(scenario == "noise")

#### Split data #### 

results_final_list <- dplyr::filter(results_final_df, measure != "synchrony", 
                                    include == "yes", part %in% c("Aboveground", "Total")) |>
  dplyr::group_by(scenario, part) |>
  dplyr::group_split()

#### Fit regression model ####

anova_df <- purrr::map_dfr(results_final_list, function(df_temp) {
  
  df_beta <- dplyr::filter(df_temp, measure == "beta")
  
  anova <- aov(value.cv ~ treatment, data = df_beta)
  
  tukey <- TukeyHSD(anova)

  # compact letter display
  cld <- multcompView::multcompLetters4(anova, tukey)
  
  tibble::tibble(treatment = names(cld$treatment$Letters), 
                 label = as.character(cld$treatment$Letters)) |> 
    dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                  .before = treatment)}) |> 
  dplyr::mutate( treatment = factor(treatment, levels = c("enrichment", "connectivity", "combined")))

lm_df <- purrr::map_dfr(results_final_list, function(i) {
  
  dplyr::group_by(i, measure) |> 
    dplyr::group_split() |> 
    purrr::map_dfr(function(j){
      
      df_temp <- dplyr::mutate(j, value.cv = log(value.cv), abiotic = log(abiotic), 
                               biotic = log(biotic))
      
      df_biotic <- dplyr::filter(df_temp, treatment == "connectivity") |>
        dplyr::mutate(value.cv = (value.cv - mean(value.cv)) / sd(value.cv), 
                      abiotic = (abiotic - mean(abiotic)) / sd(abiotic),
                      biotic = (biotic - mean(biotic)) / sd(biotic))
      
      df_abiotic <- dplyr::filter(df_temp, treatment == "enrichment") |>
        dplyr::mutate(value.cv = (value.cv - mean(value.cv)) / sd(value.cv), 
                      abiotic = (abiotic - mean(abiotic)) / sd(abiotic),
                      biotic = (biotic - mean(biotic)) / sd(biotic))
      
      lm_biotic <- lm(value.cv ~ biotic * pop_n, data = df_biotic, na.action = "na.fail")
      lm_abiotic <- lm(value.cv ~ abiotic * nutrient_input, data = df_abiotic, na.action = "na.fail")
      
      dredge_biotic <- MuMIn::dredge(lm_biotic, extra = c("R^2")) |>  
        subset(subset = 1:3) |> 
        tibble::as_tibble() |>
        dplyr::rename("continious" = "biotic", "treatment" = "pop_n", "interaction" = "biotic:pop_n")
      
      dredge_abiotic <- MuMIn::dredge(lm_abiotic, extra = c("R^2")) |> 
        subset(subset = 1:3) |> 
        tibble::as_tibble() |>
        dplyr::rename("continious" = "abiotic", "treatment" = "nutrient_input", 
                      "interaction" = "abiotic:nutrient_input")
      
      dplyr::bind_rows(connect_popn = dredge_biotic, 
                       var_enrich = dredge_abiotic, .id = "model") |> 
        tibble::add_column(part = unique(i$part), response  = unique(j$measure), 
                           .before = "(Intercept)")
      
    })}) |> 
  dplyr::mutate_if(is.numeric, round, digits = 5)
                
#### Setup ggplot ####

base_size <- 13.5

w <- 0.5

color_treatment <- c("enrichment" = "#808fe1", "connectivity" = "#97c684", "combined" = "#efc86e")
# color_treatment <- c("enrichment" = "#0f7ba2", "connectivity" = "#43b284", "combined" =  "#fab255")

gg_dummy <- data.frame(x = c(1, 2, 3), y = c(0.5, 0.5, 0.5), z = c("enrichment", "connectivity", "combined")) |> 
  dplyr::mutate(z = factor(z, levels = c("enrichment", "connectivity", "combined"))) |> 
  ggplot(aes(x = x, y = y, color = z)) + 
  geom_point() + geom_errorbar(aes(ymin = 0, ymax = 1)) +
  scale_color_manual(name = "", values = color_treatment, 
                     labels = c("enrichment" = "Variation nutrient entrichment", 
                                "connectivity" = "Fish Connectivity", 
                                "combined" = "Both")) +
  theme_classic(base_size = base_size) +
  theme(legend.position = "bottom")

#### Create ggplot: Main ####

pe_temp <- dplyr::filter(results_final_df, measure == "beta", 
                         part != "Belowground", include == "yes")
    
pe_sum <- dplyr::group_by(pe_temp, part, treatment) |> 
  dplyr::summarise(mean = mean(value.cv), lower = mean - sd(value.cv), upper = mean + sd(value.cv), 
                   .groups = "drop") |> 
  dplyr::mutate(lower = dplyr::case_when(lower < 0.0 ~ 0.0, TRUE ~ lower))
    
cv_temp <- dplyr::filter(results_final_df, measure %in% c("alpha", "gamma"),
                         part != "Belowground", include == "yes") |> 
  dplyr::select(row_id, part, treatment, measure, value.cv) |> 
  tidyr::pivot_wider(names_from = measure, values_from = value.cv)
    
gg_pe <- ggplot(data = pe_sum, aes(x = part, color = treatment, fill = treatment)) + 
  
  # adding transparent jitter
  geom_point(data = pe_temp, aes(x = part, y = value.cv, color = treatment, group = treatment), 
             alpha = 0.15, size = 0.75, shape = 19,
             position = position_jitterdodge(dodge.width = w, jitter.width = w * 0.75)) +
  
  # adding mean +- sd
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "grey25", position = position_dodge(w), 
                width = 0.0, linewidth = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(w), 
                width = 0.0, linewidth = 1.0) +
  
  geom_point(aes(y = mean), color = "grey25", position = position_dodge(w), size = 2.5) +
  geom_point(aes(y = mean), position = position_dodge(w), size = 1.0) +
  
  geom_text(data = anova_df, position = position_dodge(w), color = "black",
            aes(x = part, y = max(pe_sum$upper) * 1.1, label = label, group = treatment)) +
  annotate("text", x = 0.5, y = 30, label = "a)") +
  
  # adding hline at 1
  geom_hline(yintercept = 1, color = "black", linetype = 2) +
  
  # change scales
  scale_color_manual(name = "", values = color_treatment) +
  scale_y_continuous(limits = c(1, ceiling(max(pe_temp$value.cv) / 10) * 10), 
                     breaks = seq(1, ceiling(max(pe_temp$value.cv)/10) * 10,
                                  length.out = 5)) + 

  # adding box
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +

  # change theme
  labs(x = "", y = expression(paste("Portfolio effect ", italic(CV), italic(beta)))) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "none", axis.line = element_blank())
  
xy_lims <- list(Aboveground = c(0, 1.0), Total = c(0, 0.1)) 
    
gg_cv <- purrr::map(c(Aboveground = "Aboveground", Total = "Total"), function(part_i) {
    
  label_temp <- ifelse(test = part_i == "Aboveground", yes = "b)", no = "c)")
  
  df_temp <- dplyr::filter(cv_temp, part == part_i) 
  
  ggplot() + 
      
    # adding points
    geom_point(data = dplyr::filter(df_temp, treatment != "combined"),
               aes(x = gamma, y = alpha, color = treatment), 
               alpha = 0.15, shape = 19, size = 2.0) +
    
    geom_point(data = dplyr::filter(df_temp, treatment == "combined"),
               aes(x = gamma, y = alpha, color = treatment), 
               alpha = 0.15, shape = 19, size = 2.0) +
    
    # adding PE = 1 line
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2) +
    
    # adding text
    annotate("text", x = 0.0, y = xy_lims[names(xy_lims) == part_i][[1]][[2]], 
             label = label_temp) +
    
    # # change scales
    # coord_equal(ratio = 1) +
    scale_x_continuous(limits = xy_lims[names(xy_lims) == part_i][[1]]) +
    scale_y_continuous(limits = xy_lims[names(xy_lims) == part_i][[1]]) +
    coord_equal(ratio = 1) +
    scale_color_manual(name = "", values = color_treatment) + 
        
    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
        
    # change theme
    # labs(x = expression(paste("Meta-ecosystem scale ", italic(cv[gamma]))), 
    #      y = expression(paste("Local scale ", italic(cv[alpha])))) +
    theme_classic(base_size = base_size) + 
    theme(legend.position = "none", axis.line = element_blank(),
          strip.background = element_blank(), strip.text = element_text(hjust = 0), 
          axis.title = element_blank())
      
})
    
gg_cv <- cowplot::plot_grid(plotlist = gg_cv, ncol = 1) |>
  cowplot::ggdraw(xlim = c(-0.025, 1.0), ylim = c(-0.025, 1.0)) + 
  cowplot::draw_label(expression(paste("Meta-ecosystem scale ", italic(CV), italic(gamma))), 
                      x = 0.5, y = -0.01, angle = 0, size = base_size) + 
  cowplot::draw_label(expression(paste("Local scale ", italic(CV), italic(alpha))),
                      x = 0.01, y = 0.5, angle = 90, size = base_size)

gg_tukey <- cowplot::plot_grid(cowplot::plot_grid(gg_pe, gg_cv, ncol = 2), 
                               cowplot::get_legend(gg_dummy), nrow = 2, ncol = 1,
                               rel_heights = c(0.9, 0.1))

#### Save ggplot #### 

overwrite <- FALSE

if (overwrite) readr::write_csv2(x = lm_df, file = "04_Figures/Table-Q1.csv")

suppoRt::save_ggplot(plot = gg_tukey, filename = "Figure-2.pdf",
                     path = "04_Figures/", width = width, height = height * 0.45,
                     units = units, dpi = dpi, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_tukey, filename = "Figure-2.jpg",
                     path = "04_Figures/", width = width, height = height * 0.45,
                     units = units, dpi = dpi, overwrite = overwrite)
