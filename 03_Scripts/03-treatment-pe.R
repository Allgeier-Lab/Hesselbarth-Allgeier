##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-cv.R")
source("01_Functions/import-mortality.R")

#### Load/wrangle simulated data ####

results_noise_df <- import_cv(path = "02_Data/result-noise.rds")

#### Load mortality data ####

mortality_noise_df <- import_mortality(path = "02_Data/result-noise.rds") |> 
  dplyr::filter(pop_n != 0) |>
  dplyr::group_by(row_id, pop_n, nutrient_input) |>
  dplyr::summarise(died_total = mean(died_total), .groups = "drop")

#### Filter abundances/mortality #### 

results_final_df <- dplyr::left_join(x = results_noise_df, y = mortality_noise_df, 
                                     by = c("row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(died_total > threshold_mort ~ "no", TRUE ~ "yes")) |> 
  dplyr::filter(include == "yes") |> 
  dplyr::mutate(treatment = dplyr::case_when(abiotic != 0.0 & biotic == 0.0 ~ "enrichment", 
                                             abiotic == 0.0 & biotic != 0.0 ~ "connectivity", 
                                             abiotic != 0.0 & biotic != 0.0 ~ "combined")) |>
  dplyr::mutate(part = factor(part, levels = c("ag_production", "bg_production", "ttl_production"),
                              labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                         "ttl_production" = "Total")),
                treatment = factor(treatment, levels = c("enrichment", "connectivity", "combined")))

#### Split data #### 

results_final_list <- dplyr::filter(results_final_df, measure != "synchrony", treatment != "combined",
                                    part %in% c("Aboveground", "Total")) |>
  dplyr::group_by(part) |>
  dplyr::group_split()

#### Fit regression model ####

relimp_df <- purrr::map_dfr(results_final_list, function(i) {
  
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
      
      relimp_biotic <- lm(lm_biotic) |> 
        relaimpo::boot.relimp(type = "lmg", b = 500, level = 0.95, fixed = FALSE) |>
        relaimpo::booteval.relimp(bty = "basic")
      
      relimp_abiotic <- lm(lm_abiotic) |> 
        relaimpo::boot.relimp(type = "lmg", b = 500, level = 0.95, fixed = FALSE) |>
        relaimpo::booteval.relimp(bty = "basic")
      
      
    dplyr::bind_rows(
        tibble::tibble(
          beta = c(relimp_biotic@namen[-1], "residual"), mean = c(relimp_biotic@lmg, 1 - sum(relimp_biotic@lmg)),
          lower = c(relimp_biotic@lmg.lower, NA), higher = c(relimp_biotic@lmg.upper, NA)) |> 
          tibble::add_column(model = "biotic", part = unique(i$part), response  = unique(j$measure), 
                             .before = "beta"),
        
        tibble::tibble(
          beta = c(relimp_abiotic@namen[-1], "residual"), mean = c(relimp_abiotic@lmg, 1 - sum(relimp_abiotic@lmg)),
          lower = c(relimp_abiotic@lmg.lower, NA), higher = c(relimp_abiotic@lmg.upper, NA)) |> 
          tibble::add_column(model = "abiotic", part = unique(i$part), response  = unique(j$measure), 
                             .before = "beta")
    ) |> 
      dplyr::mutate(dplyr::across(mean:higher, ~ .x * 100), 
                    lower = dplyr::case_when(lower < 0 ~ 0.0, TRUE ~ lower))})}) |> 
  dplyr::mutate(beta = factor(beta, levels = c("abiotic", "nutrient_input", "abiotic:nutrient_input", 
                                               "biotic", "pop_n", "biotic:pop_n", "residual"), 
                              labels = c("Spatial variation", "Nutrient enrichment", "Variation:Enrichment", 
                                         "Connectivity", "Population size", "Connectivity:Population", "Residuals")))

#### Setup ggplot ####

base_size <- 13.5

w <- 1.0

arrow_width <- 0.3
arrow_size <- 1.0

color_beta <- c("Nutrient enrichment" = "#74a9cf", "Spatial variation" = "#045a8d", "Variation:Enrichment" = "#f1eef6",
                "Population size" = "#fd8d3c", "Connectivity" = "#a63603", "Connectivity:Population" = "#fdd0a2",
                "Residuals" = "grey")

# color_treatment <- c("Nutrient enrichment" = "#BD92BF", "Fish dynamics" = "#749F6F")
color_treatment <- c("Nutrient enrichment" = "#a6bddb", "Fish dynamics" = "#fdae6b")

#### Create ggplot: Main ####

pe_values <- dplyr::filter(results_final_df, measure == "beta", treatment != "combined",
                           part != "Belowground") |> 
  dplyr::mutate(treatment = factor(treatment, levels = c("enrichment", "connectivity"), 
                                   labels = c("Nutrient enrichment", "Fish dynamics")), 
                part = factor(part, levels = c("Aboveground", "Total"), labels = c("AG PP", "TTL PP")))
    
pe_sum <- dplyr::group_by(pe_values, part, treatment) |> 
  dplyr::summarise(mean = mean(value.cv), lower = mean - sd(value.cv), upper = mean + sd(value.cv), 
                   .groups = "drop") |> 
  dplyr::mutate(lower = dplyr::case_when(lower < 0.0 ~ 0.0, TRUE ~ lower))

cv_temp <- dplyr::filter(results_final_df, measure %in% c("alpha", "gamma"),
                         part != "Belowground", treatment != "combined") |> 
  dplyr::select(row_id, part, treatment, measure, value.cv) |> 
  tidyr::pivot_wider(names_from = measure, values_from = value.cv) |> 
  dplyr::mutate(treatment = factor(treatment, levels = c("enrichment", "connectivity"), 
                                   labels = c("Nutrient enrichment", "Fish dynamics")))
    
gg_pe <- ggplot(data = pe_sum, aes(x = part, color = treatment, group = treatment)) + 
  
  # adding transparent jitter
  geom_point(data = pe_values, aes(x = part, y = value.cv, color = treatment),
             alpha = 0.25, size = 1.5, shape = 19,
             position = position_jitterdodge(dodge.width = w * 0.8, jitter.width = w * 0.75)) +

  # adding mean +- sd
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "black", position = position_dodge(w * 0.8),
                width = 0.0, linewidth = 2.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(w * 0.8),
                width = 0.0, linewidth = 1.0) +
  
  # adding shade around means
  geom_point(aes(y = mean), color = "black", position = position_dodge(w * 0.8), size = 3.5) +
  geom_point(aes(y = mean), position = position_dodge(w * 0.8), size = 2.5) +
  
  # change scales
  scale_color_manual(name = "", values = color_treatment) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_y_continuous(limits = c(0, ceiling(max(pe_values$value.cv)))) +

  # adding box
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +

  # change theme
  guides(color = guide_legend(nrow = 2)) +
  labs(x = "", y = "Portfolio effect (PE)") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = c(0.2, 0.95), legend.text = element_text(size = base_size * 0.85), 
        legend.background = element_rect(colour = NA, fill = NA), axis.line = element_blank())
  
# set ylims based on AG vs TTL
xy_lims <- list(Aboveground = c(0, 1.0), Total = c(0, 0.1)) 
    
gg_cv <- purrr::map(c(Aboveground = "Aboveground", Total = "Total"), function(part_i) {
    
  x_lab <- expression(paste("AG PP - Meta-ecosystem scale (", italic(CV), italic(gamma), ")"))
  y_lab <-  expression(paste("AG PP - Local scale (", italic(CV), italic(alpha), ")"))
  
  # set location of arrows for AG
  # light blue (alpha)
  arrow_1_x_lo <- 0.0
  arrow_1_x_hi <- 0.0
  arrow_1_y_lo <- 0.95
  arrow_1_y_hi <- 0.65
  arrow_1 <- arrow(length = unit(arrow_width, "cm"))
  
  # dark blue (alpha)
  arrow_2_x_lo <- 0.075
  arrow_2_x_hi <- 0.075
  arrow_2_y_lo <- 0.65
  arrow_2_y_hi <- 0.95
  arrow_2 <- arrow(length = unit(arrow_width, "cm"))
  
  # light orange (alpha)
  arrow_3_x_lo <- 0.2
  arrow_3_x_hi <- 0.2
  arrow_3_y_lo <- 0.95
  arrow_3_y_hi <- 0.65
  arrow_3 <- arrow(length = unit(arrow_width, "cm"))
  
  # dark orange (alpha)
  arrow_4_x_lo <- 0.275
  arrow_4_x_hi <- 0.275
  arrow_4_y_lo <- 0.95
  arrow_4_y_hi <- 0.65
  arrow_4 <- arrow(length = unit(arrow_width, "cm"))
  
  # light blue (gamma)
  arrow_5_x_lo <- 0.95
  arrow_5_x_hi <- 0.65
  arrow_5_y_lo <- 0.0
  arrow_5_y_hi <- 0.0
  arrow_5 <- arrow(length = unit(arrow_width, "cm"))
  
  # dark blue (gamma)
  arrow_6_x_lo <- 0.65
  arrow_6_x_hi <- 0.95
  arrow_6_y_lo <- 0.075
  arrow_6_y_hi <- 0.075
  arrow_6 <- arrow(length = unit(arrow_width, "cm"))
  
  # light orange (gamma)
  arrow_7_x_lo <- 0.95
  arrow_7_x_hi <- 0.65
  arrow_7_y_lo <- 0.2
  arrow_7_y_hi <- 0.2
  arrow_7 <- arrow(length = unit(arrow_width, "cm"))
  
  # dark orange (gamma)
  arrow_8_x_lo <- 0.9
  arrow_8_x_hi <- 0.7
  arrow_8_y_lo <- 0.275
  arrow_8_y_hi <- 0.275
  arrow_8 <- NULL
  
  # TTL is on a different scale
  if (part_i == "Total") {
    
    x_lab <- expression(paste("TTL PP - Meta-ecosystem scale (", italic(CV), italic(gamma), ")"))
    y_lab <-  expression(paste("TTL PP - Local scale (", italic(CV), italic(alpha), ")"))
    
    # set location of arrows for AG
    # light blue (alpha)
    arrow_1_x_lo <- 0.0
    arrow_1_x_hi <- 0.0
    arrow_1_y_lo <- 0.095
    arrow_1_y_hi <- 0.065
    arrow_1 <- arrow(length = unit(arrow_width, "cm"))
    
    # dark blue (alpha)
    arrow_2_x_lo <- 0.0075
    arrow_2_x_hi <- 0.0075
    arrow_2_y_lo <- 0.065
    arrow_2_y_hi <- 0.095
    arrow_2 <- arrow(length = unit(arrow_width, "cm"))
    
    # light orange (alpha)
    arrow_3_x_lo <- 0.02
    arrow_3_x_hi <- 0.02
    arrow_3_y_lo <- 0.065
    arrow_3_y_hi <- 0.095
    arrow_3 <- arrow(length = unit(arrow_width, "cm"))
    
    # dark orange (alpha)
    arrow_4_x_lo <- 0.0275
    arrow_4_x_hi <- 0.0275
    arrow_4_y_lo <- 0.095
    arrow_4_y_hi <- 0.065
    arrow_4 <- arrow(length = unit(arrow_width, "cm"))
    
    # light blue (gamma)
    arrow_5_x_lo <- 0.09
    arrow_5_x_hi <- 0.07
    arrow_5_y_lo <- 0.0
    arrow_5_y_hi <- 0.0
    arrow_5 <- NULL
    
    # dark blue (gamma)
    arrow_6_x_lo <- 0.065
    arrow_6_x_hi <- 0.095
    arrow_6_y_lo <- 0.0075
    arrow_6_y_hi <- 0.0075
    arrow_6 <- arrow(length = unit(arrow_width, "cm"))
    
    # light orange (gamma)
    arrow_7_x_lo <- 0.095
    arrow_7_x_hi <- 0.065
    arrow_7_y_lo <- 0.02
    arrow_7_y_hi <- 0.02
    arrow_7 <- arrow(length = unit(arrow_width, "cm"))
    
    # dark orange (gamma)
    arrow_8_x_lo <- 0.09
    arrow_8_x_hi <- 0.07
    arrow_8_y_lo <- 0.0275
    arrow_8_y_hi <- 0.0275
    arrow_8 <- NULL
    
  #   arrow_x_2 <- 0.0075
  #   arrow_x_3 <- 0.02
  #   arrow_x_4 <- 0.0275
  #   
  #   arrow_y_1 <- 0.065
  #   arrow_y_2 <- 0.095
  #   
  #   arrow_y_11 <- 0.095
  #   arrow_y_22 <- 0.065
  #   
  }
  
  df_temp <- dplyr::filter(cv_temp, part == part_i)
  
  ggplot() + 
      
    # adding points
    geom_point(data = dplyr::filter(df_temp, treatment != "combined"),
               aes(x = gamma, y = alpha, color = treatment), 
               alpha = 0.5, shape = 19, size = 2.5) +
    
    geom_point(data = dplyr::filter(df_temp, treatment == "combined"),
               aes(x = gamma, y = alpha, color = treatment), 
               alpha = 0.5, shape = 19, size = 2.5) +
    
    # adding PE = 1 line
    geom_abline(slope = 1, intercept = 0, color = "grey", linetype = 2) +
    
    # add arrows y direction
    geom_segment(aes(x = arrow_1_x_lo, y = arrow_1_y_lo, xend = arrow_1_x_hi, yend = arrow_1_y_hi), 
                 color = color_beta[[1]], arrow = arrow_1, linewidth = arrow_size) +
    geom_segment(aes(x = arrow_2_x_lo, y = arrow_2_y_lo, xend = arrow_2_x_hi, yend = arrow_2_y_hi), 
                 color = color_beta[[2]], arrow = arrow_2, linewidth = arrow_size) +
    geom_segment(aes(x = arrow_3_x_lo, y = arrow_3_y_lo, xend = arrow_3_x_hi, yend = arrow_3_y_hi), 
                 color = color_beta[[4]], arrow = arrow_3, linewidth = arrow_size) +
    geom_segment(aes(x = arrow_4_x_lo, y = arrow_4_y_lo, xend = arrow_4_x_hi, yend = arrow_4_y_hi), 
                 color = color_beta[[5]], arrow = arrow_4, linewidth = arrow_size) +
    
    # add arrows x direction
    geom_segment(aes(x = arrow_5_x_lo, y = arrow_5_y_lo, xend = arrow_5_x_hi, yend = arrow_5_y_hi), 
                 color = color_beta[[1]], arrow = arrow_5, linewidth = arrow_size) +
    geom_segment(aes(x = arrow_6_x_lo, y = arrow_6_y_lo, xend = arrow_6_x_hi, yend = arrow_6_y_hi), 
                 color = color_beta[[2]], arrow = arrow_6, linewidth = arrow_size) +
    geom_segment(aes(x = arrow_7_x_lo, y = arrow_7_y_lo, xend = arrow_7_x_hi, yend = arrow_7_y_hi), 
                 color = color_beta[[4]], arrow = arrow_7, linewidth = arrow_size) +
    geom_segment(aes(x = arrow_8_x_lo, y = arrow_8_y_lo, xend = arrow_8_x_hi, yend = arrow_8_y_hi), 
                 color = color_beta[[5]], arrow = arrow_8, linewidth = arrow_size) +
    
    # change scales
    scale_x_continuous(limits = xy_lims[names(xy_lims) == part_i][[1]]) +
    scale_y_continuous(limits = xy_lims[names(xy_lims) == part_i][[1]]) +
    coord_equal(ratio = 1) +
    scale_color_manual(name = "", values = color_treatment) + 
    scale_fill_manual(values = color_beta) +
        
    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
        
    # change theme
    labs(x = x_lab, y = y_lab) +
    theme_classic(base_size = base_size) + 
    theme(legend.position = "none", axis.line = element_blank(),
          strip.background = element_blank(), strip.text = element_text(hjust = 0), 
          axis.title = element_text(size = 11))

})

# get maximum value
max_imp <- dplyr::filter(relimp_df, beta != "Residuals", response == "beta") |> 
  dplyr::pull(higher) |> 
  max()

gg_relimp <- dplyr::filter(relimp_df, beta != "Residuals", response == "beta") |>
  ggplot() + 
  
  # add errorbars?
  geom_errorbar(aes(x = beta, ymin = lower, ymax = higher, color = beta, group = part),
                position = position_dodge(w * 0.75), width = 0, linewidth = 0.5) +
  
  # adding coloumns with relative importance
  geom_col_pattern(aes(x = beta, y = mean, fill = beta, pattern = part), 
                   pattern_spacing = 0.02, pattern_density = 0.01, pattern_color = "grey",
                   position = position_dodge(w * 0.75), color = NA, width = w * 0.70) + 
  
  geom_text(data = data.frame(model = c("abiotic", "biotic"), label = c("Nutrient enrichment", "Fish dynamics"),
                              x = c(1, 1), y = c(ceiling(max_imp / 5) * 5,  
                              ceiling(max_imp / 5) * 5)), aes(x = x, y = y, label = label),
            nudge_x = c(0.4, 0.2), size = 3.5) +
  
  # facet wrap by model type
  facet_wrap(. ~ model, nrow = 1, labeller = labeller(model = c("abiotic" = "Variation nutrient entrichment",  
                                                                "biotic" = "Fish Connectivity")), 
             scales = "free_x") +
  
  # adding box
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  
  # set scales
  scale_fill_manual(name = "", values = color_beta) +
  scale_color_manual(name = "", values = color_beta) +
  scale_pattern_manual(name = "", values = c("Aboveground"  = "none", "Total" = "stripe"), 
                       labels = c("Aboveground" = "AG PP", "Total" = "TTL PP")) +
  scale_x_discrete(labels = c("Spatial variation" = "Enrich. var.", "Nutrient enrichment" = "Enrichment",
                              "Variation:Enrichment" = "Interaction", 
                              "Population size" = "Population",
                              "Connectivity:Population" = "Interaction")) +
  scale_y_continuous(limits = c(0, ceiling(max_imp / 5) * 5), 
                     breaks = seq(0, floor(max_imp / 10) * 10, 10)) +
  
  # set themes
  guides(fill = "none", color = "none", 
         pattern = guide_legend(override.aes = list(fill = "grey", color = "black", 
                                                    pattern_color = "black"))) +
  # guides(shape = guide_legend(override.aes = list(size = 5))) +
  labs(y = "Relative importance [%]", x = "") + 
  theme_classic(base_size = base_size) + 
  theme(legend.position = c(0.925, 0.5), legend.background = element_blank(),
        legend.text = element_text(size = 8.5), legend.key.size = unit(0.4, "cm"),
        axis.line = element_blank(), axis.text.x = element_text(size = 8), axis.title.y = element_text(hjust = 0.2),
        strip.text = element_blank(), strip.background = element_blank())

gg_combined <- cowplot::plot_grid(gg_pe, gg_cv$Aboveground, gg_relimp, gg_cv$Total,
                                  labels = c("A)", "C)", "B)", "D)"), label_fontface = "plain",
                                  hjust = c(-0.5, -1, -0.5, -0.75), rel_widths = c(0.6, 0.4, 0.6, 0.4))

#### Save ggplot #### 

suppoRt::save_ggplot(plot = gg_combined, filename = "Figure-2.png",
                     path = "04_Figures/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_combined, filename = "Figure-2.pdf",
                     path = "04_Figures/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = FALSE)
