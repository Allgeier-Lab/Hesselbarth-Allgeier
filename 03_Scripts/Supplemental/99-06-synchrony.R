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
source("01_Functions/plot-lm.R")

base_size <- 13.5
alpha_pts <- 1
line_width <- 0.5

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
  dplyr::mutate(treatment = dplyr::case_when(abiotic != 0.0 & biotic == 0.0 ~ "subsidies", 
                                             abiotic == 0.0 & biotic != 0.0 ~ "connectivity", 
                                             abiotic != 0.0 & biotic != 0.0 ~ "combined")) |>
  dplyr::mutate(part = factor(part, levels = c("ag_production", "bg_production", "ttl_production"),
                              labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                         "ttl_production" = "Total")),
                treatment = factor(treatment, levels = c("subsidies", "connectivity", "combined")))

#### Pre-process data ####

synchrony_pe_df <- dplyr::filter(results_final_df, measure %in% c("beta", "synchrony"), 
                                 part != "Belowground") |> 
  tidyr::pivot_wider(names_from = measure, values_from = value.cv) |> 
  dplyr::mutate(treatment = factor(treatment, levels = c("subsidies", "connectivity", "combined"), 
                                   labels = c("Nutrient enrichment", "Fish dynamics", "Combined"))) |> 
  dplyr::select(-c(value.mn, value.sd, value.prod))

#### ggplot Synchrony vs PE ####

color_treatment <- c("Nutrient enrichment" = "#a6bddb", "Fish dynamics" = "#fdae6b", "Combined" = "#ED679A")

synchrony_pe_list <- dplyr::filter(synchrony_pe_df) |> 
  dplyr::group_by(treatment) |> 
  dplyr::group_split()

# create plot for AG and Total
gg_syn_pe_list <- purrr::map(synchrony_pe_list, function(df_temp) {
  
  axis_element <- NULL
  plot_margin <- unit(c(0.5, 0.25, -0.5, -0.25), "cm") # t r b l
  foo <- function(x) sprintf("%.2f", x)
  
  if (unique(df_temp$treatment) %in% c("Nutrient enrichment", "Fish dynamics")) axis_element <- element_blank()
  if (unique(df_temp$treatment) %in% c("Fish dynamics", "Combined")) foo <- function(x) sprintf("%.1f", x)
  if (unique(df_temp$treatment)  %in% c("Fish dynamics", "Combined")) plot_margin <- unit(c(0.1, 0.25, -0.5, -0.25), "cm")
  
  ggplot(data = df_temp, aes(x = synchrony, y = beta, color = treatment)) + 
    
    # adding geoms
    geom_hline(yintercept = 1, linetype = 2, color = "lightgrey") +
    geom_point(pch = 1, alpha = alpha_pts) +
    geom_smooth(se = FALSE, span = 0.1, linewidth = line_width, linetype = 2) +
    
    # facet wrap by treatment (rows) and part (cols)
    facet_wrap(. ~ part, scales = "fixed") +
      
    # change scales
    scale_color_manual(name = "", values = color_treatment) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(labels = foo) +
    labs(x = "", y = "") + 
      
    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
    # change themes
    theme_classic(base_size = base_size) + 
    theme(legend.position = "none", axis.line = element_blank(), axis.text.x = axis_element,
          axis.ticks.x = axis_element,
          strip.background = element_blank(), strip.text = element_blank(), 
          plot.margin = plot_margin)
})

# combine plots
gg_syn_pe <- cowplot::plot_grid(plotlist = gg_syn_pe_list, nrow = 3)

gg_syn_pe <- cowplot::ggdraw(gg_syn_pe, xlim = c(-0.05, 1), ylim = c(-0.05, 1.025)) +
  cowplot::draw_label("Synchrony", x = 0.525, y = -0.025, size = base_size) + 
  cowplot::draw_label("Portfolio effect (PE)", x = -0.025, y = 0.5, size = base_size, angle = 90) +
  cowplot::draw_label("Aboveground", x = 0.27, y = 1.0, size = base_size * 0.65) + 
  cowplot::draw_label("Total", x = 0.8, y = 1.0, size = base_size * 0.65) + 
  cowplot::draw_label("Nutr. enrich.", x = 0.0, y = 0.8, size = base_size * 0.65, angle = 90) +
  cowplot::draw_label("Fish dynamics", x = 0.0, y = 0.5, size = base_size * 0.65, angle = 90) + 
  cowplot::draw_label("Combined", x = 0.0, y = 0.175, size = base_size * 0.65, angle = 90)

#### ggplot Treatment vs Synchrony ####

color_pop <- c("8" = "#888888", "16" = "#117733", "32" = "#88CCEE", 
               "64" = "#DDCC77", "128" = "#CC6677")

# create dummy plot go grab legend
gg_dummy <- data.frame(pop_n = c(8, 16, 32, 64, 128), x = 1:5, y = 1) |> 
  ggplot(aes(x = x, y = y, color = factor(pop_n))) + 
  geom_point() + geom_line() +
  scale_color_manual(name = "Pop. size", values = color_pop) + 
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

gg_legend <- cowplot::get_legend(gg_dummy)

treat_synchrony_list <- dplyr::filter(synchrony_pe_df, treatment == "Combined") |> 
  dplyr::group_by(part) |> 
  dplyr::group_split()

gg_trt_syn_list <- purrr::map(treat_synchrony_list, function(df_temp) {
  
  axis_element <- NULL
  plot_margin <- unit(c(0.5, 0.25, 0.25, 0.25), "cm")
  # foo <- function(x) sprintf("%.2f", x)
  
  if (unique(df_temp$part) == "Total") {
    axis_element <- element_blank()
    plot_margin <- unit(c(0.5, 0.25, 0.25, -0.25), "cm")
  }

  dplyr::select(df_temp, row_id, part, pop_n, nutrient_input, biotic, abiotic, synchrony) |>
    tidyr::pivot_longer(-c(row_id, part, pop_n, nutrient_input, synchrony)) |> 
    ggplot(aes(x = value, y = synchrony, color = pop_n)) + 
    
    # adding geoms
    geom_point(pch = 1, alpha = alpha_pts) + 
    geom_smooth(se = FALSE, linewidth = line_width, span = 0.75) +
    
    facet_wrap(. ~ nutrient_input + name, nrow = 5, ncol = 2, scales = "fixed") +
    
    # change scales
    scale_color_manual(name = "", values = color_pop) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1.0)) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1.0)) +
    
    labs(x = "", y = "") +
    
    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    
    # change themes
    theme_classic(base_size = base_size) + 
    theme(legend.position = "none", axis.line = element_blank(), axis.text.y = axis_element,
          axis.ticks.y = axis_element,
          strip.background = element_blank(), strip.text = element_blank(), 
          plot.margin = plot_margin)
  
})

gg_trt_syn <- cowplot::plot_grid(plotlist = gg_trt_syn_list, ncol = 2, byrow = FALSE) |> 
  cowplot::ggdraw(xlim = c(-0.01, 1), ylim = c(0, 1.0)) +
  cowplot::draw_label("Enrich. var. / Connectivity", x = 0.5, y = 0.028, size = base_size) + 
  cowplot::draw_label("Synchrony", x = 0.005, y = 0.55, size = base_size, angle = 90) +
  cowplot::draw_label("Aboveground", x = 0.28, y = 0.985, size = base_size * 0.85) +
  cowplot::draw_label("Total", x = 0.755, y = 0.985, size = base_size* 0.85) +
  cowplot::draw_label("Enrich. var.", x = 0.17, y = 0.985, size = base_size * 0.65) +
  cowplot::draw_label("Connectivity", x = 0.39, y = 0.985, size = base_size * 0.65) +
  cowplot::draw_label("Enrich. var.", x = 0.635, y = 0.985, size = base_size * 0.65) + 
  cowplot::draw_label("Connectivity", x = 0.875, y = 0.985, size = base_size * 0.65) + 
  cowplot::draw_label("Low", x = 0.025, y = 0.825, size = base_size * 0.65, angle = 90) +
  cowplot::draw_label("Medium", x = 0.025, y = 0.535, size = base_size * 0.65, angle = 90) +
  cowplot::draw_label("High", x = 0.025, y = 0.235, size = base_size * 0.65, angle = 90)

gg_trt_syn <- cowplot::plot_grid(gg_trt_syn, gg_legend, rel_heights = c(0.95, 0.05), nrow = 2)

### Save total figure ####

suppoRt::save_ggplot(plot = gg_syn_pe, filename = "Figure-S6.png",
                     path = "04_Figures/Supplemental/", width = width, height = height * 2 / 3,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_trt_syn, filename = "Figure-S7.png",
                     path = "04_Figures/Supplemental/", width = width, height = height * 1 / 2,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_syn_pe, filename = "Figure-S6.pdf",
                     path = "04_Figures/Supplemental/", width = width, height = height * 2 / 3,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_trt_syn, filename = "Figure-S7.pdf",
                     path = "04_Figures/Supplemental/", width = width, height = height * 1 / 2,
                     units = units, dpi = dpi, overwrite = FALSE)
