##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create tile plot of biotic/abiotic vs. cv

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-cv.R")

#### Load/wrangle simulated data ####

n <- 5

results_phase_df <- import_cv(path = paste0("02_Data/result-phase-", n, ".rds"))

results_noise_df <- import_cv(path = paste0("02_Data/result-noise-", n, ".rds"))

results_combinded_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, .id = "scenario") |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))

#### Setup globals ####

size_base <- 7.5

#### Tile plot of CV (fill) and biotic and abiotic variability ####

gg_abiotic_list <- purrr::map(c("phase", "noise"), function(scenario_i) {
  
  gg_part <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
    
    results_temp_df <- dplyr::filter(results_combinded_df, scenario == scenario_i, part == part_i, 
                                     measure %in% c("alha", "gamma")) |> 
      dplyr::mutate(biotic = cut(biotic, breaks = seq(0.0, 1, 0.2),labels = FALSE), 
                    abiotic = cut(abiotic, breaks = seq(0.0, 1, 0.2), labels = FALSE)) |> 
      dplyr::group_by(measure, pop_n) |> 
      dplyr::mutate(value.cv = scales::rescale(value.cv)) |> 
      dplyr::group_by(pop_n, biotic, abiotic) |> 
      dplyr::summarise(value.cv = mean(value.cv), .groups = "drop")
    
    if (part_i == "ag_production") label_part <-  "Aboveground"
    if (part_i == "bg_production") label_part <-  "Belowground"
    if (part_i == "ttl_production") label_part <- "Total      "
    
    p <- ggplot(data = results_temp_df, aes(x = biotic, y = abiotic, fill = value.cv)) +
      
      # adding geoms
      geom_tile() + 
      
      # facet panel
      facet_wrap(. ~ pop_n, ncol = 5, labeller = labeller(pop_n = function(x) paste("Pop. size", x))) +
      
      # setting scales
      scale_fill_gradientn(colors = MetBrewer::met.brewer("Demuth", n = 255, type = "continuous")) +
      scale_x_continuous(breaks = c(1, 3, 5), labels = c("low", "medium", "high")) +
      scale_y_continuous(breaks = c(1, 3, 5), labels = c("low", "medium", "high")) +
      coord_fixed(ratio = 1) + 
     
      # change theme
      # labs(title = paste0("Population size: ", pop_i)) +
      theme_classic(base_size = size_base) +
      theme(strip.background = element_blank(), strip.text = element_text(hjust = 0), 
            axis.title = element_blank(), axis.line = element_blank(), 
            legend.position = "none", 
            panel.border = element_rect(size = 0.5, fill = NA), 
            plot.margin = unit(c(t = 0.0, r = 1.0, b = 0.0, l = 1.0), "mm"))
    
    cowplot::ggdraw(p, xlim = c(0.0, 1.05)) + 
      cowplot::draw_label(label = label_part, x = 1.0, y = 0.65,
                          vjust = -0.5, angle = 270, size = size_base * 0.85)
    
  })
  
  gg_combined <- cowplot::plot_grid(plotlist = gg_part, nrow = 3)
  
  cowplot::ggdraw(gg_combined, xlim = c(-0.05, 1.0), ylim = c(-0.05, 1.0)) + 
    cowplot::draw_label("Diversiy consumer behavior", x = 0.5, y = 0, vjust = 0.5, angle = 0, size = size_base) + 
    cowplot::draw_label("Variability abiotic subsidies", x = 0.0, y = 0.5, vjust = -0.5, angle = 90, size = size_base)
    
})

names(gg_abiotic_list) <- c("phase", "noise")

#### Save figures ####

# suppoRt::save_ggplot(plot = gg_abiotic_list$phase, filename = paste0("var-cv-tiles-phase-", n, extension),
#                      path = "04_Figures/Appendix/", width = width, height = height,
#                      units = units, dpi = dpi, overwrite = FALSE)
# 
# suppoRt::save_ggplot(plot = gg_abiotic_list$noise, filename = paste0("var-cv-tiles-noise-", n, extension),
#                      path = "04_Figures/Appendix/", width = width, height = height,
#                      units = units, dpi = dpi, overwrite = FALSE)
