##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose:

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/import_data.R")

#### Load/wrangle simulated data ####

n <- 5

results_phase_df <- import_data(path = paste0("02_Data/result-phase-", n, ".rds"))

results_noise_df <- import_data(path = paste0("02_Data/result-noise-", n, ".rds"))

#### Setup globals ####

size_base <- 7.5

#### Tile plot of CV (fill) and biotic and abiotic variability ####

gg_abiotic_list <- purrr::map(list(phase = results_phase_df, noise = results_noise_df), function(results_temp_df) {
  
  gg_parts_list <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
    
    purrr::map(c("alpha", "gamma"), function(measure_i) {
      
      list_gg_pop <- purrr::map(c(8, 16, 32, 64, 128), function(pop_i) {
        
        dplyr::filter(results_temp_df, part == part_i, measure == measure_i, pop_n == pop_i) %>% 
          dplyr::mutate(biotic = cut(biotic, breaks = seq(0.0, 1, 0.2),labels = FALSE), 
                        abiotic = cut(abiotic, breaks = seq(0.0, 1, 0.2), labels = FALSE)) %>% 
          dplyr::group_by(biotic, abiotic) %>%
          dplyr::summarise(value.cv = mean(value.cv), .groups = "drop") %>%
          ggplot(aes(x = biotic, y = abiotic, fill = value.cv)) +
          geom_tile() + 
          scale_fill_gradientn(colors = MetBrewer::met.brewer("Demuth", n = 255, type = "continuous")) +
          scale_x_continuous(breaks = c(1, 3, 5), labels = c("low", "medium", "high")) +
          scale_y_continuous(breaks = c(1, 3, 5), labels = c("low", "medium", "high")) +
          labs(title = paste0("Population size: ", pop_i)) +
          theme_classic(base_size = size_base) +
          theme(legend.position = "none", plot.title = element_text(size = 6.5),
                axis.title = element_blank(), axis.line = element_blank(), 
                panel.border = element_rect(size = 0.5, fill = NA))
        
      })
      
      cowplot::plot_grid(plotlist = list_gg_pop, nrow = 2, ncol = 3) +
        cowplot::draw_label(label = "Biotic", x = 0.5, y = 0,
                            vjust = 0.5, angle = 0, size = size_base) +
        cowplot::draw_label(label = "Abiotic", x = 0, y = 0.5,
                            vjust = 0.0, angle = 90, size = size_base) +
        theme(plot.margin = unit(c(t = 0.5, r = 0.5, b = 0.5, l = 0.5), "cm"))
      
    })}) %>% 
    purrr::flatten()
  
  cowplot::plot_grid(plotlist = gg_parts_list, nrow = 3, ncol = 2, labels = "auto", 
                     label_fontface = "italic") 

})

#### Save figures ####

suppoRt::save_ggplot(plot = gg_abiotic_list$phase, filename = paste0("var-cv-tiles-phase-", n, extension),
                     path = "04_Figures/Appendix/", width = width, height = height,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_abiotic_list$noise, filename = paste0("var-cv-tiles-noise-", n, extension),
                     path = "04_Figures/Appendix/", width = width, height = height,
                     units = units, dpi = dpi, overwrite = FALSE)
