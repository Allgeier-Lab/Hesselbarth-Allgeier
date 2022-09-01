##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create total result figure including regression parameters and relative importance

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/import_data.R")

#### Load/wrangle simulated data ####

n <- 5

results_phase_df <- import_data(path = paste0("02_Data/result-phase-", n, ".rds"))

results_noise_df <- import_data(path = paste0("02_Data/result-noise-", n, ".rds"))

#### Setup ggplot ####

colors_pop <- c("8" = "#0069aa", "16" = "#00992a", "32" = "#662e8e", "64" = "#ff771c", "128" = "#f32222")

size_point <- 1.0
size_line <- 0.75
size_base <- 10.0

alpha <- 0.25

gg_dummy <- ggplot(data = results_phase_df, aes(x = biotic, y = value.cv, color = pop_n)) +
  geom_point(shape = 1, alpha = alpha) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = size_line) +
  scale_color_manual(name = "Population size", values = colors_pop) +
  theme_classic(base_size = size_base) +
  theme(legend.position = "bottom")

#### Create ggplot ####

gg_abiotic_list <- purrr::map(list(phase = results_phase_df, noise = results_noise_df), function(results_temp_df) {
  
  gg_scatter_list <- dplyr::filter(results_temp_df, part %in% c("ag_production", "bg_production", "ttl_production"),
                                   measure %in% c("alpha", "gamma")) %>%
    dplyr::group_by(part, measure) %>%
    dplyr::group_split() %>%
    purrr::map(function(temp_df) {
      
      temp_df_long <- tidyr::pivot_longer(temp_df, cols = c(biotic, abiotic),
                                          names_to = "variability", values_to = "value.var")
      
      # init ggplot
      ggplot(data = temp_df_long, aes(x = log(value.var), y = log(value.cv), color = pop_n)) +
        
        # adding geoms
        geom_point(shape = 1, alpha = alpha, size = size_point) +
        geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = size_line) +
        
        # facet wrapping
        facet_wrap(. ~ variability, ncol = 2, scales = "fixed",
                   labeller = labeller(variability = c(biotic = ifelse(test = unique(temp_df$part) == "ag_production",
                                                                       yes = "Biotic", no = ""),
                                                       abiotic = ifelse(test = unique(temp_df$part) == "ag_production",
                                                                        yes = "Abiotic", no = "")))) +
        
        # set scales
        scale_color_manual(name = "Population size", values = colors_pop) +
        scale_y_continuous(limits = log(range(temp_df_long$value.cv)),
                           breaks = seq(log(min(temp_df_long$value.cv)), log(max(temp_df_long$value.cv)),
                                        length.out = 5),
                           labels = function(x) round(exp(x), 2)) +
        scale_x_continuous(limits = log(range(temp_df_long$value.var)),
                           breaks = seq(log(min(temp_df_long$value.var)), log(max(temp_df_long$value.var)),
                                        length.out = 5),
                           labels = function(x) round(exp(x), 2)) +
        
        # labels and themes
        labs(x = "", y = "") +
        theme_classic(base_size = 10.0) +
        theme(legend.position = "none", plot.title = element_text(size = 8.0),
              strip.background = element_blank(), strip.text = element_text(hjust = 0),
              plot.margin = margin(t = 0.0, r = 5.5, b = 0.0, l = 5.5, unit = "pt"))
      
    })
  
  gg_combined <- cowplot::plot_grid(plotlist = gg_scatter_list, nrow = 3, ncol = 2, byrow = TRUE,
                                    labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
                                    label_fontface = "plain", label_size = 12.0)
  
  gg_combined <- cowplot::ggdraw(gg_combined, xlim = c(-0.015, 1.0)) +
    cowplot::draw_label("Variability", x = 0.5, y = 0, vjust = -0.5, angle = 0, size = size_base) +
    cowplot::draw_label("Coeffiecent of variation", x = 0.0, y = 0.5, vjust = 0.0, angle = 90, size = size_base)
  
  cowplot::plot_grid(gg_combined, cowplot::get_legend(gg_dummy),
                     nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))
  
})

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_abiotic_list$phase, filename = paste0("Figure-A5", extension),
                     path = "04_Figures/Appendix", width = height, height = width * 0.65,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_abiotic_list$noise, filename = paste0("Figure-A6", extension),
                     path = "04_Figures/Appendix", width = height, height = width * 0.65,
                     units = units, dpi = dpi, overwrite = FALSE)
