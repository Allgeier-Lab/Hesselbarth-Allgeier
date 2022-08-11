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
source("05_Various/import_data.R")

extension <- ".pdf"

amplitude <- "095"

#### Load/wrangle simulated data ####

file_path <- paste0("02_Data/05-variability-noise-", amplitude, ".rds")

df_results <- import_data(path = file_path)

#### Setup ggplot ####

colors_pop <- c("8" = "#447861", "16" = "#13315f", "32" = "#59386c", "64" = "#b24422")

size_point <- 1.0
size_line <- 0.75
size_base <- 10.0

alpha <- 0.1

#### Create ggplot ####

gg_scatter_list <- dplyr::filter(df_results, measure %in% c("alpha", "gamma")) %>%
  dplyr::group_by(part, measure) %>%
  dplyr::group_split() %>%
  purrr::map(function(df_temp) {
    
    df_temp_long <- tidyr::pivot_longer(df_temp, cols = c(move_meta_sd, noise_sd), 
                                        names_to = "variability", values_to = "value.var")
  
    # init ggplot
    ggplot(data = df_temp_long, aes(x = log(value.var), y = log(value.cv), color = pop_n)) + 
    
    # adding geoms
    geom_point(shape = 1, alpha = alpha, size = size_point) + 
    geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = size_line) +
    
    # facet wrapping
    facet_wrap(. ~ variability, ncol = 2, scales = "fixed",
               labeller = labeller(variability = c(move_meta_sd = ifelse(test = unique(df_temp$part) == "ag_production", 
                                                                         yes = "Biotic", no = ""), 
                                                   noise_sd = ifelse(test = unique(df_temp$part) == "ag_production", 
                                                                     yes = "Abiotic", no = "")))) +
    
    # set scales
    scale_color_manual(name = "Population size", values = colors_pop) +
    scale_y_continuous(limits = log(range(df_temp_long$value.cv)), 
                       breaks = seq(log(min(df_temp_long$value.cv)), log(max(df_temp_long$value.cv)), 
                                    length.out = 5), 
                       labels = function(x) round(exp(x), 2)) +
    scale_x_continuous(limits = log(range(df_temp_long$value.var)), 
                       breaks = seq(log(min(df_temp_long$value.var)), log(max(df_temp_long$value.var)), 
                                    length.out = 5), 
                       labels = function(x) round(exp(x), 2)) +
      
    # labels and themes
    labs(x = "", y = "") +
    theme_classic(base_size = 10.0) + 
    theme(legend.position = "none", plot.title = element_text(size = 8.0), 
          strip.background = element_blank(), strip.text = element_text(hjust = 0),
          plot.margin = margin(t = 0.0, r = 5.5, b = 0.0, l = 5.5, unit = "pt"))
  
})

gg_dummy <- ggplot(data = df_results, aes(x = move_meta_sd, y = value.cv, color = pop_n)) + 
  geom_point(shape = 1, alpha = alpha) + 
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = size_line) +
  scale_color_manual(name = "Population size", values = colors_pop) +
  theme_classic(base_size = size_base) + 
  theme(legend.position = "bottom")

gg_combined <- cowplot::plot_grid(plotlist = gg_scatter_list, nrow = 3, ncol = 2, byrow = TRUE, 
                                  labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
                                  label_fontface = "plain", label_size = 12.0)

gg_combined <- cowplot::ggdraw(gg_combined, xlim = c(-0.015, 1.0)) +
  cowplot::draw_label("Variability", x = 0.5, y = 0, vjust = -0.5, angle = 0, size = size_base) + 
  cowplot::draw_label("Coeffiecent of variation", x = 0.0, y = 0.5, vjust = 0.0, angle = 90, size = size_base)

gg_combined <- cowplot::plot_grid(gg_combined, cowplot::get_legend(gg_dummy),
                                  nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_combined, filename = paste0("09-noise-", amplitude, extension),
                     path = "04_Figures/", width = height, height = width * 0.65,
                     units = units, dpi = dpi, overwrite = FALSE)
