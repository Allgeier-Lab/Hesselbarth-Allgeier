##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Results of connectivity experiment for movement and phase variabilities

#### Load setup ####

source("05_Various/setup.R")

extension <- ".pdf"

#### Load/wrangle simulated data ####

amplitude <- "095"

file_path <- paste0("02_Data/05-variability-noise-", amplitude, ".rds")

df_cv_prod <- readr::read_rds(file_path) %>% 
  purrr::map_dfr(function(j) {
    dplyr::left_join(x = j$cv, y = j$production, by = c("part", "pop_n", "move_meta_sd", "noise_sd"), 
                     suffix = c(".cv", ".prod")) %>% 
      dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))
  }) %>% 
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass",
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                pop_n = factor(as.numeric(pop_n), ordered = TRUE)) %>% 
  dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), 
                measure %in% c("alpha", "gamma", "beta")) %>% 
  tibble::tibble()

#### Setup globals ####

base_size <- 8.0

text_size <- 2.0

colors_pop <- c("8" = "#ebcc60", "16" = "#459b75", "32" = "#374a98", "64" = "#913f98")

#### Density ####

gg_density <- dplyr::filter(df_cv_prod, measure %in% c("alpha", "gamma")) %>% 
  ggplot(aes(x = value.cv, y = ..density.., fill = pop_n, color = pop_n)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(. ~ part + measure, scales = "free", ncol = 2, nrow = 3, 
             labeller = labeller(part = c(ag_production = "Aboveground", bg_production = "Belowground", 
                                          ttl_production = "Total production"), 
                                 measure = c(alpha = "Local scale", gamma = "Meta-ecosystem scale"))) +
  scale_fill_manual(name = "Population size", values = colors_pop) + 
  scale_color_manual(name = "Population size", values = colors_pop) + 
  labs(x = "Coeffienct of variation", y = "Density") +
  theme_classic(base_size = base_size) +
  theme(legend.position = "bottom", axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0))

#### Variability vs. CV (log-linear) ####

list_gg_parts <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
  
  list_gg_measures <- purrr::map(c("alpha", "gamma"), function(measure_i) {
    
    list_gg_pop <- purrr::map(c(8, 16, 32, 64), function(pop_i) {
      
      dplyr::filter(df_cv_prod, part == part_i, measure == measure_i, pop_n == pop_i) %>% 
        dplyr::mutate(move_meta_sd = cut(move_meta_sd, breaks = seq(0, 1, 0.1)), 
                      noise_sd = cut(noise_sd, breaks = seq(0, 1, 0.1))) %>% 
        dplyr::group_by(move_meta_sd, noise_sd) %>%
        dplyr::summarise(value.cv = mean(value.cv), .groups = "drop") %>%
        ggplot(aes(x = noise_sd, y = move_meta_sd, fill = value.cv)) +
        geom_tile() + 
        scale_fill_gradientn(colors = MetBrewer::met.brewer("Demuth", n = 255, type = "continuous")) +
        # scale_x_discrete(breaks = c(0.1, 0.5, 1.0)) + scale_y_continuous(breaks = c(0.1, 0.5, 1.0)) +
        labs(title = paste0("Population size: ", pop_i), 
             x = expression(italic(noise_sd)), y = expression(italic(move_meta_sd))) +
        theme_classic(base_size = base_size) +
        theme(legend.position = "none", axis.title = element_text(size = 6.5),
              axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
              plot.title = element_text(size = 6.5))
    })
    
    cowplot::plot_grid(plotlist = list_gg_pop, ncol = 2, nrow = 2)
    
  })
  
  cowplot::plot_grid(plotlist = list_gg_measures, nrow = 1, ncol = 2)
  
})

gg_move_cv_overall <- cowplot::plot_grid(plotlist = list_gg_parts, nrow = 3, ncol = 1, 
                                         labels = c("ag)", "bg)", "ttl)"), 
                                         label_fontface = "italic", label_size = 10)

suppoRt::save_ggplot(plot = gg_move_cv_overall, filename = paste0("06-noise-cv-", amplitude, extension),
                     path = "04_Figures/", width = width, height = height,
                     units = units, dpi = dpi, overwrite = overwrite)

