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

file_path <- paste0("02_Data/05-variability-phase-", amplitude, ".rds")

df_results <- readr::read_rds(file_path) %>% 
  purrr::map_dfr(function(j) {
    dplyr::left_join(x = j$cv, y = j$production, by = c("part", "pop_n", "move_meta_sd", "phase_sd"), 
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

size_base <- 8.0

colors_pop <- c("8" = "#ebcc60", "16" = "#459b75", "32" = "#374a98", "64" = "#913f98")

#### Density ####

gg_density <- dplyr::filter(df_results, measure %in% c("alpha", "gamma")) %>% 
  ggplot(aes(x = value.cv, y = ..density.., fill = pop_n, color = pop_n)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(. ~ part + measure, scales = "free", ncol = 2, nrow = 3, 
             labeller = labeller(part = c(ag_production = "Aboveground", bg_production = "Belowground", 
                                          ttl_production = "Total production"), 
                                 measure = c(alpha = "Local scale", gamma = "Meta-ecosystem scale"))) +
  scale_fill_manual(name = "Population size", values = colors_pop) + 
  scale_color_manual(name = "Population size", values = colors_pop) + 
  labs(x = "Coeffienct of variation", y = "Density") +
  theme_classic(base_size = size_base) +
  theme(legend.position = "bottom", axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0))

#### Variability vs. CV (log-linear) ####

list_gg_parts <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
  
  purrr::map(c("alpha", "gamma"), function(measure_i) {
    
    list_gg_pop <- purrr::map(c(8, 16, 32, 64), function(pop_i) {
      
      dplyr::filter(df_results, part == part_i, measure == measure_i, pop_n == pop_i) %>% 
        dplyr::mutate(move_meta_sd = cut(move_meta_sd, breaks = seq(0.1, 1, 0.1),labels = FALSE), 
                      phase_sd = cut(phase_sd, breaks = seq(0.1, 1, 0.1), labels = FALSE)) %>% 
        dplyr::group_by(move_meta_sd, phase_sd) %>%
        dplyr::summarise(value.cv = mean(value.cv), .groups = "drop") %>%
        ggplot(aes(x = phase_sd, y = move_meta_sd, fill = value.cv)) +
        geom_tile() + 
        scale_fill_gradientn(colors = MetBrewer::met.brewer("Demuth", n = 255, type = "continuous")) +
        scale_x_continuous(breaks = c(1, 5, 9), labels = c(0.1, 0.5, 1.0)) +
        scale_y_continuous(breaks = c(1, 5, 9), labels = c(0.1, 0.5, 1.0)) +
        labs(title = paste0("Population size: ", pop_i)) +
        theme_classic(base_size = size_base) +
        theme(legend.position = "none", plot.title = element_text(size = 6.5),
              axis.title = element_blank(), axis.line = element_blank(), 
              panel.border = element_rect(size = 0.5, fill = NA))
      
    })
    
    cowplot::plot_grid(plotlist = list_gg_pop, ncol = 2, nrow = 2) + 
      cowplot::draw_label(label = expression(italic(phase_sd)), x = 0.5, y = 0, 
                          vjust = 0.5, angle = 0, size = size_base) +
      cowplot::draw_label(label = expression(italic(move_meta_sd)), x = 0, y = 0.5, 
                          vjust = 0.0, angle = 90, size = size_base) + 
      theme(plot.margin = unit(c(t = 0.5, r = 0.5, b = 0.5, l = 0.5), "cm"))
    
  })}) %>% 
  purrr::flatten()

gg_move_cv_overall <- cowplot::plot_grid(plotlist = list_gg_parts, nrow = 3, ncol = 2, 
                                         labels = "auto", label_fontface = "italic")

suppoRt::save_ggplot(plot = gg_move_cv_overall, filename = paste0("07-phase-cv-", amplitude, extension),
                     path = "04_Figures/", width = width, height = height,
                     units = units, dpi = dpi, overwrite = overwrite)

