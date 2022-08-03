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

#### Fit regression model ####

df_regression <- dplyr::filter(df_results, measure %in% c("alpha", "gamma")) %>% 
  dplyr::group_by(part, measure, pop_n) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(df_temp) {
    
    df_temp$value.cv <- log(df_temp$value.cv) - mean(log(df_temp$value.cv))
    
    lm_temp <- lm(value.cv ~ log(move_meta_sd) + log(phase_sd), data = df_temp)
    
    broom::tidy(lm_temp) %>% 
      dplyr::mutate(part = unique(df_temp$part), measure = unique(df_temp$measure), 
                    pop_n = unique(df_temp$pop_n), .before = term)
  }) 

#### Setup globals ####

color_palette <- c("Intercept" = "#ed968b", "log(move_meta_sd)" = "#88a0dc", "log(phase_sd)" = "#f9d14a")

size_point <- 7.5

#### Create ggplot ###

gg_list <- purrr::map(c("alpha", "gamma"), function(measure_i) {
  
  df_temp <- dplyr::filter(df_regression, measure == measure_i) %>% 
    dplyr::select(-c(std.error, statistic)) %>% 
    dplyr::mutate(term = dplyr::case_when(term == "(Intercept)" ~ "Intercept", TRUE ~ term), 
                  p.value = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                             p.value < 0.05 ~  "*", p.value >= 0.05 ~ ""))
  
  w <- 0.5
  
  ggplot(data = df_temp) + 
    
    # zero line
    geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
    
    # Lines
    geom_line(aes(x = pop_n, y = estimate, group = term, color = term),
              alpha = 0.5, position = position_dodge(width = w)) +
    
    # Points
    geom_point(aes(x = pop_n, y = estimate, colour = term),
               size = size_point, shape = 19, position = position_dodge(width = w)) +
    
    # Text
    geom_text(aes(x = pop_n, y = estimate, label = p.value, group = term),
              color = "white", position = position_dodge(width = w)) +
    
    # Facets
    facet_wrap(. ~ part, nrow = 3, scales = "free", 
               labeller = labeller(part = c("ag_production" = ifelse(measure_i == "alpha", yes = "Aboveground", no = ""), 
                                            "bg_production" = ifelse(measure_i == "alpha", yes = "Belowground", no = ""), 
                                            "ttl_production" = ifelse(measure_i == "alpha", yes = "Total", no = "")))) +
    
    # Stuff
    scale_color_manual(name = "Scale", values = color_palette) +
    scale_y_continuous(limits = range(df_regression$estimate), 
                       breaks = seq(min(df_regression$estimate), max(df_regression$estimate), length.out = 4), 
                       labels = function(x) round(x, digits = 2)) + 
    labs(title = ifelse(measure_i == "alpha", yes = "Local scale", no = "Meta-ecosystem scale"), 
         x = "Population size", y = ifelse(measure_i == "alpha", yes = "Parameter value", no = "")) +
    theme_classic(base_size = 10.0) + 
    theme(legend.position = "none", plot.title = element_text(size = 8.0), 
          strip.background = element_blank(), strip.text = element_text(hjust = 0))
})

gg_dummy <- dplyr::select(df_regression, -c(std.error, statistic)) %>% 
  dplyr::mutate(term = dplyr::case_when(term == "(Intercept)" ~ "Intercept", TRUE ~ term)) %>% 
  ggplot(aes(x = pop_n, y = estimate, color = term)) + 
  geom_point(size = size_point / 2) +
  scale_color_manual(name = "Scale", values = color_palette, 
                     labels = c("Intercept" = "Intercept", 
                                "log(move_meta_sd)" = "Biotic variability", 
                                "log(phase_sd)" = "Abiotic variability")) +
  theme(legend.position = "bottom")

gg_lm_coef <- cowplot::plot_grid(cowplot::plot_grid(plotlist = gg_list, ncol = 2), 
                                 cowplot::get_legend(gg_dummy), rel_heights = c(0.95, 0.05), 
                                 nrow = 2)

suppoRt::save_ggplot(plot = gg_lm_coef, filename = paste0("07-phase-coef-", amplitude, extension),
                     path = "04_Figures/", width = height, height = width,
                     units = units, dpi = dpi, overwrite = overwrite)
