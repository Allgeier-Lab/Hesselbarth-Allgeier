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

file_path <- paste0("02_Data/05-variability-noise-", amplitude, ".rds")

df_results <- readr::read_rds(file_path) %>% 
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

#### Fit regression model ####

df_regression <- dplyr::filter(df_results, measure %in% c("alpha", "gamma")) %>% 
  dplyr::group_by(part, measure, pop_n) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(df_temp) {
    
    df_temp$value.cv <- log(df_temp$value.cv)
    
    df_temp$value.cv <- df_temp$value.cv - mean(df_temp$value.cv)
    
    lm_temp <- lm(value.cv ~ move_meta_sd + noise_sd, data = df_temp)

    broom::tidy(lm_temp) %>% 
      dplyr::mutate(part = unique(df_temp$part), measure = unique(df_temp$measure), 
                    pop_n = unique(df_temp$pop_n), .before = term)
  }) 

#### Setup globals ####

color_palette <- c(alpha = "#642f70", gamma = "#007054")

size_point <- 7.5

#### Create ggplot ###

gg_list <- map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
  
  df_temp <- dplyr::filter(df_regression, part == part_i) %>% 
    dplyr::select(-c(std.error, statistic)) %>% 
    dplyr::mutate(term = dplyr::case_when(term == "(Intercept)" ~ "Intercept", TRUE ~ term), 
                  p.value = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                             p.value < 0.05 ~  "*", p.value >= 0.05 ~ ""))
  
  strip_text <- element_blank()
  
  if (part_i == "ag_production") strip_text <- element_text(hjust = 0)
   
  w <- 0.5
  
  ggplot(data = df_temp) + 
    # Lines
    geom_line(aes(x = pop_n, y = estimate, group = measure, color = measure),
              alpha = 0.5, position = position_dodge(width = w)) +
    
    # Points
    geom_point(aes(x = pop_n, y = estimate, colour = measure),
               size = size_point, shape = 19, position = position_dodge(width = w)) +

    # Text
    geom_text(aes(x = pop_n, y = estimate, label = p.value, group = measure), 
              color = "white", position = position_dodge(width = w)) +
    
    # Facets
    facet_wrap(. ~ term, ncol = 3, scales = "free_y", 
               labeller = labeller(term = c(Intercept = "Intercept", move_meta_sd = "Connectivity", 
                                            noise_sd = "Abiotic variability (noise)"))) +
    
    # Stuff
    scale_color_manual(name = "Scale", values = color_palette) +
    scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4),
                       labels = function(x) round(x, digits = 2), 
                       limits = function(x) c(min(x) - abs(min(x)) * 0.05, max(x) + max(abs(x)) * 0.05)) +
    labs(x = ifelse(test = part_i == "ttl_production", yes = "Population size", no = ""), 
         y = ifelse(test = part_i == "bg_production", yes = "Parameter value", no = ""),
         title = ifelse(test = part_i == "ag_production", yes = "Aboveground", 
                        no = ifelse(test = part_i == "bg_production", yes = "Belowground", no = "Total"))) +
    theme_classic(base_size = 10.0) + 
    theme(legend.position = "none", plot.title = element_text(size = 8.0), 
          strip.background = element_blank(), strip.text = strip_text)
  })

gg_dummy <- dplyr::select(df_regression, -c(std.error, statistic)) %>% 
  dplyr::mutate(term = dplyr::case_when(term == "(Intercept)" ~ "Intercept", TRUE ~ term)) %>% 
  ggplot(aes(x = pop_n, y = estimate, color = measure)) + 
  geom_point(size = size_point / 2) +
  scale_color_manual(name = "Scale", values = color_palette, 
                     labels = c(alpha = "Local scale", gamma = "Meta-ecosystem scale")) +
  theme(legend.position = "bottom")

gg_lm_coef <- cowplot::plot_grid(cowplot::plot_grid(plotlist = gg_list, ncol = 1), 
                                 cowplot::get_legend(gg_dummy), rel_heights = c(0.95, 0.05), 
                                 nrow = 2)

suppoRt::save_ggplot(plot = gg_lm_coef, filename = paste0("07-noise-coef-", amplitude, extension),
                     path = "04_Figures/", width = height, height = width,
                     units = units, dpi = dpi, overwrite = overwrite)
