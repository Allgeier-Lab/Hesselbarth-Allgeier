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

df_cv_prod <- readr::read_rds(file_path) %>% 
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

color_palette <- c(Intercept = "#e96052", Connectivity = "#70bdd6", Phase = "#19436d")

#### Fit regression model ####

df_regression <- dplyr::filter(df_cv_prod, measure %in% c("alpha", "gamma")) %>% 
  dplyr::group_by(part, measure, pop_n) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(i) {
    
    lm_temp <- lm(log(value.cv) ~ move_meta_sd + phase_sd, data = i)
    
    broom::tidy(lm_temp) %>% 
      dplyr::mutate(part = unique(i$part), measure = unique(i$measure), pop_n = unique(i$pop_n), 
                    .before = term)
  }) 

#### Create ggplot ###

gg_list <- map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
  
  df_temp <- dplyr::filter(df_regression, part == part_i) %>% 
    dplyr::select(-c(std.error, statistic)) %>% 
    dplyr::mutate(term = dplyr::case_when(term == "(Intercept)" ~ "Intercept", 
                                          TRUE ~ term), 
                  p.value = dplyr::case_when(p.value < 0.001 ~ "***", 
                                             p.value < 0.01 ~ "**", 
                                             p.value < 0.05 ~  "*",
                                             p.value >= 0.05 ~ "")) %>% 
    tidyr::pivot_wider(names_from = term, values_from = c(estimate, p.value))
  
  labels_y_right <- round(seq(from = min(c(df_temp$estimate_move_meta_sd, 
                                           df_temp$estimate_phase_sd)), 
                              to = max(c(df_temp$estimate_move_meta_sd, 
                                         df_temp$estimate_phase_sd)), 
                              length.out = 4), digits = 2)
  
  scales_y <- seq(from = min(df_temp$estimate_Intercept), to = max(df_temp$estimate_Intercept), 
                  length.out = 4)
  
  df_temp$estimate_move_meta_sd <- scales::rescale(df_temp$estimate_move_meta_sd, 
                                                   to = range(df_temp$estimate_Intercept))
  
  df_temp$estimate_phase_sd <- scales::rescale(df_temp$estimate_phase_sd, 
                                               to = range(df_temp$estimate_Intercept))
  
  w <- 0.25
  
  ggplot(data = df_temp) + 
    # Lines
    geom_line(aes(x = pop_n, y = estimate_Intercept, group = measure, linetype = measure, color = "Intercept"),
              alpha = 0.5, position = position_dodge(width = w)) +
    geom_line(aes(x = pop_n, y = estimate_move_meta_sd, group = measure, linetype = measure, color = "Connectivity"),
              alpha = 0.5, position = position_dodge(width = w)) +
    geom_line(aes(x = pop_n, y = estimate_phase_sd, group = measure, linetype = measure, color = "Phase"),
              alpha = 0.5, position = position_dodge(width = w)) +
    
    # Points
    geom_point(aes(x = pop_n, y = estimate_Intercept, shape = measure, color = "Intercept"), 
               size = 3.0, position = position_dodge(width = w)) +
    geom_point(aes(x = pop_n, y = estimate_move_meta_sd, shape = measure, color = "Connectivity"), 
               size = 3.0, position = position_dodge(width = w)) +
    geom_point(aes(x = pop_n, y = estimate_phase_sd, shape = measure, color = "Phase"), 
               size = 3.0, position = position_dodge(width = w)) +
    
    # Text
    geom_text(aes(x = pop_n, y = estimate_Intercept, label = p.value_Intercept, 
                  group = measure, color = "Intercept"), 
              position = position_dodge(width = w)) +
    geom_text(aes(x = pop_n, y = estimate_move_meta_sd, label = p.value_move_meta_sd, 
                  group = measure, color = "Connectivity"), 
              position = position_dodge(width = w)) +
    geom_text(aes(x = pop_n, y = estimate_phase_sd, label = p.value_phase_sd, 
                  group = measure, color = "Phase"), 
              position = position_dodge(width = w)) +
    
    # Stuff
    scale_color_manual(name = "Regression model", values = color_palette) +
    scale_shape_manual(name = "Scale", values = c(1, 2)) +
    scale_linetype_manual(name = "Scale", values = c(1, 2)) +
    scale_y_continuous(name = ifelse(test = part_i == "bg_production", yes = "Intercept", no = ""), 
                       breaks = scales_y, labels = function(x) round(x, digits = 2),
                       sec.axis = sec_axis(name = ifelse(test = part_i == "bg_production", yes = "Slope", no = ""),
                                           trans = ~ . * 1, breaks = scales_y, labels = labels_y_right)) +
    labs(x = ifelse(test = part_i == "ttl_production", yes = "Population size", no = ""), 
         title = ifelse(test = part_i == "ag_production", yes = "Aboveground", 
                        no = ifelse(test = part_i == "bg_production", 
                                    yes = "Belowground", no = "Total"))) +
    theme_classic(base_size = 10.0) + 
    theme(legend.position = "none",
          axis.title.y.left = element_text(color = color_palette[["Intercept"]]),
          axis.text.y.left = element_text(color = color_palette[["Intercept"]]), 
          axis.ticks.y.left = element_line(color =  color_palette[["Intercept"]]),
          axis.line.y.left = element_line(color =  color_palette[["Intercept"]]), 
          axis.title.y.right = element_text(color = color_palette["Connectivity"]), 
          axis.text.y.right = element_text(color = color_palette["Connectivity"]),
          axis.ticks.y.right = element_line(color =  color_palette[["Connectivity"]]),
          axis.line.y.right = element_line(color =  color_palette[["Connectivity"]]), 
          plot.title = element_text(size = 8.0))
})

gg_lm_coef <- cowplot::plot_grid(plotlist = gg_list, ncol = 1)

gg_dummy <- dplyr::select(df_regression, -c(std.error, statistic)) %>% 
  dplyr::mutate(term = dplyr::case_when(term == "(Intercept)" ~ "Intercept", 
                                        TRUE ~ term)) %>% 
  tidyr::pivot_wider(names_from = term, values_from = c(estimate, p.value)) %>% 
  ggplot(aes(x = pop_n)) + 
  geom_point(aes(y = estimate_Intercept, color = "Intercept", 
                 shape = factor(measure, labels = c("Local", "Meta-ecosystem")))) +
  geom_point(aes(y = estimate_move_meta_sd, color = "Connectivity", 
                 shape = factor(measure, labels = c("Local", "Meta-ecosystem")))) + 
  geom_point(aes(y = estimate_phase_sd, color = "Phase", 
                 shape = factor(measure, labels = c("Local", "Meta-ecosystem")))) + 
  scale_shape_manual(name = "Scale", values = c(1, 2)) +
  scale_color_manual(name = "Coeffiecent", values = color_palette) +
  theme(legend.position = "bottom")

gg_lm_coef <- cowplot::plot_grid(gg_lm_coef, cowplot::get_legend(gg_dummy), rel_heights = c(0.95, 0.05), 
                                 nrow = 2)

suppoRt::save_ggplot(plot = gg_lm_coef, filename = paste0("07-phase-coef-", amplitude, extension),
                     path = "04_Figures/", width = width, height = height,
                     units = units, dpi = dpi, overwrite = overwrite)
