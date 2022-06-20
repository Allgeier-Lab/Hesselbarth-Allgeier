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
    
    pred <- coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * i$move_meta_sd + coef(lm_temp)[[3]] * i$phase_sd
    
    r.squared <- cor(log(i$value.cv), pred)
    
    p.value <- summary(lm_temp)[["coefficients"]][2, 4]
    
    tibble::tibble(
      part = unique(i$part), measure = unique(i$measure), pop_n = unique(i$pop_n), 
      int = lm_temp$coefficients[[1]], coef_biotic = lm_temp$coefficients[[2]],
      coef_abiotic = lm_temp$coefficients[[3]], r.squared = r.squared, p.value = p.value)
  })

#### Create ggplot ###

gg_list <- map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
  
  df_temp <- dplyr::filter(df_regression, part == part_i)
  
  labels_y_right <- round(seq(from = min(c(df_temp$coef_biotic, df_temp$coef_abiotic)), 
                              to = max(c(df_temp$coef_biotic, df_temp$coef_abiotic)), 
                              length.out = 4), digits = 2)
  
  scales_y <- seq(from = min(df_temp$int), to = max(df_temp$int), length.out = 4)
  
  df_temp$coef_biotic <- scales::rescale(df_temp$coef_biotic, to = range(df_temp$int))
  
  df_temp$coef_abiotic <- scales::rescale(df_temp$coef_abiotic, to = range(df_temp$int))
  
  ggplot(data = df_temp) + 
    geom_point(aes(x = pop_n, y = int, shape = measure, color = "Intercept"), size = 3.0) +
    geom_line(aes(x = pop_n, y = int, group = measure, linetype = measure, color = "Intercept"), alpha = 0.5) +
    geom_point(aes(x = pop_n, y = coef_biotic, shape = measure, color = "Connectivity"), size = 3.0) +
    geom_line(aes(x = pop_n, y = coef_biotic, group = measure, linetype = measure, color = "Connectivity"), alpha = 0.5) +
    geom_point(aes(x = pop_n, y = coef_abiotic, shape = measure, color = "Phase"), size = 3.0) +
    geom_line(aes(x = pop_n, y = coef_abiotic, group = measure, linetype = measure, color = "Phase"), alpha = 0.5) +
    scale_color_manual(name = "Regression model", values = color_palette) +
    scale_shape_manual(name = "Scale", values = c(19, 2)) +
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

gg_dummy <- ggplot(data = df_regression, aes(x = pop_n)) + 
  geom_point(aes(y = int, color = "Intercept", shape = factor(measure, labels = c("Local", "Meta-ecosystem")))) +
  geom_point(aes(y = coef_biotic, color = "Connectivity", shape = factor(measure, labels = c("Local", "Meta-ecosystem")))) + 
  geom_point(aes(y = coef_abiotic, color = "Phase", shape = factor(measure, labels = c("Local", "Meta-ecosystem")))) + 
  scale_shape_manual(name = "Scale", values = c(19, 2)) +
  scale_color_manual(name = "Coeffiecent", values = color_palette) +
  theme(legend.position = "bottom")

gg_lm_coef <- cowplot::plot_grid(gg_lm_coef, cowplot::get_legend(gg_dummy), rel_heights = c(0.95, 0.05), 
                                 nrow = 2)

suppoRt::save_ggplot(plot = gg_lm_coef, filename = paste0("07-phase-coef-", amplitude, extension),
                     path = "04_Figures/", width = width, height = height,
                     units = units, dpi = dpi, overwrite = overwrite)

