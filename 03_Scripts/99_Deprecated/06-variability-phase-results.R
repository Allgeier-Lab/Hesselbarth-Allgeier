##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Results of connectivity experiment for phase variability

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/convert_label.R")

extension <- ".pdf"

#### Load/wrangle simulated data ####

amplitude <- "095"
frequency <- "50"

file_path <- paste0("02_Data/05-variability-phase-", amplitude, "-", frequency, ".rds")

df_cv_prod <- readr::read_rds(file_path) %>% 
  purrr::map_dfr(function(j) {
    dplyr::left_join(x = j$cv, y = j$production, by = c("part", "pop_n", "phase_sd"), 
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

#### CV Density ####

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

#### Phase vs. Synchrony ####

gg_phase_synch <- dplyr::filter(df_cv_prod, measure == "synchrony") %>% 
  ggplot(aes(x = phase_sd, y = value.cv, color = pop_n)) + 
  scale_color_manual(name = "Population size", values = colors_pop) + 
  geom_point(shape = 19, alpha = 0.5) +
  geom_smooth(formula = "y ~ x", se = FALSE, method = "loess") +
  labs(x = expression(italic(phase_sd)), y = "Synchrony") +
  theme_classic(base_size = base_size) +
  theme(legend.position = "bottom")

#### Phase vs. CV (log-linear) ####

list_gg_parts <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
  
  if (part_i == "ag_production") {
    
    labeller_measure <- c(alpha = "Local scale", "gamma" = "Meta-ecosystem scale")
    
  } else {
    
    labeller_measure = c(alpha = "", gamma = "")
    
  }
  
  if (part_i == "ttl_production") {
    
    x_text <- NULL
    
  } else {
    
    x_text <- element_blank()
    
  }
  
  df_regression <- dplyr::filter(df_cv_prod, part == part_i, measure %in% c("alpha", "gamma")) %>% 
    dplyr::group_by(measure, pop_n) %>% 
    dplyr::group_split() %>% 
    purrr::map_dfr(function(i) {
      
      lm_temp <- lm(log(value.cv) ~ phase_sd, data = i)
      
      pred <- coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * i$phase_sd
      
      r.squared <- cor(log(i$value.cv), pred)
      
      p.value <- summary(lm_temp)[["coefficients"]][2, 4]
      
      tibble::tibble(measure = unique(i$measure), pop_n = unique(i$pop_n),
                     phase_sd = i$phase_sd, value.cv = log(i$value.cv), value.cv.pred = pred, 
                     int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                     r.squared = r.squared, p.value = p.value, 
                     x_label = mean(phase_sd), # min(phase_sd) + (max(phase_sd) - min(phase_sd)) * 0.05,
                     y_label = mean(value.cv.pred)) # min(value.cv.pred) + (max(value.cv.pred) - min(value.cv.pred)) * 0.05)
    })
  
  df_label <- dplyr::select(df_regression, -c(phase_sd, value.cv, value.cv.pred)) %>%
    dplyr::group_by(measure, pop_n) %>%
    dplyr::distinct() %>%
    dplyr::mutate(dir = dplyr::case_when(coef < 0.0 ~ "-", coef > 0.0 ~ ""),
                  int = round(int, digits = 2), coef = convert_label(abs(coef), digits = 2),
                  r.squared = round(x = r.squared, digits = 2),
                  p.value = dplyr::case_when(p.value <  0.001 ~  "'***'",  p.value <  0.01 ~  "'**'",
                                             p.value <  0.05 ~  "'*'", p.value >  0.05 ~  "n.s."),
                  label = paste0("R^2==", dir, r.squared, "*';'~", p.value),
                  label = dplyr::case_when(p.value == "n.s." ~ "n.s.", TRUE ~ label))
  
  ggplot(data = df_regression) +
    geom_point(aes(x = phase_sd, y = value.cv, color = pop_n), shape = 19, alpha = 0.25) +
    geom_line(aes(x = phase_sd, y = value.cv.pred, color = pop_n)) +
    geom_label(data = df_label, aes(x = x_label, y = y_label, label = label, group = pop_n),
               parse = TRUE, size = text_size) +
    facet_wrap(. ~ measure, nrow = 1, labeller = labeller(measure = labeller_measure), 
               scales = "free") + 
    scale_color_manual(name = "Population size", values = colors_pop) +
    scale_x_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 5),
                       limits = function(x) c(min(x), max(x)), labels = function(x) round(exp(x), 2)) +
    scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 5),
                       limits = function(x) c(min(x), max(x)), labels = function(x) round(exp(x), 2)) +
    labs(x = ifelse(test = part_i == "ttl_production", yes = expression(italic(phase_sd)), no = ""), 
         y = ifelse(test = part_i == "bg_production", yes = "Coeffiecent of variation", no = "")) +
    theme_classic(base_size = base_size) + 
    theme(legend.position = "none", axis.title = element_text(size = 6.5),
          axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
          strip.background = element_blank(), axis.text.x = x_text)
  
})

gg_move_cv_overall <- cowplot::plot_grid(plotlist = list_gg_parts, nrow = 3, ncol = 1, 
                                         labels = c("ag)", "bg)", "ttl)"), 
                                         label_fontface = "italic", label_size = 10)

gg_dummy <- dplyr::filter(df_cv_prod, part == "ag_production", measure == "alpha") %>%
  ggplot(aes(x = value.move, y = value.cv, color = pop_n)) +
  geom_point(shape = 19) + geom_line() +
  scale_color_manual(name = "Population size", values = colors_pop) +
  # guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = base_size) +
  theme(legend.position = "bottom", legend.spacing.x = unit(1.0, "mm"), legend.spacing.y = unit(-1.0, "mm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 8))

gg_legend <- cowplot::get_legend(gg_dummy)

gg_move_cv_overall <- cowplot::plot_grid(gg_move_cv_overall, gg_legend, ncol = 1,
                                         rel_heights = c(0.95, 0.05))

suppoRt::save_ggplot(plot = gg_move_cv_overall, filename = paste0("06-phase-cv-", amplitude, "-", extension),
                     path = "04_Figures/", width = width, height = height / 2,
                     units = units, dpi = dpi, overwrite = overwrite)

#### CV vs. PP (log-log) ####

df_regression <- dplyr::group_by(df_cv_prod, part, measure, pop_n) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(i) {
    
    lm_temp <- lm(log(value.prod) ~ log(value.cv), data = i)
    
    pred <- coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * log(i$value.cv)
    
    r.squared <- cor(log(i$value.prod), pred)
    
    p.value <- summary(lm_temp)[["coefficients"]][2, 4]
    
    tibble::tibble(part = unique(i$part), measure = unique(i$measure), pop_n = unique(i$pop_n), 
                   value.cv = log(i$value.cv), value.prod = log(i$value.prod), value.prod.pred = pred,
                   int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                   r.squared = r.squared, p.value = p.value)
    
  })

# df_label <- dplyr::select(df_regression, -c(value.cv, value.prod.pred)) %>% 
#   dplyr::group_by(part, measure, pop_n) %>% 
#   dplyr::distinct() %>% 
#   dplyr::ungroup() %>% 
#   dplyr::mutate(dir = dplyr::case_when(coef < 0.0 ~ "-", coef > 0.0 ~ ""),
#                 int = round(int, digits = 2), coef = convert_label(coef, digits = 2),
#                 r.squared = round(x = r.squared, digits = 2), 
#                 p.value = dplyr::case_when(p.value <  0.001 ~  "'***'", p.value <  0.01 ~  "'**'",
#                                            p.value <  0.05 ~  "'*'", p.value >  0.05 ~  "n.s."), 
#                 label = paste0("R^2==", dir, r.squared, "*';'~", p.value),
#                 label = dplyr::case_when(p.value == "n.s." ~ "'---'", TRUE ~ label))

list_gg_parts <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
  
  if (part_i == "ag_production") {
    
    labeller_measure <- c("alpha" = "Local scale", "gamma" = "Meta-ecosystem scale",
                          "beta" = "Portfolio effect")
    
  } else {
    
    labeller_measure <- c("alpha" = "", "gamma" = "", "beta" = "")
    
  }
  
  labeller_part <- c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                     "ttl_production" = "Total")
  
  dplyr::filter(df_regression, part == part_i) %>% 
    ggplot() +
    geom_point(aes(x = value.cv, y = value.prod, color = pop_n), shape = 19, alpha = 0.1) +
    geom_line(data = dplyr::filter(df_regression, part == part_i),
              aes(x = value.cv, y = value.prod.pred, color = pop_n)) +
    # geom_text(data = dplyr::filter(df_label, part == part_i),
    #           aes(x = x, y = y, color = pop_n, label = label), parse = TRUE, size = text_size) +
    facet_wrap(. ~ measure + part, nrow = 3, ncol = 3, scales = "free_x",
               labeller = labeller(measure = labeller_measure, part = labeller_part)) +
    scale_color_manual(name = "Population size", values = colors_pop) +
    scale_x_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 5), 
                       limits = function(x) c(min(x), max(x)), labels = function(x) round(exp(x), 2)) +
    scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 5), 
                       limits = function(x) c(min(x), max(x)), labels = function(x) round(exp(x), 2)) +
    guides(color = guide_legend(override.aes = list(shape = 19, alpha = 1.0), 
                                nrow = 2, byrow = TRUE)) +
    labs(x = ifelse(test = part_i == "ttl_production", yes = "exp(log(Coeffiecent of variation))", no = ""),
         y = ifelse(test = part_i == "bg_production", yes = "exp(log(Total primary production))", no = "")) +
    theme_classic(base_size = base_size) + 
    theme(legend.position = "none", legend.spacing.x = unit(1.0, "mm"), legend.spacing.y = unit(-1.0, "mm"),
          legend.text = element_text(size = 6), legend.title = element_text(size = 8),
          axis.title = element_text(size = 6.5),
          axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA),
          strip.text = element_text(hjust = 0.0), strip.background = element_blank())
})

gg_cv_prod_overall <- cowplot::plot_grid(plotlist = list_gg_parts, nrow = 3)

gg_dummy <- dplyr::filter(df_cv_prod, part == "ag_production", measure == "alpha") %>% 
  ggplot(aes(x = value.move, y = value.cv, color = pop_n)) +
  geom_point(shape = 19) + geom_line() +
  scale_color_manual(name = "Population size", values = colors_pop) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", legend.spacing.x = unit(1.0, "mm"), legend.spacing.y = unit(-1.0, "mm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 8))

gg_legend <- cowplot::get_legend(gg_dummy)

gg_cv_prod_overall <- cowplot::plot_grid(gg_cv_prod_overall, gg_legend, ncol = 1, 
                                         rel_heights = c(0.95, 0.05))

suppoRt::save_ggplot(plot = gg_cv_prod_overall, filename = paste0("06-phase-cv-prod-", amplitude, "-", 
                                                                  frequency, extension),
                     path = "04_Figures/", width = height, height = width,
                     units = units, dpi = dpi, overwrite = overwrite)
