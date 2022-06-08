##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Results of connectivity experiment for low amplitude treatment 

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/convert_label.R")

extension <- ".pdf"

#### Load/wrangle simulated data ####

amplitude <- "005"
frequency <- "50"

file_path <- paste0("02_Data/05-move-variability-", amplitude, "-", frequency, ".rds")

df_cv_prod <- readr::read_rds(file_path) %>% 
  purrr::map_dfr(function(j) {
    dplyr::left_join(x = j$cv, y = j$production, by = c("part", "pop_n", "move_meta_mean", 
                                                        "move_meta_sd"), 
                     suffix = c(".cv", ".prod")) %>% 
      dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))
    }) %>% 
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass",
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                pop_n = factor(as.numeric(pop_n), ordered = TRUE),
                value.cv.log = log(value.cv), value.prod.log = log(value.prod), 
                value.move.log = log(value.move)) %>% 
  dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), 
                measure %in% c("alpha", "gamma", "beta")) %>% 
  tibble::tibble()
 
#### Setup globals ####

base_size <- 8.0

text_size <- 2.0

colors_pop <- c("8" = "#ebcc60", "16" = "#459b75", "32" = "#374a98", "64" = "#913f98")

# #### Total production ####

gg_pop_prod_agbg <- dplyr::filter(df_cv_prod, part %in% c("ag_production", "bg_production"),
                            measure == "alpha") %>%
  ggplot(aes(x = pop_n, y = value.prod, color = pop_n)) +
  # geom_jitter(shape = 19, alpha = 0.05) +
  geom_boxplot(position = position_dodge(width = 0.5), fill = NA, outlier.color = NA) +
  facet_wrap(. ~ part, scales = "free_y", nrow = 2, ncol = 1,
             labeller = labeller(part = c(ag_production = "Aboveground",
                                          bg_production = "Belowground "))) +
  scale_color_manual(values = colors_pop) +
  scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), 
                                              length.out = 4), 
                     limits = function(x) c(min(x), max(x)), 
                     labels = function(x) format(x, digits = 2, scientific = TRUE)) +
  labs(x = "Population size", y = "Total primary production") +
  theme_classic(base_size = base_size) +
  theme(legend.position = "none", axis.line = element_blank(),
        panel.border = element_rect(size = 0.5,fill = NA),
        strip.text = element_text(hjust = 0), strip.background = element_blank(), 
        plot.margin = margin(t = 5.5, r = 0.0, b = 5.5, l = 5.5, unit = "pt")) # t r b l

gg_pop_prod_ttl <- dplyr::filter(df_cv_prod, part == "ttl_production", measure == "alpha") %>%
  ggplot(aes(x = pop_n, y = value.prod, color = pop_n)) +
  # geom_jitter(shape = 19, alpha = 0.05) +
  geom_boxplot(position = position_dodge(width = 0.5), fill = NA, outlier.color = NA) +
  facet_wrap(. ~ part, scales = "free_y", nrow = 1, ncol = 1,
             labeller = labeller(part = c(ttl_production = "Total"))) +
  scale_color_manual(values = colors_pop) +
  scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), 
                                              length.out = 4), 
                     limits = function(x) c(min(x), max(x)), 
                     labels = function(x) format(x, digits = 2, scientific = TRUE)) +
  labs(x = "Population size", y = "") +
  theme_classic(base_size = base_size) +
  theme(legend.position = "none", axis.line = element_blank(),
        panel.border = element_rect(size = 0.5,fill = NA),
        strip.text = element_text(hjust = 0), strip.background = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0.0, unit = "pt")) # t r b l

gg_pop_prod_overall <- cowplot::plot_grid(gg_pop_prod_agbg, gg_pop_prod_ttl, ncol = 2)

suppoRt::save_ggplot(plot = gg_pop_prod_overall, filename = paste0("06-pop-prod-", amplitude, "-", 
                                                                   frequency, extension),
                     path = "04_Figures/", width = width, height = height / 2,
                     units = units, dpi = dpi, overwrite = overwrite)

#### Mobility vs. CV (log-linear) ####

list_gg_parts <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
  
  list_gg_measures <- purrr::map(c("alpha", "gamma"), function(measure_i) {
    
    gg_move_move <- dplyr::filter(df_cv_prod, part == part_i, measure == measure_i) %>% 
      dplyr::group_by(pop_n, move_meta_mean, move_meta_sd) %>%
      dplyr::summarise(value.cv = mean(value.cv), .groups = "drop") %>%
      dplyr::group_by(pop_n) %>% 
      dplyr::mutate(value.cv = scales::rescale(value.cv, to = c(0, 1))) %>% 
      dplyr::ungroup() %>% 
      ggplot(aes(x = move_meta_mean, y = move_meta_sd, fill = value.cv)) +
      geom_tile() +
      facet_wrap(. ~ pop_n, ncol = 4, labeller = labeller(pop_n = c("8" = "Populations size: 8", 
                                                                    "16" = "Populations size: 16", 
                                                                    "32" = "Populations size: 32", 
                                                                    "64" = "Populations size: 64"))) +
      # scale_fill_viridis_c(option = "A", na.value = "grey") +
      scale_fill_gradientn(colors = MetBrewer::met.brewer("Demuth", n = 255, type = "continuous")) +
      scale_x_continuous(breaks = c(0.1, 0.5, 1.0)) + scale_y_continuous(breaks = c(0.1, 0.5, 1.0)) +
      labs(x = expression(italic(move_meta_mean)), 
           y = ifelse(test = measure_i == "alpha", yes = expression(italic(move_meta_sd)), no = "")) +
      coord_equal() +
      theme_classic(base_size = base_size) +
      theme(legend.position = "none",  axis.title = element_text(size = 6.5),
            axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
            strip.background = element_blank(), strip.text = element_text(hjust = 0), 
            plot.margin = margin(t = 5.5, r = 5.5, b = 0.0, l = 5.5, unit = "pt")) # t r b l
    
    df_regression <- dplyr::filter(df_cv_prod, part == part_i, measure == measure_i) %>% 
      dplyr::group_by(pop_n) %>% 
      dplyr::group_split() %>% 
      purrr::map_dfr(function(i) {
        
        lm_temp <- lm(value.cv.log ~ value.move, data = i)
        
        pred <- coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * i$value.move
        
        r.squared <- cor(i$value.cv.log, pred)
        
        p.value <- summary(lm_temp)[["coefficients"]][2, 4]
        
        tibble::tibble(pop_n = unique(i$pop_n), value.move = i$value.move, 
                       value.cv.log.pred = pred, 
                       int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                       r.squared = r.squared, p.value = p.value)
        
    })
    
    df_label <- dplyr::select(df_regression, -value.move, -value.cv.log.pred) %>%
      dplyr::group_by(pop_n) %>%
      dplyr::distinct() %>%
      dplyr::mutate(dir = dplyr::case_when(coef < 0.0 ~ "-", coef > 0.0 ~ "+"),
                    int = round(int, digits = 2), coef = convert_label(abs(coef), digits = 2),
                    r.squared = round(x = r.squared, digits = 2),
                    p.value = dplyr::case_when(p.value <  0.001 ~  "'***'",  p.value <  0.01 ~  "'**'",
                                               p.value <  0.05 ~  "'*'", p.value >  0.05 ~  "n.s."),
                    # label = paste0("italic(y)==", int, dir, coef, "*italic(x)*';'~R^2==", r.squared, "*';'~", p.value),
                    label = paste0("R^2==", dir, r.squared, "*';'~", p.value),
                    label = dplyr::case_when(p.value == "n.s." ~ "'---'", TRUE ~ label)) %>% 
      dplyr::bind_cols(dplyr::filter(df_cv_prod, part == part_i, measure == measure_i) %>% 
                         dplyr::summarise(x = min(value.move) + (max(value.move) - min(value.move)) * c(0.65, 0.9, 0.65, 0.9), 
                                          y = min(value.cv.log) + (max(value.cv.log) - min(value.cv.log)) * c(0.95, 0.95, 0.85, 0.85)))
    
    gg_move_cv <- dplyr::filter(df_cv_prod, part == part_i, measure == measure_i) %>% 
      ggplot(aes(x = value.move, y = value.cv.log, color = pop_n)) +
      geom_point(shape = 1, alpha = 0.1) +
      geom_line(data = df_regression, aes(y = value.cv.log.pred)) +
      geom_text(data = df_label, aes(x = x, y = y, label = label), parse = TRUE, size = text_size) +
      scale_color_manual(name = "Population size", values = colors_pop) +
      scale_x_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 5),
                         limits = function(x) c(min(x), max(x)), labels = function(x) round(x, 1)) +
      scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 5),
                         limits = function(x) c(min(x), max(x)), labels = function(x) round(x, 1)) +
      labs(x = "Mean cross-ecosystem movement", 
           y = ifelse(measure_i == "alpha", yes = "log(Coefficient of variation)", no = "")) +
      theme_classic(base_size = base_size) + 
      theme(legend.position = "none", axis.title = element_text(size = 6.5),
            axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA),
            plot.margin = margin(t = 0.0, r = 5.5, b = 5.5, l = 5.5, unit = "pt")) # t r b l)
    
    cowplot::plot_grid(gg_move_move, gg_move_cv, nrow = 2, ncol = 1)

  })
  
  labels <- ifelse(test = c(part_i == "ag_production", part_i == "ag_production"), 
                   yes = c("Local scale", "Meta-ecosystem scale"), no = c("", ""))
  
  cowplot::plot_grid(plotlist = list_gg_measures,  nrow = 1, ncol = 2, 
                     labels = labels, hjust = c(-2.65, -1.0), 
                     label_fontface = "italic", label_size = 10)
  
})

gg_move_cv_overall <- cowplot::plot_grid(plotlist = list_gg_parts, nrow = 3, ncol = 1, 
                                         labels = c("ag)", "bg)", "ttl)"), 
                                         label_fontface = "italic", label_size = 10)

gg_dummy <- dplyr::filter(df_cv_prod, part == "ag_production", measure == "alpha") %>% 
  ggplot(aes(x = value.move, y = value.cv, color = pop_n)) +
  geom_point(shape = 19) + geom_line() +
  scale_color_manual(name = "Population size", values = colors_pop) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", legend.spacing.x = unit(1.0, "mm"), legend.spacing.y = unit(-1.0, "mm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 8))

gg_legend <- cowplot::get_legend(gg_dummy)

gg_move_cv_overall <- cowplot::plot_grid(gg_move_cv_overall, gg_legend, ncol = 1, 
                                         rel_heights = c(0.95, 0.05))

suppoRt::save_ggplot(plot = gg_move_cv_overall, filename = paste0("06-connect-cv-", amplitude, "-", 
                                                                  frequency, extension),
                     path = "04_Figures/", width = width, height = height,
                     units = units, dpi = dpi, overwrite = overwrite)

#### CV vs. PP (log-log) ####

df_regression <- dplyr::group_by(df_cv_prod, part, measure, pop_n) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(i) {
    
    lm_temp <- lm(value.prod.log ~ value.cv.log, data = i)
    
    pred <- coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * i$value.cv.log
    
    r.squared <- cor(i$value.prod.log, pred)
    
    p.value <- summary(lm_temp)[["coefficients"]][2, 4]
    
    tibble::tibble(part = unique(i$part), measure = unique(i$measure), pop_n = unique(i$pop_n), 
                   value.cv.log = i$value.cv.log, value.prod.log.pred = pred,
                   int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                   r.squared = r.squared, p.value = p.value)
    
  })

df_label <- dplyr::select(df_regression, -value.cv.log, -value.prod.log.pred) %>% 
  dplyr::group_by(part, measure, pop_n) %>% 
  dplyr::distinct() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(dir = dplyr::case_when(coef < 0.0 ~ "-", coef > 0.0 ~ ""),
                int = round(int, digits = 2), coef = convert_label(coef, digits = 2),
                r.squared = round(x = r.squared, digits = 2), 
                p.value = dplyr::case_when(p.value <  0.001 ~  "'***'", p.value <  0.01 ~  "'**'",
                                           p.value <  0.05 ~  "'*'", p.value >  0.05 ~  "n.s."), 
                label = paste0("R^2==", dir, r.squared, "*';'~", p.value),
                label = dplyr::case_when(p.value == "n.s." ~ "'---'", TRUE ~ label)) %>% 
  dplyr::left_join(x = ., y = dplyr::group_by(df_cv_prod, part, measure) %>% 
                     dplyr::summarise(pop_n = factor(c(8, 16, 32, 64), ordered = TRUE), 
                                      x = min(value.cv.log) + 
                                        (max(value.cv.log) - min(value.cv.log)) * c(0.7, 0.9, 0.7, 0.9), 
                                      y = min(value.prod.log) + 
                                        (max(value.prod.log) - min(value.prod.log)) * c(0.9, 0.9, 0.8, 0.8), 
                                      .groups = "drop"), 
                   by = c("part", "measure", "pop_n"))

list_gg_parts <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
  
  
  if (part_i == "ag_production") {
    
    labeller_measure <- c("alpha" = "Local scale", "gamma" = "Meta-ecosystem scale",
                          "beta" = "Portfolio effect")
    
  } else {
    
    labeller_measure <- c("alpha" = "", "gamma" = "", "beta" = "")
    
  }
  
  labeller_part <- c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                     "ttl_production" = "Total")
  
  dplyr::filter(df_cv_prod, part == part_i) %>% 
    ggplot(aes(x = value.cv.log, y = value.prod.log, color = pop_n)) +
    geom_point(shape = 1, alpha = 0.1) +
    geom_line(data = dplyr::filter(df_regression, part == part_i), 
              aes(x = value.cv.log, y = value.prod.log.pred, color = pop_n)) +
    geom_text(data = dplyr::filter(df_label, part == part_i), 
              aes(x = x, y = y, color = pop_n, label = label), parse = TRUE, size = text_size) +
    facet_wrap(. ~ measure + part, nrow = 3, ncol = 3, scales = "free_x",
               labeller = labeller(measure = labeller_measure, part = labeller_part)) +
    scale_color_manual(name = "Population size", values = colors_pop) +
    scale_x_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 5), 
                       limits = function(x) c(min(x), max(x)), labels = function(x) round(x, 2)) +
    scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 5), 
                       limits = function(x) c(min(x), max(x)), labels = function(x) round(x, 2)) +
    guides(color = guide_legend(override.aes = list(shape = 19, alpha = 1.0), 
                                nrow = 2, byrow = TRUE)) +
    labs(x = ifelse(test = part_i == "ttl_production", yes = "log(Coeffiecent of variation)", no = ""),
         y = ifelse(test = part_i == "bg_production", yes = "log(Total primary production)", no = "")) +
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

suppoRt::save_ggplot(plot = gg_cv_prod_overall, filename = paste0("06-cv-prod-", amplitude, "-", 
                                                                  frequency, extension),
                     path = "04_Figures/", width = height, height = width,
                     units = units, dpi = dpi, overwrite = overwrite)

#### Mobility vs PE #### 

df_regression <- dplyr::filter(df_cv_prod, measure == "beta") %>% 
  dplyr::group_by(part, pop_n) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(i) {
    
    lm_temp <- lm(value.cv.log ~ value.move, data = i)
    
    pred <- coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * i$value.move
    
    r.squared <- cor(i$value.cv.log, pred)
    
    p.value <- summary(lm_temp)[["coefficients"]][2, 4]
    
    tibble::tibble(part = unique(i$part), pop_n = unique(i$pop_n), 
                   value.move = i$value.move, value.cv.log.pred = pred,
                   int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                   r.squared = r.squared, p.value = p.value)
    
  })

df_label <- dplyr::select(df_regression, -value.move, -value.cv.log.pred) %>% 
  dplyr::group_by(part, pop_n) %>% 
  dplyr::distinct() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(dir = dplyr::case_when(coef < 0.0 ~ "-", coef > 0.0 ~ "+"),
                int = round(int, digits = 2), coef = convert_label(coef, digits = 2),
                r.squared = round(x = r.squared, digits = 2), 
                p.value = dplyr::case_when(p.value <  0.001 ~  "'***'", p.value <  0.01 ~  "'**'",
                                           p.value <  0.05 ~  "'*'", p.value >  0.05 ~  "n.s."), 
                label = paste0("R^2==", dir, r.squared, "*';'~", p.value),
                label = dplyr::case_when(p.value == "n.s." ~ "'---'", TRUE ~ label)) %>% 
  dplyr::left_join(x = ., y = dplyr::filter(df_cv_prod, measure == "beta") %>% 
                     dplyr::group_by(part) %>% 
                     dplyr::summarise(pop_n = factor(c(8, 16, 32, 64), ordered = TRUE), 
                                      x = min(value.move) + 
                                        (max(value.move) - min(value.move)) * 0.95, 
                                      y = min(value.cv.log) + 
                                        (max(value.cv.log) - min(value.cv.log)) * c(1.0, 0.9, 0.8, 0.7), 
                                      .groups = "drop"), 
                   by = c("part", "pop_n"))

gg_move_pe_overall <- dplyr::filter(df_cv_prod, measure == "beta") %>% 
  ggplot(aes(x = value.move, y = value.cv.log, color = pop_n)) +
  geom_point(shape = 1, alpha = 0.1) +
  geom_line(data = df_regression, aes(x = value.move, y = value.cv.log.pred, color = pop_n)) +
  geom_text(data = df_label, aes(x = x, y = y, color = pop_n, label = label), parse = TRUE, size = 2.0) +
  facet_grid(rows = dplyr::vars(part), scales = "free_y",
             labeller = labeller(part = c("ag_production" = "Aboveground", "bg_production" = "Belowground",
                                          "ttl_production" = "Total"))) +
  scale_color_manual(name = "Population size", values = colors_pop) +
  scale_x_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 5), 
                     limits = function(x) c(min(x), max(x)), labels = function(x) round(x, 2)) +
  scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 5), 
                     limits = function(x) c(min(x), max(x)), labels = function(x) round(x, 2)) +
  guides(color = guide_legend(override.aes = list(shape = 19, alpha = 1.0), 
                              ncol = 1, byrow = TRUE)) +
  labs(x = "Mean cross-ecosystem movement", y = "Portfolio effect") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "right", legend.spacing.y = unit(-1, "mm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 8),
        axis.title = element_text(size = 6.5),
        strip.text = element_text(hjust = 0.0), strip.background = element_blank(), 
        axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA))

suppoRt::save_ggplot(plot = gg_move_pe_overall, filename = paste0("06-connect-pe-", amplitude, "-", 
                                                                  frequency, extension),
                     path = "04_Figures/", width = width, height = height * 0.40,
                     units = units, dpi = dpi, overwrite = overwrite)
