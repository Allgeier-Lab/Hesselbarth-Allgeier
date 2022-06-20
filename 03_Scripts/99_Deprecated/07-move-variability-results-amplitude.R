##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Results of simulation experiment for low/high amplitude treatment

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/convert_label.R")

extension <- ".pdf"

colors_pop <- c("8" = "#ebcc60", "16" = "#459b75", "32" = "#374a98", "64" = "#913f98")

#### Load/wrangle simulated data ####

file_names <- list.files(path = "02_Data/", pattern = "05-move-variability-*", 
                         full.names = TRUE)

names(file_names) <- stringr::str_sub(file_names, start = -10, end = -5)

df_cv_prod <- purrr::map_dfr(file_names, function(i) {
  readr::read_rds(i) %>% 
    purrr::map_dfr(function(j) {
      dplyr::left_join(x = j$cv, y = j$production, by = c("part", "pop_n", "move_meta_mean", 
                                                          "move_meta_sd"), 
                       suffix = c(".cv", ".prod")) %>% 
        dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))
    })
}, .id = "experiment") %>%
  tidyr::separate(col = experiment, sep = "-", into = c("amplitude", "cycles")) %>% 
  dplyr::mutate(amplitude = factor(amplitude, levels = c("005", "095"), labels = c("5%", "95%")), 
                cycles = factor(cycles, levels = c("50", "12"), labels = c("50", "12")),
                part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass",
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                pop_n = factor(as.numeric(pop_n), ordered = TRUE),
                value.cv.log = log(value.cv), value.prod.log = log(value.prod), 
                value.move.log = log(value.move)) %>% 
  dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), 
                measure %in% c("alpha", "gamma", "beta")) %>% 
  tibble::tibble()

df_regression <- dplyr::group_by(df_cv_prod, amplitude, cycles, part, measure, pop_n) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(i) {
    
    lm_temp <- lm(value.prod.log ~ value.cv.log, data = i)
    
    pred <- coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * i$value.cv.log
    
    r.squared <- cor(i$value.prod.log, pred)
    
    p.value <- summary(lm_temp)[["coefficients"]][2, 4]
    
    tibble::tibble(cycles = unique(i$cycles), amplitude = unique(i$amplitude),
                   part = unique(i$part), measure = unique(i$measure), pop_n = unique(i$pop_n), 
                   value.cv.log = i$value.cv.log, value.prod.log.pred = pred,
                   int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                   r.squared = r.squared, p.value = p.value)
  }) %>% 
  dplyr::mutate(dir = dplyr::case_when(coef < 0.0 ~ "-", coef > 0.0 ~ ""), 
                r.squared = round(x = r.squared, digits = 2),
                p.value = dplyr::case_when(p.value <  0.001 ~  "'***'", p.value <  0.01 ~  "'**'",
                                           p.value <  0.05 ~  "'*'", p.value >  0.05 ~  "n.s."),
                label = paste0("R^2==", dir, r.squared, "*';'~", p.value), 
                label = dplyr::case_when(p.value == "n.s." ~ "", TRUE ~ label))

gg_coef <- dplyr::select(df_regression, -c(value.cv.log, value.prod.log.pred)) %>%
  dplyr::distinct() %>%
  # dplyr::filter(cycles == 50) %>%
  dplyr::mutate(measure = factor(measure, levels = c("beta", "gamma", "alpha"))) %>% 
  ggplot(aes(x = pop_n, y = measure, fill = coef)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = label), parse = TRUE, size = 1.5) +
  scale_fill_gradient2(name = "Regression\ncoeffiecent", low = "#e96052", mid = "#ffe6b8", high = "#19436d",
                       midpoint = 0.0, limits = function(x) c(-max(abs(x)), max(abs(x)))) +
  scale_y_discrete(labels = c(alpha = "Local\nscale", gamma = "Meta-ecosystem\nscale", beta = "Portfolio\neffect")) +
  coord_equal() +
  labs(x = "Population size", y = "Scale") +
  facet_wrap(. ~ part + cycles + amplitude, ncol = 4, nrow = 3, 
             labeller = labeller(part = c(ag_production = "Aboveground", bg_production = "Belowground", 
                                          ttl_production = "Total"), 
                                 cycles = c("50" = "Yearly cycles", "12" = "4 year cycles"),
                                 amplitude = c("5%" = "Low amplitude", "95%" = "High amplitude"))) +  
  theme_classic(base_size = 8.0) +
  theme(legend.position = "right", legend.text = element_text(size = 6), legend.title = element_text(size = 8),
        axis.title = element_text(size = 6.5), axis.line = element_blank(), 
        panel.border = element_rect(size = 0.5, fill = NA), strip.text = element_text(hjust = 0.0), 
        strip.background = element_blank())

suppoRt::save_ggplot(plot = gg_coef, filename = paste0("07-pop-scale-coeff", extension),
                     path = "04_Figures/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = overwrite)  
