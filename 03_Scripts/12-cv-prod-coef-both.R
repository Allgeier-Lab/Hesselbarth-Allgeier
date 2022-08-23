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

#### Load/wrangle simulated data ####

n <- 5

df_phase <- import_data(path = paste0("02_Data/05-variability-phase-", n, ".rds"))

df_noise <- import_data(path = paste0("02_Data/05-variability-noise-", n, ".rds"))

df_total <- dplyr::bind_rows(phase = df_phase, noise = df_noise, .id = "scenario")

#### Fit regression model ####

df_regression <- dplyr::filter(df_total, part %in% c("ag_production", "bg_production", "ttl_production"), 
                               measure %in% c("alpha", "gamma")) %>%
  dplyr::group_by(scenario, part, measure, pop_n) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(df_temp) {
    
    df_temp$value.prod <- log(df_temp$value.prod) - mean(log(df_temp$value.prod))
    
    df_temp$value.cv <- log(df_temp$value.cv)
    df_temp$biotic <- log(df_temp$biotic)
    df_temp$abiotic <- log(df_temp$abiotic)
    
    lm_temp_cv <- lm(value.prod ~ value.cv, data = df_temp)
    lm_temp_sep <- lm(value.prod ~ biotic + abiotic, data = df_temp)
    
    purrr::map_dfr(list(cv = lm_temp_cv, separated = lm_temp_sep), 
                   broom::tidy, .id = "lm") %>%
      dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                    measure = unique(df_temp$measure), pop_n = unique(df_temp$pop_n), .before = term)
  }) %>%
  dplyr::mutate(p.value = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                                           p.value < 0.05 ~  "*", p.value >= 0.05 ~ ""),
                direction = dplyr::case_when(term != "(Intercept)" & estimate < 0.0 ~ "negative", 
                                             term != "(Intercept)" & estimate > 0.0 ~ "positive", 
                                             term == "(Intercept)" ~ as.character(NA)))

#### Setup ggplot ####

color_scenario <- c("phase" = "#1e466e", "noise" = "#e76254")

color_term <- c("biotic" = "#41b282", "abiotic" = "#007aa1")

size_point <- 5
size_line <- 0.75
size_text <- 2.5
size_base <- 10.0

alpha <- 0.5

width_pos <- 0.5

#### Create ggplot cv model ####

gg_coef_cv <- ggplot(data = dplyr::filter(df_regression, term != "(Intercept)", lm == "cv"), 
                     aes(group = scenario)) + 
  
  # adding geoms
  geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
  geom_linerange(aes(x = pop_n, ymin = 0.0, ymax = estimate, color = scenario),
                 position = position_dodge(width = width_pos), size = size_line) +
  geom_point(aes(x = pop_n, y = estimate, fill = scenario, color = scenario, shape = direction),
             position = position_dodge(width = width_pos), size = size_point) +
  geom_text(aes(x = pop_n, y = estimate, label = p.value), position = position_dodge(width = width_pos), 
            vjust = 0.75, size = size_text, color = "white") +
  
  # facet gridding
  facet_grid(rows = dplyr::vars(part), cols = dplyr::vars(measure), scales = "fixed", 
             labeller = labeller(measure = c("alpha" = "Local scale", "gamma" = "Meta-ecosystem scale"), 
                                 part = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                          "ttl_production" = "Total"))) + 
  
  # set scales and labs
  scale_color_manual(name = "", values = color_scenario,
                     labels = c("phase" = "Abiotic phase", "noise" = "Abiotic noise")) +
  scale_fill_manual(name = "", values = color_scenario,
                    labels = c("phase" = "Abiotic phase", "noise" = "Abiotic noise")) +
  scale_shape_manual(name = "", values = c("negative" = 25, "positive" = 24)) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(df_regression$pop_n))) +
  scale_y_continuous(limits = function(x) c(min(x), max(x)), 
                     breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4), 
                     labels = function(x) round(x, 2)) +

  # setup theme
  labs(x = "Population size", y = "Regression parameter estimate") +
  guides(color = guide_legend(order = 1, override.aes = list(shape = 17)), 
         fill = "none", shape = guide_legend(order = 2)) +
  theme_classic(base_size = size_base) + 
  theme(legend.position = "bottom", plot.title = element_text(size = 8.0), 
        panel.border = element_rect(size = 0.5, fill = NA),
        strip.background = element_blank(), strip.text = element_text(hjust = 0))

#### Create ggplot separated model ####

gg_list <- purrr::map(c("alpha", "gamma"), function(measure_i) {
  
  purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
    
    df_temp <- dplyr::filter(df_regression, term != "(Intercept)", lm == "separated", 
                             part == part_i, measure == measure_i)
    
    ggplot(data = df_temp, aes(group = term)) + 
      
      # adding geoms
      geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
      geom_linerange(aes(x = pop_n, ymin = 0.0, ymax = estimate, color = term),
                     position = position_dodge(width = width_pos), size = size_line) +
      geom_point(aes(x = pop_n, y = estimate, fill = term, color = term, shape = direction),
                 position = position_dodge(width = width_pos), size = size_point * 0.75) +
      geom_text(aes(x = pop_n, y = estimate, label = p.value), position = position_dodge(width = width_pos), 
                vjust = 0.75, size = size_text, color = "white") +
      
      # facet gridding
      facet_wrap(. ~ scenario, ncol = 2, labeller = labeller(term = c("biotic" = "Biotic", "abiotic" = "Abiotic"), 
                                                             scenario = c("noise" = "Noise scenario", "phase" = "Phase scenario"))) +
      
      # set scales and labs
      # scale_color_manual(name = "", values = color_scenario,
      #                    labels = c("phase" = "Abiotic phase", "noise" = "Abiotic noise")) +
      # scale_fill_manual(name = "", values = color_scenario,
      #                   labels = c("phase" = "Abiotic phase", "noise" = "Abiotic noise")) +
      scale_color_manual(name = "", values = color_term) +
      scale_fill_manual(name = "", values = color_term) +
      scale_shape_manual(name = "", values = c("negative" = 25, "positive" = 24)) +
      coord_flip() +
      scale_x_discrete(limits = rev(levels(df_regression$pop_n))) +
      scale_y_continuous(limits = function(x) c(min(x), max(x)), 
                         breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4), 
                         labels = function(x) round(x, 2)) +
      
      # setup theme
      labs(x = "", y = "") +
      guides(color = guide_legend(order = 1, override.aes = list(shape = 17)), 
             fill = "none", shape = guide_legend(order = 2)) +
      theme_classic(base_size = size_base) + 
      theme(legend.position = "none", plot.title = element_text(size = 8.0), 
            panel.border = element_rect(size = 0.5, fill = NA),
            strip.background = element_blank(), strip.text = element_text(hjust = 0))
  })
}) %>% purrr::flatten()

gg_coef_sep <- cowplot::plot_grid(plotlist = gg_list, nrow = 3, ncol = 2, byrow = FALSE,
                                  labels = c("a)", "b)", "c)", "d)", "e)", "f)"), 
                                  label_fontface = "plain", label_size = 12.0)

gg_coef_sep <- cowplot::ggdraw(gg_coef_sep, xlim = c(-0.015, 1.0)) + 
  cowplot::draw_label("Regression parameter estimate", x = 0.5, y = 0, vjust = -0.5, angle = 0, size = size_base) + 
  cowplot::draw_label("Population size", x = 0.0, y = 0.5, vjust = 0.0, 
                      angle = 90, size = size_base)

gg_dummy <- ggplot(data = data.frame(pop_n = 1, value = 1, term = c("biotic", "abiotic"), 
                                     direction = c("positive", "negative")), 
                   aes(x = pop_n, y = value, color = term, shape = direction)) +
  geom_point(size = size_point) + 
  scale_shape_manual(name = "", values = c("negative" = 25, "positive" = 24)) +
  scale_color_manual(name = "", values = color_term, 
                     labels = c("biotic" = "Biotic", "abiotic" = "Abiotic")) + 
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme(legend.position = "bottom")

gg_coef_sep <- cowplot::plot_grid(gg_coef_sep, cowplot::get_legend(gg_dummy),
                                  nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))

#### Save plot ####

suppoRt::save_ggplot(plot = gg_coef_cv, filename = paste0("12-coef-cv-", n, extension),
                     path = "04_Figures/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_coef_sep, filename = paste0("12-coef-sep-", n, extension),
                     path = "04_Figures/", width = width, height = height,
                     units = units, dpi = dpi, overwrite = FALSE)
