##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Variation partitioning noise variability

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/import_data.R")

#### Load/wrangle simulated data ####

n <- 5

results_phase_df <- import_data(path = paste0("02_Data/result-phase-", n, ".rds"))

results_noise_df <- import_data(path = paste0("02_Data/result-noise-", n, ".rds"))

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario")

#### Fit regression model ####

regression_df <- dplyr::filter(results_combined_df, part %in% c("ag_production", "bg_production", 
                                                                "ttl_production"), 
                               measure %in% c("alpha", "gamma")) %>%
  dplyr::group_by(scenario, part, measure, pop_n) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(temp_df) {
  
    temp_df_stand <- dplyr::select(temp_df, value.prod, value.cv) %>% 
      dplyr::mutate(across(.fns = function(x) log(x))) %>%
      dplyr::mutate(across(.fns = function(x) (x - mean(x)) / sd(x)))
    
    lm_temp <- lm(value.prod ~ value.cv, data = temp_df_stand)
    
    broom::tidy(lm_temp) %>%
      dplyr::mutate(scenario = unique(temp_df$scenario), part = unique(temp_df$part), 
                    measure = unique(temp_df$measure), pop_n = unique(temp_df$pop_n), .before = term) %>% 
      dplyr::mutate(p.value = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                                               p.value < 0.05 ~  "*", p.value >= 0.05 ~ "n.s."),
                    r2 = summary(lm_temp)$adj.r.squared,
                    direction = dplyr::case_when(term != "(Intercept)" & estimate < 0.0 ~ "increase", 
                                                 term != "(Intercept)" & estimate > 0.0 ~ "decrease", 
                                                 term == "(Intercept)" ~ as.character(NA)))
})

#### Setup ggplot ####

colors_pop <- c("8" = "#64be35", "16" = "#ac4288", "32" = "#ffe16a", "64" = "#ff939b", "128" = "#fe2a2d")

size_point <- 1.0
size_line <- 0.75
size_base <- 10.0

alpha <- 0.25

gg_dummy <- ggplot(data = dplyr::filter(results_combined_df, measure != "beta"), 
                   aes(x = value.cv, y = value.prod, color = pop_n)) + 
  geom_point(shape = 1, alpha = alpha) + 
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = size_line) +
  scale_color_manual(name = "Population size", values = colors_pop) +
  theme_classic(base_size = size_base) + 
  theme(legend.position = "bottom")

#### Create ggplot ####

gg_abiotic_list <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  gg_scatter_list <- dplyr::filter(results_combined_df, scenario == scenario_i, 
                                   measure %in% c("alpha", "gamma")) %>%
    dplyr::group_by(part, measure) %>%
    dplyr::group_split() %>%
    purrr::map(function(temp_df) {
      
      temp_df_stand <- dplyr::select(temp_df, value.prod, value.cv, pop_n) %>% 
        dplyr::mutate(across(-pop_n, function(x) log(x))) %>%
        dplyr::mutate(across(-pop_n, function(x) (x - mean(x)) / sd(x)))
      
      # init ggplot
      ggplot(data = temp_df_stand, aes(x = value.cv, y = value.prod, color = pop_n)) + 
        
        # adding geoms
        geom_point(shape = 1, alpha = alpha, size = size_point) + 
        geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = size_line) +
        
        # facet wrap 
        facet_wrap(. ~ pop_n, scales = "free", nrow = 2, ncol = 3) + 
        
        # set scales
        scale_color_manual(name = "Population size", values = colors_pop) +
        scale_x_continuous(limits = function(x) range(x), labels = function(x) round(x, 2),
                           breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4)) +
        scale_y_continuous(limits = function(x) range(x), labels = function(x) round(x, 2),
                           breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4)) +
        
        # labels and themes
        labs(x = "", y = "") +
        theme_classic(base_size = 10.0) + 
        theme(legend.position = "none", plot.title = element_text(size = 8.0), 
              strip.background = element_blank(), strip.text = element_blank(),
              plot.margin = margin(t = 5.5, r = 5.5, b = 0.0, l = 5.5, unit = "pt"))
      
    })
  
  gg_combined <- cowplot::plot_grid(plotlist = gg_scatter_list, nrow = 3, ncol = 2, byrow = TRUE, 
                                    labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
                                    label_fontface = "plain", label_size = 12.0)
  
  gg_combined <- cowplot::ggdraw(gg_combined, xlim = c(-0.015, 1.0)) +
    cowplot::draw_label("Coeffiecent of variation", x = 0.5, y = 0, vjust = -0.5, angle = 0, size = size_base) + 
    cowplot::draw_label("Primary production", x = 0.0, y = 0.5, vjust = 0.0, angle = 90, size = size_base)
  
  cowplot::plot_grid(gg_combined, cowplot::get_legend(gg_dummy),
                     nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))
  
})

#### Save plot ####

suppoRt::save_ggplot(plot = gg_abiotic_list$phase, filename = paste0("cv-vs-production-phase", extension),
                     path = "04_Figures/Appendix", width = height, height = width * 0.75,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_abiotic_list$noise, filename = paste0("cv-vs-production-noise", extension),
                     path = "04_Figures/Appendix", width = height, height = width * 0.75,
                     units = units, dpi = dpi, overwrite = FALSE)
