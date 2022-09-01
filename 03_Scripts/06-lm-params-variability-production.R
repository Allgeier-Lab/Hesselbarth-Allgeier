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
source("01_Functions/import_data.R")

#### Load/wrangle simulated data ####

n <- 5

results_phase_df <- import_data(path = paste0("02_Data/result-phase-", n, ".rds"))

results_noise_df <- import_data(path = paste0("02_Data/result-noise-", n, ".rds"))

#### Fit regression model ####

regression_df <- purrr::map_dfr(list(phase = results_phase_df, noise = results_noise_df), function(results_temp_df) {
  
  dplyr::filter(results_temp_df, part %in% c("ag_production", "bg_production", "ttl_production"), 
                measure %in% c("alpha", "gamma")) %>%
    dplyr::group_by(part, measure, pop_n) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(temp_df) {
      
      temp_df_stand <- dplyr::select(temp_df, value.prod, value.cv, biotic, abiotic) %>% 
        dplyr::mutate(across(.fns = function(x) log(x))) %>%
        dplyr::mutate(across(.fns = function(x) (x - mean(x)) / sd(x)))
      
      purrr::map_dfr(c(cv = "value.cv", separated = "biotic * abiotic"), y = temp_df_stand, function(x, y) {
        
        lm_temp <- lm(as.formula(paste0("value.prod ~ ", paste0(x))), data = y) 
        
        broom::tidy(lm_temp) %>% 
          dplyr::mutate(part = unique(temp_df$part), measure = unique(temp_df$measure), pop_n = unique(temp_df$pop_n), 
                        .before = term) %>% 
          dplyr::mutate(p.value = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", 
                                                   p.value >= 0.05 ~ ""), 
                        r2 = summary(lm_temp)$adj.r.squared) %>% 
          dplyr::filter(term != "(Intercept)")
      }, .id = "explanatory")
    }) %>% 
    dplyr::mutate(p.value = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                                             p.value < 0.05 ~  "*", p.value >= 0.05 ~ ""),
                  direction = dplyr::case_when(term != "(Intercept)" & estimate < 0.0 ~ "negative", 
                                               term != "(Intercept)" & estimate > 0.0 ~ "positive", 
                                               term == "(Intercept)" ~ as.character(NA)), 
                  term = dplyr::case_when(term == "value.cv" ~ "cv", term != "value.cv" ~ term))
}, .id = "scenario")

#### Setup ggplot ####

color_term <- c("biotic" = "#00af73", "abiotic" = "#006d9a", "biotic:abiotic" = "#ffaa3a", 
                "cv" = "#ff3b18")

size_point <- 5
size_line <- 0.75
size_text <- 2.5
size_base <- 10.0

alpha <- 0.5

width_pos <- 0.65

#### Create ggplot model parameters ####

gg_coef_scenario <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  ggplot(data = dplyr::filter(regression_df, scenario == scenario_i, term != "(Intercept)")) + 
    
    # adding geoms
    geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
    geom_linerange(aes(x = pop_n, ymin = 0.0, ymax = estimate, color = term),
                   position = position_dodge(width = width_pos), size = size_line) +
    geom_point(aes(x = pop_n, y = estimate, fill = term, color = term, shape = direction),
               position = position_dodge(width = width_pos), size = size_point * 0.75) +
    geom_text(aes(x = pop_n, y = estimate, label = p.value, group = term), position = position_dodge(width = width_pos), 
              vjust = 0.75, size = size_text, color = "white") +
    
    # facet gridding
    facet_grid(rows = dplyr::vars(part), cols = dplyr::vars(measure), 
               labeller = labeller(part = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                            "ttl_production" = "Total"), 
                                   measure = c("alpha" = "Local", "gamma" = "Meta-ecosystem"))) +
    
    # set scales and labs
    scale_color_manual(name = "", values = color_term) +
    scale_fill_manual(name = "", values = color_term) +
    scale_shape_manual(name = "", values = c("negative" = 25, "positive" = 24)) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(regression_df$pop_n))) +
    scale_y_continuous(limits = function(x) range(x), labels = function(x) round(x, 2),
                       breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4)) +
    
    # setup theme
    labs(x = "", y = "") +
    guides(color = guide_legend(order = 1, override.aes = list(shape = 17)), 
           fill = "none", shape = guide_legend(order = 2)) +
    theme_classic(base_size = size_base) + 
    theme(legend.position = "bottom", plot.title = element_text(size = 8.0), 
          panel.border = element_rect(size = 0.5, fill = NA),
          strip.background = element_blank(), strip.text = element_text(hjust = 0))
})

#### Save plot ####

suppoRt::save_ggplot(plot = gg_coef_scenario$phase, filename = paste0("Figure-4", extension),
                     path = "04_Figures/", width = width, height = height * 0.85,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_coef_scenario$noise, filename = paste0("Figure-5", extension),
                     path = "04_Figures/", width = width, height = height * 0.85,
                     units = units, dpi = dpi, overwrite = FALSE)
