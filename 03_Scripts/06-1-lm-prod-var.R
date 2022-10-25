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

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") %>% 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))

#### Split data #### 

results_combined_list <- dplyr::filter(results_combined_df, 
                                       part %in% c("ag_production", "bg_production", "ttl_production"), 
                                       measure == "gamma") %>%
  dplyr::group_by(scenario, part, pop_n) %>%
  dplyr::group_split()

#### Fit regression model ####

regression_df <- purrr::map_dfr(results_combined_list, function(temp_df) {
  
  temp_df_stand <- dplyr::select(temp_df, value.prod, value.cv, biotic, abiotic) %>% 
    dplyr::mutate(across(.fns = function(x) log(x))) %>%
    dplyr::mutate(across(.fns = function(x) (x - mean(x)) / sd(x)))
  
  purrr::map_dfr(c(cv = "value.cv", separated = "biotic * abiotic"), y = temp_df_stand, function(x, y) {
    
    lm_temp <- lm(as.formula(paste0("value.prod ~ ", paste0(x))), data = y) 
    
    broom::tidy(lm_temp) %>% 
      dplyr::mutate(scenario = unique(temp_df$scenario), part = unique(temp_df$part), 
                    pop_n = unique(temp_df$pop_n),
                    r2 = summary(lm_temp)$adj.r.squared, .before = term) %>% 
      dplyr::filter(term != "(Intercept)")}, .id = "explanatory")}) %>% 
  dplyr::mutate(explanatory = factor(explanatory, levels = c("cv", "separated")), 
                term = dplyr::case_when(term == "value.cv" ~ "cv", term != "value.cv" ~ term),
                term = factor(term, levels = c("cv", "biotic", "abiotic", "biotic:abiotic")),
                p.value = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                                           p.value < 0.05 ~  "*", p.value >= 0.05 ~ ""),
                direction = dplyr::case_when(estimate < 0.0 ~ "decrease", 
                                             estimate > 0.0 ~ "increase"))
#### Relative importance R2 ####

importance_df <- purrr::map_dfr(results_combined_list, function(temp_df) {
  
  message("\r> ", unique(temp_df$scenario), "; ", unique(temp_df$part), "; ",
          unique(temp_df$measure), "; ", unique(temp_df$pop_n), "\t\t", appendLF = FALSE)
  
  temp_df_stand <- dplyr::select(temp_df, value.prod, biotic, abiotic) %>% 
    dplyr::mutate(across(.fns = function(x) log(x))) %>%
    dplyr::mutate(across(.fns = function(x) (x - mean(x)) / sd(x)))
  
  lm_temp <- lm(value.prod ~ biotic * abiotic, data = temp_df_stand) 
  
  rel_r2 <- relaimpo::boot.relimp(lm_temp, type = "lmg", b = 100, level = 0.95, fixed = FALSE) %>%
    relaimpo::booteval.relimp(bty = "basic")
 
  tibble::tibble(
    scenario = unique(temp_df$scenario), part = unique(temp_df$part), 
    measure = unique(temp_df$measure), pop_n = unique(temp_df$pop_n),
    beta = c(rel_r2@namen[-1], "residual"), mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)), 
    lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA)
  )}) %>% 
  dplyr::mutate(beta = factor(beta, levels = c("residual", "biotic:abiotic", "abiotic", "biotic")))

#### Setup ggplot ####

# color_term <- c("biotic" = "#00af73", "abiotic" = "#006d9a", "biotic:abiotic" = "#ffaa3a", 
#                 "cv" = "#ff3b18")

color_term <- c("biotic" = "#00af73", "abiotic" = "#006d9a", "biotic:abiotic" = "#ffaa3a", 
                "residuals" = "grey")

size_point <- 5
size_line <- 0.75
size_text <- 2.5
size_base <- 10.0

alpha <- 0.5

width_pos <- 0.65

#### Create ggplot model parameters ####

gg_coef_scenario <- ggplot(data = dplyr::filter(regression_df, term %in% c("biotic", "abiotic",
                                                                           "biotic:abiotic"))) + 
  
  # adding geoms
  geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
  geom_linerange(aes(x = pop_n, ymin = 0.0, ymax = estimate, color = term),
                 position = position_dodge(width = width_pos), size = size_line) +
  geom_point(aes(x = pop_n, y = estimate, fill = term, color = term, shape = direction),
             position = position_dodge(width = width_pos), size = size_point * 0.75) +
  geom_text(aes(x = pop_n, y = estimate, label = p.value, group = term), position = position_dodge(width = width_pos), 
            vjust = 0.75, size = size_text, color = "white") +
  
  # facet gridding
  facet_grid(rows = dplyr::vars(part), cols = dplyr::vars(scenario),
             labeller = labeller(part = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                          "ttl_production" = "Total"), 
                                 scenario = c("phase" = "Phase scenario","noise" = "Noise scenario"))) +
  
  # set scales and labs
  scale_color_manual(name = "", values = color_term[-4], 
                     labels = c("Consumer behavior", "Abiotic subsidies", "Interaction Behavior:Subsidies", "Residuals")) +
  scale_fill_manual(name = "", values = color_term[-4], 
                    labels = c("Consumer behavior", "Abiotic subsidies", "Interaction Behavior:Subsidies", "Residuals")) +
  scale_shape_manual(name = "", values = c("decrease" = 25, "increase" = 24)) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(regression_df$pop_n))) +
  scale_y_continuous(limits = function(x) range(x), labels = function(x) round(x, 2),
                     breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4)) +
  
  # setup theme
  labs(x = "Population size", y = "Parameter estimate") +
  guides(color = guide_legend(order = 1, override.aes = list(shape = 17)), 
         fill = "none", shape = guide_legend(order = 2)) +
  theme_classic(base_size = size_base) + 
  theme(legend.position = "bottom", plot.title = element_text(size = 8.0), 
        axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA),
        strip.background = element_blank(), strip.text = element_text(hjust = 0))

#### Create ggplot relative importance ####

gg_relimp_scenario <- ggplot(data = importance_df) + 
    
  # relative importance bars
  geom_col(aes(x = pop_n, y = mean * 100, fill = beta)) + 
  
  # set scales
  scale_fill_manual(name = "", values = color_term) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) + 
  
  # facet gridding
  facet_grid(rows = dplyr::vars(part), cols = dplyr::vars(scenario),
             labeller = labeller(part = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                          "ttl_production" = "Total"), 
                                 scenario = c("phase" = "Phase scenario","noise" = "Noise scenario"))) +
    
  # setup theme
  labs(x = "Population size", y = "Parameter estimate") +
  theme_classic(base_size = size_base) + 
  theme(legend.position = "bottom",
        axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA),
        strip.background = element_blank(), strip.text = element_text(hjust = 0), 
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt"))

#### Save plot ####

overwrite <- FALSE

# Coef

suppoRt::save_ggplot(plot = gg_coef_scenario, filename = paste0("Figure-4", extension),
                     path = "04_Figures/", width = width, height = height * 0.75,
                     units = units, dpi = dpi, overwrite = overwrite)

# Rel importance

# suppoRt::save_ggplot(plot = gg_relimp_scenario, filename = paste0("Figure-A5", extension),
#                      path = "04_Figures/Appendix/", width = width, height = height * 0.85,
#                      units = units, dpi = dpi, overwrite = overwrite)
