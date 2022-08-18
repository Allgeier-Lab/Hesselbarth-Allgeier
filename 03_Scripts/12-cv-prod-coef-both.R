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

amplitude <- "095"

df_phase <- import_data(path = paste0("02_Data/05-variability-phase-", amplitude, ".rds"))

df_noise <- import_data(path = paste0("02_Data/05-variability-noise-", amplitude, ".rds"))

df_total <- dplyr::bind_rows(phase = df_phase, noise = df_noise, .id = "scenario")

#### Fit regression model ####

df_regression <- dplyr::filter(df_total, part %in% c("ag_production", "bg_production", "ttl_production"), 
                               measure %in% c("alpha", "gamma")) %>%
  dplyr::group_by(scenario, part, measure, pop_n) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(df_temp) {
    
    df_temp$value.prod <- log(df_temp$value.prod) - mean(log(df_temp$value.prod))
    
    lm_temp <- lm(value.prod ~ log(value.cv), data = df_temp)
    
    broom::tidy(lm_temp) %>%
      dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                    measure = unique(df_temp$measure), pop_n = unique(df_temp$pop_n), .before = term) %>% 
      dplyr::mutate(term = dplyr::case_when(term == "(Intercept)" ~ "Intercept", TRUE ~ term), 
                    p.value = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                                               p.value < 0.05 ~  "*", p.value >= 0.05 ~ ""),
                    direction = dplyr::case_when(term != "Intercept" & estimate < 0.0 ~ "negative", 
                                                 term != "Intercept" & estimate > 0.0 ~ "positive", 
                                                 term == "Intercept" ~ as.character(NA)))
    
  })

#### Setup ggplot ####

color_scenario <- c("phase" = "#1e466e", "noise" = "#e76254")

size_point <- 5
size_line <- 0.75
size_text <- 2.5
size_base <- 10.0

alpha <- 0.5

width_pos <- 0.5

#### Create ggplot ####

# ypos_text <- dplyr::filter(df_regression, term != "Intercept") %>% 
#   dplyr::pull(estimate) %>% 
#   max() %>% 
#   magrittr::multiply_by(1.15)

gg_coef <- ggplot(data = dplyr::filter(df_regression, term != "Intercept"), 
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

#### Save plot ####

suppoRt::save_ggplot(plot = gg_coef, filename = paste0("12-both-", amplitude, extension),
                     path = "04_Figures/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = FALSE)
