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

extension <- ".pdf"

amplitude <- "095"

#### Load/wrangle simulated data ####

file_path_phase <- paste0("02_Data/05-variability-phase-", amplitude, ".rds")

file_path_noise <- paste0("02_Data/05-variability-noise-", amplitude, ".rds")

df_phase <- import_data(path = file_path_phase)

df_noise <- import_data(path = file_path_noise)

df_total <- dplyr::bind_rows(phase = df_phase, noise = df_noise, .id = "scenario") %>% 
  dplyr::filter(measure %in% c("alpha", "gamma")) %>% 
  dplyr::group_by(scenario, part, measure, pop_n) %>%
  dplyr::group_split() %>% 
  purrr::map_dfr(function(i) {
    
    density_temp <- dplyr::pull(i, value.cv) %>% 
      density()
    
    data.frame(scenario = unique(i$scenario), part = unique(i$part), 
               measure = unique(i$measure), pop_n = unique(i$pop_n), 
               cv = density_temp$x, density = density_temp$y)
    }) %>% 
  tidyr::unite("color_id", pop_n, scenario, remove = FALSE) %>% 
  dplyr::mutate(color_id = factor(color_id, 
                                  levels = c("8_phase", "8_noise", "16_phase", "16_noise",
                                             "32_phase", "32_noise", "64_phase", "64_noise"), 
                                  labels = c("8 Indiv.; Phase", "8 Indiv.; Noise",
                                             "16 Indiv.; Phase", "16 Indiv.; Noise",
                                             "32 Indiv.; Phase", "32 Indiv.; Noise",
                                             "64 Indiv.; Phase", "64 Indiv.; Noise")), 
                scenario = factor(scenario, levels = c("phase", "noise")))

#### Calculate means #### 

df_total_sum <- dplyr::group_by(df_total, scenario, part, measure, pop_n) %>% 
  dplyr::summarise(cv_mn = mean(cv), cv_sd = sd(cv), .groups = "drop")

dplyr::filter(df_total) %>% 
  dplyr::group_by(scenario, part, measure, pop_n) %>% 
  dplyr::summarise(mean = mean(cv), min = min(cv), max = max(cv), groups = "drop") %>% 
  dplyr::mutate(range = max - min) %>% 
  dplyr::group_by(scenario, part, measure) %>% 
  dplyr::summarise(range = mean(range)) %>% 
  dplyr::mutate(range = round(range, 3)) %>% 
  dplyr::filter(scenario == "phase", measure == "gamma")

#### Setup ggplot ####

colors_pop <- c("8 Indiv.; Phase" = "#447861", "8 Indiv.; Noise" = "#7caf5c",
                "16 Indiv.; Phase" = "#13315f", "16 Indiv.; Noise" = "#4457a5",
                "32 Indiv.; Phase" = "#59386c", "32 Indiv.; Noise" = "#b1a1cc",
                "64 Indiv.; Phase" = "#b24422", "64 Indiv.; Noise" = "#c44d76")

size_base <- 10

#### Create ggplot ####

gg_cv_densities <- ggplot(df_total, aes(x = cv, y = density)) +
  geom_line(aes(color = color_id, linetype = scenario)) + 
  geom_area(aes(fill = color_id), alpha = 0.25) +
  facet_wrap(. ~ part + measure, nrow = 3, ncol = 2,  scales = "free",
             labeller = labeller(part = c(ag_production = "Aboveground", 
                                          bg_production = "Belowground", 
                                          ttl_production = "Total production"), 
                                 measure = c(alpha = "Local scale", 
                                             gamma = "Meta-ecosystem scale"))) +
  scale_x_continuous(limits = function(x) c(min(x), max(x)), 
                     breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 5), 
                     labels = function(x) round(x, 2)) +
  scale_fill_manual(name = "", values = colors_pop) +
  scale_color_manual(name = "", values = colors_pop) + 
  scale_linetype_manual(name = "", values = c(phase = 1, noise = 2), 
                        labels = c("phase" = "Abiotic phase", noise = "Abiotic noise")) +
  guides(linetype = guide_legend(nrow = 2, order = 1)) +
  labs(x = "Coefficient of variation", y = "Density") +
  theme_classic(base_size = size_base) +
  theme(legend.position = "bottom",
        axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0.0))

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_cv_densities, filename = paste0("06-cv-densities-", amplitude, extension),
                     path = "04_Figures/", width = width, height = height * 0.85,
                     units = units, dpi = dpi, overwrite = FALSE)
