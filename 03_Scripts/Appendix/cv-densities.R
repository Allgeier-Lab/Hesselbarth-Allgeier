##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Densities of CV values

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/import_data.R")

#### Load data####

n <- 5

phase_df <- import_data(path = paste0("02_Data/result-phase-", n, ".rds"))

noise_df <- import_data(path =  paste0("02_Data/result-noise-", n, ".rds"))

#### Wrangle data ####

combined_df <- dplyr::bind_rows(phase = phase_df, noise = noise_df, .id = "scenario") |> 
  dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"),
                measure %in% c("alpha", "gamma")) |> 
  dplyr::group_by(scenario, part, measure, pop_n) |> 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(i) {
    
    density_temp <- dplyr::pull(i, value.cv) |> 
      density()
    
    data.frame(scenario = unique(i$scenario), part = unique(i$part), 
               measure = unique(i$measure), pop_n = unique(i$pop_n), 
               cv = density_temp$x, density = density_temp$y)
    
  }) |> 
  tidyr::unite("color_id", pop_n, scenario, remove = FALSE) |> 
  dplyr::mutate(color_id = factor(color_id, 
                                  levels = c("8_phase", "8_noise", "16_phase", "16_noise",
                                             "32_phase", "32_noise", "64_phase", "64_noise"), 
                                  labels = c("8 Indiv.; Phase", "8 Indiv.; Noise",
                                             "16 Indiv.; Phase", "16 Indiv.; Noise",
                                             "32 Indiv.; Phase", "32 Indiv.; Noise",
                                             "64 Indiv.; Phase", "64 Indiv.; Noise")), 
                scenario = factor(scenario, levels = c("phase", "noise")))

#### Calculate means #### 

combined_df_sum <- dplyr::group_by(combined_df, scenario, part, measure, pop_n) |> 
  dplyr::summarise(cv_mn = mean(cv), cv_sd = sd(cv), .groups = "drop")

#### Setup ggplot ####

colors_pop <- c("8 Indiv.; Phase" = "#92c6de", "8 Indiv.; Noise" = "#0069aa",
                "16 Indiv.; Phase" = "#9cde81", "16 Indiv.; Noise" = "#00992a",
                "32 Indiv.; Phase" = "#c9a6cf", "32 Indiv.; Noise" = "#662e8e",
                "64 Indiv.; Phase" = "#ffba68", "64 Indiv.; Noise" = "#ff771c",
                "128 Indiv.; Phase" = "#ff908f", "128 Indiv.; Noise" = "#f32222")

size_base <- 10

#### Create ggplot ####

gg_cv_densities <- ggplot(combined_df, aes(x = cv, y = density)) +
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
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.0), 
        axis.line = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
        legend.position = "bottom",)

#### Save ggplot ####

# suppoRt::save_ggplot(plot = gg_cv_densities, filename = paste0("Figure-A2", extension),
#                      path = "04_Figures/Appendix/", width = width, height = height * 0.85,
#                      units = units, dpi = dpi, overwrite = FALSE)
