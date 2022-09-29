##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Figure of nutrient input curves

#### Load setup ####

source("05_Various/setup.R")

#### Adapt parameters ####

amplitude_mn <- 0.95

frequency <- years

#### Stable values #### 

stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
                                              ag_biomass = starting_values_list$ag_biomass,
                                              parameters = parameters_list)

starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool

starting_values_list$detritus_pool <- stable_values_list$detritus_pool

#### Simulate nutrient inputs variability ####

variability <- 0.65

n <- 5

input_values_phase <- meta.arrR::simulate_nutrient_sine(n = n, max_i = max_i, frequency = frequency, 
                                                        input_mn = nutrient_input, amplitude_mn = amplitude_mn, 
                                                        phase_sd = variability) %>% 
  meta.arrR::get_input_df(gamma = FALSE, long = TRUE)

input_values_noise <- meta.arrR::simulate_nutrient_noise(n = n, max_i = max_i, frequency = frequency, 
                                                         input_mn = nutrient_input, amplitude_mn = amplitude_mn, 
                                                         noise_sd = variability) %>% 
  
  meta.arrR::get_input_df(gamma = FALSE, long = TRUE)

input_values_combined <- dplyr::bind_rows(phase = input_values_phase, noise = input_values_noise, 
                                          .id = "scenario") %>% 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise"), 
                                  labels = c("Phase scenario", "Noise scenario")))

#### Setup plots ####

filter_factor <- 10

breaks_x <- seq(from = 0, to = max_i / filter_factor, by = max_i / years)
# labels_x <- paste(seq(from = 0, to = years / filter_factor, by = 1), "years")

line_size <- 0.25

line_color <- MetBrewer::met.brewer(name = "Java", n = n, type = "discrete")

base_size <- 7.5

extension <- ".png"

#### Create ggplot ####

gg_input_curves <- dplyr::filter(input_values_combined, timestep < max_i / filter_factor) %>% 
  
  ggplot(aes(x = timestep, y = value, color = meta)) +
  
  # adding geoms
  geom_hline(yintercept = nutrient_input, linetype = 2, color = "grey") +
  geom_line(size = line_size) + 
  
  # facets
  facet_wrap(. ~ scenario, ncol = 2) +
  
  
  # set scales
  scale_x_continuous(breaks = breaks_x) + 
  scale_y_continuous(breaks = c(nutrient_input * 0.05, nutrient_input, nutrient_input * 1.95),
                     labels = c("5%", "Mean", "195%")) + 
  scale_color_manual(values = line_color) + 
  
  # theming
  labs(x = "Time", y = "Abiotic nutrient subsidies") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "none", axis.text = element_blank(),
        strip.background = element_blank(), strip.text = element_text(hjust = 0))


#### Save result ####

suppoRt::save_ggplot(plot = gg_input_curves, filename = paste0("Figure-1-nutr", extension),
                     path = "04_Figures/", width = 100, height = 35,
                     units = units, dpi = dpi, overwrite = T)

# ##--------------------------------------------##
# ##    Author: Maximilian H.K. Hesselbarth     ##
# ##    Coastal Ecology and Conservation Lab    ##
# ##    University of Michigan                  ##
# ##    mhessel@umich.edu                       ##
# ##    www.github.com/mhesselbarth             ##
# ##--------------------------------------------##
# 
# # Purpose: Figure of nutrient input curves
# 
# #### Load setup ####
# 
# source("05_Various/setup.R")
# 
# #### Adapt parameters ####
# 
# amplitude_mn <- 0.95
# 
# frequency <- years
# 
# #### Stable values #### 
# 
# stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
#                                               ag_biomass = starting_values_list$ag_biomass,
#                                               parameters = parameters_list)
# 
# starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool
# 
# starting_values_list$detritus_pool <- stable_values_list$detritus_pool
# 
# #### Simulate nutrient inputs variability ####
# 
# variability <- 0.65
# 
# n <- 5
# 
# input_values_phase <- meta.arrR::simulate_nutrient_sine(n = n, max_i = max_i, frequency = frequency, 
#                                                         input_mn = nutrient_input, amplitude_mn = amplitude_mn, 
#                                                         phase_sd = variability) %>% 
#   meta.arrR::get_input_df(gamma = FALSE, long = TRUE)
# 
# input_values_noise <- meta.arrR::simulate_nutrient_noise(n = n, max_i = max_i, frequency = frequency, 
#                                                          input_mn = nutrient_input, amplitude_mn = amplitude_mn, 
#                                                          noise_sd = variability) %>% 
#   
#   meta.arrR::get_input_df(gamma = FALSE, long = TRUE)
# 
# input_values_combined <- dplyr::bind_rows(phase = input_values_phase, noise = input_values_noise, 
#                                           .id = "scenario") %>% 
#   dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise"), 
#                                   labels = c("Phase scenario", "Noise scenario")))
# 
# #### Setup plots ####
# 
# filter_factor <- 10
# 
# breaks_x <- seq(from = 0, to = max_i / filter_factor, by = max_i / years)
# 
# labels_x <- paste(seq(from = 0, to = years / filter_factor, by = 1), "years")
# 
# label_size <- 7.5
# 
# line_size <- 0.25
# 
# line_color <- MetBrewer::met.brewer(name = "Java", n = n, type = "discrete")
# 
# base_size <- 5.5
# 
# extension <- ".png"
# 
# #### Create ggplot ####
# 
# gg_input_curves <- dplyr::filter(input_values_combined, timestep < max_i / filter_factor) %>% 
#   
#   ggplot(aes(x = timestep, y = value, color = meta)) +
#   
#   # adding geoms
#   geom_hline(yintercept = nutrient_input, linetype = 2, color = "grey") +
#   # geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
#   geom_line(size = line_size) + 
#   
#   # facets
#   facet_wrap(. ~ scenario, ncol = 2) +
#   
#   
#   # set scales
#   scale_x_continuous(breaks = breaks_x, labels = labels_x) + 
#   scale_y_continuous(breaks = c(nutrient_input * 0.05, nutrient_input,nutrient_input * 1.95),
#                      labels = c("5%", "Mean", "195%")) + 
#   scale_color_manual(values = line_color) + 
#   
#   # theming
#   labs(x = "Time", y = "Abiotic nutrient subsidies") +
#   theme_classic(base_size = base_size) + 
#   theme(legend.position = "none", axis.text = element_blank(),
#         strip.background = element_blank(), strip.text = element_text(hjust = 0))
# 
# 
# #### Save result ####
# 
# suppoRt::save_ggplot(plot = gg_input_curves, filename = paste0("Figure-1", extension),
#                      path = "04_Figures/", width = 125, height = 35,
#                      units = units, dpi = dpi, overwrite = T)
# 
# # ##--------------------------------------------##
# # ##    Author: Maximilian H.K. Hesselbarth     ##
# # ##    Coastal Ecology and Conservation Lab    ##
# # ##    University of Michigan                  ##
# # ##    mhessel@umich.edu                       ##
# # ##    www.github.com/mhesselbarth             ##
# # ##--------------------------------------------##
# # 
# # # Purpose: Figure of nutrient input curves
# # 
# # #### Load setup ####
# # 
# # source("05_Various/setup.R")
# # 
# # #### Adapt parameters ####
# # 
# # amplitude_mn <- 0.95
# # 
# # frequency <- years
# # 
# # #### Stable values #### 
# # 
# # stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
# #                                               ag_biomass = starting_values_list$ag_biomass,
# #                                               parameters = parameters_list)
# # 
# # starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool
# # 
# # starting_values_list$detritus_pool <- stable_values_list$detritus_pool
# # 
# # #### Simulate nutrient inputs variability ####
# # 
# # variability <- 0.5
# # 
# # n <- 5
# # 
# # # simulate nutrient input
# # input_values_null <- meta.arrR::simulate_nutrient_sine(n = n, max_i = max_i, frequency = frequency, 
# #                                                        input_mn = nutrient_input, amplitude_mn = amplitude_mn)
# # 
# # input_values_sine <- meta.arrR::simulate_nutrient_sine(n = n, max_i = max_i, frequency = frequency, 
# #                                                        input_mn = nutrient_input, amplitude_mn = amplitude_mn, 
# #                                                        phase_sd = variability)
# # 
# # input_values_noise <- meta.arrR::simulate_nutrient_noise(n = n, max_i = max_i, frequency = frequency, 
# #                                                          input_mn = nutrient_input, amplitude_mn = amplitude_mn, 
# #                                                          noise_sd = variability)
# # 
# # #### Individuals plots ####
# # 
# # filter_factor <- 10
# # 
# # breaks_x <- seq(from = 0, to = max_i / filter_factor, by = max_i / years)
# # 
# # labels_x <- paste(seq(from = 0, to = years / filter_factor, by = 1), "years")
# # 
# # label_size <- 7.5
# # 
# # line_size <- 0.25
# # 
# # line_color <- MetBrewer::met.brewer(name = "Java", n = n, type = "discrete")
# # 
# # base_size <- 7.5
# # 
# # gg_input_null <- meta.arrR::get_input_df(input_values_null, gamma = FALSE, long = TRUE) %>% 
# #   dplyr::filter(timestep < max_i / filter_factor) %>% 
# #   ggplot(aes(x = timestep, y = value)) +
# #   geom_hline(yintercept = nutrient_input, linetype = 2, color = "grey") +
# #   # geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
# #   geom_line(size = line_size) + 
# #   # annotate(geom = "text", x = max_i / filter_factor / 2, y = nutrient_input * 1.1, 
# #   #          label = "Mean nutrient input", color = "darkgrey") +
# #   scale_x_continuous(breaks = breaks_x, labels = labels_x) + 
# #   scale_y_continuous(breaks = c(nutrient_input * 0.05, nutrient_input,nutrient_input * 1.95),
# #                      labels = c("5%", "Mean", "195%")) + 
# #   # scale_color_manual(values = line_color) + 
# #   labs(x = "", y = "") +
# #   theme_classic(base_size = base_size) + 
# #   theme(legend.position = "none", plot.margin = unit(c(t = 5.5, r = 5.5, b = 0.0, l = 5.5), "pt"), 
# #         axis.text = element_text(color = "black"))
# # 
# # gg_input_sine <- meta.arrR::get_input_df(input_values_sine, gamma = FALSE, long = TRUE) %>% 
# #   dplyr::filter(timestep < max_i / filter_factor) %>% 
# #   ggplot(aes(x = timestep, y = value, color = meta)) +
# #   geom_hline(yintercept = nutrient_input, linetype = 2, color = "grey") +
# #   # geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
# #   geom_line(size = line_size) + 
# #   scale_x_continuous(breaks = breaks_x, labels = labels_x) + 
# #   scale_y_continuous(breaks = c(nutrient_input * 0.05, nutrient_input,nutrient_input * 1.95),
# #                      labels = c("5%", "Mean", "195%")) + 
# #   scale_color_manual(values = line_color) + 
# #   labs(x = "", y = "") +
# #   theme_classic(base_size = base_size) + 
# #   theme(legend.position = "none", plot.margin = unit(c(t = 0.0, r = 1.5, b = 5.5, l = 5.5), "pt"), 
# #         axis.text = element_text(color = "black"))
# # 
# # gg_input_noise <- meta.arrR::get_input_df(input_values_noise, gamma = FALSE, long = TRUE) %>% 
# #   dplyr::filter(timestep < max_i / filter_factor) %>% 
# #   ggplot(aes(x = timestep, y = value, color = meta)) +
# #   geom_hline(yintercept = nutrient_input, linetype = 2, color = "grey") +
# #   # geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
# #   geom_line(size = line_size) + 
# #   scale_x_continuous(breaks = breaks_x, labels = labels_x) + 
# #   scale_y_continuous(breaks = c(nutrient_input * 0.05, nutrient_input,nutrient_input * 1.95),
# #                      labels = c("5%", "Mean", "195%")) + 
# #   scale_color_manual(values = line_color) + 
# #   labs(x = "", y = "") +
# #   theme_classic(base_size = base_size) + 
# #   theme(legend.position = "none", plot.margin = unit(c(t = 0.0, r = 5.5, b = 5.5, l = 1.5), "pt"), 
# #         axis.text = element_text(color = "black"))
# # 
# # #### Combine overall figure #### 
# # 
# # gg_input_var <- cowplot::plot_grid(gg_input_sine, gg_input_noise, ncol = 2, 
# #                                    labels = c("b)", "c)"), label_size = label_size)
# # 
# # gg_input_overall <- cowplot::plot_grid(gg_input_null, gg_input_var, nrow = 2, 
# #                                        labels = c("a)", ""), label_size = label_size)
# # 
# # #### Save result ####
# # 
# # suppoRt::save_ggplot(plot = gg_input_overall, filename = paste0("Figure-1", extension),
# #                      path = "04_Figures/", width = width, height = height * 0.35,
# #                      units = units, dpi = dpi, overwrite = FALSE)
