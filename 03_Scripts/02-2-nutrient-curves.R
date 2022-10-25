##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create figure of nutrient curves for increasing variability

#### Load setup ####

source("05_Various/setup.R")

#### Adapt parameters ####

amplitude_mn <- 0.95

frequency <- years

variability <- 0.75

n <- 5

years_filter <- 45

#### Stable values #### 

stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
                                              ag_biomass = starting_values_list$ag_biomass,
                                              parameters = parameters_list)

starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool

starting_values_list$detritus_pool <- stable_values_list$detritus_pool

#### Simulate nutrient inputs variability ####

input_values_noise <- meta.arrR::simulate_nutrient_noise(n = n, max_i = max_i, frequency = frequency, 
                                                         input_mn = nutrient_input, amplitude_mn = amplitude_mn, 
                                                         noise_sd = variability) |>  
  
  meta.arrR::get_input_df(gamma = FALSE, long = TRUE) |> 
  dplyr::filter(timestep > (max_i / years) * years_filter)

input_values_meta <- dplyr::group_by(input_values_noise, timestep) |> 
  dplyr::summarise(value = sum(value)) |> 
  dplyr::filter(timestep > (max_i / years) * years_filter)

#### Setup plots ####

breaks_x <- seq(from = (max_i / years) * years_filter, to = max_i, by = max_i / years)

line_size <- 0.25

line_color <- MetBrewer::met.brewer(name = "Java", n = n, type = "discrete")

base_size <- 7.5

extension <- ".png"

#### Create ggplot local ####

gg_local <- ggplot(data = input_values_noise, aes(x = timestep, y = value, color = meta)) +

  # adding geoms
  geom_hline(yintercept = nutrient_input, linetype = 2, size = line_size, color = "grey") +
  geom_line(size = line_size) +

  # set scales
  scale_x_continuous(breaks = breaks_x) +
  scale_y_continuous(limits = c(nutrient_input * 0.05, nutrient_input * 1.95),
                     breaks = c(nutrient_input * 0.05, nutrient_input, nutrient_input * 1.95),
                     labels = c("5%", "Mean", "195%")) +
  scale_color_manual(values = line_color) +

  # theming
  labs(x = "", y = "Abiotic nutrient subsidies") +
  theme_classic(base_size = base_size) +
  theme(legend.position = "none", axis.text = element_blank())

gg_meta <- ggplot(data = input_values_meta, aes(x = timestep, y = value)) +
  
  # adding geoms
  geom_hline(yintercept = nutrient_input * n, linetype = 2, size = line_size, color = "grey") +
  geom_line(size = line_size, color = "#32b2daff") +
  
  # set scales
  scale_x_continuous(breaks = breaks_x) +
  scale_y_continuous(limits = c(nutrient_input * 0.05, nutrient_input * 1.95) * n,
                     breaks = c(nutrient_input * 0.05, nutrient_input, nutrient_input * 1.95) * n,
                     labels = c("5%", "Mean", "195%")) +
  scale_color_manual(values = line_color) +
  
  # theming
  labs(x = "", y = "") +
  theme_classic(base_size = base_size) +
  theme(legend.position = "none", axis.text = element_blank())

gg_input_curves <- cowplot::plot_grid(gg_local, gg_meta, ncol = 2)

gg_input_curves <- cowplot::ggdraw(gg_input_curves, ylim = c(-0.05, 1.0)) + 
  cowplot::draw_label("Time", x = 0.5, y = 0, angle = 0, size = base_size)

#### Save result ####

suppoRt::save_ggplot(plot = gg_input_curves, filename = paste0("Figure-1-nutr", extension),
                     path = "04_Figures/", width = 85, height = 35,
                     units = units, dpi = dpi, overwrite = FALSE)
