##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create figure of nutrient curves for increasing variability

#### Load setup ####

source("01_Functions/setup.R")

#### Adapt parameters ####

variability <- c(0.1, 0.5, 0.9)

years_filter <- 40

#### Simulate nutrient inputs variability ####

input_values_noise <- purrr::map_dfr(variability, function(i) {
  
  meta.arrR::simulate_nutrient_noise(n = n, max_i = max_i, frequency = frequency, 
                                     input_mn = 0.000121005, noise = i) |>  
    
    meta.arrR::get_input_df(gamma = FALSE, long = TRUE) |> 
    dplyr::filter(timestep > (max_i / years) * years_filter) |> 
    dplyr::mutate(meta = stringr::str_replace(meta, "meta_", "Local system "), 
                  meta = factor(meta))
}, .id = "variability") |> 
  dplyr::group_by(variability, meta) |> 
  dplyr::mutate(timestep_yr = 1:dplyr::n() * (min_per_i / (60 * 24 * 365)))

input_values_meta <- dplyr::group_by(input_values_noise, variability, timestep_yr) |> 
  dplyr::summarise(value = sum(value))

#### Setup plots ####

breaks_x <- seq(from = 0, to = 10, by = 2)

breaks_y <- c(0, 0.000121005, 0.000121005 * 2)

line_size <- 0.5

color_meta <- MetBrewer::met.brewer(name = "Java", n = n, type = "discrete")

names(color_meta) <- levels(input_values_noise$meta) |> as.character()

color_meta <- c(color_meta, "Meta-ecosystem" = "black")

base_size <- 13.5

#### Create ggplot local ####

gg_combined <- ggplot() + 
  
  # adding geoms with nutrients per time
  geom_line(data = input_values_noise, aes(x = timestep_yr, y = value, color = meta), 
            linetype = 2) +
  geom_line(data = input_values_meta, aes(x = timestep_yr, y = (value / 5), color = "Meta-ecosystem"), 
            linetype = 1) +

  # facet by variability
  facet_wrap(. ~ variability, nrow = 3, 
             labeller = labeller(variability = c("1" = "Low variability", "2" = "Medium variability", 
                                                 "3" = "High variability"))) +
  
  # adding box
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +

  # set scales
  scale_color_manual(name = "", values = color_meta) +
  scale_x_continuous(breaks = breaks_x) +
  scale_y_continuous(breaks = breaks_y, limits = range(breaks_y), labels = function(x) round(x, 5)) +
  
  # themes
  # guides(color = "none") +
  labs(x = "Time [years]", y = expression(paste("Abiotic nutrient subsidies [", g~h^-2~m^-2, "]"))) +
  theme_classic(base_size = base_size) +
  theme(legend.position = "bottom", axis.line = element_blank(),
        strip.background = element_blank(), strip.text = element_text(hjust = 0))

#### Save result ####

suppoRt::save_ggplot(plot = gg_combined, filename = "Figure-S1.png",
                     path = "04_Figures/Supplemental/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = FALSE)
