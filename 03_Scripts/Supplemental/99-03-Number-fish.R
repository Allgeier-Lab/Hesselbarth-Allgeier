##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose:

#### Load setup ####

source("01_Functions/setup.R")

#### Adapt parameters ####

# mem_per_cpu <- "7G"
# time <- "05:00:00" # hh:mm:ss

#### Stable values #### 

stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
                                              ag_biomass = starting_values_list$ag_biomass,
                                              parameters = parameters_list)

starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool

starting_values_list$detritus_pool <- stable_values_list$detritus_pool

#### Run model ####

starting_values_list$pop_n <- 32

# update move meta_sd parameters
parameters_list$move_meta_sd <- 0.5

# setup metaecosystems
metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                       starting_values = starting_values_list, parameters = parameters_list,
                                       dimensions = dimensions, grain = grain, 
                                       use_log = FALSE, verbose = FALSE)

# simulate nutrient input
input_temp <- meta.arrR::simulate_nutrient_noise(n = n, max_i = max_i, frequency = frequency, 
                                                 input_mn = 0.000121005, noise = 0.5, 
                                                 verbose = FALSE)

result_temp <- purrr::map_dfr(c(0.1, 0.5, 0.9), function(i) {
  
  # update move meta_sd parameters
  parameters_list$move_meta_sd <- i
  
  meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                 nutrients_input = input_temp, movement = "behav", 
                                 torus_diffusion = TRUE, max_i = max_i, min_per_i = min_per_i,
                                 seagrass_each = seagrass_each, save_each = save_each)  |> 
    meta.arrR::filter_meta(filter = c((max_i / years) * years_filter, max_i),
                           reset = TRUE, verbose = FALSE) |> 
    meta.arrR::get_abundance() |> 
    dplyr::mutate(meta = paste0("Local system ", meta), meta = factor(meta))}, .id = "connectivity") |> 
  dplyr::mutate(timestep_yr = timestep * 120 / 60 / 24 / 365 - years_filter)

#### Create ggplot ####

breaks_x <- seq(from = 0, to = 10, by = 2)

color_meta <- MetBrewer::met.brewer(name = "Java", n = n, type = "discrete")

names(color_meta) <- levels(result_temp$meta) |> as.character()

base_size <- 13.5

gg_abundance <- ggplot(data = result_temp) + 
  
  # adding geoms with nutrients per time
  geom_line(aes(x = timestep_yr, y = abundance, color = meta)) +
  
  # facet by variability
  facet_wrap(. ~ connectivity, nrow = 3,
             labeller = labeller(connectivity = c("1" = "Low connectivity", "2" = "Medium connectivity",
                                                 "3" = "High connectivity"))) +
  
  # adding box
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  
  # set scales
  scale_color_manual(name = "", values = color_meta) +
  scale_x_continuous(breaks = breaks_x) +
  
  # themes
  # guides(color = "none") +
  labs(x = "Time [years]", y = "Local fish population size") +
  theme_classic(base_size = base_size) +
  theme(legend.position = "bottom", axis.line = element_blank(),
        strip.background = element_blank(), strip.text = element_text(hjust = 0)) 

#### Save result ####

suppoRt::save_ggplot(plot = gg_abundance, filename = "Figure-S3.png",
                     path = "04_Figures/Supplemental/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = FALSE)
  
