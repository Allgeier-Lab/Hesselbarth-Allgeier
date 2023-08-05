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

#### Stable values #### 

stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
                                              ag_biomass = starting_values_list$ag_biomass,
                                              parameters = parameters_list)

starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool

starting_values_list$detritus_pool <- stable_values_list$detritus_pool

#### Run simulation model ####

starting_values_list$pop_n <- 64

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

# run model
result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                              nutrients_input = input_temp, movement = "behav", 
                                              torus_diffusion = TRUE, max_i = max_i, min_per_i = min_per_i,
                                              seagrass_each = seagrass_each, save_each = save_each)  |> 
  meta.arrR::filter_meta(filter = c((max_i / years) * years_filter, max_i),
                         reset = TRUE, verbose = FALSE)

#### Create ggplot ####

base_size <- 13.5

gg_biomass <- purrr::map_dfr(result_temp$seafloor, function(i) dplyr::filter(i, timestep == max(timestep)), 
               .id = "meta") |> 
  dplyr::mutate(meta = stringr::str_replace(meta, pattern = "Meta_", replacement = "Local ecosystem")) |> 
  dplyr::select(meta, x, y, ag_biomass) |> 
  ggplot() +
  
  # adding raster geom
  geom_raster(aes(x = x, y = y, fill = ag_biomass)) + 
  
  # faceet wrap singlular meta ecosystems
  facet_wrap(. ~ meta, nrow = 2, ncol = 3) + 
  
  # scale fill
  scale_fill_gradientn(name = expression(paste("AG biomass [", gDW~d^-1~m^-2, "]")),
                       colors = rev(MetBrewer::met.brewer(name = "Hokusai1", n = 255))) +
  coord_fixed()+
  
  # set themes
  theme_classic(base_size = base_size) +
  theme(legend.position = c(0.85, 0.25), axis.text = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(),
        strip.background = element_blank(), strip.text = element_text(hjust = 0)) 

#### Save result ####

suppoRt::save_ggplot(plot = gg_biomass, filename = "Figure-S4.png",
                     path = "04_Figures/Supplemental/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = FALSE)

