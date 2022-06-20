##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Simulation experiment of movement parameters and CV/Production

#### Load setup ####

source("05_Various/setup.R")

#### Adapt parameters ####

#### Stable values #### 

stable_values <- arrR::get_req_nutrients(bg_biomass = starting_list$bg_biomass,
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

starting_list$nutrients_pool <- stable_values$nutrients_pool

starting_list$detritus_pool <- stable_values$detritus_pool

# setup metaecosystems
metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                       starting_values = starting_list, parameters = parameters_list,
                                       dimensions = dimensions, grain = grain, 
                                       use_log = FALSE)

# simulate nutrient input
input_temp <- meta.arrR::simulate_nutrient_noise(n = n, max_i = max_i, frequency = years, 
                                                 input_mn = nutrient_input, amplitude_mn = 0.95)

# run model
result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                              nutrients_input = input_temp, movement = "behav",
                                              max_i = max_i, min_per_i = min_per_i,
                                              seagrass_each = seagrass_each, save_each = save_each)

# filter model run
result_filt <- meta.arrR::filter_meta(x = result_temp, filter = c((max_i / years) * years_filter, max_i), 
                                      reset = TRUE, verbose = FALSE)

result_sum <- meta.arrR::summarize_meta(result = result_filt, biomass = FALSE)[[2]]

x_i <- dplyr::select(result_sum, meta, timestep, ag_production) %>% 
  tidyr::pivot_wider(names_from = meta, values_from = ag_production, 
                     names_prefix = "meta_") %>% 
  dplyr::slice(-1) %>% 
  dplyr::select(-timestep) %>% 
  as.matrix()

x_m <- apply(x_i, MARGIN = 1, FUN = sum)

sd_i <- apply(x_i, MARGIN = 2, FUN = sd)
mn_i <- apply(x_i, MARGIN = 2, FUN = mean)

alpha_average <- mean((sd_i / mn_i))

alpha_wang <- sum(apply(x_i, MARGIN = 2, FUN = sd)) / mean(x_m)

alpha_wilcox <- mn_i / mn_m

x_i <- matrix(data = runif(n = 27000), ncol = 3)

# cv = sd / mean

alpha_a <- mean(apply(x_i, MARGIN = 2, FUN = sd) / apply(x_i, MARGIN = 2, FUN = mean))

x_m <- apply(x_i, MARGIN = 1, FUN = sum)


alpha_b <- sum(apply(x_i, MARGIN = 2, FUN = sd)) / mean(x_m)

alpha_c <- sum((apply(x_i, MARGIN = 2, FUN = mean) / mean(x_m)) * 
                 (apply(x_i, MARGIN = 2, FUN = sd) / apply(x_i, MARGIN = 2, FUN = mean)))


alpha_a
alpha_b
alpha_c
