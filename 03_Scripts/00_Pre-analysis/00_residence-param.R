##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

#### Change parameters and starting values 

parameters_list$nutrients_loss <- 0.0 

# check if all parameters are present and meaningful
check_parameters(starting_values = starting_list, parameters = parameters_list)

#### Setup environment #### 

# create 5 reef cells in center of seafloor
reef_matrix <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0),
                      ncol = 2, byrow = TRUE)

# get stable nutrient/detritus values
stable_values <- arrR::get_stable_values(bg_biomass = starting_list$bg_biomass,
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

starting_list$nutrients_pool <- stable_values$nutrients_pool

starting_list$detritus_pool <- stable_values$detritus_pool

#### Setup experiment #### 

n <- 3

# one iterations equals 120 minutes
min_per_i <- 120

# run the model for ten years
years <- 1
max_i <- (60 * 24 * 365 * years) / min_per_i

# run seagrass once each day
days_seagrass <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days_seagrass

# create seafloor
input_seafloor <- arrR::setup_seafloor(dimensions = c(50, 50), grain = 1, 
                                       reef = reef_matrix, starting_values = starting_list)

res_par <- seq(from = 6, to = 102, by = 6)

result_moved <- purrr::map_dfr(seq_along(res_par), function(i) {
  
  parameters_list$move_residence <- res_par[[i]]
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                         starting_values = starting_list, parameters = parameters_list,
                                         dimensions = dimensions, grain = grain)

  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                                nutrients_input = 0.0, movement = "behav",
                                                max_i = max_i, min_per_i = min_per_i,
                                                seagrass_each = seagrass_each, save_each = 1)
  
  dplyr::bind_rows(result_temp$fishpop) %>% 
    dplyr::mutate(residence = dplyr::case_when(residence == 0 ~ 0, residence != 0 ~ 1)) %>% 
    dplyr::group_by(id) %>% 
    dplyr::group_split() %>% 
    purrr::map_dbl(function(j) {
      counts_temp <- dplyr::pull(j, residence) %>% rle()
      mean(counts_temp$lengths[which(counts_temp$values == 1)])
    }) %>% 
    mean(na.rm = TRUE) %>% 
    magrittr::multiply_by(min_per_i) %>% magrittr::divide_by(60) %>%
    dplyr::bind_cols(residence = result_temp$parameters$move_residence, 
                     residence_h = result_temp$parameters$move_residence * min_per_i / 24,
                     local_h = .)

})
