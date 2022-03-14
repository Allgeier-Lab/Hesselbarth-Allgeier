##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### default parameters ####
default_parameters <- list(
  
  # belowground biomass
  bg_biomass_min = 275.89,
  bg_biomass_max = 933.03,
  bg_v_max = 9.8,
  bg_k_m = 178.1,
  bg_gamma = 0.0082,
  
  # aboveground biomass
  ag_biomass_min = 8.87,
  ag_biomass_max = 193.01,
  ag_v_max = 8.1,
  ag_k_m = 12.6,
  ag_gamma = 0.0144,
  
  # seagrass
  seagrass_thres = -1/4,
  seagrass_slope = 2.0,
  seagrass_slough = 0.01,
  
  # nutrients
  nutrients_diffusion = 2/3,
  nutrients_loss = 0.01,
  
  # detritus
  detritus_mineralization = 0.01,
  detritus_diffusion = 1/3,
  detritus_fish_decomp = 0.5,
  detritus_fish_diffusion = 1/3,
  detritus_loss = 0.0,
  
  # fishpop movement
  move_mean = 10.0,
  move_var = 5.0,
  move_border = 2.0,
  move_reef = 1.0,
  move_return = 15.0,
  move_residence = 0.0,
  move_residence_var = 0.0,
  move_lambda = 0.0,
  
  # fishpop dimensions
  pop_a = 0.0121,
  pop_b = 3.161,
  pop_k = 0.2,
  pop_linf = 41.6,
  pop_n_body = 0.02999,
  
  # fishpop reserves
  pop_reserves_max = 0.25,
  pop_reserves_thres_mean = 0.05,
  pop_reserves_thres_var = 0.0,
  pop_reserves_consump = 0.1,
  
  # fishpop respiration
  resp_intercept = 0.0108,
  resp_slope = -0.2,
  resp_temp_low = 2.1,
  resp_temp_optm = 36,
  resp_temp_max = 40
)

#### default starting values ####

# ag <- default_parameters$ag_biomass_min +
#   (default_parameters$ag_biomass_max - default_parameters$ag_biomass_min) * 1/3
# 
# bg <- default_parameters$bg_biomass_min +
#   (default_parameters$bg_biomass_max - default_parameters$bg_biomass_min) * 1/3

default_starting <- list(
  
  # biomass (mean field data)
  bg_biomass = 547.01948, # bg,
  ag_biomass = 37.84915, # ag,
  
  # nutrients
  nutrients_pool = 0.0,
  detritus_pool = 0.0,
  
  # fishpop related
  pop_n = 8,
  pop_mean_size = 9.0,
  pop_var_size = 10.0
)

((default_starting$bg_biomass - default_parameters$bg_biomass_min) /
  (default_parameters$bg_biomass_max - default_parameters$bg_biomass_min)) * 100

((default_starting$ag_biomass - default_parameters$ag_biomass_min) /
    (default_parameters$ag_biomass_max - default_parameters$ag_biomass_min)) * 100

#### save results ####

overwrite <- TRUE

suppoRt::save_rds(object = default_parameters, filename = "00_default_parameters.rds", 
                  path = "02_Data/", overwrite = overwrite)

suppoRt::save_rds(object = default_starting, filename = "00_default_starting.rds", 
                  path = "02_Data/", overwrite = overwrite)
