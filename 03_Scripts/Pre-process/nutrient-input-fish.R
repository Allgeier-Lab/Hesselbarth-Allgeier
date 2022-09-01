##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Calculate average excretion of local fish populations

#### Load setup ####

source("05_Various/setup.R")

#### Change parameters and starting values ####

parameters_list$nutrients_loss <- 0.0 

#### Setup environment #### 

# create 5 reef cells in center of seafloor
reef_matrix <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0),
                      ncol = 2, byrow = TRUE)

# get stable nutrient/detritus values
stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass, 
                                       ag_biomass = starting_values_list$ag_biomass,
                                       parameters = parameters_list)

starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool

starting_values_list$detritus_pool <- stable_values_list$detritus_pool

#### Setup experiment #### 

days_save <- 5 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)
save_each <- (24 / (min_per_i / 60)) * days_save

foo <- function(itr) {
  
  # create seafloor
  input_seafloor <- arrR::setup_seafloor(dimensions = dimensions, grain = grain, 
                                         reef = reef_matrix, starting_values = starting_values_list, 
                                         verbose = FALSE)
  
  # create fishpop
  input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor, 
                                       starting_values = starting_values_list, parameters = parameters_list,
                                       use_log = use_log, verbose = FALSE)
    
  # run simulation
  result_temp <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                                      parameters = parameters_list, movement = "behav",
                                      max_i = max_i, min_per_i = min_per_i, seagrass_each = seagrass_each, 
                                      save_each = save_each, verbose = FALSE) %>% 
    arrR::filter_mdlrn(filter = c((max_i / years) * years_filter, max_i))
  
  # calculate total excretion per time step averaged over all time steps
  excretion_ttl <- dplyr::select(result_temp$fishpop, timestep, id, excretion) %>% 
    dplyr::group_by(timestep) %>% 
    dplyr::summarise(excretion = sum(excretion, na.rm = TRUE), .groups = "drop") %>% 
    dplyr::mutate(excretion_last = dplyr::lag(excretion), 
                  excretion_diff = (excretion - excretion_last) / save_each) %>% 
    dplyr::pull(excretion_diff) %>% 
    mean(na.rm = TRUE)
  
  # calculate mean consumption per individual and time step
  consumption_mn <- dplyr::select(result_temp$fishpop, timestep, id, consumption) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(consumption_last = dplyr::lag(consumption), 
                  consumption_diff = (consumption - consumption_last) / save_each) %>% 
    dplyr::pull(consumption_diff) %>% 
    mean(na.rm = TRUE)
  
  data.frame(itr = itr, excretion_ttl = excretion_ttl, consumption_mn = consumption_mn)
  
}

globals <- c("dimensions", "grain", "reef_matrix", "starting_values_list", # input_seafloor
             "parameters_list", "use_log", # input_fishpop
             "max_i", "min_per_i", "seagrass_each", "save_each", # run_simulation
             "years", "years_filter") # filter_mdlrn

input_df <- data.frame(itr = 1:iterations)

#### Submit to HPC ####

# create .sh script
sbatch_fish <- rslurm::slurm_apply(f = foo, params = input_df,
                                   global_objects = globals, jobname = "nutr_fish",
                                   nodes = nrow(input_df), cpus_per_node = 1, 
                                   slurm_options = list("account" = account, 
                                                        "partition" = "standard",
                                                        "time" = "00:30:00", ## hh:mm::ss
                                                        "mem-per-cpu" = "7G"),
                                   pkgs = c("arrR", "dplyr"), rscript_path = rscript_path,
                                   submit = FALSE)

# check results
suppoRt::rslurm_missing(sbatch_fish)

# get results as list
fish_input <- rslurm::get_slurm_out(sbatch_fish, outtype = "table")

# save results to disk
suppoRt::save_rds(object = fish_input, file = "02_Data/nutrient-input-fish.rds", 
                  overwrite = FALSE)

# delete .sh scripts
rslurm::cleanup_files(sbatch_fish)

#### Results of simulations runs ####

result_rds <- readr::read_rds("02_Data/nutrient-input-fish.rds")

# calculate nutrient input each cell so that total input equals total excretion
(nutrient_input_cell <- mean(result_rds$excretion_ttl) / prod(dimensions))
# 4.67148e-05

# calculate mean consumption per individual
mean(result_rds$consumption_mn) * 8
mean(result_rds$consumption_mn) * 128

# # calculate nutrients loss parameter: nutrients_loss = nutrients_input / nutrients_pool
# (nutrients_loss <- nutrient_input_cell / stable_values_list$nutrients_pool)
# # 0.02423493
