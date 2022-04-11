##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

#### Change parameters and starting values ####

parameters_list$nutrients_loss <- 0.0 

#### Setup environment #### 

# create 5 reef cells in center of seafloor
reef_matrix <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0),
                      ncol = 2, byrow = TRUE)

# get stable nutrient/detritus values
stable_values <- arrR::get_req_nutrients(bg_biomass = starting_list$bg_biomass, 
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

starting_list$nutrients_pool <- stable_values$nutrients_pool

starting_list$detritus_pool <- stable_values$detritus_pool

#### Setup experiment #### 
 
days_save <- 5 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)
save_each <- (24 / (min_per_i / 60)) * days_save

foo <- function(itr) {
  
  # create seafloor
  input_seafloor <- arrR::setup_seafloor(dimensions = dimensions, grain = grain, 
                                         reef = reef_matrix, starting_values = starting_list, 
                                         verbose = FALSE)
  
  input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor, 
                                       starting_values = starting_list, parameters = parameters_list,
                                       use_log = use_log, verbose = FALSE)
    
  result_temp <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                                      parameters = parameters_list, movement = "behav",
                                      max_i = max_i, min_per_i = min_per_i, seagrass_each = seagrass_each, 
                                      save_each = save_each, verbose = FALSE)
  
  dplyr::select(result_temp$fishpop, timestep, id, age, weight, excretion) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(excretion_last = dplyr::lag(excretion), 
                  excretion_diff = (excretion - excretion_last) / save_each) %>% 
    dplyr::group_by(timestep) %>% 
    dplyr::summarise(excretion_ttl = sum(excretion_diff)) %>% 
    dplyr::summarise(excretion_mn = mean(excretion_ttl, na.rm = TRUE), 
                     excretion_sd = sd(excretion_ttl, na.rm = TRUE)) %>% 
    dplyr::mutate(i = itr)
  
}

globals <- c("dimensions", "grain", "reef_matrix", "starting_list", 
             "parameters_list", "use_log",
             "max_i", "min_per_i", "seagrass_each", "save_each")

input_df <- data.frame(itr = 1:50)

#### Submit to HPC ####

# create .sh script
sbatch_fish <- rslurm::slurm_apply(f = foo, params = input_df,
                                   global_objects = globals, jobname = "nutr_fish",
                                   nodes = nrow(input_df), cpus_per_node = 1, 
                                   slurm_options = list("account" = account, 
                                                        "partition" = "standard",
                                                        "time" = "00:30:00", ## hh:mm::ss
                                                        "mem-per-cpu" = "7G"),
                                   pkgs = c("arrR", "dplyr"),
                                   rscript_path = rscript_path, sh_template = sh_template, 
                                   submit = FALSE)

# check results
suppoRt::rslurm_missing(sbatch_fish)

# get results as list
fish_input <- rslurm::get_slurm_out(sbatch_fish, outtype = "table")

# save results to disk
suppoRt::save_rds(object = fish_input, file = "02_Data/00_nutr_input_fish.rds", 
                  overwrite = overwrite)

# delete .sh scripts
rslurm::cleanup_files(sbatch_fish)

#### Results of simulations runs ####

result <- readr::read_rds("02_Data/00_nutr_input_fish.rds")

# nutrients_loss = nutrients_input / nutrients_pool
# 0.001936472 / 0.001936472 = 0.02344663886736086597096
(mean(result$excretion_mn) / prod(dimensions)) / stable_values$nutrients_pool

#### Run model with paramerization ####

# create seafloor
input_seafloor <- arrR::setup_seafloor(dimensions = dimensions, grain = grain, 
                                       reef = reef_matrix, starting_values = starting_list)

input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor, starting_values = starting_list, 
                                     parameters = parameters_list, use_log = use_log)

result <- purrr::map2(c(0.0, nutrient_input), c(0.0, 0.02344663886736086597096), function(i, j) {

  input_nutrients <- meta.arrR::sim_nutr_input(n = 1, max_i = max_i, input_mn = i, 
                                               freq_mn = freq_mn, amplitude_mod = 0.05)
  
  parameters_list$nutrients_loss <- j
  
  arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                       nutrients_input = input_nutrients$values$meta_1$input, 
                       parameters = parameters_list, movement = "behav", max_i = max_i, 
                       min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each)
  
})

plot(result[[1]], summarize = FALSE)
plot(result[[2]], summarize = FALSE)
