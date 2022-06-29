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

list_parameters$nutrients_loss <- 0.0 

#### Setup environment #### 

# create 5 reef cells in center of seafloor
matrix_reef <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0),
                      ncol = 2, byrow = TRUE)

# get stable nutrient/detritus values
stable_values <- arrR::get_req_nutrients(bg_biomass = list_starting$bg_biomass, 
                                         ag_biomass = list_starting$ag_biomass,
                                         parameters = list_parameters)

list_starting$nutrients_pool <- stable_values$nutrients_pool

list_starting$detritus_pool <- stable_values$detritus_pool

#### Setup experiment #### 
 
days_save <- 5 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)
save_each <- (24 / (min_per_i / 60)) * days_save

foo <- function(itr) {
  
  # create seafloor
  input_seafloor <- arrR::setup_seafloor(dimensions = dimensions, grain = grain, 
                                         reef = matrix_reef, starting_values = list_starting, 
                                         verbose = FALSE)
  
  input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor, 
                                       starting_values = list_starting, parameters = list_parameters,
                                       use_log = use_log, verbose = FALSE)
    
  result_temp <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                                      parameters = list_parameters, movement = "behav",
                                      max_i = max_i, min_per_i = min_per_i, seagrass_each = seagrass_each, 
                                      save_each = save_each, verbose = FALSE)
  
  dplyr::select(result_temp$fishpop, timestep, id, age, weight, excretion) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(excretion_last = dplyr::lag(excretion), 
                  excretion_diff = (excretion - excretion_last) / save_each) %>%
    dplyr::ungroup() %>% 
    dplyr::summarise(excretion_mn = mean(excretion_diff, na.rm = TRUE), 
                     excretion_sd = sd(excretion_diff, na.rm = TRUE)) %>% 
    dplyr::mutate(i = itr)
  
}

globals <- c("dimensions", "grain", "matrix_reef", "list_starting", 
             "list_parameters", "use_log",
             "max_i", "min_per_i", "seagrass_each", "save_each")

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
                                   pkgs = c("arrR", "dplyr"),
                                   rscript_path = rscript_path, sh_template = sh_template, 
                                   submit = FALSE)

# check results
suppoRt::rslurm_missing(sbatch_fish)

# get results as list
fish_input <- rslurm::get_slurm_out(sbatch_fish, outtype = "table")

# save results to disk
suppoRt::save_rds(object = fish_input, file = "02_Data/01-nutr-input-fish.rds", 
                  overwrite = overwrite)

# delete .sh scripts
rslurm::cleanup_files(sbatch_fish)

#### Results of simulations runs ####

result <- readr::read_rds("02_Data/01-nutr-input-fish.rds")

# calculate mean total input of fish population each timestep
total_excretion <- mean(result$excretion_mn) * list_starting$pop_n

# calculate nutrient input each cell so that total input equals total excretion
(nutrient_input_cell <- total_excretion / prod(dimensions))
# 4.584905e-05

# calculate nutrients loss parameter: nutrients_loss = nutrients_input / nutrients_pool
(nutrients_loss <- nutrient_input_cell / stable_values$nutrients_pool)
# 0.02367659

#### Run model with paramerization ####

# create seafloor
input_seafloor <- arrR::setup_seafloor(dimensions = dimensions, grain = grain, 
                                       reef = matrix_reef, starting_values = list_starting)

input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor, starting_values = list_starting, 
                                     parameters = list_parameters, use_log = use_log)

result <- purrr::map2(c(0.0, nutrient_input_cell), c(0.0, nutrients_loss), function(i, j) {

  input_nutrients <- meta.arrR::simulate_nutr_input(n = 1, max_i = max_i, input_mn = i, 
                                                    frequency = years, amplitude_mod = 0.05)
  
  list_parameters$nutrients_loss <- j
  
  arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                       nutrients_input = input_nutrients$values$meta_1$input, 
                       parameters = list_parameters, movement = "behav", max_i = max_i, 
                       min_per_i = min_per_i, seagrass_each = seagrass_each, 
                       save_each = save_each, verbose = TRUE)}) %>% 
  purrr::set_names(c("no_input", "input"))

summarize <- FALSE

plot(result$no_input, summarize = summarize)

plot(result$input, summarize = summarize)
