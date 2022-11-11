##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Calculate average excretion of local fish populations

#### Load setup ####

source("01_Functions/setup.R")

#### Change parameters and starting values ####



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

foo <- function(itr, pop_n) {
  
  starting_values_list$pop_n <- pop_n
  
  # create seafloor
  input_seafloor <- arrR::setup_seafloor(dimensions = dimensions, grain = grain, 
                                         reef = reef_matrix, starting_values = starting_values_list, 
                                         verbose = FALSE)
  
  # create fishpop
  input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor, 
                                       starting_values = starting_values_list, parameters = parameters_list,
                                       use_log = FALSE, verbose = FALSE)
    
  # run simulation
  result_temp <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                                      parameters = parameters_list, movement = "behav",
                                      max_i = max_i, min_per_i = min_per_i, seagrass_each = seagrass_each, 
                                      save_each = save_each, verbose = FALSE) |> 
    arrR::filter_mdlrn(filter = c((max_i / years) * years_filter, max_i), reset = TRUE)
  
  # calculate total excretion per time step averaged over all time steps
  excretion_ttl <- dplyr::select(result_temp$fishpop, timestep, id, excretion) |> 
    dplyr::group_by(timestep) %>% 
    dplyr::summarise(excretion = sum(excretion, na.rm = TRUE), .groups = "drop") |> 
    dplyr::mutate(excretion_last = dplyr::lag(excretion), 
                  excretion_diff = (excretion - excretion_last) / save_each) |> 
    dplyr::pull(excretion_diff) %>% 
    mean(na.rm = TRUE)
  
  data.frame(itr = itr, pop_n = pop_n, excretion_ttl = excretion_ttl)
  
}

globals <- c("dimensions", "grain", "reef_matrix", "starting_values_list", # input_seafloor
             "parameters_list", # input_fishpop
             "max_i", "min_per_i", "seagrass_each", "save_each", # run_simulation
             "years", "years_filter") # filter_mdlrn

input_df <- data.frame(pop_n = rep(x = c(8, median(8:128), 128), each = iterations)) |> 
  dplyr::mutate(itr = 1:dplyr::n(), .before = "pop_n")

#### Submit to HPC ####

# create .sh script
sbatch_fish <- rslurm::slurm_apply(f = foo, params = input_df,
                                   global_objects = globals, jobname = "nutr_fish",
                                   nodes = nrow(input_df), cpus_per_node = 1, 
                                   slurm_options = list("account" = account, 
                                                        "partition" = "standard",
                                                        "time" = "01:00:00", ## hh:mm::ss
                                                        "mem-per-cpu" = "7G"),
                                   pkgs = c("arrR", "dplyr"), rscript_path = rscript_path,
                                   submit = FALSE)

#### Collect results ####

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

ggplot(data = result_rds, aes(x = factor(pop_n, ordered = TRUE), y = excretion_ttl)) + 
  geom_boxplot() + 
  geom_jitter(alpha = 0.25) + 
  labs(x = "Population size", y = "Total excretion") +
  scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4), 
                     labels = function(x) round(x, 3)) +
  theme_classic()

nutrient_input_cell <- dplyr::group_by(result_rds, pop_n) |> 
  dplyr::summarise(excretion_ttl = mean(excretion_ttl), 
                   excretion_cell = excretion_ttl /  prod(dimensions)) |> 
  dplyr::mutate(level = c("low", "medium", "high"))

# save results to disk
suppoRt::save_rds(object = nutrient_input_cell, file = "02_Data/nutrient-input-cell.rds", 
                  overwrite = FALSE)

# # calculate nutrients loss parameter: nutrients_loss = nutrients_input / nutrients_pool
# (nutrients_loss <- nutrient_input_cell / stable_values_list$nutrients_pool)
