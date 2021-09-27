##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

source("05_Various/setup.R")

#### Load data ####

default_starting <- readRDS("02_Data/default_starting.rds")

default_parameters <- readRDS("02_Data/default_parameters.rds")

#### Basic parameters ####

# set min_per_i
min_per_i <- 120

# run the model for n years
years <- 25

max_i <- (60 * 24 * 365 * years) / min_per_i

# run seagrass only 1 day
days <- 1

seagrass_each <- (24 / (min_per_i / 60)) * days

# save results only every m days
days <- 25 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)

save_each <- (24 / (min_per_i / 60)) * days

max_i %% save_each

# set frequency of input peaks
freq_mn <- years / 10

# set number of repetitions
itr <- 50

# number of local metaecosystems
n <- 9

# create no reefs
reefs <- NULL

# # create 5 reef cells in center of seafloor
# reefs <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0),
#                 ncol = 2, byrow = TRUE)

# setup extent and grain
dimensions <- c(100, 100)

grain <- 1

#### Adapt parameters ####

default_starting$pop_n <- 0

default_parameters$nutrients_diffusion <- 0.0

default_parameters$detritus_diffusion <- 0.0

default_parameters$detritus_fish_diffusion <- 0.0

#### Stable values ####

stable_values <- arrR::get_stable_values(starting_values = default_starting,
                                       parameters = default_parameters)

default_starting$nutrients_pool <- stable_values$nutrients_pool

default_starting$detritus_pool <- stable_values$detritus_pool

input_mn <- stable_values$nutr_input

# setup metaecosystems
metasyst <- meta.arrR::setup_meta(n = n, max_i = max_i, 
                                  starting_values = default_starting, parameters = default_parameters, 
                                  dimensions = dimensions, grain = grain, reefs = reefs,)

#### Setup experiment ####

itr <- 50

# create variability data.frame with all combinations 
variability_experiment <- expand.grid(amplitude = c(0, 0.5, 1), 
                                      phase = c(0, 0.5, 1)) %>% 
  dplyr::mutate(amplitude_label = dplyr::case_when(amplitude == 0 ~ "Low", 
                                                   amplitude == 0.5 ~ "Medium", 
                                                   TRUE ~ "High"), 
                phase_label = dplyr::case_when(phase == 0 ~ "Low", 
                                               phase == 0.5 ~ "Medium", 
                                               TRUE ~ "High")) %>% 
  tidyr::unite("combined_label", amplitude_label, phase_label, sep = "_", remove = FALSE) %>% 
  dplyr::slice(rep(1:n(), each = itr))

#### Setup future plan ####

# login node -> cluster nodes -> core
login <- tweak(remote, workers = "greatlakes.arc-ts.umich.edu", user = "mhessel")

sbatch <- tweak(batchtools_slurm, template = "future_slurm.tmpl",
                resources = list(job_name = "run-meta",
                                 log_file = "run-meta.log",
                                 walltime = "01:30:00", # walltime <hh:mm:ss>
                                 mem_cpu  = "12G")) # memory per core in mb

plan(list(
  login,
  sbatch,
  sequential
))

# create globals
globals_meta <- list(n = n, max_i = max_i, variability_experiment = variability_experiment,
                     input_mn = input_mn, freq_mn = freq_mn, metasyst = metasyst, 
                     default_parameters = default_parameters, min_per_i = min_per_i, 
                     seagrass_each = seagrass_each, save_each = save_each)

#### Run model on HPC ####

metarn_nofish %<-% future.apply::future_lapply(1:nrow(variability_experiment), FUN = function(i) {
  
  result %<-% {
    
    # simulate input 
    input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                            variability = as.numeric(variability_experiment[i, 1:2]),
                                            input_mn = input_mn, freq_mn = freq_mn, 
                                            verbose = FALSE)
    
    # run model
    result_temp <- meta.arrR::run_meta(metasyst = metasyst, nutr_input = input_temp,
                                       parameters = default_parameters,
                                       max_i = max_i, min_per_i = min_per_i, 
                                       seagrass_each = seagrass_each,
                                       save_each = save_each, verbose = FALSE)
    
    # create filename
    file_name <- paste0("/home/mhessel/results/meta-rn_nofish_", i, ".rds")

    # save result explicit in folder
    saveRDS(object = result_temp, file = file_name)

    # only return string
    file_name
    
  }
}, future.globals = globals_meta, future.seed = TRUE)

#### Save data ####

overwrite <- FALSE

# get data from HPC and save to ~/Downloads/results

# get id of all results
result_id <- list.files(path = "~/Downloads/results",
                               full.names = TRUE, pattern = "^meta-rn_nofish_*") %>%
  stringr::str_sort(numeric = TRUE) %>% 
  map_int(function(i) stringr::str_extract(i, pattern = "[0-9]+") %>% as.integer)

# check if all ids are there
which(!1:nrow(variability_experiment) %in% result_id)

suppoRt::save_rds(object = variability_experiment, filename = "variability_experiment.rds", 
                  path = "02_Data/", overwrite = overwrite)
