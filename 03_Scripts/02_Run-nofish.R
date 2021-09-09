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
days <- 25

save_each <- (24 / (min_per_i / 60)) * days

max_i %% save_each

# set frequency of input peaks
freq_mn <- 5

# set number of repetitions
itr <- 50

# number of local metaecosystems
n <- 7

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
metasyst <- setup_meta(n = n, max_i = max_i, dimensions = dimensions, grain = grain, reefs = reefs,
                       starting_values = default_starting, parameters = default_parameters)

#### Setup experiment ####

# define parts to sample CV
part <- c("bg_biomass", "ag_biomass", "bg_production", "ag_production")

# set low and high variability treatment levels
variability_lo <- 0.0

variability_hi <- 0.5

# combine full grid
variability_full <- expand.grid(ampl = c(variability_lo, variability_hi), 
                                freq = c(variability_lo, variability_hi)) %>% 
  dplyr::bind_cols(label = factor(x = c("lo_lo", "hi_lo", "lo_hi", "hi_hi"), 
                                  levels = c("lo_lo", "hi_lo", "lo_hi", "hi_hi"), 
                                  labels = c("Low amplitude - Low frequency",
                                             "High amplitude - Low frequency",
                                             "Low amplitude - High frequency",
                                             "High amplitude - High frequency")))

# create experiment data.frame
variability_exp <- dplyr::bind_cols(
  ampl = rep(variability_full$ampl, each = itr),
  freq = rep(variability_full$freq, each = itr), 
  label = rep(variability_full$label, each = itr), 
  itr = rep(x = 1:itr, times = nrow(variability_full))
)

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
globals_meta <- list(n = n, max_i = max_i, variability_exp = variability_exp,
                     input_mn = input_mn, freq_mn = freq_mn, metasyst = metasyst, 
                     default_parameters = default_parameters, min_per_i = min_per_i, 
                     seagrass_each = seagrass_each, save_each = save_each, part = part)

#### Run model on HPC ####

sampled_cv_nofish %<-% future.apply::future_lapply(1:nrow(variability_exp), FUN = function(i) {
  
  result %<-% {
    
    # simulate input 
    input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                            variability = as.numeric(variability_exp[i, 1:2]),
                                            input_mn = input_mn, freq_mn = freq_mn)
    
    # run model
    result_temp <- meta.arrR::run_meta(metasyst = metasyst, nutr_input = input_temp,
                                       parameters = default_parameters,
                                       max_i = max_i, min_per_i = min_per_i, 
                                       seagrass_each = seagrass_each,
                                       save_each = save_each, verbose = FALSE)
    
    # sample CV for increasing scale
    cv_temp <- purrr::map_dfr(part, function(j) meta.arrR::sample_cv_gamma(result_temp, what = j)) 
    
    # add information of run
    cv_temp <- dplyr::bind_cols(label = variability_exp[i, 3], itr = variability_exp[i, 4], 
                                part = rep(part, each = n), cv_temp)
    
    # reshape to long format
    cv_temp <- tidyr::pivot_longer(cv_temp, -c(label, itr, part, n), names_to = "stat")
    
    # create filename
    file_name <- paste0("/home/mhessel/results/meta-rn_nofish_", i, ".rds")

    # save result explicit in folder
    saveRDS(object = cv_temp, file = file_name)
    
    # only return string
    file_name
    
  }
}, future.globals = globals_meta, future.seed = TRUE)

#### Save data ####

overwrite <- FALSE

# get data from HPC and save to ~/Downloads/results
sampled_cv_nofish <- list.files(path = "~/Downloads/results",
                               full.names = TRUE, pattern = "^meta-rn_nofish_*") %>%
  stringr::str_sort(numeric = TRUE) %>% 
  purrr::map_dfr(readRDS) %>% 
  dplyr::mutate(n = factor(n, ordered = TRUE))

# check if all runs are present (input treatments x local systems x ag/bg x stat)
nrow(variability_exp) * n * length(part) * 4 == nrow(sampled_cv_nofish)

suppoRt::save_rds(object = variability_exp, filename = "variability_exp.rds", 
                  path = "02_Data/", overwrite = overwrite)

suppoRt::save_rds(object = sampled_cv_nofish, filename = "sampled_cv_nofish.rds", 
                  path = "02_Data/", overwrite = overwrite)
