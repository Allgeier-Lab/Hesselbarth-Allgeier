##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

#### Adapt parameters ####

default_parameters$nutrients_diffusion <- 0.0
default_parameters$detritus_diffusion <- 0.0
default_parameters$detritus_fish_diffusion <- 0.0

default_starting$bg_biomass <- default_parameters$bg_biomass_max
default_starting$ag_biomass <- default_parameters$ag_biomass_max

default_starting$pop_n <- 0

# save results only every m days
days_save <- 125 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)

save_each <- (24 / (min_per_i / 60)) * days_save

#### Stable values ####

stable_values <- arrR::get_stable_values(starting_values = default_starting,
                                         parameters = default_parameters)

default_starting$nutrients_pool <- stable_values$nutrients_pool

default_starting$detritus_pool <- stable_values$detritus_pool

input_mn <- stable_values$nutr_input

#### Setup experiment ####

itr <- 1000

# create variability data.frame with all combinations 
sim_experiment <- data.frame(amplitude = runif(n = itr, min = 0, max = 1), 
                             phase = runif(n = itr, min = 0, max = 1))

#### Create function ####

# create globals
globals <- list(n = n, max_i = max_i, input_mn = input_mn, freq_mn = freq_mn, 
                default_starting = default_starting, default_parameters = default_parameters, 
                dimensions = dimensions, grain = grain,
                min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each) 

foo <- function(amplitude, phase) {
  
  # simulate data #
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i, 
                                         starting_values = globals$default_starting, 
                                         parameters = globals$default_parameters, 
                                         dimensions = globals$dimensions, grain = globals$grain, 
                                         reefs = NULL, verbose = FALSE)
  
  # simulate input 
  input_temp <- meta.arrR::sim_nutr_input(n = globals$n, max_i = globals$max_i,
                                          variability = c(amplitude, phase),
                                          input_mn = globals$input_mn, freq_mn = globals$freq_mn, 
                                          verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutr_input = input_temp,
                                     parameters = globals$default_parameters,
                                     max_i = globals$max_i, min_per_i = globals$min_per_i, 
                                     seagrass_each = globals$seagrass_each,
                                     save_each = globals$save_each, verbose = FALSE)
  
  # filter data #
  
  # filter input
  input_temp <- meta.arrR::filter_meta(x = input_temp,
                                       filter = seq(from = globals$max_i / 2, 
                                                    to = globals$max_i,
                                                    by = globals$save_each), 
                                       verbose = FALSE)

  # fiter result
  result_temp <- meta.arrR::filter_meta(x = result_temp, filter = c(globals$max_i / 2, 
                                                                    globals$max_i))

  # calc cv #
    
  # sample variability for input
  cv_temp_in <- meta.arrR::calc_variability(x = result_temp$nutr_input, 
                                            verbose = FALSE)
  
  # sample variability for output and biomass and production
  cv_temp_out <- purrr::map_dfr(c("biomass", "production", "turnover"), function(i) {
    meta.arrR::calc_variability(x = result_temp, what = i, verbose = FALSE)
  })
  
  # create data.frame #
  
  # combine to df
  cbind(amplitude = amplitude, phase = phase, rbind(cv_temp_in, cv_temp_out))
  
}

#### Submit to HPC #### 

nofish_sbatch <- rslurm::slurm_apply(f = foo, params = sim_experiment, 
                                     global_objects = "globals", jobname = "nofish_cv",
                                     nodes = nrow(sim_experiment), cpus_per_node = 1, 
                                     slurm_options = list("account" = "jeallg1", 
                                                          "partition" = "standard",
                                                          "time" = "01:00:00", ## hh:mm::ss
                                                          "mem-per-cpu" = "5G"
                                                          ),
                                     pkgs = c("meta.arrR", "purrr"),
                                     rscript_path = rscript_path, sh_template = sh_template, 
                                     submit = FALSE)

#### Collect results ####

nofish_result <- rslurm::get_slurm_out(nofish_sbatch, outtype = "table")

rslurm::cleanup_files(nofish_sbatch)

#### Save data ####

suppoRt::save_rds(object = nofish_result, filename = "04_nofish_cv.rds", 
                  path = "02_Data/", overwrite = overwrite)

#### Load data ####

nofish_result <- readRDS(file = "02_Data/nofish_cv.rds")

#### Pre-process data ####

#### Create ggplot ####

#### Save ggplot ####
