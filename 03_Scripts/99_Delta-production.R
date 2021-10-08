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
years <- 50

max_i <- (60 * 24 * 365 * years) / min_per_i

# run seagrass only 1 day
days <- 1

seagrass_each <- (24 / (min_per_i / 60)) * days

# save results only every m days
days <- 25 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)

save_each <- (24 / (min_per_i / 60)) * days

max_i %% save_each

# set frequency of input peaks
freq_mn <- years * seq(from = 1/4, to = 1, by = 1/4)

# set number of repetitions
itr <- 50

# number of local metaecosystems
n <- 1

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

#### Setup experiment ####

itr <- 1000

# create variability data.frame with all combinations 
amp_freq_experiment <- data.frame(amplitude_mod = runif(n = itr, min = 0, max = 1), 
                                  freq_mn = sample(x = freq_mn, size = itr, replace = TRUE))

#### Create function ####

# create globals
globals <- list(n = n, max_i = max_i, default_starting = default_starting, 
                default_parameters = default_parameters, dimensions = dimensions, 
                grain = grain, input_mn = stable_values$nutr_input,
                min_per_i = min_per_i, seagrass_each = seagrass_each, 
                save_each = save_each) 

foo <- function(amplitude_mod, freq_mn) {
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i, 
                                         starting_values = globals$default_starting, 
                                         parameters = globals$default_parameters, 
                                         dimensions = globals$dimensions, grain = globals$grain, 
                                         reefs = NULL, verbose = FALSE)
  
  # simulate input 
  input_temp <- meta.arrR::sim_nutr_input(n = globals$n, max_i = globals$max_i,
                                          input_mn = globals$input_mn, freq_mn = freq_mn, 
                                          amplitude_mod = amplitude_mod,
                                          verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutr_input = input_temp,
                                     parameters = globals$default_parameters,
                                     max_i = globals$max_i, min_per_i = globals$min_per_i, 
                                     seagrass_each = globals$seagrass_each,
                                     save_each = globals$save_each, verbose = FALSE)
  
  # calculate production
  prod <- dplyr::summarize(dplyr::group_by(get_meta_production(result = result_temp, 
                                                               lag = TRUE, turnover = FALSE),
                                           meta, part), 
                           value_min = min(value, na.rm = TRUE), 
                           value_max = max(value, na.rm = TRUE), .groups = "drop")
  
  delta_in <- max(input_temp$values$Meta_1) - min(input_temp$values$Meta_1)
  
  # delta_in_rel <- delta_in / max(input_temp$values$Meta_1) 
  
  delta_out <- dplyr::mutate(prod, delta_out = value_max - value_min)
  
  data.frame(freq = freq_mn, amplitude = amplitude_mod, 
             delta_in = delta_in, delta_out[, c('part', "delta_out")])

}

#### Submit to HPC #### 

deltaprod_sbatch <- rslurm::slurm_apply(f = foo, params = amp_freq_experiment, 
                                        global_objects = "globals", jobname = "delta_prod",
                                        nodes = nrow(amp_freq_experiment), cpus_per_node = 1, 
                                        slurm_options = list("account" = "jeallg1", 
                                                             "partition" = "standard",
                                                             "time" = "00:30:00", ## hh:mm::ss
                                                             "mem-per-cpu" = "5G"),
                                        pkgs = c("meta.arrR", "dplyr"),
                                        rscript_path = rscript_path, sh_template = sh_template, 
                                        submit = FALSE)

#### Collect results ####

deltaprod_result <- rslurm::get_slurm_out(deltaprod_sbatch, outtype = "table")

suppoRt::save_rds(object = deltaprod_result, filename = "nofish_delta-prod.rds", 
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(deltaprod_sbatch)

#### Load data ####

deltaprod_result <- readRDS(file = "02_Data/nofish_delta-prod.rds")

#### Pre-process data ####

#### Create ggplot ####

gg_delta_prod <- ggplot(data = deltaprod_result, 
                        aes(x = delta_in, y = delta_out, col = factor(freq))) +
  geom_point(pch = 1) +
  geom_smooth(size = 0.25) + 
  facet_wrap(. ~ part, scales = "free_y") +
  scale_color_viridis_d(name = "Frequency", option = "A") +
  labs(x = "Delta input", y = "Delta production") + 
  theme_classic() + 
  theme(legend.position = "bottom")

#### Save ggplot ####
suppoRt::save_ggplot(plot = gg_delta_prod, filename = "gg_nofish_delta-production.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
