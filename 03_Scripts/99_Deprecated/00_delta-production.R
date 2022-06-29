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

n <- 1 

list_starting$pop_n <- 0

list_parameters$nutrients_diffusion <- 0.0
list_parameters$detritus_diffusion <- 0.0
list_parameters$detritus_fish_diffusion <- 0.0

#### Stable values ####

stable_values <- arrR::get_stable_values(bg_biomass = list_starting$bg_biomass,
                                         ag_biomass = list_starting$ag_biomass,
                                         parameters = list_parameters)

list_starting$nutrients_pool <- stable_values$nutrients_pool

list_starting$detritus_pool <- stable_values$detritus_pool

#### Setup experiment ####

freq_mn <- years * seq(from = 1/4, to = 1, by = 1/4)

itr <- 500

# create variability data.frame with all combinations
amp_freq_experiment <- data.frame(amplitude_mod = runif(n = itr, min = 0, max = 1),
                                  freq_mn = sample(x = freq_mn, size = itr, replace = TRUE))

#### Create HPC function ####

# create globals
globals <- list(n = n, max_i = max_i, list_starting = list_starting, 
                list_parameters = list_parameters, dimensions = dimensions, 
                grain = grain, input_mn = stable_values$nutr_input,
                min_per_i = min_per_i, seagrass_each = seagrass_each, 
                save_each = save_each) 

foo <- function(amplitude_mod, freq_mn) {
  
  # simulate data #
  
  # create vector to classify cycles
  cycle_temp <- unique(c(seq(from = 0, to = globals$max_i, by = globals$max_i / freq_mn), 
                         globals$max_i))
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i, 
                                         starting_values = globals$list_starting, 
                                         parameters = globals$list_parameters, 
                                         dimensions = globals$dimensions, grain = globals$grain, 
                                         reef = NULL, verbose = FALSE)
  
  # simulate input 
  input_temp <- meta.arrR::sim_nutr_input(n = globals$n, max_i = globals$max_i,
                                          input_mn = globals$input_mn, freq_mn = freq_mn, 
                                          amplitude_mod = amplitude_mod,
                                          verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutr_input = input_temp,
                                     parameters = globals$list_parameters,
                                     max_i = globals$max_i, min_per_i = globals$min_per_i, 
                                     seagrass_each = globals$seagrass_each,
                                     save_each = globals$save_each, verbose = FALSE)
  
  # calc delta input #
  
  # filter input
  input_temp <- meta.arrR::filter_meta(x = input_temp,
                                       filter = seq(from = globals$max_i / 2, to = globals$max_i,
                                                    by = globals$save_each), verbose = FALSE)
  
  # calculate difference for cycles
  delta_in <- dplyr::mutate(input_temp$values$Meta_1,
                            cycle = cut(timestep, breaks = cycle_temp, labels = FALSE,
                                        include.lowest = TRUE))
  
  # calculate delta difference per cycle
  delta_in <- dplyr::summarise(dplyr::group_by(delta_in, cycle),
                               input_min = min(input, na.rm = TRUE),
                               input_max = max(input, na.rm = TRUE),
                               delta = input_max - input_min, .groups = "drop")
  
  # calculate mean per part
  delta_in <- mean(delta_in$delta)
  
  # calc delta output #
  
  result_temp <- meta.arrR::filter_meta(x = result_temp, 
                                        filter = c(globals$max_i / 2, globals$max_i))
  
  # calculate production
  prod_temp <- dplyr::mutate(get_meta_production(result = result_temp, lag = TRUE, 
                                                 turnover = FALSE), 
                             cycle = cut(timestep, breaks = cycle_temp, labels = FALSE, 
                                         include.lowest = TRUE))
  
  # remove NA case due to lag = TRUE
  prod_temp <- prod_temp[complete.cases(prod_temp), ]
  
  # calc total prod
  prod_temp <- dplyr::mutate(tidyr::pivot_wider(prod_temp, names_from = part, 
                                                values_from = value), 
                             ttl_production = ag_production + bg_production)
  
  # reshape to long
  prod_temp <- tidyr::pivot_longer(prod_temp, -c(meta, timestep, cycle), 
                                   names_to = "part", values_to = "value")
  
  # calculate delta difference per cycle
  delta_out <- dplyr::summarise(dplyr::group_by(prod_temp, part, cycle), 
                                value_min = min(value, na.rm = TRUE), 
                                value_max = max(value, na.rm = TRUE),
                                delta = value_max - value_min, .groups = "drop")
  
  # calculate mean per part
  delta_out <- dplyr::summarise(dplyr::group_by(delta_out, part), delta_out = mean(delta))
  
  # combine data #
  
  data.frame(freq = freq_mn, amplitude = amplitude_mod, 
             delta_in = delta_in, delta_out[, c('part', "delta_out")])

}

#### Submit to HPC #### 

deltaprod_sbatch <- rslurm::slurm_apply(f = foo, params = amp_freq_experiment, 
                                        global_objects = "globals", jobname = "deltaprod",
                                        nodes = nrow(amp_freq_experiment), cpus_per_node = 1, 
                                        slurm_options = list("account" = account, 
                                                             "partition" = "standard",
                                                             "time" = "00:30:00", ## hh:mm::ss
                                                             "mem-per-cpu" = "5G"),
                                        pkgs = c("meta.arrR", "dplyr", "tidyr"),
                                        rscript_path = rscript_path, sh_template = sh_template, 
                                        submit = FALSE)

#### Collect results ####

deltaprod_result <- rslurm::get_slurm_out(deltaprod_sbatch, outtype = "table")

suppoRt::save_rds(object = deltaprod_result, filename = "02_nofish_delta-prod.rds", 
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(deltaprod_sbatch)

#### Load data ####

deltaprod_result <- readRDS(file = "02_Data/02_nofish_delta-prod.rds")

#### Pre-process data ####

#### Create ggplot ####

gg_delta_prod <- ggplot(data = deltaprod_result, 
                        aes(x = delta_in, y = delta_out, col = factor(freq))) +
  geom_point(pch = 19) +
  facet_wrap(. ~ part, scales = "free_y", nrow = 3) +
  scale_color_viridis_d(name = "Frequency", option = "C") +
  labs(x = "Delta input", y = "Delta production") + 
  theme_classic() + 
  theme(legend.position = "bottom")

#### Save ggplot ####
suppoRt::save_ggplot(plot = gg_delta_prod, filename = "02_gg_nofish_delta-production.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
