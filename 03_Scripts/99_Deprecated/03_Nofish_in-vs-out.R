##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

source("05_Various/setup.R")

#### Load data ####

list_starting <- readRDS("02_Data/list_starting.rds")

list_parameters <- readRDS("02_Data/list_parameters.rds")

#### Adapt parameters ####

list_starting$pop_n <- 0

list_parameters$nutrients_diffusion <- 0.0

list_parameters$detritus_diffusion <- 0.0

list_parameters$detritus_fish_diffusion <- 0.0

#### Stable values ####

list_stable <- arrR::get_stable_values(bg_biomass = list_starting$bg_biomass,
                                         ag_biomass = list_starting$ag_biomass,
                                         parameters = list_parameters)

list_starting$nutrients_pool <- list_stable$nutrients_pool

list_starting$detritus_pool <- list_stable$detritus_pool

input_mn <- list_stable$nutr_input

#### Setup experiment ####

itr <- 50

# create variability data.frame with all combinations 
variability_experiment <- expand.grid(amplitude = c(0, 0.5, 1), 
                                      phase = c(0, 0.5, 1)) %>% 
  dplyr::slice(rep(1:n(), each = itr))

#### Create function ####

# create globals
globals <- list(n = n, max_i = max_i, input_mn = input_mn, freq_mn = freq_mn, 
                list_starting = list_starting, list_parameters = list_parameters, 
                dimensions = dimensions, grain = grain,
                min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each) 

foo <- function(amplitude, phase) {
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i, 
                                         starting_values = globals$list_starting, 
                                         parameters = globals$list_parameters, 
                                         dimensions = globals$dimensions, grain = globals$grain, 
                                         reef = NULL, verbose = FALSE)
  
  # simulate input 
  input_temp <- meta.arrR::sim_nutr_input(n = globals$n, max_i = globals$max_i,
                                          variability = c(amplitude, phase),
                                          input_mn = globals$input_mn, freq_mn = globals$freq_mn, 
                                          verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutr_input = input_temp,
                                     parameters = globals$list_parameters,
                                     max_i = globals$max_i, min_per_i = globals$min_per_i, 
                                     seagrass_each = globals$seagrass_each,
                                     save_each = globals$save_each, verbose = FALSE)
  
  # sample variability for input
  cv_temp_in <- meta.arrR::calc_variability(x = result_temp$nutr_input, 
                                            verbose = FALSE)
  
  # sample variability for output and biomass and production
  cv_temp_out <- purrr::map_dfr(c("biomass", "production", "turnover"), function(i) {
    meta.arrR::calc_variability(x = result_temp, what = i, verbose = FALSE)
  })
  
  out_gamma <- cv_temp_out[which(cv_temp_out$measure == "gamma"), c(1, 3)]
  
  out_synchrony <- cv_temp_out[which(cv_temp_out$measure == "synchrony"), c(1, 3)]
  
  # combine to df
  data.frame(amplitude = amplitude, phase = phase, part = out_gamma[, "part"], 
             in_gamma = cv_temp_in[3, "value"], in_synchrony = cv_temp_in[4, "value"], 
             out_gamma = out_gamma[, "value"], out_synchrony = out_synchrony[, "value"])

}

#### Submit to HPC #### 

nofish_sbatch <- rslurm::slurm_apply(f = foo, params = variability_experiment, 
                                     global_objects = "globals", jobname = "nofish_inout",
                                     nodes = nrow(variability_experiment), cpus_per_node = 1, 
                                     slurm_options = list("account" = account, 
                                                          "partition" = "standard",
                                                          "time" = "01:00:00", ## hh:mm::ss
                                                          "mem-per-cpu" = "3G"),
                                     pkgs = c("meta.arrR", "purrr"),
                                     rscript_path = rscript_path, sh_template = sh_template, 
                                     submit = FALSE)

#### Collect results ####

nofish_result <- rslurm::get_slurm_out(nofish_sbatch, outtype = "table")

rslurm::cleanup_files(nofish_sbatch)

#### Save data ####

suppoRt::save_rds(object = nofish_result, filename = "nofish_in-vs-out.rds", 
                  path = "02_Data/", overwrite = overwrite)

#### Load data ####

nofish_result <- readRDS(file = "02_Data/nofish_in-vs-out.rds")

#### Pre-process data ####

nofish_result <- dplyr::mutate(nofish_result,
                               part = factor(part, levels = c("bg_biomass", "ag_biomass",
                                                              "bg_production", "ag_production", 
                                                              "bg_turnover", "ag_turnover")))

# reshape to long gamma syncro

gg_inout <- ggplot(data = nofish_result) + 
  geom_point(aes(x = in_synchrony, y = out_synchrony), pch = 1) + 
  geom_smooth(aes(x = in_synchrony, y = out_synchrony)) + 
  facet_wrap(. ~ part, scales = "free_y", ncol = 2) +
  labs(x = "gamma variability input", y = "gamma variability output") +
  theme_classic()

suppoRt::save_ggplot(plot = gg_inout, filename = "gg_nofish_in-vs-out_synchrony.png", 
                     path = "04_Figures/", width = width, height = height, dpi = dpi, 
                     units = units, overwrite = overwrite)

