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
days <- 125 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)

save_each <- (24 / (min_per_i / 60)) * days

max_i %% save_each

# set frequency of input peaks
freq_mn <- years * 1/4

# set number of repetitions
itr <- 50

# number of local metaecosystems
n <- 9

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

#### Setup experiment ####

itr <- 50

# create variability data.frame with all combinations 
variability_experiment <- expand.grid(amplitude = c(0, 0.5, 1), 
                                      phase = c(0, 0.5, 1)) %>% 
  dplyr::slice(rep(1:n(), each = itr))

#### Create function ####

# create globals
globals <- list(n = n, max_i = max_i, input_mn = input_mn, freq_mn = freq_mn, 
                default_starting = default_starting, default_parameters = default_parameters, 
                dimensions = dimensions, grain = grain,
                min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each) 

foo <- function(amplitude, phase) {
  
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
  
  # sample variability for input
  cv_temp_in <- meta.arrR::calc_variability(x = result_temp$nutr_input, 
                                            verbose = FALSE)
  
  # sample variability for output and biomass and production
  cv_temp_out <- purrr::map_dfr(c("biomass", "production", "turnover"), function(i) {
    meta.arrR::calc_variability(x = result_temp, what = i, verbose = FALSE)
  })
  
  # combine to df
  cbind(amplitude = amplitude, phase = phase, rbind(cv_temp_in, cv_temp_out))
}

#### Submit to HPC #### 

nofish_sbatch <- rslurm::slurm_apply(f = foo, params = variability_experiment, 
                                     global_objects = "globals", jobname = "nofish_calc_vari",
                                     nodes = nrow(variability_experiment), cpus_per_node = 1, 
                                     slurm_options = list("account" = "jeallg1", 
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

suppoRt::save_rds(object = nofish_result, filename = "nofish_sample-variability.rds", 
                  path = "02_Data/", overwrite = overwrite)

#### Load data ####

nofish_result <- readRDS(file = "02_Data/nofish_calc-variability.rds")

#### Pre-process data ####

nofish_result <- dplyr::mutate(nofish_result,
                               amplitude_label = dplyr::case_when(amplitude == 0 ~ "Low",
                                                                  amplitude == 0.5 ~ "Medium",
                                                                  TRUE ~ "High"),
                               phase_label = dplyr::case_when(phase == 0 ~ "Low",
                                                              phase == 0.5 ~ "Medium",
                                                              TRUE ~ "High")) %>%
  tidyr::unite("input", amplitude_label, phase_label, sep = "_", remove = FALSE) %>% 
  dplyr::mutate(input = factor(input, levels = input, 
                               labels = paste0(amplitude_label, " Amplitude\n// ", 
                                               phase_label, " Phase")), 
                measure = factor(measure, levels = c("alpha", "beta", "gamma", "synchrony")))

#### Create ggplot ####

parts <- list(c("bg_biomass", "ag_biomass"), c("bg_production", "ag_production"), 
              c("bg_turnover", "ag_turnover"))

gg_variability <- purrr::map(parts, function(i) { 
  dplyr::filter(nofish_result, part %in% i) %>% 
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2, col = "grey") +
    geom_boxplot(aes(x = input, y = value, fill = part)) +
    facet_wrap(. ~ measure, scales = "free_x", ncol = 2) +
    scale_fill_manual(name = "", values = c("#ED5B66", "#60BAE4")) +
    labs(x = "Input classification variability", y = "Output variability") +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "bottom")
  
})

#### Save ggplot ####

parts <- c("biomas", "production", "turnover")

# save all output figures
purrr::walk(seq_along(parts), function(i) {
  
  part_temp <- stringr::str_replace(string = parts[i], pattern = "_", replacement = "-")
  
  name_temp <- paste0("gg_nofish_calc-variability_", part_temp, ".png")
  
  suppoRt::save_ggplot(plot = gg_variability[[i]], filename = name_temp, 
                       path = "04_Figures/", width = height, height = width, dpi = dpi, 
                       units = units, overwrite = overwrite)
})
