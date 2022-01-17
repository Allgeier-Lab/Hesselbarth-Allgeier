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

n <- 1 

default_starting$pop_n <- 0

#### Stable values ####

stable_values <- arrR::get_stable_values(bg_biomass = default_starting$bg_biomass,
                                         ag_biomass = default_starting$ag_biomass,
                                         parameters = default_parameters)

default_starting$nutrients_pool <- stable_values$nutrients_pool

default_starting$detritus_pool <- stable_values$detritus_pool

input_mn <- stable_values$nutr_input

#### Setup experiment ####

# create variability data.frame with all combinations 
sim_experiment <- expand.grid(amplitude_mod = seq(from = 0, to = 1, by = 0.1), 
                              freq_mn = seq(from = years * 1/4, to = years, by = 3.75))

#### Create HPC function ####

# create globals
globals <- list(n = n, max_i = max_i, default_starting = default_starting, 
                default_parameters = default_parameters, dimensions = dimensions, 
                grain = grain, input_mn = input_mn,
                min_per_i = min_per_i, seagrass_each = seagrass_each, 
                save_each = save_each) 

foo <- function(amplitude_mod, freq_mn) {
  
  # simulate data #
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i, 
                                         starting_values = globals$default_starting, 
                                         parameters = globals$default_parameters, 
                                         dimensions = globals$dimensions, grain = globals$grain, 
                                         reef = NULL, verbose = FALSE)
  
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
  
  # calc total production #
  result_temp <- meta.arrR::filter_meta(x = result_temp, 
                                        filter = c(globals$max_i / 2, globals$max_i))
  
  prod_temp <- dplyr::summarize(dplyr::filter(result_temp$seafloor$Meta_1, 
                                              timestep == max(timestep)), 
                                bg_production = sum(bg_production),
                                ag_production = sum(ag_production))
  
  # calc total production
  prod_temp <- dplyr::mutate(prod_temp, ttl_production = bg_production + ag_production)
  
  data.frame(freq = freq_mn, amplitude = amplitude_mod, prod_temp)

}

#### Submit to HPC #### 

ttlprodcont_sbatch <- rslurm::slurm_apply(f = foo, params = sim_experiment, 
                                          global_objects = "globals", jobname = "ttlprodcont",
                                          nodes = nrow(sim_experiment), cpus_per_node = 1, 
                                          slurm_options = list("account" = "jeallg1", 
                                                               "partition" = "standard",
                                                               "time" = "00:30:00", ## hh:mm::ss
                                                               "mem-per-cpu" = "5G"),
                                          pkgs = c("meta.arrR", "dplyr"),
                                          rscript_path = rscript_path, sh_template = sh_template, 
                                          submit = FALSE)

#### Collect results ####

ttlprodcont_result <- rslurm::get_slurm_out(ttlprodcont_sbatch, outtype = "table")

suppoRt::save_rds(object = ttlprodcont_result, filename = "02_nofish_ttl-prod-cont.rds", 
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(ttlprodcont_sbatch)

#### Load data #### 

ttlprodcont_result <- readRDS("02_Data/02_nofish_ttl-prod-cont.rds")

#### Pre-process data ####

ttlprodcont_result <- tidyr::pivot_longer(ttlprodcont_result, -c(freq, amplitude), 
                                          names_to = "part") %>% 
  dplyr::mutate(freq_rel = 1 - abs((freq - max(freq)) / (max(freq) - min(freq))))

#### Create ggplot ####

x_lab <- c("", "Amplitude", "")
y_lab <- c("Frequency", "", "")

unique_parts <- c("bg_production", "ag_production", "ttl_production")

gg_continuous <- purrr::map(seq_along(unique_parts), function(i) {
  
  data_temp <- dplyr::filter(ttlprodcont_result, part == unique_parts[i]) 
  
  fill_breaks <- seq(from = min(data_temp$value), 
                     to = max(data_temp$value), length.out = 3)
  
  ggplot(data = data_temp) + 
    geom_raster(aes(x = amplitude, y = freq_rel, fill = value)) + 
    scale_fill_gradientn(name = "", colours = viridis::viridis(n = 255), 
                         breaks = fill_breaks) + 
    scale_y_continuous(breaks = unique(data_temp$freq_rel),
                       labels = unique(data_temp$freq)) +
    geom_vline(xintercept = 0.5, linetype = 2, col = "lightgrey") +
    geom_hline(yintercept = 0.5, linetype = 2, col = "lightgrey") +
    scale_x_continuous(breaks = unique(data_temp$amplitude)) +
    labs(x = x_lab[i], y = y_lab[i], subtitle = unique_parts[i]) +
    coord_equal() + 
    theme_classic() + 
    theme(legend.position = "bottom", legend.key.width = unit(9, "mm"))
  
})

gg_continuous <- cowplot::plot_grid(plotlist = gg_continuous, nrow = 1)

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_continuous, filename = "02_gg_nofish_total-production-cont.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
