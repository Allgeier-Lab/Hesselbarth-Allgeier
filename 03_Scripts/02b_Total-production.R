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

stable_values <- arrR::get_stable_values(starting_values = default_starting,
                                         parameters = default_parameters)

default_starting$nutrients_pool <- stable_values$nutrients_pool

default_starting$detritus_pool <- stable_values$detritus_pool

input_mn <- stable_values$nutr_input * seq(from = 1/2, to = 3/2, by = 1/4)

#### Setup experiment ####

freq_mn <- years * seq(from = 1/4, to = 1, by = 1/4)

itr <- 1000

# create variability data.frame with all combinations 
sim_experiment <- data.frame(amplitude_mod = runif(n = itr, min = 0, max = 1), 
                             freq_mn = sample(x = freq_mn, size = itr, replace = TRUE), 
                             input_mn = sample(x = input_mn, size = itr, replace = TRUE))

#### Create function ####

# create globals
globals <- list(n = n, max_i = max_i, default_starting = default_starting, 
                default_parameters = default_parameters, dimensions = dimensions, 
                grain = grain, min_per_i = min_per_i, seagrass_each = seagrass_each, 
                save_each = save_each) 

foo <- function(amplitude_mod, freq_mn, input_mn) {
  
  # simulate data #
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i, 
                                         starting_values = globals$default_starting, 
                                         parameters = globals$default_parameters, 
                                         dimensions = globals$dimensions, grain = globals$grain, 
                                         reefs = NULL, verbose = FALSE)
  
  # simulate input 
  input_temp <- meta.arrR::sim_nutr_input(n = globals$n, max_i = globals$max_i,
                                          input_mn = input_mn, freq_mn = freq_mn, 
                                          amplitude_mod = amplitude_mod,
                                          verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutr_input = input_temp,
                                     parameters = globals$default_parameters,
                                     max_i = globals$max_i, min_per_i = globals$min_per_i, 
                                     seagrass_each = globals$seagrass_each,
                                     save_each = globals$save_each, verbose = FALSE)
  
  # calc total production #
  
  # fiter result
  result_temp <- meta.arrR::filter_meta(x = result_temp, 
                                        filter = c(globals$max_i / 2, globals$max_i))
  
  prod_temp <- dplyr::summarize(dplyr::filter(result_temp$seafloor$Meta_1, 
                                              timestep == max(timestep)), 
                                bg_production = sum(bg_production),
                                ag_production = sum(ag_production))
  
  # calc total production
  prod_temp <- dplyr::mutate(prod_temp, ttl_production = bg_production + ag_production)
  
  # combine data #
  
  data.frame(freq = freq_mn, amplitude = amplitude_mod, input_mn = input_mn,
             prod_temp)
  
}

#### Submit to HPC #### 

ttlprod_sbatch <- rslurm::slurm_apply(f = foo, params = sim_experiment, 
                                      global_objects = "globals", jobname = "ttlprod",
                                      nodes = nrow(sim_experiment), cpus_per_node = 1, 
                                      slurm_options = list("account" = "jeallg1", 
                                                           "partition" = "standard",
                                                           "time" = "00:30:00", ## hh:mm::ss
                                                           "mem-per-cpu" = "5G"),
                                      pkgs = c("meta.arrR", "dplyr"),
                                      rscript_path = rscript_path, sh_template = sh_template, 
                                      submit = FALSE)

#### Collect results ####

ttlprod_result <- rslurm::get_slurm_out(ttlprod_sbatch, outtype = "table")

suppoRt::save_rds(object = ttlprod_result, filename = "02_nofish_ttl-prod.rds", 
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(ttlprod_sbatch)

#### Load data ####

ttlprod_result <- readRDS(file = "02_Data/02_nofish_ttl-prod.rds")

#### Pre-process data ####

ttlprod_result <- tidyr::pivot_longer(ttlprod_result, -c(freq, amplitude, input_mn), 
                                      names_to = "part") %>% 
  dplyr::group_by(freq, part) %>% 
  dplyr::mutate(value_rel = (value - max(value)) / max(value) * 100) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(input_mn = factor(input_mn, ordered = TRUE, labels = c("lo", "med-lo",
                                                                        "med", "med-hi", "hi")))

#### Create ggplot ####

# legend_position <- c("12.5" = "none", "25" = "none", "37.5" = "none", "50" = "right")

gg_ttl_prod <- purrr::map(sort(unique(ttlprod_result$freq)), function(i) {
  
  dplyr::filter(ttlprod_result, freq == i) %>% 
    ggplot(aes(x = amplitude, y = value, col = factor(input_mn))) +
    geom_point(pch = 19) + 
    geom_smooth() + 
    facet_wrap(. ~ part, scales = "free_y", nrow = 3) +
    scale_color_viridis_d(name = "Enrichment level", option = "C") +
    labs(x = "Amplitude", y = "Total production", subtitle = paste("Frequency: ", i, "cycles")) + 
    theme_classic() + 
    theme(legend.position = "none")
  
})

gg_ttl_prod <- cowplot::plot_grid(plotlist = gg_ttl_prod, nrow = 2, ncol = 2)

#### Save ggplot ####
suppoRt::save_ggplot(plot = gg_ttl_prod, filename = "02_gg_nofish_total-production.png", 
                     path = "04_Figures/", width = width, height = height, dpi = dpi, 
                     units = units, overwrite = overwrite)
