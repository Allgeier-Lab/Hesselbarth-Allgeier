##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/get_modifier.R")

#### Adapt parameters ####

n <- 1

default_starting$pop_n <- 0

default_parameters$nutrients_diffusion <- 0.0
default_parameters$detritus_diffusion <- 0.0
default_parameters$detritus_fish_diffusion <- 0.0

# default_parameters$seagrass_thres <- 1/3

#### Stable values ####

stable_values <- arrR::get_stable_values(bg_biomass = default_starting$bg_biomass,
                                         ag_biomass = default_starting$ag_biomass,
                                         parameters = default_parameters)

default_starting$nutrients_pool <- stable_values$nutrients_pool

default_starting$detritus_pool <- stable_values$detritus_pool

#### Simulate input ####

# itr <- 50

amplitude_lvls <- c(low = 0.05, medium = 0.5, high = 1)

enrichment_lvls <- c(low = 0.5, medium = 1.0, high = 1.5)

sim_experiment <- expand.grid(amplitude = amplitude_lvls, enrichment = enrichment_lvls) # %>% 
  # dplyr::slice(rep(x = 1:dplyr::n(), each = itr))

#### Setup HPC function ####

globals <- list(n = n, max_i = max_i, default_starting = default_starting, 
                default_parameters = default_parameters, dimensions = dimensions, 
                grain = grain, nutr_input = stable_values$nutr_input, freq_mn = freq_mn, 
                min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each) 

foo <- function(amplitude, enrichment) {
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i,
                                         starting_values = globals$default_starting,
                                         parameters = globals$default_parameters,
                                         dimensions = globals$dimensions, grain = globals$grain,
                                         reef = NULL, verbose = FALSE)
  
  # simulate input
  input_temp <- meta.arrR::sim_nutr_input(n = globals$n, max_i = globals$max_i,
                                          amplitude_mod = amplitude,
                                          input_mn = globals$nutr_input * enrichment,
                                          freq_mn = globals$freq_mn)
  
  # run model
  result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutr_input = input_temp,
                                     parameters = globals$default_parameters,
                                     max_i = globals$max_i, min_per_i = globals$min_per_i,
                                     seagrass_each = globals$seagrass_each,
                                     save_each = globals$save_each, verbose = FALSE)
  
  # # filter only second half of timesteps
  # result_temp <- meta.arrR::filter_meta(x = result_temp, filter = c(globals$max_i / 2,
  #                                                                   globals$max_i) ,
  #                                       verbose = FALSE)
  
  # calc CV
  cv <- purrr::map_dfr(c("biomass", "production"), function(i) {
    meta.arrR::calc_variability(x = result_temp, what = i, lag = TRUE,
                                verbose = FALSE)})
  
  # add parameter values
  cv <- cbind(amplitude = amplitude, enrichment = enrichment, cv)
  
  # only alpha needed because n = 1
  cv <- cv[cv$measure %in% c("alpha"), ]
  
  # only get last timestep
  result_temp <- meta.arrR::filter_meta(x = result_temp, filter = globals$max_i, 
                                        verbose = FALSE)
  
  # calculate total biomass and production
  biomass <- data.frame(amplitude = rep(x = amplitude, times = 4), 
                        enrichment = rep(x = enrichment, times = 4),
                        part = c("ag_biomass", "bg_biomass", 
                                 "ag_production", "bg_production"),
                        measure = rep(x = "total", times = 4),
                        value = c(sum(result_temp$seafloor$Meta_1$ag_biomass), 
                                  sum(result_temp$seafloor$Meta_1$bg_biomass), 
                                  sum(result_temp$seafloor$Meta_1$ag_production), 
                                  sum(result_temp$seafloor$Meta_1$bg_production)))
  
  # combine to one df
  rbind(cv, biomass, make.row.names = FALSE)
  
}

#### Submit to HPC #### 

variability_sbatch <- rslurm::slurm_apply(f = foo, params = sim_experiment, 
                                          global_objects = "globals", jobname = "var_enrich",
                                          nodes = nrow(sim_experiment), cpus_per_node = 1, 
                                          slurm_options = list("account" = "jeallg1", 
                                                               "partition" = "standard",
                                                               "time" = "02:00:00", ## hh:mm::ss
                                                               "mem-per-cpu" = "7G", 
                                                               "exclude" = exclude_nodes),
                                          pkgs = c("meta.arrR"),
                                          rscript_path = rscript_path, sh_template = sh_template, 
                                          submit = FALSE)

#### Collect results ####

var_enrich_result <- rslurm::get_slurm_out(variability_sbatch, outtype = "table")

row.names(var_enrich_result) <- 1:nrow(var_enrich_result)

suppoRt::save_rds(object = var_enrich_result, filename = "03_var_enrich.rds", 
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(variability_sbatch)

#### Load data ####

var_enrich_result <- readRDS("02_Data/03_var_enrich.rds")

var_enrich_result <- dplyr::mutate(var_enrich_result ,
                                   amplitude = factor(amplitude, levels = amplitude_lvls,
                                                      labels = c("Low", "Medium", "High")), 
                                   enrichment = factor(enrichment, levels = enrichment_lvls, 
                                                       labels = c("Low enrichment", "Medium enrichment", "High enrichment")),
                                   part = factor(part, levels = c("ag_biomass", "bg_biomass", 
                                                                  "ag_production", "bg_production")), 
                                   measure = factor(measure, levels = c("alpha", "total")))

#### Create ggplot CV ####

gg_cv_biomass <- dplyr::filter(var_enrich_result, part %in% c("ag_biomass", "bg_biomass"), 
                               measure == "alpha") %>%
  ggplot() + 
  geom_point(aes(x = amplitude, y = value, group = enrichment, col = enrichment)) + 
  geom_line(aes(x = amplitude, y = value, group = enrichment, col = enrichment), alpha = 0.5) +
  facet_wrap(. ~ part, scales = "fixed") +
  scale_color_viridis_d(name = "", option = "C") +
  labs(y = "CV", x = "") +
  theme_classic() + 
  theme(legend.position = "none")

gg_cv_prod <- dplyr::filter(var_enrich_result, part %in% c("ag_production", "bg_production"), 
                            measure == "alpha") %>%
  ggplot() + 
  geom_point(aes(x = amplitude, y = value, group = enrichment, col = enrichment)) + 
  geom_line(aes(x = amplitude, y = value, group = enrichment, col = enrichment), alpha = 0.5) +
  facet_wrap(. ~ part, scales = "fixed") +
  scale_color_viridis_d(name = "", option = "C") +
  labs(y = "CV", x = "Amplitude level") +
  theme_classic() + 
  theme(legend.position = "bottom")

gg_cv <- cowplot::plot_grid(gg_cv_biomass, gg_cv_prod, nrow = 2, ncol = 1)

#### Create ggplot total ####

gg_ttl_biomass <- dplyr::filter(var_enrich_result, part %in% c("ag_biomass", "bg_biomass"),
                                measure == "total") %>%
  ggplot() + 
  geom_point(aes(x = amplitude, y = value, group = enrichment, col = enrichment)) + 
  geom_line(aes(x = amplitude, y = value, group = enrichment, col = enrichment), alpha = 0.5) +
  facet_wrap(. ~ part, scales = "fixed") + 
  scale_color_viridis_d(name = "", option = "C") +
  labs(y = "Absolute value", x = "") +
  theme_classic() + 
  theme(legend.position = "none")

gg_ttl_prod <- dplyr::filter(var_enrich_result, part %in% c("ag_production", "bg_production"), 
                             measure == "total") %>%
  ggplot() + 
  geom_point(aes(x = amplitude, y = value, group = enrichment, col = enrichment)) + 
  geom_line(aes(x = amplitude, y = value, group = enrichment, col = enrichment), alpha = 0.5) +
  facet_wrap(. ~ part, scales = "fixed") + 
  scale_color_viridis_d(name = "", option = "C") +
  labs(y = "Absolute value", x = "Amplitude level") +
  theme_classic() + 
  theme(legend.position = "bottom")

gg_ttl <- cowplot::plot_grid(gg_ttl_biomass, gg_ttl_prod, nrow = 2, ncol = 1)

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_cv, filename = "03_gg_nofish_enrich_cv.png", 
                     path = "04_Figures/", width = width, height = height, dpi = dpi, 
                     units = units, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_ttl, filename = "03_gg_nofish_enrich_ttl.png", 
                     path = "04_Figures/", width = width, height = height, dpi = dpi, 
                     units = units, overwrite = overwrite)
