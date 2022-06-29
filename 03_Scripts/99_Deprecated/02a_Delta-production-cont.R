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

list_parameters$nutrients_diffusion <- 0.0
list_parameters$detritus_diffusion <- 0.0
list_parameters$detritus_fish_diffusion <- 0.0

list_starting$bg_biomass <- list_parameters$bg_biomass_max
list_starting$ag_biomass <- list_parameters$ag_biomass_max

n <- 1 

list_starting$pop_n <- 0

#### Stable values ####

stable_values <- arrR::get_stable_values(bg_biomass = list_starting$bg_biomass,
                                         ag_biomass = list_starting$ag_biomass,
                                         parameters = list_parameters)

list_starting$nutrients_pool <- stable_values$nutrients_pool

list_starting$detritus_pool <- stable_values$detritus_pool

#### Setup experiment ####

# create variability data.frame with all combinations 
sim_experiment <- expand.grid(amplitude_mod = seq(from = 0, to = 1, by = 0.1), 
                              freq_mn = seq(from = years * 1/4, to = years, by = 3.75))

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
  
  # result #
  
  data.frame(freq = freq_mn, amplitude = amplitude_mod, 
             delta_in = delta_in, delta_out[, c('part', "delta_out")])
  
}

#### Submit to HPC #### 

deltaprodcont_sbatch <- rslurm::slurm_apply(f = foo, params = sim_experiment, 
                                            global_objects = "globals", jobname = "deltaprodcont",
                                            nodes = nrow(sim_experiment), cpus_per_node = 1, 
                                            slurm_options = list("account" = account, 
                                                                 "partition" = "standard",
                                                                 "time" = "00:30:00", ## hh:mm::ss
                                                                 "mem-per-cpu" = "5G"),
                                            pkgs = c("meta.arrR", "dplyr", "tidyr"),
                                            rscript_path = rscript_path, 
                                            sh_template = sh_template, 
                                            submit = FALSE)

#### Collect results ####

deltaprodcont_result <- rslurm::get_slurm_out(deltaprodcont_sbatch, outtype = "table")

suppoRt::save_rds(object = deltaprodcont_result, filename = "02_nofish_delta-prod-cont.rds", 
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(deltaprodcont_sbatch)

#### Load data ####

deltaprodcont_result <- readRDS(file = "02_Data/02_nofish_delta-prod-cont.rds")

#### Pre-process data ####

deltaprodcont_result <- dplyr::mutate(deltaprodcont_result, 
                                      freq_rel = 1 - abs((freq - max(freq)) /
                                                           (max(freq) - min(freq))))

#### Create ggplot ####

x_lab <- c("", "Amplitude", "")
y_lab <- c("Frequency", "", "")

unique_parts <- c("bg_production", "ag_production", "ttl_production")

gg_continuous <- purrr::map(seq_along(unique_parts), function(i) {
  
  data_temp <- dplyr::filter(deltaprodcont_result, part == unique_parts[i]) 
  
  fill_breaks <- seq(from = min(data_temp$delta_out), 
                     to = max(data_temp$delta_out), length.out = 3)
  
  ggplot(data = data_temp) + 
    geom_raster(aes(x = amplitude, y = freq_rel, fill = delta_out)) + 
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

# gg_continuous <- ggplot(deltaprodcont_result) + 
#   geom_raster(aes(x = amplitude, y = freq, fill = delta_out_rel)) + 
#   facet_wrap(. ~ part, nrow = 1) +
#   scale_fill_continuous(name = "Delta out", type = "viridis") +
#   guides(fill = guide_colorbar(order = 1), size = "none") + 
#   labs(x = "Amplitude", y = "Frequency") +
#   theme_classic() + 
#   theme(legend.position = "bottom")

#### Save ggplot ####
suppoRt::save_ggplot(plot = gg_continuous, filename = "02_gg_nofish_delta-production-cont.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
