##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

#### Setup experiment ####

stable_values <- arrR::get_stable_values(bg_biomass = list_starting$bg_biomass,
                                         ag_biomass = list_starting$ag_biomass,
                                         parameters = list_parameters)

input_mn <- stable_values$nutr_input

itr <- 50

# create variability data.frame with all combinations 
variability_experiment <- expand.grid(amplitude = seq(from = 0, to = 1, by = 0.2), 
                                      phase = seq(from = 0, to = 1, by = 0.2)) %>% 
  dplyr::slice(rep(1:n(), each = itr))

#### Create HPC function ####

foo <- function(amplitude, phase) {
  
  # simulate input
  input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                          variability = c(amplitude, phase),
                                          input_mn = input_mn, freq_mn = freq_mn)
  
  # sample cv
  cv_temp <- meta.arrR::calc_variability(input_temp)
  
  # combine to vector
  cbind(amplitude = amplitude, phase = phase, cv_temp)
  
}

#### Submit to HPC #### 

variability_sbatch <- rslurm::slurm_apply(f = foo, params = variability_experiment, 
                                          global_objects = c("n", "max_i", "input_mn", "freq_mn"),
                                          jobname = "example_cont_vari",
                                          nodes = nrow(variability_experiment), cpus_per_node = 1, 
                                          slurm_options = list("account" = account, 
                                                               "partition" = "standard",
                                                               "time" = "00:05:00"), ## hh:mm::ss
                                          pkgs = "meta.arrR",
                                          rscript_path = rscript_path, sh_template = sh_template, 
                                          submit = FALSE)

#### Collect results ####

variability_result <- rslurm::get_slurm_out(variability_sbatch, outtype = "table")

rslurm::cleanup_files(variability_sbatch)

#### Save data ####

suppoRt::save_rds(object = variability_result, filename = "01_example_cont-variability.rds", 
                  path = "02_Data/", overwrite = overwrite)

#### Load data #### 

variability_result <- readRDS("02_Data/01_example_cont-variability.rds")

#### Pre-process data ####

variability_result <- dplyr::group_by(variability_result, amplitude, phase, measure) %>% 
  dplyr::summarise(mean = mean(value), sd = sd(value), .groups = "drop")

#### Create ggplot ####

gg_continuous <- purrr::map(unique(variability_result$measure), function(i) {
  
  dplyr::filter(variability_result, measure == i) %>% 
    ggplot() + 
    geom_raster(aes(x = amplitude, y = phase, fill = mean)) + 
    geom_point(aes(x = amplitude, y = phase, size = sd), pch = 1) +
    scale_fill_continuous(name = "gamma", type = "viridis") +
    guides(fill = guide_colorbar(order = 1), size = "none") + 
    labs(x = "Variability Amplitude", y = "Variability Phase", subtitle = i) +
    coord_equal() + 
    theme_classic()
  
})

gg_continuous <- cowplot::plot_grid(plotlist = gg_continuous, nrow = 2, ncol = 2)

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_continuous, filename = "01_gg_example_cont-variability.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
