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

list_stable <- get_stable_values(starting_values = list_starting, 
                                   parameters = list_parameters)

# create variability data.frame with all combinations 
variability_experiment <- expand.grid(amplitude = c(0.05, 0.5, 1), 
                                      phase = c(0.05, 0.5, 1)) %>% 
  dplyr::slice(rep(1:n(), each = iterations))

#### Create HPC function ####

foo <- function(amplitude, phase) {
  
  # simulate input
  input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                          variability = c(amplitude, phase),
                                          input_mn = list_stable$nutr_input, freq_mn = freq_mn)
  
  # sample cv
  cv_temp <- meta.arrR::calc_variability(input_temp)
  
  # combine to vector
  cbind(amplitude = amplitude, phase = phase, cv_temp)
  
}

#### Submit to HPC #### 

variability_sbatch <- rslurm::slurm_apply(f = foo, params = variability_experiment, 
                                          global_objects = c("n", "max_i", "list_stable", "freq_mn"),
                                          jobname = "example_calc_vari",
                                          nodes = nrow(variability_experiment), cpus_per_node = 1, 
                                          slurm_options = list("account" = account, 
                                                               "partition" = "standard",
                                                               "time" = "00:15:00", ## hh:mm::ss
                                                               "mem-per-cpu" = "7G"),
                                          pkgs = "meta.arrR",
                                          rscript_path = rscript_path, sh_template = sh_template, 
                                          submit = FALSE)

#### Collect results ####

variability_result <- rslurm::get_slurm_out(variability_sbatch, outtype = "table")

rslurm::cleanup_files(variability_sbatch)

#### Save data ####

suppoRt::save_rds(object = variability_result, filename = "01_example_calc-variability.rds", 
                  path = "02_Data/", overwrite = overwrite)

#### Load data ####

variability_result <- readRDS("02_Data/01_example_calc-variability.rds")

#### Pre-process data ####

variability_result <- dplyr::mutate(variability_result, 
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

gg_variability <- ggplot(data = variability_result) +
  geom_hline(yintercept = 0, linetype = 2, col = "grey") +
  geom_boxplot(aes(x = input, y = value)) +
  facet_wrap(. ~ measure, scales = "free_x", ncol = 2) +
  labs(x = "Input classification variability", y = "Output variability") +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "bottom")

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_variability, filename = "01_gg_example_calc-variability.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
