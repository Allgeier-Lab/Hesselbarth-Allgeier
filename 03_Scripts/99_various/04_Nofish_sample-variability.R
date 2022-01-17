##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

source("05_Various/setup.R")
source("01_Functions/fit_nls.R")
source("01_Functions/predict_nls.R")

#### Load data ####

default_starting <- readRDS("02_Data/default_starting.rds")

default_parameters <- readRDS("02_Data/default_parameters.rds")

#### Adapt parameters ####

default_starting$pop_n <- 0

default_parameters$nutrients_diffusion <- 0.0

default_parameters$detritus_diffusion <- 0.0

default_parameters$detritus_fish_diffusion <- 0.0

#### Stable values ####

stable_values <- arrR::get_stable_values(bg_biomass = default_starting$bg_biomass,
                                         ag_biomass = default_starting$ag_biomass,
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
                min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each, 
                itr = itr) 

foo <- function(amplitude, phase) {
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i, 
                                         starting_values = globals$default_starting, 
                                         parameters = globals$default_parameters, 
                                         dimensions = globals$dimensions, grain = globals$grain, 
                                         reef = NULL, verbose = FALSE)
  
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
  
  # sample variability for output and biomass and production
  cv_temp_in <- meta.arrR::sample_variability(x = result_temp$nutr_input, itr = globals$itr, 
                                              verbose = FALSE)
  
  # sample variability for output and biomass and production
  cv_temp_out <- purrr::map_dfr(c("biomass", "production", "turnover"), function(i) {
    meta.arrR::sample_variability(x = result_temp, what = i, itr = globals$itr, 
                                  verbose = FALSE)
  })
  
  # combine to df
  cbind(amplitude = amplitude, phase = phase, rbind(cv_temp_in, cv_temp_out))
}

#### Submit to HPC #### 

nofish_sbatch <- rslurm::slurm_apply(f = foo, params = variability_experiment, 
                                     global_objects = "globals", jobname = "nofish_sample_vari",
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

suppoRt::save_rds(object = nofish_result, filename = "nofish_variability.rds", 
                  path = "02_Data/", overwrite = overwrite)

#### Load data ####

nofish_result <- readRDS(file = "02_Data/nofish_sample-variability.rds")

#### Pre-process data ####

nofish_result <- dplyr::mutate(nofish_result, 
                               amplitude_label = dplyr::case_when(amplitude == 0 ~ "Low",
                                                                  amplitude == 0.5 ~ "Medium",
                                                                  TRUE ~ "High"),
                               phase_label = dplyr::case_when(phase == 0 ~ "Low",
                                                              phase == 0.5 ~ "Medium",
                                                              TRUE ~ "High")) %>%
  tidyr::unite("input", amplitude_label, phase_label, sep = "_", remove = FALSE) %>%
  # tidyr::pivot_longer(c(alpha, beta, gamma, synchrony), names_to = "stat") %>% 
  dplyr::mutate(part = factor(part), 
                n = factor(n, ordered = TRUE),
                input = factor(input, levels = input, 
                               labels = paste0(amplitude_label, " Amplitude // ", 
                                               phase_label, " Phase")), 
                stat = factor(stat, levels = c("alpha", "beta", "gamma", "synchrony")))

#### Fit decay model #####

# fit NLS model
variability_fit <- dplyr::mutate(nofish_result, n = as.numeric(n)) %>% 
  dplyr::group_by(part, input, stat) %>% 
  dplyr::group_split()

# create names
variability_names <- purrr::map_chr(variability_fit, function(i) {
  
  paste(unique(i$part), unique(i$input), unique(i$stat), sep = "-") 
  
})

# set names and fit model
variability_fit <- purrr::set_names(variability_fit, variability_names) %>%
  purrr::map(fit_nls)

# predict data
variability_pred <- purrr::map_dfr(variability_fit, predict_nls, x = 1:9, 
                                   .id = "input") %>% 
  tidyr::separate(col = input, sep = "-", into = c("part", "input", "stat")) %>% 
  dplyr::mutate(input = factor(input, levels = unique(input)))

# get coefficients of model
variability_coef <- dplyr::bind_rows(variability_fit, .id = "input") %>% 
  tidyr::separate(col = input, sep = "-", into = c("part", "input", "stat")) %>% 
  dplyr::mutate(input = factor(input, levels = unique(input)))


#### Create ggplot ####

size_text <- 2.75

x_text <- 8
y_text_a <- 0.9
y_text_b <- 0.8

digits_text <- 5

parts <- c("bg_biomass", "ag_biomass", "bg_production", "ag_production", 
           "bg_turnover", "ag_turnover")

gg_variability_output <- purrr::map(parts, function(i) {
  
  # filter data by part
  output_temp <- dplyr::filter(nofish_result, part %in% i, 
                               stat %in% c("gamma", "synchrony"))
  
  pred_temp <- dplyr::filter(variability_pred, part %in% i, 
                             stat %in% c("gamma", "synchrony"))
  
  coef_temp <- dplyr::filter(variability_coef, part %in% i, 
                             stat %in% c("gamma", "synchrony"))
  
  ggplot() +
    geom_hline(yintercept = 0, linetype = 2, col = "grey") +
    geom_boxplot(data = output_temp, aes(x = n, y = mean, fill = stat), alpha = 0.25) +
    geom_line(data = pred_temp, aes(x = x, y = value, col = stat)) +
    geom_label(data = dplyr::filter(coef_temp, stat == "gamma", term %in% c("n", "beta")),
               aes(x = x_text, y = y_text_a, col = "gamma",
                   label = paste(term, ":", round(estimate, digits = digits_text))),
               size = size_text, show.legend = FALSE) +
    geom_label(data = dplyr::filter(coef_temp, stat == "synchrony", term %in% c("n", "beta")),
               aes(x = x_text, y = y_text_b, col = "synchrony",
                   label = paste(term, ":", round(estimate, digits = digits_text))),
               size = size_text, show.legend = FALSE) +
    facet_wrap(. ~ input, nrow = 3, ncol = 3) +
    scale_fill_manual(name = "", values = c("#ED5B66", "#60BAE4")) +
    scale_color_manual(name = "", values = c("#ED5B66", "#60BAE4")) +
    labs(x = "# local ecosystems", y = "Variability value", title = i) +
    theme_classic() +
    theme(legend.position = "bottom")
  
})

#### Save ggplot ####

# save all output figures
purrr::walk(seq_along(parts), function(i) {
  
  part_temp <- stringr::str_replace(string = parts[i], pattern = "_", replacement = "-")
  
  name_temp <- paste0("gg_nofish_sample-variability_", part_temp, ".png")
  
  suppoRt::save_ggplot(plot = gg_variability_output[[i]], filename = name_temp, 
                       path = "04_Figures/", width = height, height = width, dpi = dpi, 
                       units = units, overwrite = overwrite)
})
