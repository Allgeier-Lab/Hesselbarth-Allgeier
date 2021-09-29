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
days <- 25 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)

save_each <- (24 / (min_per_i / 60)) * days

max_i %% save_each

# set frequency of input peaks
freq_mn <- years / 10

# set number of repetitions
itr <- 50

# number of local metaecosystems
n <- 9

# setup extent and grain
dimensions <- c(100, 100)

grain <- 1

# create vector with parts to calc variability
part <- c("bg_biomass", "ag_biomass", "bg_production", "ag_production")

# create vector for lag
lag <- c(FALSE, FALSE, TRUE, TRUE)

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
                min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each, 
                part = part, lag = lag) 

foo <- function(amplitude, phase, globals) {
  
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
  cv_temp_in <- cbind(what = "input", part = NA,
                      meta.arrR::sample_variability(x = result_temp$nutr_input, 
                                                    verbose = FALSE))
  
  # sample variability for output and biomass and production
  cv_temp_out <- purrr::map_dfr(seq_along(globals$part), function(i) {
    cbind(what = "output", part = globals$part[i],
          meta.arrR::sample_variability(x = result_temp, 
                                        what = globals$part[i], lag = globals$lag[i], 
                                        verbose = FALSE))
  })
  
  # combine to df
  cbind(amplitude = amplitude, phase = phase, rbind(cv_temp_in, cv_temp_out))
}

#### Submit to HPC #### 

nofish_sbatch <- rslurm::slurm_apply(f = foo, params = variability_experiment, 
                                     globals = globals, jobname = "nofish",
                                     nodes = nrow(variability_experiment), cpus_per_node = 1, 
                                     slurm_options = list("account" = "jeallg1", 
                                                          "partition" = "standard",
                                                          "time" = "01:00:00", ## hh:mm::ss
                                                          "mem-per-cpu" = "12G"),
                                     pkgs = c("meta.arrR", "purrr"),
                                     rscript_path = rscript_path, sh_template = sh_template, 
                                     submit = FALSE)

#### Collect results ####

nofish_result <- rslurm::get_slurm_out(nofish_sbatch, outtype = "table")

rslurm::cleanup_files(nofish_sbatch)

#### Save data ####

suppoRt::save_rds(object = nofish_result, filename = "variability_nofish.rds", 
                  path = "02_Data/", overwrite = FALSE)

#### Load data ####

nofish_result <- readRDS(file = "02_Data/variability_nofish.rds")

nofish_result <- dplyr::mutate(nofish_result,
                               amplitude_label = dplyr::case_when(amplitude == 0 ~ "Low",
                                                                  amplitude == 0.5 ~ "Medium",
                                                                  TRUE ~ "High"),
                               phase_label = dplyr::case_when(phase == 0 ~ "Low",
                                                              phase == 0.5 ~ "Medium",
                                                              TRUE ~ "High")) %>%
  tidyr::unite("input", amplitude_label, phase_label, sep = "_", remove = FALSE) %>% 
  dplyr::mutate(input = factor(input, levels = input, 
                               labels = paste0(amplitude_label, " Amplitude // ", 
                                               phase_label, " Phase")), 
                part = factor(part, levels = c("bg_biomass", "ag_biomass", 
                                               "bg_production", "ag_production")),
                n = factor(n, ordered = TRUE))

variability_output_names <- paste(rep(levels(variability_output$input), each = 4 * 2), 
                                  rep(levels(variability_output$stat), each = 4),
                                  levels(variability_output$part), sep = "-")

#### Output variability ####

# filter only input data, get required columns to reshape to long
variability_output <- dplyr::filter(nofish_result, what == "output") %>% 
  dplyr::select(-c(what, alpha, beta, amplitude, phase)) %>% 
  tidyr::pivot_longer(cols = c(gamma, synchrony), names_to = "stat") %>% 
  dplyr::mutate(stat = factor(stat, levels = c("gamma", "synchrony")))

#### Fit decay model #####

# fit NLS model
variability_output_fit <- dplyr::mutate(variability_output, n = as.numeric(n)) %>% 
  dplyr::group_by(input, stat, part,) %>% 
  dplyr::group_split() %>%
  purrr::set_names(variability_output_names) %>%
  purrr::map(fit_nls)

# predict data
variability_output_pred <- purrr::map_dfr(variability_output_fit, predict_nls, x = 1:9, 
                                          .id = "input") %>% 
  tidyr::separate(col = input, sep = "-", into = c("input", "stat", "part")) %>% 
  dplyr::mutate(input = factor(input, levels = levels(nofish_result$input)))

# get coefficients of model
variability_output_coef <- dplyr::bind_rows(variability_output_fit, .id = "input") %>% 
  tidyr::separate(col = input, sep = "-", into = c("input", "stat", "part")) %>% 
  dplyr::mutate(input = factor(input, levels = levels(nofish_result$input)))

#### Create ggplot ####

size_text <- 2.75

x_text <- 8
y_text_a <- 0.9
y_text_b <- 0.8

digits_text <- 5

gg_variability_output <- purrr::map(parts, function(i) {
  
  # filter data by part
  output_temp <- dplyr::filter(variability_output, part == i)
  
  pred_temp <- dplyr::filter(variability_output_pred, part == i)
  
  coef_temp <- dplyr::filter(variability_output_coef, part == i)
  
  ggplot() +
    geom_hline(yintercept = 0, linetype = 2, col = "grey") +
    geom_boxplot(data = output_temp, aes(x = n, y = value, fill = stat), alpha = 0.25) +
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

overwrite <- FALSE

# save all output figures
purrr::walk(seq_along(parts), function(i) {
  
  part_temp <- stringr::str_replace(string = parts[i], pattern = "_", replacement = "-")
  
  name_temp <- paste0("gg_output_variability_", part_temp, ".png")
  
  suppoRt::save_ggplot(plot = gg_variability_output[[i]], filename = name_temp, 
                       path = "04_Figures/", width = height, height = width, dpi = dpi, 
                       units = units, overwrite = overwrite)
})
