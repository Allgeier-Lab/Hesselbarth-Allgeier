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

# nothing to change here

#### Stable values #### 

list_stable <- arrR::get_req_nutrients(bg_biomass = list_starting$bg_biomass,
                                         ag_biomass = list_starting$ag_biomass,
                                         parameters = list_parameters)

list_starting$nutrients_pool <- list_stable$nutrients_pool / 2

list_starting$detritus_pool <- list_stable$detritus_pool / 2

#### Setup experiment ####

# create random sd for each iteration of amplitude levels and stochastic treatments
amplitude_sd <- runif(n = iterations * length(amplitude_levels), min = 0.0, max = 1.0)

phase_sd <- runif(n = iterations * length(amplitude_levels), min = 0.0, max = 1.0)

null_sd <- rep(x = 0, times = iterations * length(amplitude_levels))

# create data.frame with all combinations
experiment_df <- tibble::tibble(amplitude_mean = rep(x = rep(x = amplitude_levels, 
                                                             each = iterations), times = 3), 
                                amplitude_sd = c(amplitude_sd, null_sd, amplitude_sd), 
                                phase_sd = c(null_sd, phase_sd, phase_sd), 
                                stochastic = rep(x = c("amplitude", "phase", "both"), 
                                                 each = iterations * 3))

# setup metaecosystems
metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = matrix_reef,
                                       starting_values = list_starting, parameters = list_parameters,
                                       dimensions = dimensions, grain = grain, use_log = use_log, 
                                       verbose = FALSE)

# setup HPC function
foo_hpc <- function(amplitude_mean, amplitude_sd, phase_sd, stochastic) {
  
  library(dplyr)
  
  # simulate nutrient input
  input_temp <- meta.arrR::simulate_nutr_input(n = n, max_i = max_i, frequency = freq_mn, 
                                               input_mn = nutrient_input, 
                                               amplitude_mn = amplitude_mean, amplitude_sd = amplitude_sd, 
                                               phase_mn = 0.0, phase_sd = phase_sd, 
                                               verbose = FALSE)

  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = list_parameters,
                                                nutrients_input = input_temp, movement = "behav",
                                                max_i = max_i, min_per_i = min_per_i,
                                                seagrass_each = seagrass_each,
                                                save_each = save_each, verbose = FALSE)
  
  # filter model run
  result_temp <- meta.arrR::filter_meta(x = result_temp, filter = c((max_i / years) * years_filter, max_i), 
                                        reset = TRUE, verbose = FALSE)
  
  # calc cv
  cv <- meta.arrR::calc_variability(x = result_temp, lag = c(FALSE, TRUE))
  
  # calculate biomass/production
  production <- meta.arrR::summarize_meta(result = result_temp, biomass = TRUE, production = TRUE,
                                          lag = c(FALSE, FALSE)) %>% 
    purrr::map(tidyr::pivot_longer, cols = -c(meta, timestep), names_to = "part") %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(timestep == max_i) %>%
    dplyr::group_by(part) %>% 
    dplyr::summarise(alpha = mean(value), gamma = sum(value)) %>% 
    tidyr::pivot_longer(-part, names_to = "measure")
    
  # combine to result data.frame and list
  list(cv = dplyr::mutate(dplyr::bind_rows(cv), 
                          amplitude_mean = amplitude_mean, amplitude_sd = amplitude_sd, 
                          phase_mean = 0.0, phase_sd = phase_sd, 
                          move_meta_sd = 0.0, 
                          stochastic = stochastic), 
       production = dplyr::mutate(production,
                                  amplitude_mean = amplitude_mean, amplitude_sd = amplitude_sd, 
                                  phase_mean = 0.0, phase_sd = phase_sd, 
                                  move_meta_sd = 0.0, 
                                  stochastic = stochastic))
  
}

#### Submit HPC ####

globals <- c("n", "max_i", "freq_mn", "nutrient_input", "metasyst_temp", "list_parameters", 
             "min_per_i", "seagrass_each", "save_each", "years", "years_filter") 

sbatch_cv <- rslurm::slurm_apply(f = foo_hpc, params = experiment_df, 
                                 global_objects = globals, jobname = "pp_abiotic_cv",
                                 nodes = nrow(experiment_df), cpus_per_node = 1, 
                                 slurm_options = list("account" = account, 
                                                      "partition" = "standard",
                                                      "time" = "01:00:00", ## hh:mm::ss
                                                      "mem-per-cpu" = "7G", 
                                                      "exclude" = exclude_nodes),
                                 pkgs = c("dplyr", "meta.arrR"),
                                 rscript_path = rscript_path, sh_template = sh_template, 
                                 submit = FALSE)

#### Collect results #### 

suppoRt::rslurm_missing(x = sbatch_cv)

cv_result <- rslurm::get_slurm_out(sbatch_cv, outtype = "raw")

suppoRt::save_rds(object = cv_result, filename = "01-abiotic-variability.rds",
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_cv)

#### Load results CV ####

cv_result <- readr::read_rds("02_Data/01-abiotic-variability.rds") %>% 
  purrr::map_dfr(function(i) i$cv) %>% 
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass", 
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")), 
                amplitude_mean = factor(amplitude_mean, labels = c("Low mean", "Medium mean", "High mean"), 
                                        ordered = TRUE),
                stochastic = factor(stochastic, levels = c("amplitude", "phase", "both"))) %>% 
  tibble::tibble()

#### Create ggplots CV ####

# create facet labels
label_part <- c(ag_production = "ag production", bg_production = "bg production",
                ttl_production = "ttl production")

label_stochastic <- c(amplitude = "Amplitude variability", phase = "Phase variability", 
                      both = "Simultaneously variability")

# create ggplot variability vs cv
gg_cv_abiotic <- dplyr::filter(cv_result, measure %in% c("alpha", "gamma"), 
                               part %in% c("ag_production", "bg_production", "ttl_production"), 
                               stochastic %in% c("amplitude", "phase")) %>%
  dplyr::mutate(single_sd = dplyr::case_when(stochastic == "amplitude" ~ amplitude_sd, 
                                             stochastic == "phase" ~ phase_sd)) %>% 
  ggplot(aes(x = single_sd, y = value, linetype = measure, color = amplitude_mean)) + 
  geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
  geom_point(alpha = 0.25, shape = 19, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5, formula = 'y ~ x') +
  facet_grid(rows = vars(part), cols = vars(stochastic), scales = "free", 
             labeller = labeller(part = label_part, stochastic = label_stochastic)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1.0, length.out = 5), limits = c(0, 1)) +
  # scale_y_continuous(breaks = seq(from = 0, to = 0.5, by = 0.1), limits = c(0, 0.5)) +
  scale_color_manual(name = "Amplitude base level", values = c("#007c2aff", "#ffcc00ff", "#2f62a7ff")) +
  scale_linetype_manual(name = "Scale", values = c(1, 2), labels = c(expression(paste(alpha, " scale")),
                                                                     expression(paste(gamma, " scale")))) +
  labs(x = "Abiotic variability", y = expression("CV"["primary production"])) +
  guides(linetype = guide_legend(order = 1, nrow = 2, keywidth = unit(10, "mm"),
                                 override.aes = list(color = "black"), title.position = "top"),
         color = guide_legend(order = 2, nrow = 3, title.position = "top")) +
  theme_bw() + theme(legend.position = "bottom", text = element_text(family = "Georgia"), 
                     panel.grid = element_blank())

# create ggplot alpha vs gamma
gg_pe_abiotic <- dplyr::filter(cv_result, measure %in% c("alpha", "gamma"), 
                               part %in% c("ag_production", "bg_production", "ttl_production"), 
                               stochastic %in% c("amplitude", "phase")) %>% 
  dplyr::mutate(single_sd = dplyr::case_when(stochastic == "amplitude" ~ amplitude_sd, 
                                             stochastic == "phase" ~ phase_sd)) %>% 
  tidyr::pivot_wider(names_from = measure, values_from = value) %>% 
  ggplot(aes(x = alpha, y = gamma, color = amplitude_mean, alpha = single_sd)) + 
  geom_abline(slope = 1, linetype = 2, color = "grey") +
  geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
  annotate(geom = "text", x = 0.25, y = 0.3, angle = 25, 
           color = "darkgrey", label = "paste(beta, '=1')", parse = TRUE) +
  geom_point(shape = 19) +
  facet_grid(rows = vars(part), cols = vars(stochastic), scales = "free",
             labeller = labeller(part = label_part, stochastic = label_stochastic)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 0.5, by = 0.1), limits = c(0, 0.5)) +
  # scale_y_continuous(breaks = seq(from = 0, to = 0.5, by = 0.1), limits = c(0, 0.5)) +
  scale_color_manual(name = "Amplitude base level", values = c("#007c2a", "#ffcc00", "#2f62a7")) + 
  labs(x = expression(paste(alpha, " CV"["primary production"])),
       y = expression(paste(gamma, " CV"["primary production"]))) +
  guides(alpha = "none", color = guide_legend(nrow = 1, title.position = "top")) +
  theme_bw() + theme(legend.position = "bottom", text = element_text(family = "Georgia"), 
                     panel.grid = element_blank())

#### Load results prod ####

prod_result <- readr::read_rds("02_Data/01-abiotic-variability.rds") %>% 
  purrr::map_dfr(function(i) i$production) %>% 
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass", 
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma")), 
                amplitude_mean = factor(amplitude_mean, labels = c("Low mean", "Medium mean", "High mean"), 
                                        ordered = TRUE),
                stochastic = factor(stochastic, levels = c("amplitude", "phase", "both"))) %>% 
  tibble::tibble()

#### Create ggplots prod ####

gg_prod_abiotic <- dplyr::filter(prod_result, measure == "gamma", 
                                 stochastic %in% c("amplitude", "phase")) %>% 
  dplyr::mutate(single_sd = dplyr::case_when(stochastic == "amplitude" ~ amplitude_sd, 
                                             stochastic == "phase" ~ phase_sd)) %>% 
  ggplot(aes(x = single_sd, y = value, color = amplitude_mean)) + 
  geom_point(shape = 19) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5, formula = 'y ~ x') +
  facet_grid(rows = vars(part), cols = vars(stochastic), scales = "free",
             labeller = labeller(part = label_part, stochastic = label_stochastic)) + 
  scale_color_manual(name = "Amplitude base level", values = c("#007c2a", "#ffcc00", "#2f62a7")) + 
  labs(x = "Abiotic variability", y = expression(paste(gamma, " primary production"))) +
  guides(alpha = "none", color = guide_legend(nrow = 1, title.position = "top")) +
  theme_bw() + theme(legend.position = "bottom", text = element_text(family = "Georgia"), 
                     panel.grid = element_blank())
