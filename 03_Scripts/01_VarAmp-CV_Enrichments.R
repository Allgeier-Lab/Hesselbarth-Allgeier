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

default_starting$pop_n <- 0

default_parameters$nutrients_diffusion <- 0.0
default_parameters$detritus_diffusion <- 0.0
default_parameters$detritus_fish_diffusion <- 0.0

# default_parameters$seagrass_thres <- 1/3

#### Stable values ####

stable_values <- arrR::get_stable_values(starting_values = default_starting,
                                         parameters = default_parameters)

default_starting$nutrients_pool <- stable_values$nutrients_pool

default_starting$detritus_pool <- stable_values$detritus_pool

#### Simulate input ####

itr <- 50

amplitude_levels <- rep(c(low = 0.05, medium = 0.5, high = 1), each = itr)

# simulate all treatment levels
variability_input <- purrr::map(1:length(amplitude_levels), function(i) {
  
  message("\r> Amplitude level: ", i, "/", length(amplitude_levels), "\t\t",
          appendLF = FALSE)
  
  get_modifier(n = n, local = amplitude_levels[i], modifier = amplitude_levels,
               method = "sample")}) %>%
  purrr::flatten() %>%
  purrr::set_names(paste(rep(x = names(amplitude_levels), each = n + 1),
                         seq(from = 0, to = 9, by = 1), sep = "_"))

# setupt enrichment levels
enrichment_levels <- c(low = 0.75, medium = 1.0, high = 1.25)

# create name vector....
names_var_input <- rep(x = names(amplitude_levels), each = n + 1) %>% 
  paste(seq(from = 0, to = 9, by = 1), sep = "_") %>% 
  rep(each = length(enrichment_levels)) %>% 
  paste(enrichment_levels, sep = "_")

# get only 4 variability levels, repeat each n times and add enrichment col
variability_input <- names(variability_input) %>% 
  stringr::str_which(pattern = "_0|_3|_6|_9") %>%
  magrittr::extract(variability_input, .) %>% 
  purrr::map(function(i) {
    purrr::map(1:length(enrichment_levels), function(j) {
      dplyr::bind_cols(i, enrich = enrichment_levels[j])})}) %>% 
  purrr::flatten() %>% 
  purrr::set_names(names_var_input)

#### Setup HPC function ####

globals <- list(n = n, max_i = max_i, default_starting = default_starting, 
                default_parameters = default_parameters, dimensions = dimensions, 
                grain = grain, input_mn = stable_values$nutr_input, freq_mn = freq_mn,
                min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each) 

foo <- function(nutr_input) {
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i,
                                         starting_values = globals$default_starting,
                                         parameters = globals$default_parameters,
                                         dimensions = globals$dimensions, grain = globals$grain,
                                         reefs = NULL, verbose = FALSE)
  
  # simulate input
  input_temp <- meta.arrR::sim_nutr_input(n = globals$n, max_i = globals$max_i,
                                          amplitude_mod = nutr_input[, 2],
                                          phase_mod = nutr_input[, 3],
                                          input_mn = globals$input_mn * unique(nutr_input[, 4]), 
                                          freq_mn = globals$freq_mn)
  
  # run model
  result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutr_input = input_temp,
                                     parameters = globals$default_parameters,
                                     max_i = globals$max_i, min_per_i = globals$min_per_i,
                                     seagrass_each = globals$seagrass_each,
                                     save_each = globals$save_each, verbose = FALSE)
  
  # # filter only second half of timesteps
  # result_temp <- meta.arrR::filter_meta(x = result_temp, filter = c(globals$max_i / 2, 
  #                                                                   globals$max_i))
  
  # calc CV
  cv_temp <- dplyr::bind_rows(meta.arrR::calc_variability(x = result_temp, lag = c(FALSE, TRUE), 
                                                          verbose = FALSE))
  
  # add type col
  cv_temp <- dplyr::bind_cols(type = "cv", cv_temp)
  
  # get the total sum of each timestep and metaecosystem
  sum_temp <- dplyr::bind_rows(meta.arrR::summarize_meta(result = result_temp))
  
  # filter only final timestep
  sum_temp <- dplyr::filter(sum_temp, timestep == max(timestep))
  
  # calculate mean (alpha) and total (gamma) biomass values
  sum_temp <- apply(X = sum_temp[, -c(1,2)] , MARGIN = 2, 
                    FUN = function(x) c(alpha = mean(x, na.rm = TRUE), 
                                        gamma = sum(x, na.rm = TRUE)))
  
  # reshape to long format
  sum_temp <- tidyr::pivot_longer(data.frame(type = "absolute", measure = row.names(sum_temp), 
                                             sum_temp), -c(type, measure), 
                                  names_to = "part")
  
  # rearrange cols
  sum_temp <- dplyr::select(sum_temp, type, part, measure, value)
  
  # combine to one df
  dplyr::arrange(dplyr::bind_rows(cv_temp, sum_temp), type, part, measure)
  
}

#### Submit to HPC #### 

sbatch_var_cv <- rslurm::slurm_map(f = foo, x = variability_input, 
                                   global_objects = "globals", jobname = "VarAmp_CV_Enrich",
                                   nodes = length(variability_input), cpus_per_node = 1, 
                                   slurm_options = list("account" = "jeallg1", 
                                                        "partition" = "standard",
                                                        "time" = "02:00:00", ## hh:mm::ss
                                                        "mem-per-cpu" = "7G", 
                                                        "exclude" = exclude_nodes),
                                   pkgs = c("dplyr", "meta.arrR", "tidyr"),
                                   rscript_path = rscript_path, sh_template = sh_template, 
                                   submit = FALSE)

#### Collect results ####

df_var_cv <- rslurm::get_slurm_out(sbatch_var_cv, outtype = "raw") %>% 
  dplyr::bind_rows(.id = "id")

suppoRt::save_rds(object = var_prod_result, filename = "01_VarAmp-CV_Enrich.rds", 
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_var_cv)

#### Load data ####

df_var_cv <- readRDS("02_Data/01_VarAmp-CV_Enrich.rds") %>% 
  tidyr::separate(col = id, sep = "_",  into = c("amplitude", "variability", "enrichment")) %>% 
  dplyr::group_by(amplitude, variability, enrichment, type, part, measure) %>%
  dplyr::summarise(value_mn = mean(value), value_sd = sd(value), .groups = "drop") %>%
  dplyr::mutate(amplitude = factor(amplitude, levels = c("low", "medium", "high"), 
                                   labels = c("Low amplitude", "Medium amplitude", "High amplitude")),
                variability = factor(variability, ordered = TRUE),
                enrichment = factor(enrichment, levels = c(0.75, 1, 1.25), 
                                    labels = c("Low enrichment", "Medium enrichment", "High enrichment")),
                type = factor(type, levels = c("cv", "absolute")),
                part = factor(part, levels = c("ag_biomass", "ag_production",
                                               "bg_biomass", "bg_production")),
                measure = factor(measure, levels = c("alpha", "beta", "gamma", "synchrony")))

#### Create CV ggplot ####

col_palette <- c("#586F7E", "#168B98", "#ED5B66", "#DAAD4F")

part_select <- c("bg_biomass", "bg_production", "ag_biomass", "ag_production")

gg_alpha <- dplyr::filter(df_var_cv, type == "cv", measure == "alpha",
                          part %in% part_select) %>%
  ggplot() +
  geom_line(aes(x = variability, y = value_mn, col = part, group = part), alpha = 0.25) +
  geom_point(aes(x = variability, y = value_mn, col = part)) +
  geom_linerange(aes(x = variability, ymin = value_mn - value_sd, ymax = value_mn + value_sd, 
                     col = part)) +
  # geom_boxplot(aes(x = variability, y = value, fill = part)) +
  facet_wrap(. ~ enrichment + amplitude, scales = "fixed", nrow = 3) + 
  # scale_fill_manual(name = "", values = col_palette) +
  scale_color_manual(name = "", values = col_palette) +
  labs(y = expression(paste("Coefficient of variation ", alpha)), 
       x = "Variability of amplitude") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

gg_gamma <- dplyr::filter(df_var_cv, type == "cv", measure == "gamma",
                          part %in% part_select) %>%
  ggplot() +
  geom_line(aes(x = variability, y = value_mn, col = part, group = part), alpha = 0.25) +
  geom_point(aes(x = variability, y = value_mn, col = part)) +
  geom_linerange(aes(x = variability, ymin = value_mn - value_sd, ymax = value_mn + value_sd, 
                     col = part)) +
  # geom_boxplot(aes(x = variability, y = value, fill = part)) +
  facet_wrap(. ~ enrichment + amplitude, scales = "fixed", nrow = 3) + 
  # scale_fill_manual(name = "", values = col_palette) +
  scale_color_manual(name = "", values = col_palette) +
  labs(y = expression(paste("Coefficient of variation ", gamma)), 
       x = "Variability of amplitude") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

gg_beta <- dplyr::filter(df_var_cv, type == "cv", measure == "beta",
                         part %in% part_select) %>%
  ggplot() +
  geom_line(aes(x = variability, y = value_mn, col = part, group = part), alpha = 0.25) +
  geom_point(aes(x = variability, y = value_mn, col = part)) +
  geom_linerange(aes(x = variability, ymin = value_mn - value_sd, ymax = value_mn + value_sd, 
                     col = part)) +
  # geom_boxplot(aes(x = variability, y = value, fill = part)) +
  facet_wrap(. ~ enrichment + amplitude, scales = "fixed", nrow = 3) + 
  # scale_fill_manual(name = "", values = col_palette) +
  scale_color_manual(name = "", values = col_palette) +
  labs(y = expression(paste("Portfolio effect ", beta)), 
       x = "Variability of amplitude") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")
  
#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_alpha, filename = "01_VarAmp-CV_alpha_enrich.png",
                     path = "04_Figures", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_gamma, filename = "01_VarAmp-CV_gamma_enrich.png",
                     path = "04_Figures", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_beta, filename = "01_VarAmp-CV_beta_enrich.png",
                     path = "04_Figures", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
