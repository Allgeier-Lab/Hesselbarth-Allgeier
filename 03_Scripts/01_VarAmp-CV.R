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

# default_parameters$seagrass_thres <- -0.5

#### Stable values ####

stable_values <- arrR::get_stable_values(bg_biomass = default_starting$bg_biomass,
                                         ag_biomass = default_starting$ag_biomass,
                                         parameters = default_parameters)

default_starting$nutrients_pool <- stable_values$nutrients_pool

default_starting$detritus_pool <- stable_values$detritus_pool

#### Simulate input ####

# number of iterations
itr <- 50

# set enrichment and amplitude levels
enrichment_levels <- c(low = 0.75, medium = 1, high = 1.25)

amplitude_levels <- c(low = 0.05, medium = 0.5, high = 1.0)

# sample amplitude variability for all treatments
variability_input <- tidyr::expand_grid(enrichment = enrichment_levels,
                                        amplitude = amplitude_levels) %>% 
  dplyr::slice(rep(x = 1:dplyr::n(), each = itr)) %>% 
  purrr::pmap(function(enrichment, amplitude) {
    
    result_temp <- get_modifier(n = n, local = amplitude, modifier = amplitude_levels,
                                method = "sample")
    
    purrr::map(result_temp, function(i) cbind(enrichment_lvl = enrichment, amplitude_lvl = amplitude, i))}) %>% 
  purrr::flatten()

#### Setup HPC function ####

globals <- list(n = n, max_i = max_i, default_starting = default_starting, 
                default_parameters = default_parameters, dimensions = dimensions, 
                grain = grain, input_mn = stable_values$nutrients_input, freq_mn = freq_mn,
                min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each) 

foo <- function(nutr_input) {
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i,
                                         starting_values = globals$default_starting,
                                         parameters = globals$default_parameters,
                                         dimensions = globals$dimensions, grain = globals$grain,
                                         reef = NULL, verbose = FALSE)
  
  # simulate input
  input_temp <- meta.arrR::sim_nutr_input(n = globals$n, max_i = globals$max_i,
                                          seagrass_each = globals$seagrass_each,
                                          amplitude_mod = nutr_input[, 4],
                                          phase_mod = nutr_input[, 5],
                                          input_mn = globals$input_mn * unique(nutr_input[, 1]), 
                                          freq_mn = globals$freq_mn)
  
  # plot(input_temp, gamma = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutrients_input = input_temp,
                                     parameters = globals$default_parameters,
                                     max_i = globals$max_i, min_per_i = globals$min_per_i,
                                     seagrass_each = globals$seagrass_each,
                                     save_each = globals$save_each, verbose = FALSE)
  
  # plot(result_temp, summarize = TRUE)
  
  # filter only second half of timesteps
  result_temp <- meta.arrR::filter_meta(x = result_temp, filter = c(globals$max_i / 2,
                                                                    globals$max_i), 
                                        reset = TRUE, verbose = FALSE)
  
  # calc CV
  cv_temp <- dplyr::bind_rows(meta.arrR::calc_variability(x = result_temp, lag = c(FALSE, TRUE), 
                                                          verbose = FALSE))
  
  # add type col
  cv_temp <- dplyr::bind_cols(type = "cv", cv_temp)
  
  # get the total sum of each timestep and metaecosystem
  result_sum <- dplyr::bind_rows(meta.arrR::summarize_meta(result = result_temp))
  
  # filter only final timestep
  result_sum <- dplyr::filter(result_sum, timestep == max(timestep))
  
  # calculate mean (alpha) and total (gamma) biomass values
  result_sum <- apply(X = result_sum[, -c(1,2)] , MARGIN = 2, 
                      FUN = function(x) c(alpha = mean(x, na.rm = TRUE), 
                                          gamma = sum(x, na.rm = TRUE)))
  
  # reshape to long format
  result_sum <- tidyr::pivot_longer(data.frame(type = "absolute", measure = row.names(result_sum), 
                                               result_sum), -c(type, measure), 
                                    names_to = "part")
  
  # rearrange cols
  result_sum <- dplyr::select(result_sum, type, part, measure, value)
  
  # combine to one df
  dplyr::bind_cols(enrichment_lvl = unique(nutr_input[, 1]),
                   amplitude_lvl = unique(nutr_input[, 2]), 
                   n_diff = unique(nutr_input[, 3]), 
                   dplyr::arrange(dplyr::bind_rows(cv_temp, result_sum), type, part, measure))
  
}

#### Submit to HPC #### 

sbatch_var_cv <- rslurm::slurm_map(f = foo, x = variability_input, 
                                   global_objects = "globals", jobname = "VarAmp_CV",
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

# suppoRt::rslurm_missing(x = sbatch_var_cv)

df_var_cv <- rslurm::get_slurm_out(sbatch_var_cv, outtype = "table")

suppoRt::save_rds(object = df_var_cv, filename = "01_VarAmp-CV.rds", 
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_var_cv)

#### Load data ####

df_var_cv <- readRDS("02_Data/01_VarAmp-CV.rds") %>% 
  dplyr::group_by(enrichment_lvl, amplitude_lvl, n_diff, type, part, measure) %>%
  dplyr::summarise(value_mn = mean(value), value_sd = sd(value), .groups = "drop") %>% 
  dplyr::mutate(enrichment_lvl = factor(enrichment_lvl, levels = c(0.75, 1.0, 1.25), labels = c("low", "medium", "high")), 
                amplitude_lvl = factor(amplitude_lvl, levels = c(0.05, 0.5, 1.0), labels = c("low", "medium", "high")),
                n_diff = as.integer(n_diff),
                type = factor(type, levels = c("cv", "absolute")),
                part = factor(part, levels = c("ag_biomass", "ag_production", "ttl_biomass",
                                               "bg_biomass", "bg_production", "ttl_production")),
                measure = factor(measure, levels = c("alpha", "beta", "gamma", "synchrony")))

#### Create ggplot ####

# create switch for biomass or production
switch <- "combined"

parts <- list(Aboveground = c("ag_biomass", "ag_production"), 
              Belowground = c("bg_biomass", "bg_production"), 
              Total = c("ttl_biomass", "ttl_production"))

col_palette <- c("#5ABCD6", "#FAD510", "#F22301")

legend_position <- c("none", "none", "bottom")

x_axis <- c(" ", " ", "Variability of amplitude")

x_labels <- list(element_blank(), element_blank(), NULL)

labels_facet <- list(c(low = "Low enrichment", medium = "Medium enrichment", high = "High enrichment"), 
                     c(low = "", medium = "", high = ""), 
                     c(low = "", medium = "", high = ""))

rel_heights <- c(1, 1, 1.25)

#### CV ####

y_digits <- function(x) sprintf("%.2f", x)

# y_axis <- c("", expression(paste("Coefficient of variation ", alpha)), "")
y_axis <- c(" ", "Coefficient of variation (alpha)", " ")

gg_alpha_cv <- purrr::map(seq_along(parts), function(i){
  
  dplyr::filter(df_var_cv, type == "cv", measure == "alpha", part %in% parts[[i]]) %>%
    ggplot() +
    geom_point(aes(x = n_diff, y = value_mn, col = amplitude_lvl, shape = part)) +
    geom_line(aes(x = n_diff, y = value_mn, col = amplitude_lvl, linetype = part)) +
    # geom_linerange(aes(x = n_diff, ymin = value_mn - value_sd,
    #                    ymax = value_mn + value_sd, col = amplitude_lvl, group = part)) +
    facet_wrap(. ~ enrichment_lvl, ncol = 3, nrow = 1, 
               labeller = labeller(enrichment_lvl = labels_facet[[i]])) +
    scale_x_continuous(breaks = seq(from = 0, to = 9, by = 1)) +
    scale_y_continuous(labels = y_digits) +
    scale_color_manual(name = "Amplitude treatment", values = col_palette) +
    scale_shape_manual(name = "Measure", values = c(1, 19), 
                       labels = c("biomass", "production")) +
    scale_linetype_manual(name = "Measure", values = c(2, 1), 
                          labels = c("biomass", "production")) +
    labs(y = y_axis[i], x = x_axis[i], subtitle = names(parts)[i]) +
    theme_classic(base_size = base_size) + 
    theme(legend.position = legend_position[i], strip.background = element_blank(), 
          strip.text = element_text(hjust = 0), 
          axis.text.x = x_labels[[i]])
  
})

gg_alpha_cv <- cowplot::plot_grid(plotlist = gg_alpha_cv, ncol = 1, nrow = 3, 
                                  rel_heights = rel_heights)

# y_axis <- c("", expression(paste("Coefficient of variation ", gamma)), "")
y_axis <- c(" ", "Coefficient of variation (gamma)", " ")

gg_gamma_cv <- purrr::map(seq_along(parts), function(i){
  
  dplyr::filter(df_var_cv, type == "cv", measure == "gamma", part %in% parts[[i]]) %>%
    ggplot() +
    geom_point(aes(x = n_diff, y = value_mn, col = amplitude_lvl, shape = part)) +
    geom_line(aes(x = n_diff, y = value_mn, col = amplitude_lvl, linetype = part)) +
    # geom_linerange(aes(x = n_diff, ymin = value_mn - value_sd,
    #                    ymax = value_mn + value_sd, col = amplitude_lvl, group = part)) +
    facet_wrap(. ~ enrichment_lvl, ncol = 3, nrow = 1, 
               labeller = labeller(enrichment_lvl = labels_facet[[i]])) +
    scale_x_continuous(breaks = seq(from = 0, to = 9, by = 1)) +
    scale_y_continuous(labels = y_digits) +
    scale_color_manual(name = "Amplitude treatment", values = col_palette) +
    scale_shape_manual(name = "Measure", values = c(1, 19), 
                       labels = c("biomass", "production")) +
    scale_linetype_manual(name = "Measure", values = c(2, 1), 
                          labels = c("biomass", "production")) +
    labs(y = y_axis[i], x = x_axis[i], subtitle = names(parts)[i]) +
    theme_classic(base_size = base_size) + 
    theme(legend.position = legend_position[i], strip.background = element_blank(), 
          strip.text = element_text(hjust = 0), 
          axis.text.x = x_labels[[i]])
  
})

gg_gamma_cv <- cowplot::plot_grid(plotlist = gg_gamma_cv, ncol = 1, nrow = 3, 
                                  rel_heights = rel_heights)

# y_axis <- c("", expression(paste("Coefficient of variation ", beta)), "")
y_axis <- c(" ", "Coefficient of variation (beta)", " ")

gg_beta_cv <- purrr::map(seq_along(parts), function(i){
  
  dplyr::filter(df_var_cv, type == "cv", measure == "beta", part %in% parts[[i]]) %>%
    ggplot() +
    geom_point(aes(x = n_diff, y = value_mn, col = amplitude_lvl, shape = part)) +
    geom_line(aes(x = n_diff, y = value_mn, col = amplitude_lvl, linetype = part)) +
    # geom_linerange(aes(x = n_diff, ymin = value_mn - value_sd,
    #                    ymax = value_mn + value_sd, col = amplitude_lvl, group = part)) +
    facet_wrap(. ~ enrichment_lvl, ncol = 3, nrow = 1, 
               labeller = labeller(enrichment_lvl = labels_facet[[i]])) +
    scale_x_continuous(breaks = seq(from = 0, to = 9, by = 1)) +
    scale_y_continuous(labels = y_digits) +
    scale_color_manual(name = "Amplitude treatment", values = col_palette) +
    scale_shape_manual(name = "Measure", values = c(1, 19), 
                       labels = c("biomass", "production")) +
    scale_linetype_manual(name = "Measure", values = c(2, 1), 
                          labels = c("biomass", "production")) +
    labs(y = y_axis[i], x = x_axis[i], subtitle = names(parts)[i]) +
    theme_classic(base_size = base_size) + 
    theme(legend.position = legend_position[i], strip.background = element_blank(), 
          strip.text = element_text(hjust = 0), 
          axis.text.x = x_labels[[i]])
  
})

gg_beta_cv <- cowplot::plot_grid(plotlist = gg_beta_cv, ncol = 1, nrow = 3, 
                                 rel_heights = rel_heights)

#### Absolute values ####

y_digits <- function(x) sprintf("%.0f", x)

# y_axis <- c("", expression(paste("Absolute values ", alpha)), "")
y_axis <- c(" ", "Absolute values (alpha)", " ")

gg_alpha_abs <- purrr::map(seq_along(parts), function(i){
  
  dplyr::filter(df_var_cv, type == "absolute", measure == "alpha", part %in% parts[[i]]) %>%
    ggplot() +
    geom_point(aes(x = n_diff, y = value_mn, col = amplitude_lvl, shape = part)) +
    geom_line(aes(x = n_diff, y = value_mn, col = amplitude_lvl, linetype = part)) +
    # geom_linerange(aes(x = n_diff, ymin = value_mn - value_sd,
    #                    ymax = value_mn + value_sd, col = amplitude_lvl, group = part)) +
    facet_wrap(. ~ enrichment_lvl, ncol = 3, nrow = 1, 
               labeller = labeller(enrichment_lvl = labels_facet[[i]])) +
    scale_x_continuous(breaks = seq(from = 0, to = 9, by = 1)) +
    scale_y_continuous(labels = y_digits) +
    scale_color_manual(name = "Amplitude treatment", values = col_palette) +
    scale_shape_manual(name = "Measure", values = c(1, 19), 
                       labels = c("biomass", "production")) +
    scale_linetype_manual(name = "Measure", values = c(2, 1), 
                          labels = c("biomass", "production")) +
    labs(y = y_axis[i], x = x_axis[i], subtitle = names(parts)[i]) +
    theme_classic(base_size = base_size) + 
    theme(legend.position = legend_position[i], strip.background = element_blank(), 
          strip.text = element_text(hjust = 0), 
          axis.text.x = x_labels[[i]])
  
})

gg_alpha_abs <- cowplot::plot_grid(plotlist = gg_alpha_abs, ncol = 1, nrow = 3, 
                                   rel_heights = rel_heights)

# y_axis <- c("", expression(paste("Absolute values ", gamma)), "")
y_axis <- c(" ", "Absolute values (gamma)", " ")

gg_gamma_abs <- purrr::map(seq_along(parts), function(i){
  
  dplyr::filter(df_var_cv, type == "absolute", measure == "gamma", part %in% parts[[i]]) %>%
    ggplot() +
    geom_point(aes(x = n_diff, y = value_mn, col = amplitude_lvl, shape = part)) +
    geom_line(aes(x = n_diff, y = value_mn, col = amplitude_lvl, linetype = part)) +
    # geom_linerange(aes(x = n_diff, ymin = value_mn - value_sd,
    #                    ymax = value_mn + value_sd, col = amplitude_lvl, group = part)) +
    facet_wrap(. ~ enrichment_lvl, ncol = 3, nrow = 1, 
               labeller = labeller(enrichment_lvl = labels_facet[[i]])) +
    scale_x_continuous(breaks = seq(from = 0, to = 9, by = 1)) +
    scale_y_continuous(labels = y_digits) +
    scale_color_manual(name = "Amplitude treatment", values = col_palette) +
    scale_shape_manual(name = "Measure", values = c(1, 19), 
                       labels = c("biomass", "production")) +
    scale_linetype_manual(name = "Measure", values = c(2, 1), 
                          labels = c("biomass", "production")) +
    labs(y = y_axis[i], x = x_axis[i], subtitle = names(parts)[i]) +
    theme_classic(base_size = base_size) + 
    theme(legend.position = legend_position[i], strip.background = element_blank(), 
          strip.text = element_text(hjust = 0), 
          axis.text.x = x_labels[[i]])
  
})

gg_gamma_abs <- cowplot::plot_grid(plotlist = gg_gamma_abs, ncol = 1, nrow = 3, 
                                   rel_heights = rel_heights)

#### Save ggplot ####

#### CV ####
suppoRt::save_ggplot(plot = gg_alpha_cv, filename = paste0("01_VarAmp_CV_alpha_", switch, ".png"),
                     path = "04_Figures", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_gamma_cv, filename = paste0("01_VarAmp_CV_gamma_", switch, ".png"),
                     path = "04_Figures", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_beta_cv, filename = paste0("01_VarAmp_CV_beta_", switch, ".png"),
                     path = "04_Figures", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

#### Abs #### 
suppoRt::save_ggplot(plot = gg_alpha_abs, filename = paste0("01_VarAmp_abs_alpha_", switch, ".png"),
                     path = "04_Figures", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_gamma_abs, filename = paste0("01_VarAmp_abs_gamma_", switch, ".png"),
                     path = "04_Figures", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
