##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Relationship between movement parameters and mean cross-meta-ecosystem movement

#### Load setup ####

source("05_Various/setup.R")

#### Adapt parameters ####

# max_i <- max_i / 5

#### Stable values #### 

stable_values <- arrR::get_req_nutrients(bg_biomass = list_starting$bg_biomass,
                                         ag_biomass = list_starting$ag_biomass,
                                         parameters = list_parameters)

list_starting$nutrients_pool <- stable_values$nutrients_pool

list_starting$detritus_pool <- stable_values$detritus_pool

#### Setup experiment ####

iterations <- 50

# create data.frame with all combinations
experiment_df <- tibble::tibble(
  move_meta_mean = c(c(0, runif(n = iterations - 2, min = 0.0, max = 1.0), 1), 
                     rep(x = 0, times = iterations)),
  move_meta_sd = c(rep(x = 0, times = iterations), 
                   c(0, runif(n = iterations - 2, min = 0.0, max = 1.0), 1)), 
  rand = rep(x = c("mean", "sd"), each = iterations))

experiment_df <- tibble::tibble(
  move_meta_mean = c(c(0, runif(n = iterations - 2, min = 0.0, max = 1.0), 1), 
                     rep(x = 0, times = iterations)),
  move_meta_sd = c(rep(x = 0, times = iterations), 
                   c(0, runif(n = iterations - 2, min = 0.0, max = 1.0), 1)), 
  rand = rep(x = c("mean", "sd"), each = iterations))
  
# setup metaecosystems
metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = matrix_reef,
                                       starting_values = list_starting, parameters = list_parameters,
                                       dimensions = dimensions, grain = grain, use_log = use_log, 
                                       verbose = FALSE)

foo_hpc <- function(move_meta_mean, move_meta_sd, rand) {
  
  library(dplyr)
  
  # update move meta sd
  list_parameters$move_meta_mean <- move_meta_mean
  
  list_parameters$move_meta_sd <- move_meta_sd
  
  # create new attributed matrix
  attr_replace <- meta.arrR:::setup_attributes(fishpop = metasyst_temp$fishpop, 
                                               parameters = list_parameters, max_i = max_i)
  
  # replace matrix
  metasyst_temp$fishpop_attr[, 3] <- attr_replace[, 3]
  
  # simulate nutrient input h
  input_temp <-  meta.arrR::simulate_nutr_input(n = n, max_i = max_i, frequency = years, 
                                                input_mn = nutrient_input, verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = list_parameters,
                                                nutrients_input = input_temp, movement = "behav",
                                                max_i = max_i, min_per_i = min_per_i,
                                                seagrass_each = seagrass_each,
                                                save_each = save_each, verbose = FALSE)
  
  # get moved counts
  dplyr::bind_rows(result_temp$fishpop) %>% 
    dplyr::filter(timestep == max_i) %>% 
    dplyr::left_join(y = as.data.frame(metasyst_temp$fishpop_attr), by = "id") %>% 
    dplyr::select(id, length, weight, moved, move_prob) %>% 
    dplyr::left_join(y = dplyr::select(dplyr::bind_rows(metasyst_temp$fishpop), id, length, weight), 
                     by = "id", suffix = c(".end", ".start")) %>% 
    dplyr::mutate(rand = rand, move_meta_mean = move_meta_mean, move_meta_sd = move_meta_sd,
                  .before = "id")
  
  # range <- meta.arrR::get_abundance(result_temp) %>% 
  #   dplyr::group_by(timestep) %>% 
  #   dplyr::summarise(lo = quantile(abundance, probs = 0.25), hi = quantile(abundance, probs = 0.75),
  #                    .groups = "drop") %>% 
  #   dplyr::filter(timestep != 0) %>%  
  #   tidyr::pivot_longer(-timestep) %>% 
  #   dplyr::group_by(name) %>% 
  #   dplyr::summarise(mean = mean(value), sd = sd(value)) %>% 
  #   dplyr::mutate(name = factor(name, levels = c("lo", "hi")), 
  #                 move_meta_sd = move_sd, .before = "name")
  # 
  # return(list(moved = moved, range = range))
  
}

#### Submit HPC

globals <- c("list_parameters", "metasyst_temp", "max_i", "n", "years", "nutrient_input", 
             "min_per_i", "seagrass_each", "save_each")

sbatch_move_meta <- rslurm::slurm_apply(f = foo_hpc, params = experiment_df, 
                                        global_objects = globals, jobname = "move_meta_param",
                                        nodes = nrow(experiment_df), cpus_per_node = 1, 
                                        slurm_options = list("account" = account, 
                                                             "partition" = "standard",
                                                             "time" = "01:00:00", ## hh:mm::ss
                                                             "mem-per-cpu" = "7G"),
                                        pkgs = c("dplyr", "meta.arrR"),
                                        rscript_path = rscript_path, sh_template = sh_template, 
                                        submit = FALSE)

#### Collect results #### 

suppoRt::rslurm_missing(x = sbatch_move_meta)

move_meta_result <- rslurm::get_slurm_out(sbatch_move_meta, outtype = "table")

suppoRt::save_rds(object = move_meta_result, filename = "04-move-meta-param.rds",
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_move_meta)

#### Load results ####

move_meta_result <- readr::read_rds("02_Data/04-move-meta-param.rds")

#### Create figure moved ####

move_meta_result_sum <- dplyr::mutate(move_meta_result, 
                                      param_rand = dplyr::case_when(rand == "mean" ~ move_meta_mean, 
                                                                    rand == "sd" ~ move_meta_sd)) %>% 
  dplyr::group_by(rand, param_rand) %>% 
  dplyr::summarise(mean = mean(moved), sd = sd(moved)) %>% 
  dplyr::mutate(lo = dplyr::case_when(mean - sd > 0 ~ mean - sd,
                                      mean - sd < 0 ~ 0),
                hi = mean + sd)
  
gg_probs_moved <- ggplot(data = move_meta_result_sum, 
                         aes(x = param_rand, y = mean, color = rand)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") + 
  geom_errorbar(aes(ymin = lo, ymax = hi), alpha = 0.35) +
  geom_point(shape = 19) +
  geom_smooth(se = FALSE, linetype = 2, size = 0.5) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(name = "Parameter", values = c("#d0413d", "#0c214e")) + 
  labs(x = expression(paste("Parameter ", italic(move_meta_mean), " / ", italic(move_meta_sd))),
                 y = "Mean cross-system movement") +
  theme_classic() 

suppoRt::save_ggplot(plot = gg_probs_moved, filename = "04-param-connect.pdf",
                     path = "04_Figures/", width = width, height = height / 2,
                     units = units, dpi = dpi, overwrite = overwrite)
