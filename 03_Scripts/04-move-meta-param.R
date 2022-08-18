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

#### Stable values #### 

list_stable <- arrR::get_req_nutrients(bg_biomass = list_starting$bg_biomass,
                                       ag_biomass = list_starting$ag_biomass,
                                       parameters = list_parameters)

list_starting$nutrients_pool <- list_stable$nutrients_pool

list_starting$detritus_pool <- list_stable$detritus_pool

#### Setup experiment ####

reps <- 50

matrix_lhs <- lhs::improvedLHS(n = reps, k = 1, dup = 2)

matrix_lhs[, 1] <- qunif(matrix_lhs[, 1], 0.1, 1.0) 

df_experiment <- tibble::tibble(move_meta_sd = matrix_lhs[, 1]) %>% 
  dplyr::slice(rep(x = 1:dplyr::n(), times = 10))

# setup metaecosystems
metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = matrix_reef,
                                       starting_values = list_starting, parameters = list_parameters,
                                       dimensions = dimensions, grain = grain, use_log = use_log, 
                                       verbose = FALSE)

foo_hpc <- function(move_meta_sd) {
  
  list_parameters$move_meta_sd <- move_meta_sd
  
  # create new attributed matrix
  attr_replace <- meta.arrR:::setup_attributes(fishpop = metasyst_temp$fishpop, 
                                               parameters = list_parameters, max_i = max_i)
  
  # replace matrix
  metasyst_temp$fishpop_attr[, 3] <- attr_replace[, 3]
  
  # simulate nutrient input h
  input_temp <- meta.arrR::simulate_nutrient_sine(n = n, max_i = max_i, frequency = years, 
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
    dplyr::mutate(move_meta_sd = move_meta_sd,.before = "id")
  
}

#### Submit HPC ####

globals <- c("list_parameters", "metasyst_temp", "max_i", "n", "years", "nutrient_input", 
             "min_per_i", "seagrass_each", "save_each")

sbatch_move_meta <- rslurm::slurm_apply(f = foo_hpc, params = df_experiment, 
                                        global_objects = globals, jobname = "move_meta_param",
                                        nodes = nrow(df_experiment), cpus_per_node = 1, 
                                        slurm_options = list("account" = account, 
                                                             "partition" = "standard",
                                                             "time" = "01:00:00", ## hh:mm::ss
                                                             "mem-per-cpu" = "7G"),
                                        pkgs = c("dplyr", "meta.arrR"),
                                        rscript_path = rscript_path, submit = FALSE)

#### Collect results #### 

suppoRt::rslurm_missing(x = sbatch_move_meta)

move_meta_result <- rslurm::get_slurm_out(sbatch_move_meta, outtype = "table")

suppoRt::save_rds(object = move_meta_result, filename = "04-move-meta-param.rds",
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_move_meta)

#### Load results ####

move_meta_result <- readr::read_rds("02_Data/04-move-meta-param.rds")

#### Summarize results ####

move_meta_result_sum <- dplyr::mutate(move_meta_result) %>% 
  dplyr::group_by(move_meta_sd) %>% 
  dplyr::summarise(mean = mean(moved), sd = sd(moved), 
                   quant_lo = quantile(moved, probs = 0.35), 
                   quant_hi = quantile(moved, probs = 0.65), .groups = "drop") %>% 
  dplyr::mutate(lo = dplyr::case_when(mean - sd > 0 ~ mean - sd,
                                      mean - sd < 0 ~ 0), hi = mean + sd)

#### Create figure moved ####

gg_probs_moved <- ggplot(data = move_meta_result_sum, aes(x = move_meta_sd, y = mean)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") + 
  # geom_ribbon(aes(ymin = quant_lo, ymax = quant_hi), alpha = 0.2, fill = "#fff178") +
  geom_smooth(se = FALSE, linetype = 1, size = 0.75, color = "#32b2da") +
  geom_errorbar(aes(ymin = lo, ymax = hi)) +
  geom_point(shape = 19) +
  scale_x_continuous(limits = c(0.1, 1), breaks = seq(from = 0.1, to = 1.0, by = 0.1)) +
  scale_color_manual(name = "Parameter", values = c("#d0413d", "#0c214e")) + 
  labs(x = "Biotic connectivity variability", y = "Mean cross-system movement") +
  theme_classic(base_size = 10.0) 

suppoRt::save_ggplot(plot = gg_probs_moved, filename = paste0("04-param-connect", extension),
                     path = "04_Figures/", width = width, height = height * 1/3,
                     units = units, dpi = dpi, overwrite = FALSE)
