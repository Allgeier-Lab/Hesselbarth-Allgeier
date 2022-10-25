##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Relationship between biotic variability and mean cross-meta-ecosystem movement

#### Load setup ####

source("05_Various/setup.R")

n <- 5

#### Adapt parameters ####

parameters_list$nutrients_loss <- 0.0

#### Stable values #### 

stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
                                       ag_biomass = starting_values_list$ag_biomass,
                                       parameters = parameters_list)

starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool

starting_values_list$detritus_pool <- stable_values_list$detritus_pool

#### Setup experiment ####

reps <- 50

matrix_lhs <- lhs::improvedLHS(n = reps, k = 1, dup = 2)

matrix_lhs[, 1] <- qunif(matrix_lhs[, 1], 0.1, 1.0) 

experiment_df <- tibble::tibble(move_meta_sd = matrix_lhs[, 1]) %>% 
  dplyr::slice(rep(x = 1:dplyr::n(), times = 10))

# setup metaecosystems
metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                       starting_values = starting_values_list, parameters = parameters_list,
                                       dimensions = dimensions, grain = grain, use_log = use_log, 
                                       verbose = FALSE)

foo_hpc <- function(move_meta_sd) {
  
  parameters_list$move_meta_sd <- move_meta_sd
  
  # create new attributed matrix
  attr_replace <- meta.arrR:::setup_attributes(fishpop = metasyst_temp$fishpop, 
                                               parameters = parameters_list, max_i = max_i)
  
  # replace matrix
  metasyst_temp$fishpop_attr[, 3] <- attr_replace[, 3]
  
  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                                nutrients_input = 0.0, movement = "behav",
                                                max_i = max_i, min_per_i = min_per_i,
                                                seagrass_each = seagrass_each,
                                                save_each = save_each, verbose = FALSE) %>% 
    meta.arrR::filter_meta(filter = c((max_i / years) * years_filter, max_i), 
                           reset = TRUE, verbose = FALSE)
  
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

globals <- c("metasyst_temp", "parameters_list", "max_i", # setup_attributes
             "min_per_i", "seagrass_each", "save_each", # run_simulation_meta
             "years", "years_filter") # filter_meta

sbatch_move_meta <- rslurm::slurm_apply(f = foo_hpc, params = experiment_df, 
                                        global_objects = globals, jobname = "move_meta",
                                        nodes = nrow(experiment_df), cpus_per_node = 1, 
                                        slurm_options = list("account" = account, 
                                                             "partition" = "standard",
                                                             "time" = "01:00:00", ## hh:mm::ss
                                                             "mem-per-cpu" = "7G"),
                                        pkgs = c("dplyr", "meta.arrR"),
                                        rscript_path = rscript_path, submit = FALSE)

#### Collect results #### 

suppoRt::rslurm_missing(x = sbatch_move_meta)

move_meta_result <- rslurm::get_slurm_out(sbatch_move_meta, outtype = "table")

suppoRt::save_rds(object = move_meta_result, filename = "biotic-vs-connectivity.rds",
                  path = "02_Data/", overwrite = FALSE)

rslurm::cleanup_files(sbatch_move_meta)

#### Load results ####

move_meta_result <- readr::read_rds("02_Data/biotic-vs-connectivity.rds")

#### Summarize results ####

move_meta_result_sum <- dplyr::mutate(move_meta_result) %>% 
  dplyr::group_by(move_meta_sd) %>% 
  dplyr::summarise(mean = mean(moved), sd = sd(moved), 
                   quant_lo = quantile(moved, probs = 0.35), 
                   quant_hi = quantile(moved, probs = 0.65), .groups = "drop") %>% 
  dplyr::mutate(lo = dplyr::case_when(mean - sd > 0 ~ mean - sd,
                                      mean - sd < 0 ~ 0), hi = mean + sd)

#### Setup ggplot ####

breaks_x <- seq(from = 0, to = max_i / filter_factor, by = max_i / years)
# labels_x <- paste(seq(from = 0, to = years / filter_factor, by = 1), "years")

line_size <- 0.25

point_size <- 0.5

base_size <- 7.5

#### Create figure moved ####

# gg_probs_moved <- ggplot(data = move_meta_result_sum, aes(x = move_meta_sd, y = mean)) +
#   geom_hline(yintercept = 0, linetype = 2, color = "grey") + 
#   # geom_ribbon(aes(ymin = quant_lo, ymax = quant_hi), alpha = 0.2, fill = "#fff178") +
#   geom_smooth(se = FALSE, linetype = 1, size = 0.75, color = "#32b2da") +
#   geom_errorbar(aes(ymin = lo, ymax = hi)) +
#   geom_point(shape = 19) +
#   scale_x_continuous(limits = c(0.1, 1), breaks = seq(from = 0.1, to = 1.0, by = 0.1)) +
#   scale_color_manual(name = "Parameter", values = c("#d0413d", "#0c214e")) + 
#   labs(x = "Biotic variability", y = "Mean cross-system movement (connectivity)") +
#   theme_classic(base_size = 10.0)

gg_probs_moved <- ggplot(data = move_meta_result_sum, aes(x = move_meta_sd, y = mean)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), size = 0.25, color = "#1e466e", alpha = 0.75) +
  geom_smooth(se = FALSE, linetype = 1, size = line_size, color = "#1e466e") + # #32b2da
  geom_point(shape = 19, size = point_size, color = "#1e466e") +
  scale_x_continuous(limits = c(0.1, 1), breaks = seq(from = 0.1, to = 1.0, by = 0.1)) +
  scale_color_manual(name = "Parameter", values = c("#d0413d", "#0c214e")) + 
  labs(x = "Diversity consumer behavior", y = "Connectivity") +
  theme_classic(base_size = base_size) +
  theme(axis.text = element_blank())

#### Save ggplot ####

# suppoRt::save_ggplot(plot = gg_probs_moved, filename = paste0("Figure-A1", extension),
#                      path = "04_Figures/Appendix", width = width, height = height * 1/3,
#                      units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_probs_moved, filename = "Figure-1-connect.png",
                     path = "04_Figures/", width = 85 / 2, height = 35,
                     units = units, dpi = dpi, overwrite = FALSE)
