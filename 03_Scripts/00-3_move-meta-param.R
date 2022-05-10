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

# nothing to change

#### Stable values #### 

stable_values <- arrR::get_req_nutrients(bg_biomass = starting_list$bg_biomass,
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

starting_list$nutrients_pool <- stable_values$nutrients_pool

starting_list$detritus_pool <- stable_values$detritus_pool

#### Setup experiment ####

iterations <- 50

# create data.frame with all combinations
experiment_df <- tibble::tibble(move_sd = c(0, runif(n = iterations - 2, min = 0.0, max = 1.0), 1))

# setup metaecosystems
metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                       starting_values = starting_list, parameters = parameters_list,
                                       dimensions = dimensions, grain = grain, use_log = use_log, 
                                       verbose = FALSE)

foo_hpc <- function(move_sd) {
  
  library(dplyr)
  
  # update move meta sd
  parameters_list$move_meta_sd <- move_sd
  
  # create new attributed matrix
  attr_replace <- meta.arrR:::setup_attributes(fishpop = metasyst_temp$fishpop, parameters = parameters_list, 
                                               max_i = max_i)
  
  # replace matrix
  metasyst_temp$fishpop_attr <- attr_replace
  
  # simulate nutrient input h
  input_temp <-  meta.arrR::simulate_nutr_input(n = n, max_i = max_i, frequency = freq_mn, 
                                                input_mn = nutrient_input, verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
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
    dplyr::mutate(move_meta_sd = move_sd, .before = "id")
  
}

#### Submit HPC

globals <- c("parameters_list", "metasyst_temp", "max_i", "n", "freq_mn", "nutrient_input", 
             "min_per_i", "seagrass_each", "save_each")

sbatch_move_meta <- rslurm::slurm_apply(f = foo_hpc, params = experiment_df, 
                                        global_objects = globals, jobname = "move_meta_sd",
                                        nodes = nrow(experiment_df), cpus_per_node = 1, 
                                        slurm_options = list("account" = account, 
                                                             "partition" = "standard",
                                                             "time" = "01:30:00", ## hh:mm::ss
                                                             "mem-per-cpu" = "7G"),
                                        pkgs = c("dplyr", "tidyr", "meta.arrR"),
                                        rscript_path = rscript_path, sh_template = sh_template, 
                                        submit = FALSE)

#### Collect results #### 

suppoRt::rslurm_missing(x = sbatch_move_meta)

move_meta_result <- rslurm::get_slurm_out(sbatch_move_meta, outtype = "table")

suppoRt::save_rds(object = move_meta_result, filename = "00-move_meta.rds",
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_move_meta)

#### Load results ####

move_meta_result <- readr::read_rds("02_Data/00-move_meta.rds")

move_meta_result_sum <- dplyr::group_by(move_meta_result, move_meta_sd) %>% 
  dplyr::summarise(mean = mean(moved), sd = sd(moved)) %>% 
  dplyr::mutate(lo = dplyr::case_when(mean - sd > 0 ~ mean - sd,
                                      mean - sd < 0 ~ 0),
                hi = mean + sd)
  
#### Create figures ####

gg_probs_moved <- ggplot(data = move_meta_result_sum, aes(x = move_meta_sd, y = mean)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") + 
  geom_smooth(se = FALSE, linetype = 2, color = "#2755a8", size = 1) +
  geom_point(shape = 19) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = expression(paste(sigma, " Parameter ", italic(move_meta_sd))), y = "Mean cross-system movement") +
  theme_classic() 
