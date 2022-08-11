##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Simulation experiment of movement parameters and CV/Production

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

reps <- 500

matrix_lhs <- lhs::improvedLHS(n = reps, k = 2, dup = 2)

matrix_lhs[, 1] <- qunif(matrix_lhs[, 1], 0.1, 1.0) 

matrix_lhs[, 2] <- qunif(matrix_lhs[, 2], 0.1, 1.0) 

table(cut(matrix_lhs[, 1], breaks = seq(0.1, 1, 0.1)),
      cut(matrix_lhs[, 2], breaks = seq(0.1, 1, 0.1)))

df_experiment <- tibble::tibble(move_meta_sd = matrix_lhs[, 1], 
                                phase_sd = matrix_lhs[, 2]) %>% 
  dplyr::slice(rep(x = 1:dplyr::n(), times = 4)) %>% 
  dplyr::mutate(pop_n =  rep(x = c(8, 16, 32, 64), each = reps))

amplitude_mn <- 0.95

frequency <- years

#### Init HPC function ####

foo_hpc <- function(pop_n, move_meta_sd, phase_sd) {
  
  list_starting$pop_n <- pop_n
  
  # update move meta_sd parameters
  list_parameters$move_meta_sd <- move_meta_sd
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = matrix_reef,
                                         starting_values = list_starting, parameters = list_parameters,
                                         dimensions = dimensions, grain = grain, 
                                         use_log = FALSE, verbose = FALSE)
  
  # simulate nutrient input
  input_temp <- meta.arrR::simulate_nutrient_sine(n = n, max_i = max_i, frequency = frequency, 
                                                  input_mn = nutrient_input, amplitude_mn = amplitude_mn,
                                                  phase_sd = phase_sd, verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = list_parameters,
                                                nutrients_input = input_temp, movement = "behav",
                                                max_i = max_i, min_per_i = min_per_i,
                                                seagrass_each = seagrass_each, save_each = save_each, 
                                                verbose = FALSE)
  
  # filter model run
  result_temp <- meta.arrR::filter_meta(x = result_temp, filter = c((max_i / years) * years_filter, max_i), 
                                        reset = TRUE, verbose = FALSE)
  
  # get moved counts
  moved <- dplyr::bind_rows(result_temp$fishpop) %>% 
    dplyr::filter(timestep == max_i) %>% 
    dplyr::left_join(y = as.data.frame(metasyst_temp$fishpop_attr), by = "id") %>% 
    dplyr::select(id, moved, move_prob) %>% 
    dplyr::arrange(id)
  
  # calc cv
  cv <- meta.arrR::calc_variability(x = result_temp, lag = c(FALSE, TRUE))
  
  # calculate biomass/production
  prod <- meta.arrR::summarize_meta(result = result_temp, biomass = TRUE, production = TRUE, 
                                    lag = c(FALSE, FALSE)) %>% 
    purrr::map(function(i) { 
      dplyr::filter(i, timestep == max_i) %>% 
        tidyr::pivot_longer(-c(meta, timestep), names_to = "part") %>% 
        dplyr::group_by(part) %>% 
        dplyr::summarise(alpha = mean(value), gamma = sum(value)) %>% 
        tidyr::pivot_longer(-part, names_to = "measure", values_to = "value")
      
    })
  
  # combine to result data.frame and list
  list(moved = dplyr::mutate(moved, pop_n = pop_n, move_meta_sd = move_meta_sd, phase_sd = phase_sd), 
       cv = dplyr::mutate(dplyr::bind_rows(cv), pop_n = pop_n, move_meta_sd = move_meta_sd, phase_sd = phase_sd), 
       prod = dplyr::mutate(dplyr::bind_rows(prod), pop_n = pop_n, move_meta_sd = move_meta_sd, phase_sd = phase_sd))
  
}

#### Submit HPC

globals <- c("n", "max_i", "matrix_reef", "list_starting", "list_parameters", "dimensions", "grain", # setup_meta
             "years", "frequency", "nutrient_input", "amplitude_mn", # simulate_nutr_input
             "min_per_i", "seagrass_each", "save_each", # run_simulation_meta
             "years_filter") # filter_meta 

sbatch_cv <- rslurm::slurm_apply(f = foo_hpc, params = df_experiment, 
                                 global_objects = globals, jobname = "phase_sd",
                                 nodes = nrow(df_experiment), cpus_per_node = 1, 
                                 slurm_options = list("account" = account, 
                                                      "partition" = "standard",
                                                      "time" = "02:00:00", ## hh:mm::ss
                                                      "mem-per-cpu" = "7G",
                                                      "exclude" = exclude_nodes),
                                 pkgs = c("arrR", "dplyr", "meta.arrR", "purrr", "tidyr"),
                                 rscript_path = rscript_path, submit = FALSE)

#### Collect results #### 

suppoRt::rslurm_missing(x = sbatch_cv)

cv_result <- rslurm::get_slurm_out(sbatch_cv, outtype = "raw")

suppoRt::save_rds(object = cv_result, 
                  filename = paste0("05-variability-phase-", 
                                    stringr::str_remove(amplitude_mn, pattern = "\\."), 
                                    ".rds"),
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_cv)
