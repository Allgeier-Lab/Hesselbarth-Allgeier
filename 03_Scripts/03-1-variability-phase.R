##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Run simulation experiment for phase variability

#### Load setup ####

source("01_Functions/setup.R")

#### Adapt parameters ####

mem_per_cpu <- "7G" # "7G" "15G"
time <- "05:00:00" # hh:mm::ss # "02:00:00" "05:00:00"

#### Stable values #### 

stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
                                              ag_biomass = starting_values_list$ag_biomass,
                                              parameters = parameters_list)

starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool

starting_values_list$detritus_pool <- stable_values_list$detritus_pool

#### Setup experiment ####

experiment_df <- readRDS("02_Data/experiment-parameters.rds")

nutrient_input_cell <- readRDS("02_Data/nutrient-input-cell.rds")

experiment_full_df <- dplyr::slice(experiment_df, rep(x = 1:dplyr::n(), times = 3)) |> 
  dplyr::mutate(nutrient_input = rep(x = nutrient_input_cell$excretion_cell, 
                                     each = nrow(experiment_df)))

#### Init HPC function ####

foo_hpc <- function(pop_n, biotic, abiotic, nutrient_input) {
  
  starting_values_list$pop_n <- pop_n
  
  # update move meta_sd parameters
  parameters_list$move_meta_sd <- biotic
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                         starting_values = starting_values_list, parameters = parameters_list,
                                         dimensions = dimensions, grain = grain, 
                                         use_log = FALSE, verbose = FALSE)
  
  # simulate nutrient input
  input_temp <- meta.arrR::simulate_nutrient_sine(n = n, max_i = max_i, frequency = frequency, 
                                                  input_mn = nutrient_input, amplitude_mn = amplitude_mn,
                                                  phase_sd = abiotic, verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                                nutrients_input = input_temp, movement = "behav", 
                                                torus_diffusion = FALSE, max_i = max_i, min_per_i = min_per_i,
                                                seagrass_each = seagrass_each, save_each = save_each, 
                                                verbose = FALSE) |> 
    meta.arrR::filter_meta(filter = c((max_i / years) * years_filter, max_i), 
                           reset = TRUE, verbose = FALSE)
  
  # get moved counts
  moved <- dplyr::bind_rows(result_temp$fishpop) |> 
    dplyr::filter(timestep == max_i) |> 
    dplyr::left_join(y = as.data.frame(metasyst_temp$fishpop_attr), by = "id") |> 
    dplyr::select(id, moved, move_prob) |> 
    dplyr::arrange(id)
  
  # calc cv
  cv <- meta.arrR::calc_variability(x = result_temp, lag = c(FALSE, TRUE))
  
  # calculate biomass/production
  prod <- meta.arrR::summarize_meta(result = result_temp, biomass = TRUE, production = TRUE, 
                                    lag = c(FALSE, FALSE)) |> 
    purrr::map(function(i) { 
      dplyr::filter(i, timestep == max_i) |>
        tidyr::pivot_longer(-c(meta, timestep), names_to = "part") |> 
        dplyr::group_by(part) |> 
        dplyr::summarise(alpha = mean(value), gamma = sum(value)) |> 
        tidyr::pivot_longer(-part, names_to = "measure", values_to = "value")
      
    })
  
  # calculate lagged production over time
  prod_time <- meta.arrR::summarize_meta(result = result_temp, biomass = FALSE, production = TRUE, 
                                         lag = c(FALSE, TRUE))[["production"]]
  
  # combine to result data.frame and list
  list(fishpop_init = dplyr::mutate(dplyr::bind_rows(metasyst_temp$fishpop), pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input), 
       moved = dplyr::mutate(moved, pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input), 
       cv = dplyr::mutate(dplyr::bind_rows(cv), pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input), 
       prod = dplyr::mutate(dplyr::bind_rows(prod), pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input), 
       prod_time = dplyr::mutate(prod_time, pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input))
  
}

#### Submit HPC

globals <- c("n", "max_i", "reef_matrix", "starting_values_list", "parameters_list", "dimensions", "grain", # setup_meta
             "frequency", "amplitude_mn", # simulate_nutr_input
             "min_per_i", "seagrass_each", "save_each", # run_simulation_meta
             "years", "years_filter") # filter_meta 

sbatch_cv <- rslurm::slurm_apply(f = foo_hpc, params = experiment_full_df, 
                                 global_objects = globals, jobname = "phase_sd",
                                 nodes = nrow(experiment_full_df), cpus_per_node = 1, 
                                 slurm_options = list("account" = account, 
                                                      "partition" = "standard",
                                                      "time" = time,
                                                      "mem-per-cpu" = mem_per_cpu),
                                 pkgs = c("arrR", "dplyr", "meta.arrR", "purrr", "tidyr"),
                                 rscript_path = rscript_path, submit = FALSE)

#### Collect results #### 

suppoRt::rslurm_missing(x = sbatch_cv)

cv_result <- rslurm::get_slurm_out(sbatch_cv, outtype = "raw")

suppoRt::save_rds(object = cv_result, path = "02_Data/", 
                  overwrite = FALSE, filename = "result-phase.rds")

rslurm::cleanup_files(sbatch_cv)
