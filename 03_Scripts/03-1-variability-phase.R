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

mem_per_cpu <- "7G"
time <- "05:00:00" # hh:mm:ss

#### Stable values #### 

stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
                                              ag_biomass = starting_values_list$ag_biomass,
                                              parameters = parameters_list)

starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool

starting_values_list$detritus_pool <- stable_values_list$detritus_pool

#### Setup experiment ####

experiment_df <- readRDS("02_Data/experiment-parameters.rds")

nutrient_input_fish <- readRDS("02_Data/nutrient-input-fish.rds") |> 
  dplyr::group_by(pop_n) |> 
  dplyr::summarise(excretion_cell = mean(excretion_cell)) |> 
  dplyr::mutate(pop_n = factor(pop_n, ordered = TRUE))

# nutrient_levels <- nutrient_input_fish$excretion_cell
nutrient_levels <- nutrient_input_fish$excretion_cell[1] * seq(from = 1, to = 3, length.out = 3)

experiment_full_df <- dplyr::slice(experiment_df, rep(x = 1:dplyr::n(), times = length(nutrient_levels))) |> 
  dplyr::mutate(nutrient_input = rep(x = nutrient_levels, each = nrow(experiment_df)))

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
                                                  input_mn = nutrient_input, noise = abiotic,
                                                  verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                                nutrients_input = input_temp, movement = "behav", 
                                                torus_diffusion = TRUE, max_i = max_i, min_per_i = min_per_i,
                                                seagrass_each = seagrass_each, save_each = save_each, 
                                                verbose = FALSE) |> 
    meta.arrR::filter_meta(filter = c((max_i / years) * years_filter, max_i),
                           reset = TRUE, verbose = FALSE)
  
  # filter local system for cells only near reef
  result_temp_near <- result_temp
  
  for (i in 1:n) {
    
    result_temp_near$seafloor[[i]] <- dplyr::mutate(result_temp_near$seafloor[[i]], 
                                                    distance = sqrt((0 - x)^2 + (0 - y)^2)) |>
      dplyr::filter(distance <= 12.5)
    
  }
  
  # get moved counts
  connectivity <- dplyr::bind_rows(result_temp$fishpop) |> 
    dplyr::filter(timestep == max_i) |> 
    dplyr::left_join(y = as.data.frame(metasyst_temp$fishpop_attr), by = "id") |> 
    dplyr::select(id, moved, move_prob) |> 
    dplyr::arrange(id)
  
  # calc mortality
  mortality <- dplyr::bind_rows(result_temp$fishpop) |> 
    dplyr::filter(timestep == max_i) |> 
    dplyr::mutate(died_total = died_consumption + died_background) |> 
    dplyr::select(id, died_consumption, died_background, died_total) |> 
    dplyr::arrange(id)

  abundance <- meta.arrR::get_abundance(result_temp) |> 
    dplyr::group_by(meta) |> 
    dplyr::summarise(min = min(abundance), mean = mean(abundance), max = max(abundance))
  
  # calculate biomass/production
  prod_cumulative <- meta.arrR::summarize_meta(result = result_temp, biomass = FALSE, production = TRUE, 
                                               lag = c(FALSE, FALSE))[["production"]] |> 
    dplyr::filter(timestep == max_i) |>
    tidyr::pivot_longer(-c(meta, timestep), names_to = "part") |> 
    dplyr::group_by(part) |> 
    dplyr::summarise(alpha = mean(value), gamma = sum(value),  .groups = "drop") |> 
    tidyr::pivot_longer(-part, names_to = "measure", values_to = "value")
  
  prod_cumulative_near <- meta.arrR::summarize_meta(result = result_temp_near, biomass = FALSE, production = TRUE, 
                                                    lag = c(FALSE, FALSE))[["production"]] |> 
    dplyr::filter(timestep == max_i) |>
    tidyr::pivot_longer(-c(meta, timestep), names_to = "part") |> 
    dplyr::group_by(part) |> 
    dplyr::summarise(alpha = mean(value), gamma = sum(value),  .groups = "drop") |> 
    tidyr::pivot_longer(-part, names_to = "measure", values_to = "value")
  
  # calc cv
  cv <- meta.arrR::calc_variability(x = result_temp, biomass = FALSE, production = TRUE,
                                    lag = c(FALSE, TRUE))[["production"]]
  
  cv_near <- meta.arrR::calc_variability(x = result_temp_near, biomass = FALSE, production = TRUE,
                                         lag = c(FALSE, TRUE))[["production"]]
  
  # combine to result data.frame and list
  list(fishpop_init = dplyr::mutate(dplyr::bind_rows(metasyst_temp$fishpop), pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input), 
       connectivity = dplyr::mutate(connectivity, pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input), 
       mortality = dplyr::mutate(mortality, pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input), 
       abundance = dplyr::mutate(abundance, pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input), 
       prod_cumulative = dplyr::mutate(prod_cumulative, pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input), 
       prod_cumulative_near = dplyr::mutate(prod_cumulative_near, pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input), 
       cv = dplyr::mutate(dplyr::bind_rows(cv), pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input),
       cv_near = dplyr::mutate(dplyr::bind_rows(cv_near), pop_n = pop_n, biotic = biotic, abiotic = abiotic, nutrient_input = nutrient_input))
  
}

#### Submit HPC

globals <- c("n", "max_i", "reef_matrix", "starting_values_list", "parameters_list", "dimensions", "grain", # setup_meta
             "frequency", # simulate_nutr_input
             "min_per_i", "seagrass_each", "save_each", # run_simulation_meta
             "years", "years_filter") # filter_meta 

sbatch_phase <- rslurm::slurm_apply(f = foo_hpc, params = experiment_full_df, 
                                    global_objects = globals, jobname = "phase_sd",
                                    nodes = nrow(experiment_full_df), cpus_per_node = 1, 
                                    slurm_options = list("account" = account, 
                                                         "partition" = "standard",
                                                         "time" = time,
                                                         "mem-per-cpu" = mem_per_cpu),
                                    pkgs = c("arrR", "dplyr", "meta.arrR", "purrr", "tidyr"),
                                    rscript_path = rscript_path, submit = FALSE)

#### Collect results #### 

suppoRt::rslurm_missing(x = sbatch_phase)

cv_phase <- rslurm::get_slurm_out(sbatch_phase, outtype = "raw")

suppoRt::save_rds(object = cv_phase, path = "02_Data/", 
                  overwrite = FALSE, filename = "result-phase.rds")

rslurm::cleanup_files(sbatch_phase)
