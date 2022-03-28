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

parameters_list$move_residence <- NULL
parameters_list$move_residence_var <- NULL
parameters_list$move_lambda <- NULL

# check if all parameters are present and meaningful
check_parameters(starting_values = starting_list, parameters = parameters_list)

#### setup_inputs ####
# create 5 reef cells in center of seafloor
reef_matrix <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0),
                      ncol = 2, byrow = TRUE)

# get stable nutrient/detritus values
stable_values <- arrR::get_stable_values(bg_biomass = starting_list$bg_biomass,
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

starting_list$nutrients_pool <- stable_values$nutrients_pool

starting_list$detritus_pool <- stable_values$detritus_pool

# create seafloor
input_seafloor <- arrR::setup_seafloor(dimensions = c(50, 50), grain = 1, 
                                       reef = reef_matrix, starting_values = starting_list)

# create fishpop
input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor, 
                                     starting_values = starting_list, 
                                     parameters = parameters_list, 
                                     use_log = use_log)

#### run_sim ####

# one iterations equals 120 minutes
min_per_i <- 120

# run the model for ten years
years <- 50
max_i <- (60 * 24 * 365 * years) / min_per_i

# run seagrass once each day
days <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days

# save results only every 365 days
days <- 365 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)
save_each <- (24 / (min_per_i / 60)) * days

df_experiment <- expand.grid(enrichment_levels = enrichment_levels, 
                             amplitude_mod = amplitude_mod) %>% 
  dplyr::slice(rep(1:dplyr::n(), each = iterations))

#### Setup function ####

foo <- function(enrichment_levels, amplitude_mod) {
  
  # simulate input
  input_nutrients <- meta.arrR::sim_nutr_input(n = 1, max_i = max_i,
                                               amplitude_mod = amplitude_mod, phase_mod = 0,
                                               input_mn = stable_values$nutrients_input * enrichment_levels,
                                               freq_mn = freq_mn, verbose = FALSE)
  
  # create seafloor
  input_seafloor <- arrR::setup_seafloor(dimensions = c(50, 50), grain = 1, 
                                         reef = reef_matrix, starting_values = starting_list)
  
  # create fishpop
  input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor, 
                                       starting_values = starting_list, 
                                       parameters = parameters_list,
                                       use_log = use_log)
  
  # create list with no fishpop and fishpop
  input_list <- list(NULL, input_fishpop)
  
  # run model
  result_temp <- purrr::map(input_list, function(i) {
    
    run_simulation(seafloor = input_seafloor, fishpop = i,
                   nutrients_input = input_nutrients$values$meta_1$input,
                   parameters = parameters_list, movement = "behav",
                   max_i = max_i, min_per_i = min_per_i,
                   seagrass_each = seagrass_each, save_each = save_each, 
                   verbose = FALSE)
    
  })
  
  # loop through all model runs
  result_temp <- purrr::map_dfr(result_temp, function(j) {
    
    # filter only last timestep
    seafloor_temp <- dplyr::filter(j$seafloor, timestep == max_i, reef == 0)
    
    # calculate mean and sum biomass and production
    seafloor_temp <- dplyr::summarise(seafloor_temp, 
                                      dplyr::across(c(ag_biomass, bg_biomass, ag_production, bg_production),
                                                    function(x) c(mean = mean(x), sum = sum(x))))
    
    # add function column
    seafloor_temp <- tibble::add_column(seafloor_temp, fun = c("mean", "sum"), 
                                        .before = "ag_biomass") 
    
    # reshape to longer format
    seafloor_temp <- tidyr::pivot_longer(seafloor_temp, -fun) 
    
    # add treatment levels
    dplyr::mutate(seafloor_temp, enrichment = enrichment_levels, amplitude = amplitude_mod, 
                  fishpop = ifelse(test = nrow(j$fishpop) == 0, yes = "nofish", no = "local"))
    
  })
  
  # reshape wider using nofish and local
  result_temp <- tidyr::pivot_wider(result_temp, names_from = fishpop, values_from = value)
  
  return(result_temp)
  
}

#### Submit to HPC ####

globals <- c("stable_values", "freq_mn", "reef_matrix", "starting_list", "parameters_list", 
             "use_log", "max_i", "min_per_i", "seagrass_each", "save_each")

# create .sh script
sbatch_pp <- rslurm::slurm_apply(f = foo, params = df_experiment,
                                 global_objects = globals, jobname = "PP_nofish_local",
                                 nodes = nrow(df_experiment), cpus_per_node = 1, 
                                 slurm_options = list("account" = account, 
                                                      "partition" = "standard",
                                                      "time" = "00:30:00", ## hh:mm::ss
                                                      "mem-per-cpu" = "7G"),
                                 pkgs = c("arrR", "purrr", "dplyr", "tidyr", "tibble"),
                                 rscript_path = rscript_path, sh_template = sh_template, 
                                 submit = FALSE)

# check results
suppoRt::rslurm_missing(sbatch_pp)

# get results as list
pp_nofish_local <- rslurm::get_slurm_out(sbatch_pp, outtype = "table")

# save results to disk
suppoRt::save_rds(object = pp_nofish_local, file = "02_Data/00_pp_nofish_local.rds", 
                  overwrite = overwrite)

# delete .sh scripts
rslurm::cleanup_files(sbatch_pp)

#### Results treatments  ####

pp_nofish_local <- readr::read_rds("02_Data/00_pp_nofish_local.rds") %>% 
  dplyr::mutate(reldiff = ((local - nofish) / nofish) * 100) %>% 
  dplyr::mutate(fun = factor(fun), name = factor(name), 
                enrichment = factor(enrichment, ordered = TRUE),
                amplitude = factor(amplitude, ordered = TRUE))

dplyr::filter(pp_nofish_local, fun == "sum") %>% 
  ggplot() +
  geom_boxplot(aes(x = enrichment, y = reldiff, fill = amplitude),
               position = position_dodge(), alpha = 0.75) +
  geom_point(aes(x = enrichment, y = reldiff, fill = amplitude, col = amplitude), 
             position = position_jitterdodge(), alpha = 0.25) + 
  facet_wrap(. ~ name, scales = "free_y") + 
  labs(x = "Enrichment treatment", y = expression(paste(Delta, "no fish, local fish [%]"))) +
  scale_fill_viridis_d(name = "Amplitude treatment") +
  theme_classic() + theme(legend.position = "bottom")

dplyr::group_by(pp_nofish_local, fun, name, enrichment, amplitude) %>%
  dplyr::summarise(reldiff_mn = mean(reldiff), reldiff_sd = sd(reldiff),
                   .groups = "drop") %>%
  dplyr::arrange(desc(fun), name, enrichment, amplitude) %>% 
  dplyr::filter(fun == "sum")
