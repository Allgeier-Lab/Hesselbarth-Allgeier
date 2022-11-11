##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Calculate time period individuals are within different behavior stages

#### Load setup ####

source("01_Functions/setup.R")

#### Change parameters and starting values ####

# Nothing to change

#### Setup environment #### 

# create 5 reef cells in center of seafloor
reef_matrix <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0),
                      ncol = 2, byrow = TRUE)

# get stable nutrient/detritus values
stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
                                         ag_biomass = starting_values_list$ag_biomass,
                                         parameters = parameters_list)

starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool

starting_values_list$detritus_pool <- stable_values_list$detritus_pool

#### Setup experiment #### 

# one iterations equals 120 minutes
min_per_i <- 120

# run the model for ten years
years <- 1
max_i <- (60 * 24 * 365 * years) / min_per_i

# run seagrass once each day
days_seagrass <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days_seagrass

values <- c(0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2)

input_df <- expand.grid(pop_reserves_max = values, pop_reserves_consump = values, 
                        pop_reserves_thres_mean = values) |> 
  dplyr::slice(rep(1:dplyr::n(), each = 5))  

foo <- function(pop_reserves_max, pop_reserves_consump, pop_reserves_thres_mean) {
  
  parameters_list$pop_reserves_max <- pop_reserves_max
  
  parameters_list$pop_reserves_consump <- pop_reserves_consump
  
  parameters_list$pop_reserves_thres_mean <- pop_reserves_thres_mean
  
  # create seafloor
  input_seafloor <- arrR::setup_seafloor(dimensions = dimensions, reef = reef_matrix, 
                                         starting_values = starting_values_list, verbose = FALSE)
  
  # create fishpop
  input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor, starting_values = starting_values_list, 
                                       parameters = parameters_list, use_log = use_log, 
                                       verbose = FALSE)
  
  result <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                                 nutrients_input = nutrient_input, parameters = parameters_list, 
                                 movement = "behav", max_i = max_i, min_per_i = min_per_i,
                                 seagrass_each = seagrass_each, save_each = 1, verbose = FALSE)

  behavior <- purrr::map_dfr(1:nrow(input_fishpop), function(i) {

    counts_temp <- dplyr::mutate(result$fishpop,
                                 behavior = dplyr::case_when(behavior %in% c(1, 2) ~ "shelter",
                                                             behavior == 3 ~ "forage"))

    counts_temp <- dplyr::filter(counts_temp, id == i, timestep > 0)

    counts_temp <- dplyr::pull(counts_temp, behavior)

    counts_temp <- rle(counts_temp)

    data.frame(pop_reserves_max = pop_reserves_max, pop_reserves_consump = pop_reserves_consump,
               pop_reserves_thres_mean = pop_reserves_thres_mean, id = i,
               shelter = mean(counts_temp$lengths[which(counts_temp$values == "shelter")]),
               forage = mean(counts_temp$lengths[which(counts_temp$values == "forage")]))
    
  })
  
  return(behavior)

}

#### Submit to HPC model #### 

globals <- c("dimensions", "reef_matrix", "starting_values_list", 
             "parameters_list", "use_log",
             "nutrient_input", "max_i", "min_per_i", "seagrass_each")

# create .sh script
sbatch_behav <- rslurm::slurm_apply(f = foo, params = input_df,
                                    global_objects = globals, jobname = "move_behav",
                                    nodes = nrow(input_df), cpus_per_node = 1, 
                                    slurm_options = list("account" = account, 
                                                         "partition" = "standard",
                                                         "time" = "00:30:00", ## hh:mm::ss
                                                         "mem-per-cpu" = "7G"),
                                    pkgs = c("arrR", "purrr", "dplyr"),
                                    rscript_path = rscript_path, submit = FALSE)

# check results
suppoRt::rslurm_missing(sbatch_behav)

# get results as list
behavior_states <- rslurm::get_slurm_out(sbatch_behav, outtype = "table")

# save results to disk
suppoRt::save_rds(object = behavior_states, file = "02_Data/behavior-states.rds", 
                  overwrite = overwrite)

# delete .sh scripts
rslurm::cleanup_files(sbatch_behav)

#### Results parameter space ####

behavior_states <- readr::read_rds("02_Data/behavior-states.rds")

cutoff <- 3.5

(behavior_states_mn <- dplyr::group_by(behavior_states, pop_reserves_max, pop_reserves_consump, 
                                      pop_reserves_thres_mean) |> 
  dplyr::summarise(shelter = mean(shelter), forage = mean(forage), .groups = "drop") |> 
  dplyr::mutate(shelter = shelter * min_per_i / 60, forage = forage * min_per_i / 60, 
                diff_shelter = abs(12 - shelter),  diff_forage = abs(12 - forage)) |> 
  dplyr::filter(diff_shelter <= cutoff, diff_forage <= cutoff))

#### Run model for longer time with parametrization ####

parameters_list$pop_reserves_max <- 0.01

parameters_list$pop_reserves_consump <- 0.20

parameters_list$pop_reserves_thres_mean <- 0.01

# create seafloor
input_seafloor <- arrR::setup_seafloor(dimensions = dimensions, reef = reef_matrix, 
                                       starting_values = starting_values_list)

# create fishpop
input_fishpop <- arrR::setup_fishpop(seafloor = input_seafloor, starting_values = starting_values_list, 
                                     parameters = parameters_list, use_log = use_log)

# one iterations equals 120 minutes
min_per_i <- 120

# run the model for n years
years <- 100
max_i <- (60 * 24 * 365 * years) / min_per_i

# run seagrass once each day
days_seagrass <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days_seagrass

# # save results only every n days
days_save <- 365 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)
save_each <- (24 / (min_per_i / 60)) * days_save

(result <- arrR::run_simulation(seafloor = input_seafloor, fishpop = input_fishpop,
                                nutrients_input = nutrient_input, parameters = parameters_list, 
                                movement = "behav", max_i = max_i, min_per_i = min_per_i,
                                seagrass_each = seagrass_each, save_each = save_each))

plot(result)
plot(result, summarize = TRUE)
# plot(result, what = "fishpop", verbose = FALSE)

dplyr::filter(result$seafloor, timestep == max_i) |> 
  dplyr::select(x, y, consumption, excretion) |> 
  tidyr::pivot_longer(-c(x, y)) |> 
  dplyr::mutate(value = log(value)) |> 
  ggplot(aes(x = x, y = y)) +
  geom_raster(aes(fill = value)) + 
  facet_wrap(. ~ name) + 
  scale_fill_gradientn(name = "nutr", colours = c("#368AC0", "#F4B5BD", "#EC747F"), 
                       na.value = "#9B964A") +
  coord_fixed(ratio = 1) + 
  theme_classic()

dplyr::filter(result$seafloor, timestep == max_i) |> 
  dplyr::select(x, y, consumption, excretion) |> 
  dplyr::mutate(net = excretion - consumption) |> 
  dplyr::mutate(net = log(net + abs(min(net)))) |> 
  ggplot(aes(x = x, y = y)) +
  geom_raster(aes(fill = net)) + 
  scale_fill_gradientn(name = "net(nutr)", colours = c("#368AC0", "#F4B5BD", "#EC747F"), 
                       na.value = "#9B964A") +
  coord_fixed(ratio = 1) + 
  theme_classic()

dplyr::select(result$fishpop, id, x, y, behavior) |> 
  ggplot(aes(x = x, y = y)) + 
  geom_point(aes(col = factor(behavior)), shape = 1 ) + 
  geom_polygon(data = data.frame(x = c(result$extent[1], result$extent[2], 
                                       result$extent[2], result$extent[1]), 
                                 y = c(result$extent[3], result$extent[3], 
                                       result$extent[4], result$extent[4])), 
               aes(x = x, y = y), fill = "blue", col = NA, alpha = 0.15) +
  facet_wrap(. ~ behavior) + 
  coord_fixed(ratio = 1) + 
  theme_classic() + 
  theme(legend.position = "none")

table(result$fishpop$behavior)
