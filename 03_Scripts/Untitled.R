library(arrR)
library(meta.arrR)

source("01_Functions/setup.R")

# starting_values_list <- readRDS("02_Data/starting-values.rds")
# parameters_list <- readRDS("02_Data/parameter-values.rds")
nutrient_input_cell <- readRDS("02_Data/nutrient-input-cell.rds")

starting_values_list$pop_n <- 8 # 8, 16, 32, 64, 128 

# update move meta_sd parameters
parameters_list$move_meta_sd <- 0.5 # biotic
abiotic <- 0.5
nutrient_input <- nutrient_input_cell$excretion_cell |> max()

# parameters_list$nutrients_loss <- 0.01
# parameters_list$detritus_loss <- 0.01

# parameters_list$nutrients_diffusion <- 0.5
parameters_list$detritus_diffusion <- 0.01
parameters_list$detritus_fish_diffusion <- 0.01

# parameters_list$seagrass_slough <- 0.01
# parameters_list$detritus_mineralization <- 0.01

amplitude_mn <- 0.95

stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
                                              ag_biomass = starting_values_list$ag_biomass,
                                              parameters = parameters_list)

starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool

starting_values_list$detritus_pool <- stable_values_list$detritus_pool

# create seafloor
input_seafloor <- setup_seafloor(dimensions = dimensions, grain = grain, 
                                 reef = reef_matrix, starting_values = starting_values_list)

# create fishpop
input_fishpop <- setup_fishpop(seafloor = input_seafloor, starting_values = starting_values_list, 
                               parameters = parameters_list)

# simulate nutrient input
input_temp <- meta.arrR::simulate_nutrient_noise(n = 1, max_i = max_i, frequency = frequency, 
                                                 input_mn = nutrient_input, amplitude_mn = amplitude_mn,
                                                 noise_sd = abiotic, verbose = FALSE)

# plot(input_temp$values$meta_1$timestep, input_temp$values$meta_1$input, type = "l")

(result <- run_simulation(seafloor = input_seafloor, fishpop = input_fishpop, input_temp$values$meta_1$input,
                         parameters = parameters_list, movement = "behav", torus_diffusion = FALSE,
                         max_i = max_i, min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each))

plot(result)
plot_production(result)

plot(result, summarize = TRUE)
plot_production(result, summarize = TRUE, lag = TRUE)

plot(result, summarize = TRUE, what = "fishpop")

result$fishpop |> dplyr::group_by(id) |> dplyr::summarise(died_consumption = max(died_consumption)) |> dplyr::arrange(-died_consumption)
result$fishpop |> dplyr::group_by(id) |> dplyr::summarise(age = mean(age)) |> dplyr::arrange(-age)

