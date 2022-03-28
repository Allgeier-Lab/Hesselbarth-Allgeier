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
                                     parameters = parameters_list)

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

# simulate input
nutrients_input <- meta.arrR::sim_nutr_input(n = 1, max_i = max_i,
                                             amplitude_mod = 0.05, phase_mod = 0,
                                             input_mn = stable_values$nutrients_input,
                                             freq_mn = freq_mn, verbose = FALSE)

df_experiment <- expand.grid(amplitude_mod = c(0.05, 0.5, 1.0), 
                             fishpop = c(1, 2)) %>%
  dplyr::mutate(move = dplyr::case_when(fishpop == 1 ~ "rand", 
                                        fishpop == 2 ~ "behav"))

input_list <- list(NULL, input_fishpop)

model_runs <- purrr::map(1:nrow(df_experiment), function(i) {
  
  # simulate input
  nutrients_input <- meta.arrR::sim_nutr_input(n = 1, max_i = max_i,
                                               amplitude_mod = df_experiment[i, 1], phase_mod = 0,
                                               input_mn = stable_values$nutrients_input,
                                               freq_mn = freq_mn, verbose = FALSE)
  
  run_simulation(seafloor = input_seafloor, fishpop = input_list[[df_experiment[i, 2]]],
                 nutrients_input = nutrients_input$values$meta_1$input,
                 parameters = parameters_list, movement = df_experiment[i, 3],
                 max_i = max_i, min_per_i = min_per_i,
                 seagrass_each = seagrass_each, save_each = save_each)})

seagrass_diff <- purrr::map_dfr(seq_along(model_runs), function(i) {
  
  dplyr::filter(model_runs[[i]]$seafloor, timestep == max_i, reef == 0) %>%
    dplyr::summarise(dplyr::across(c(ag_biomass, bg_biomass, ag_production, bg_production),
                                   function(x) c(mean = mean(x), sum = sum(x)))) %>%
    tibble::add_column(fun = c("mean", "sum"), .before = "ag_biomass") %>%
    tidyr::pivot_longer(-fun) %>% 
    dplyr::mutate(amplitude = df_experiment[i, 1], 
                  fishpop = ifelse(test = df_experiment[i, 2] == 1, yes = "nofish", no = "immob"))}) %>%
  tidyr::pivot_wider(names_from = fishpop, values_from = value) %>% 
  dplyr::mutate(value.diff = ((immob - nofish) / immob) * 100)

purrr::map_dfr(seq_along(model_runs), function(i) { 
  arrR::get_production(model_runs[[i]]) %>% 
  dplyr::mutate(amplitude = df_experiment[i, 1], 
                fishpop = ifelse(test = df_experiment[i, 2] == 1, 
                                 yes = "nofish", no = "immob"))
  }) %>%
  tidyr::pivot_longer(-c(timestep, amplitude, fishpop), names_to = "part") %>% 
  ggplot() + 
  geom_line(aes(x = timestep, y = value, col = fishpop)) + 
  facet_wrap(. ~ part * amplitude, scales = "free_y") + 
  labs(x = "Timesteps", y = expression(paste(Delta, "PP (ttl environ)"))) +
  theme_classic()
