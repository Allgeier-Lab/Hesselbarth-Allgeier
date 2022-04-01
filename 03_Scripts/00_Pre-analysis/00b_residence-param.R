##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

Rcpp::sourceCpp("01_Functions/random_move.cpp")

#### Change parameters and starting values ####

# parameters_list$nutrients_loss <- 0.0 

# # check if all parameters are present and meaningful
# check_parameters(starting_values = starting_list, parameters = parameters_list)

#### Setup environment #### 

# # create 5 reef cells in center of seafloor
# reef_matrix <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0),
#                       ncol = 2, byrow = TRUE)
# 
# # get stable nutrient/detritus values
# stable_values <- arrR::get_stable_values(bg_biomass = starting_list$bg_biomass,
#                                          ag_biomass = starting_list$ag_biomass,
#                                          parameters = parameters_list)
# 
# starting_list$nutrients_pool <- stable_values$nutrients_pool
# 
# starting_list$detritus_pool <- stable_values$detritus_pool

#### Setup experiment #### 

# n <- 3
# 
# # one iterations equals 120 minutes
# min_per_i <- 120
# 
# # run the model for ten years
# years <- 1
# max_i <- (60 * 24 * 365 * years) / min_per_i
# 
# # run seagrass once each day
# days_seagrass <- 1
# seagrass_each <- (24 / (min_per_i / 60)) * days_seagrass
# 
# # create seafloor
# input_seafloor <- arrR::setup_seafloor(dimensions = c(50, 50), grain = 1, 
#                                        reef = reef_matrix, starting_values = starting_list)
# 
# res_par <- seq(from = 6, to = 102, by = 6)
# 
# result_moved <- purrr::map_dfr(seq_along(res_par), function(i) {
#   
#   parameters_list$move_residence <- res_par[[i]]
#   
#   # setup metaecosystems
#   metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
#                                          starting_values = starting_list, parameters = parameters_list,
#                                          dimensions = dimensions, grain = grain)
# 
#   # run model
#   result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
#                                                 nutrients_input = 0.0, movement = "behav",
#                                                 max_i = max_i, min_per_i = min_per_i,
#                                                 seagrass_each = seagrass_each, save_each = 1)
#   
#   dplyr::bind_rows(result_temp$fishpop) %>% 
#     dplyr::mutate(residence = dplyr::case_when(residence == 0 ~ 0, residence != 0 ~ 1)) %>% 
#     dplyr::group_by(id) %>% 
#     dplyr::group_split() %>% 
#     purrr::map_dbl(function(j) {
#       counts_temp <- dplyr::pull(j, residence) %>% rle()
#       mean(counts_temp$lengths[which(counts_temp$values == 1)])
#     }) %>% 
#     mean(na.rm = TRUE) %>% 
#     magrittr::multiply_by(min_per_i) %>% magrittr::divide_by(60) %>%
#     dplyr::bind_cols(residence = result_temp$parameters$move_residence, 
#                      residence_h = result_temp$parameters$move_residence * min_per_i / 24,
#                      local_h = .)
# 
# })

#### Get distribution of residence ####

residence_parameters <- seq(from = 6, to = 102, by = 6)

reps <- 5000000

result_df <- purrr::map_dfr(seq_along(residence_parameters), function(i) {

  result_temp <- random_move(n = reps, move_residence = residence_parameters[i])
  
  data.frame(param = rep(x = residence_parameters[i], times = reps),
             moved = result_temp)}) %>%
  dplyr::mutate(param_h = param * min_per_i / 60, moved_h = moved * min_per_i / 60)

result_df_sum <- dplyr::group_by(result_df, param, param_h) %>% 
  dplyr::summarise(moved_h_mn = mean(moved_h), moved_h_sd = sd(moved_h))

param_filter <- 96

mean_move <- dplyr::filter(result_df_sum, param_h == param_filter) %>% dplyr::pull(moved_h_mn)

sd_move <- dplyr::filter(result_df_sum, param_h == param_filter) %>% dplyr::pull(moved_h_sd)

max_move <- dplyr::filter(result_df, param_h == param_filter) %>% dplyr::pull(moved_h) %>% max()

gg_move_probs <- ggplot(data = dplyr::filter(result_df, param_h == param_filter), 
                        aes(x = moved_h, y = ..density..)) + 
  geom_histogram(fill = "#85D4E3", col = "black", alpha = 1, binwidth = 2) +
  geom_density(fill = NA, col = "black", bw = 2) +
  geom_vline(xintercept = mean_move, col = "#EF6567", linetype = 2, size = 1) +
  annotate(geom = "segment", x = mean_move,  xend = mean_move - sd_move, y = 0.0475, yend = 0.0475,
           arrow = ggplot2::arrow(length = unit(0.25, "cm")), col = "#EF6567") +
  annotate(geom = "segment", x = mean_move,  xend = mean_move + sd_move, y = 0.0475, yend = 0.0475,
           arrow = ggplot2::arrow(length = unit(0.25, "cm")), col = "#EF6567") +
  annotate(geom = "label", x = mean_move, y = 0.05, label = "mean ± sd", col = "#EF6567") +
  annotate(geom = "text", x = 40, y = 0.05, label = paste0("Theoretical max. residence = ", param_filter, "h")) +
  annotate(geom = "text", x = 40, y = 0.0475, label = paste0("Mean residence = ", round(mean_move, 2), "h")) +
  annotate(geom = "text", x = 40, y = 0.045, label = paste0("Actual max. residence = ", max_move, 2, "h")) +
  scale_x_continuous(breaks = seq(from = 0, to = max_move, by = 4), limits = c(0, max_move)) +
  labs(x = "Moved after x hours", y = "Density") +
  theme_classic()
