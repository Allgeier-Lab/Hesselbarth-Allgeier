##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

# Rcpp::sourceCpp("01_Functions/random_move.cpp")

#### Increase sd ####

residence_sd_parameters <- seq(from = 0, to = 1, by = 0.25)

fishpop_attr <- purrr::map_dfr(1:iterations, function(i) {
  
  purrr::map_dfr(residence_sd_parameters, function(j) {
    
    message("\r> Progress: i=", i, "/", iterations, appendLF = FALSE)
  
    sd_temp <- parameters_list$move_residence_mean * j
    
    parameters_list$move_residence_sd <- sd_temp
    
    input_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                        starting_values = starting_list, parameters = parameters_list,
                                        dimensions = dimensions, use_log = use_log, 
                                        verbose = FALSE)
    
    dplyr::bind_cols(itr = i, mean = parameters_list$move_residence_mean, sd = round(sd_temp, 2),
                     input_temp$fishpop_attr[, c(1,3)])})}) %>% 
  dplyr::mutate(sd_h = sd * min_per_i / 60,
                residence_h = residence * min_per_i / 60, 
                sd_h = factor(sd_h, ordered = TRUE))
  
dplyr::filter(fishpop_attr, sd_h != 0) %>% 
  ggplot(aes(x = residence_h, color = sd_h, fill = sd_h)) +
  geom_vline(xintercept = parameters_list$move_residence_mean * 120 / 60,
             color = "black", linetype = 1) + 
  geom_density(alpha = 0.15, bw = 1.0) +
  annotate(geom = "text", x = 22.5, y = 0.095, 
           label = "sd residence time=0") + 
  scale_color_manual(name = "sd residence time",
                     values = c("#a82203", "#208cc0", "#f1af3a", "#637b31")) + 
  scale_fill_manual(name = "sd residence time",
                    values = c("#a82203", "#208cc0", "#f1af3a", "#637b31")) + 
  scale_x_continuous(breaks = seq(from = 0, to = 60, by = 5)) +
  labs(x = "Individual moves after x hours", y = "Density") +
  theme_classic() + theme(legend.position = c(0.9, 0.5))

# #### Get distribution of residence ####
# 
# residence_parameters <- seq(from = 6, to = 102, by = 6)
# 
# reps <- 1000000
# 
# result_df <- purrr::map_dfr(seq_along(residence_parameters), function(i) {
#   
#   message("\n\n> Progress: ", i, "/", length(residence_parameters))
# 
#   result_temp <- random_move(n = reps, move_residence = residence_parameters[i])
#   
#   data.frame(param = rep(x = residence_parameters[i], times = reps),
#              moved = result_temp)}) %>%
#   dplyr::mutate(param_h = param * min_per_i / 60, moved_h = moved * min_per_i / 60)
# 
# result_df_sum <- dplyr::group_by(result_df, param, param_h) %>% 
#   dplyr::summarise(moved_h_mn = mean(moved_h), moved_h_sd = sd(moved_h))
# 
# param_filter <- 96
# 
# mean_move <- dplyr::filter(result_df_sum, param_h == param_filter) %>% dplyr::pull(moved_h_mn)
# 
# sd_move <- dplyr::filter(result_df_sum, param_h == param_filter) %>% dplyr::pull(moved_h_sd)
# 
# max_move <- dplyr::filter(result_df, param_h == param_filter) %>% dplyr::pull(moved_h) %>% max()
# 
# gg_move_probs <- ggplot(data = dplyr::filter(result_df, param_h == param_filter), 
#                         aes(x = moved_h, y = ..density..)) + 
#   geom_histogram(fill = "#85D4E3", col = "black", alpha = 1, binwidth = 2) +
#   geom_density(fill = NA, col = "black", bw = 2) +
#   geom_vline(xintercept = mean_move, col = "#EF6567", linetype = 2, size = 1) +
#   annotate(geom = "segment", x = mean_move,  xend = mean_move - sd_move, y = 0.0475, yend = 0.0475,
#            arrow = ggplot2::arrow(length = unit(0.25, "cm")), col = "#EF6567") +
#   annotate(geom = "segment", x = mean_move,  xend = mean_move + sd_move, y = 0.0475, yend = 0.0475,
#            arrow = ggplot2::arrow(length = unit(0.25, "cm")), col = "#EF6567") +
#   annotate(geom = "label", x = mean_move, y = 0.05, label = "mean ± sd", col = "#EF6567") +
#   annotate(geom = "text", x = 40, y = 0.05, label = paste0("Theoretical max. residence = ", param_filter, "h")) +
#   annotate(geom = "text", x = 40, y = 0.0475, label = paste0("Mean residence = ", round(mean_move, 2), "h")) +
#   annotate(geom = "text", x = 40, y = 0.045, label = paste0("Actual max. residence = ", max_move, "h")) +
#   scale_x_continuous(breaks = seq(from = 0, to = max_move, by = 4), limits = c(0, max_move)) +
#   labs(x = "Moved after x hours", y = "Density") +
#   theme_classic()

# #### Save ggplot ####
# 
# suppoRt::save_ggplot(plot = gg_move_probs, filename = "00_residence-param.png",
#                      path = "04_Figures", width = height, height = width, dpi = dpi, 
#                      units = units, overwrite = overwrite)
