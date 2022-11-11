##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Plot dimensions of inital fish populations

#### Load setup ####

source("01_Functions/setup.R")

#### Simulate example populations ####

# create seafloor
input_seafloor <- arrR::setup_seafloor(dimensions = dimensions, grain = grain, 
                                       reef = reef_matrix, starting_values = starting_values_list, 
                                       verbose = FALSE)

fishpop_dim_df <- purrr::map_dfr(c(8, 16, 32, 64, 128), function(i) {
  
  purrr::map_dfr(1:iterations, function(j) {
    
    starting_values_list$pop_n <- i
    
    # create fishpop
    arrR::setup_fishpop(seafloor = input_seafloor, 
                        starting_values = starting_values_list, parameters = parameters_list,
                        use_log = FALSE, verbose = FALSE) |> 
      dplyr::summarise(length_sum = sum(length),
                       weight_sum = sum(weight)) |>
      dplyr::mutate(pop_n = i, itr = j, .before = "length_sum")
    
  })
})

#### Create ggplot ####

gg_fishpop_dim <- tidyr::pivot_longer(fishpop_dim_df, cols = c(length_sum, weight_sum)) |> 
  dplyr::mutate(pop_n = factor(pop_n, ordered = TRUE), 
                name = factor(name, levels = c("length_sum", "weight_sum"))) |> 
  ggplot(aes(x = pop_n, y = value)) + 
  geom_boxplot() +
  geom_jitter(alpha = 0.25) +
  facet_wrap(. ~ name, scales = "free_y", labeller = labeller(name = c("length_sum" = "Population length", 
                                                                       "weight_sum" = "Population weight"))) + 
  scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4), 
                     labels = function(x) round(x, 2)) +
  labs(x = "Population size", y = "Population sum") +
  theme_classic()

#### Save ggplot ####

# suppoRt::save_ggplot(plot = gg_fishpop_dim, filename = "gg_fishpop_dim.pdf", path = "04_Figures/Appendix/", 
#                      width = width, height = height * 0.65, units = units, dpi = dpi, 
#                      overwrite = FALSE)
