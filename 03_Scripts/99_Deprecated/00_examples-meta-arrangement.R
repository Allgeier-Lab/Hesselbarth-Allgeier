##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/get_modifier.R")

#### Adapt parameters ####



#### Create reef and coords object #### 

seafloor_xy_rand <- cbind(id = 1:n, x = runif(n = 9, min = -1, max = 1), 
                          y = runif(n = 9, min = -1, max = 1))

seafloor_xy_clust <- cbind(id = 1:n, x = c(runif(n = 8, min = -1, max = 0), 
                                           runif(n = 1, min = 0.75, max = 1)), 
                           y = c(runif(n = 8, min = -1, max = 0), 
                                 runif(n = 1, min = 0.75, max = 1)))

seafloor_xy_reg <- cbind(id = 1:n, expand.grid(seq(from = -0.75, to = 0.75, length.out = 3), 
                                               seq(from = -0.75, to = 0.75, length.out = 3)))

names(seafloor_xy_reg)[2:3] <- c("x", "y")

seafloor_xy_full <- list(rand = seafloor_xy_rand, reg = seafloor_xy_reg, clust = seafloor_xy_clust)

#### Stable values ####

list_stable <- arrR::get_stable_values(bg_biomass = list_starting$bg_biomass,
                                         ag_biomass = list_starting$ag_biomass,
                                         parameters = list_parameters)

list_starting$nutrients_pool <- list_stable$nutrients_pool

list_starting$detritus_pool <- list_stable$detritus_pool

#### Different ecosystem placements ####

metasyst_full <- purrr::map(seafloor_xy_full, function(i) 
  meta.arrR::setup_meta(n = n, max_i = max_i, seafloor_xy = i,
                        starting_values = list_starting, parameters = list_parameters, reef = matrix_reef,
                        dimensions = dimensions, grain = grain, verbose = FALSE))

purrr::map(metasyst_full, plot, viridis_option = "B", base_size = 12)

#### Residence time variability ####

metassyst <- meta.arrR::setup_meta(n = n, max_i = max_i, seafloor_xy = NULL,
                                   starting_values = list_starting, parameters = list_parameters, reef = matrix_reef,
                                   dimensions = dimensions, grain = grain, verbose = FALSE)

variability_vec <- c(0.1, 0.25, 0.35)

residence_var <- purrr::map_dfr(variability_vec, function(i) {
  
  # save results only every m days
  days <- 15 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)
  move_each <- (24 / (min_per_i / 60)) * days
  
  data.frame(id = 1:100, residence = rnorm(n = 100, mean = move_each, sd = move_each * i))}, 
  .id = "variability") %>%
  dplyr::mutate(variability = dplyr::case_when(variability == 1 ~ variability_vec[1], 
                                            variability == 2 ~ variability_vec[2], 
                                            variability == 3 ~ variability_vec[3]), 
                variability = factor(variability, labels = c("low", "medium", "high")), 
                residence = ((residence * min_per_i) / 60) / 24)

ggplot(data = residence_var, aes(x = residence)) + 
  geom_histogram(aes(y = ..count.., fill = variability, col = variability), 
                 position = "identity", binwidth = 1, alpha = 0.25) + 
  geom_density(aes(y = ..count.., fill = variability, col = variability), alpha = 0.5) +
  geom_vline(xintercept = 15, linetype = 2, size = 1) +
  scale_fill_viridis_d(name = "Variability") + scale_color_viridis_d(name = "Variability") +
  labs(x = "Residence time [days]", y = "# Indiividuals") + 
  theme_classic()

#### Distance decay ####

variability_vec <- c(0.01, 0.25, 0.5)

dist <- seq(from = 0, to = 1, by = 0.01)

decay_var <- purrr::map_dfr(variability_vec, function(i) {
  
  lambda <- runif(n = 25, min = 1 * (1 - i), max = 1 * (1 + i))
  
  purrr::map_dfr(lambda, function(j) data.frame(dist = dist, p = exp(-dist * j)), .id = "id")}, 
  .id = "variability") %>% 
  dplyr::mutate(variability = dplyr::case_when(variability == 1 ~ variability_vec[1], 
                                               variability == 2 ~ variability_vec[2], 
                                               variability == 3 ~ variability_vec[3]), 
                variability = factor(variability, labels = c("Low variability", "Medium variability", "High variability")))

ggplot(data = decay_var) + 
  geom_line(aes(x = dist, y = p, group = id)) + 
  facet_wrap(. ~ variability) +
  scale_y_continuous(limits = c(0, 1)) + 
  # scale_fill_viridis_d(name = "Variability") + scale_color_viridis_d(name = "Variability") +
  labs(x = expression(paste("Distance between ", italic("i,j"), " [unitless]")),
       y = expression(paste("Probability reaching ", italic("j")))) + 
  theme_classic()
