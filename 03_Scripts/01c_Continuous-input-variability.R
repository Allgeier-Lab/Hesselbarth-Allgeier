##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

source("05_Various/setup.R")

#### Basic parameters ####

# number of local metaecosystems
n <- 9

# set min_per_i
min_per_i <- 120

# run the model for n years
years <- 25

max_i <- (60 * 24 * 365 * years) / min_per_i

# setup nutrient input to be maximum 10% of inital nutrients
input_mn <- meta.arrR::meta.arrR_starting_values$nutrients_pool * 0.1

freq_mn <- years / 10

#### Setup experiment ####

itr <- 50

# create variability data.frame with all combinations 
variability_experiment <- expand.grid(amplitude = seq(from = 0, to = 1, by = 0.1), 
                                      phase = seq(from = 0, to = 1, by = 0.1)) %>% 
  dplyr::slice(rep(1:n(), each = itr))

#### Variability of input ####

# sample variability for different treatment lvls and increasing scale
variability_df <- purrr::map_dfr(1:nrow(variability_experiment), function(i) {
  
    # print progress
    message("\r> Variability: ", i, " / ", nrow(variability_experiment), "\t\t", 
            appendLF = FALSE)
    
    # simulate input
    input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                            variability = as.numeric(variability_experiment[i, 1:2]),
                                            input_mn = input_mn, freq_mn = freq_mn)
    
    # sample cv
    cv_temp <- meta.arrR::calc_variability(input_temp)
    
    # create data.frame
    dplyr::bind_cols(amplitude = variability_experiment[i, 1], phase = variability_experiment[i, 2], 
                     gamma = cv_temp$gamma)}) %>% 
  dplyr::group_by(amplitude, phase) %>% 
  dplyr::summarise(mean = mean(gamma), sd = sd(gamma), .groups = "drop")

#### Create ggplot ####

gg_input_continuous <- ggplot(data = variability_df) + 
  geom_raster(aes(x = amplitude, y = phase, fill = mean)) + 
  geom_point(aes(x = amplitude, y = phase, size = sd), pch = 1) +
  scale_fill_continuous(name = "gamma", type = "viridis") +
  scale_size_continuous(name = "SD") + 
  guides(fill = guide_colorbar(order = 1), size = guide_legend(order = 0)) + 
  labs(x = "Variability Amplitude", y = "Variability Phase") +
  coord_equal() + 
  theme_classic() + theme(legend.position = "bottom", legend.box = "vertical", 
                          legend.key.width = unit(1.5, 'cm'))

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_input_continuous, filename = "gg_example-gg_input_continuous.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = FALSE)
