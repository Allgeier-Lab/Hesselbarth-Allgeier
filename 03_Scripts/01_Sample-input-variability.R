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
n <- 7

# set min_per_i
min_per_i <- 120

# run the model for n years
years <- 25

max_i <- (60 * 24 * 365 * years) / min_per_i

# setup nutrient input to be maximum 10% of inital nutrients
input_mn <- meta.arrR::meta.arrR_starting_values$nutrients_pool * 0.1

# set frequency of input peaks
freq_mn <- 5

# set low and high variability treatment levels
variability_lo <- 0.0

variability_hi <- 0.5

# combine full grid
variability_full <- expand.grid(ampl = c(variability_lo, variability_hi), 
                                freq = c(variability_lo, variability_hi)) %>% 
  dplyr::bind_cols(label = factor(x = c("lo_lo", "hi_lo", "lo_hi", "hi_hi"), 
                                  levels = c("lo_lo", "hi_lo", "lo_hi", "hi_hi"), 
                                  labels = c("Low amplitude - Low frequency",
                                             "High amplitude - Low frequency",
                                             "Low amplitude - High frequency",
                                             "High amplitude - High frequency")))

# set number of repetitions
itr <- 50

#### Variability for increasing scale ####

# sample variability for different treatment lvls and increasing scale
input_sampled <- purrr::map_dfr(1:nrow(variability_full), function(i) {
  
  # repeat for itr repetitions
  purrr::map_dfr(1:itr, function(j) {
    
    # print progress
    message("\r> Variability: ", i, "/", nrow(variability_full), 
            " || Iterations: ", j, "/", itr, "\t\t", appendLF = FALSE)
    
    # simulate input
    input_temp <- sim_nutr_input(n = n, max_i = max_i,
                                 variability = as.numeric(variability_full[i, 1:2]),
                                 input_mn = input_mn, freq_mn = freq_mn)
    
    # sample cv
    cv_temp <- sample_cv_gamma(input_temp)
    
    # create data.frame
    dplyr::bind_cols(label = variability_full[i, 3], itr = j, cv_temp)})
  
  }) %>%
  dplyr::mutate(n = factor(n, ordered = TRUE)) %>%
  tidyr::pivot_longer(-c(label, itr, n), names_to = "stat")

#### Create ggplot ####

gg_input_sampled <- dplyr::filter(input_sampled, stat %in% c("gamma", "synchrony")) %>%
  ggplot() +
  geom_boxplot(aes(x = n, y = value), fill = "grey") +
  geom_hline(yintercept = 0, linetype = 2, col = "grey") +
  facet_wrap(. ~ stat + label,  ncol = 4, nrow = 2, scales = "fixed") +
  labs(x = "# local ecosystems", y = "Variability value") +
  theme_classic() +
  theme(legend.position = "bottom")

#### Save ggplot ####

overwrite <- FALSE

suppoRt::save_rds(object = input_sampled, filename = "sampled_cv_input.rds", 
                  path = "02_Data/", overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_input_sampled, filename = "gg_input_sampled.pdf", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
