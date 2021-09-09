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
  dplyr::bind_cols(label = c("lo_lo", "hi_lo", "lo_hi", "hi_hi"))

#### Simulate input ####

# simulate all treatment levels
input_full <- purrr::map(1:nrow(variability_full), function(i) {
  
  sim_nutr_input(n = n, max_i = max_i, variability = as.numeric(variability_full[i, 1:2]),
                 input_mn = input_mn, freq_mn = freq_mn)
  
})

# name list
names(input_full) <- variability_full$label

### Flatten data ####

# create data.frame with all values
input_df <- purrr::map_dfr(input_full, function(i) {
  
  purrr::map_dfr(seq_along(i$values), function(j) {
    
    data.frame(id = j, timestep = 1:max_i, values = i$values[[j]])
    
    })
  }, .id = "input") %>%
  dplyr::mutate(input = factor(input, levels = c("lo_lo", "hi_lo", "lo_hi", "hi_hi"),
                               labels = c("Low amplitude - Low frequency",
                                          "High amplitude - Low frequency",
                                          "Low amplitude - High frequency",
                                          "High amplitude - High frequency")),
                id = factor(id))

# calculate gamma value
input_df_regional <- dplyr::group_by(input_df, input, timestep) %>%
  dplyr::summarise(values = sum(values))

#### Create ggplot ####

gg_input_local <- ggplot(data = input_df) +
  geom_line(aes(x = timestep, y = values, col = id)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  facet_wrap(. ~ input) +
  scale_color_viridis_d(name = "Meta-ecosystems") +
  labs(x = "Time step", y = "Nutrient input") +
  guides(color = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")

gg_input_regional <- ggplot(data = input_df_regional) +
  geom_line(aes(x = timestep, y = values, col = "Gamma")) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  facet_wrap(. ~ input) +
  scale_color_manual(name = "Meta-ecosystems", values = "black") +
  labs(x = "Time step", y = "Nutrient input") +
  guides(color = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")

#### Save ggplot ####

overwrite <- FALSE

suppoRt::save_ggplot(plot = gg_input_local, filename = "gg_input_local.pdf", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_input_regional, filename = "gg_input_regional.pdf", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
