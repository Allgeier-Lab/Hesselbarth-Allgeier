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

itr <- 5

# create variability data.frame with all combinations 
variability_experiment <- expand.grid(amplitude = c(0, 0.5, 1), 
                                      phase = c(0, 0.5, 1)) %>% 
  dplyr::mutate(amplitude_label = dplyr::case_when(amplitude == 0 ~ "Low", 
                                                   amplitude == 0.5 ~ "Medium", 
                                                   TRUE ~ "High"), 
                phase_label = dplyr::case_when(phase == 0 ~ "Low", 
                                               phase == 0.5 ~ "Medium", 
                                               TRUE ~ "High")) %>% 
  tidyr::unite("combined_label", amplitude_label, phase_label, sep = "_", remove = FALSE) %>% 
  dplyr::slice(rep(1:n(), each = itr))

#### Simulate input ####

# simulate all treatment levels
input_full <- purrr::map(1:nrow(variability_experiment), function(i) {
  
  meta.arrR::sim_nutr_input(n = n, max_i = max_i, 
                            variability = as.numeric(variability_experiment[i, 1:2]),
                            input_mn = input_mn, freq_mn = freq_mn)
  
})

# name list
names(input_full) <- variability_experiment$combined_label

### Flatten data ####

# create data.frame with all values
input_df <- purrr::map_dfr(input_full, meta.arrR::get_input_df, 
                           gamma = FALSE, long = TRUE, .id = "Input") %>% 
  dplyr::group_by(Input, Timestep, Meta) %>% 
  dplyr::summarise(mean = mean(Value), sd = sd(Value), 
                   .groups = "drop") %>% 
  dplyr::mutate(Input = factor(Input, levels = variability_full$combined_label, 
                               labels = paste0(variability_full$amplitude_label, " Amplitude // ", 
                                               variability_full$phase_label, " Phase")))

# calculate gamma value
input_df_gamma <- dplyr::group_by(input_df, Input, Timestep) %>%
  dplyr::summarise(Value = sum(Value))

#### Create ggplot ####

gg_input_local <- ggplot(data = input_df) +
  geom_line(aes(x = Timestep, y = mean, col = Meta)) +
  geom_ribbon(aes(x = Timestep, ymin = mean - sd, ymax = mean + sd, col = Meta)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  facet_wrap(. ~ Input) +
  coord_cartesian(xlim = c(0, max_i / 25 * 3)) +
  scale_color_viridis_d(name = "Meta-ecosystems") +
  labs(x = "Time step", y = "Nutrient input") +
  guides(color = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")

gg_input_regional <- ggplot(data = input_df_gamma) +
  geom_line(aes(x = Timestep, y = Value, col = "Gamma")) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  facet_wrap(. ~ Input) +
  coord_cartesian(xlim = c(0, max_i / 25 * 3)) +
  scale_color_manual(name = "Meta-ecosystems", values = "black") +
  labs(x = "Time step", y = "Nutrient input") +
  guides(color = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")

#### Save ggplot ####

overwrite <- FALSE

suppoRt::save_ggplot(plot = gg_input_local, filename = "gg_example-input_local.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_input_regional, filename = "gg_example-input_regional.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
