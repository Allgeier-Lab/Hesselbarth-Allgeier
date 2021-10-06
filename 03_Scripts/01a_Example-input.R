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

freq_mn <- years * 1/4

# create variability data.frame with all combinations 
variability_input <- expand.grid(amplitude = c(0, 0.5, 1), 
                                phase = c(0, 0.5, 1)) %>% 
  dplyr::mutate(amplitude_label = dplyr::case_when(amplitude == 0 ~ "Low", 
                                                   amplitude == 0.5 ~ "Medium", 
                                                   TRUE ~ "High"), 
                phase_label = dplyr::case_when(phase == 0 ~ "Low", 
                                               phase == 0.5 ~ "Medium", 
                                               TRUE ~ "High")) %>% 
  tidyr::unite("combined_label", amplitude_label, phase_label, sep = "_", remove = FALSE)

#### Simulate input ####

itr <- 50

# simulate all treatment levels
alpha_df <- purrr::map(1:nrow(variability_input), function(i) {
  
  message("\r> Variability: ", i, "/", nrow(variability_input), "\t\t", 
          appendLF = FALSE)
  
  meta.arrR::sim_nutr_input(n = n, max_i = max_i, 
                            variability = as.numeric(variability_input[i, 1:2]),
                            input_mn = input_mn, freq_mn = freq_mn) %>% 
    meta.arrR::get_input_df(gamma = FALSE, long = TRUE)}) %>% 
  purrr::set_names(variability_input$combined_label) %>% 
  dplyr::bind_rows(.id = "Input") %>% 
  dplyr::mutate(Input = factor(Input, levels = variability_input$combined_label, 
                               labels = paste0(variability_input$amplitude_label, " Amplitude // ", 
                                               variability_input$phase_label, " Phase")))

# simulate all treatment levels
gamma_df <- purrr::map(1:nrow(variability_input), function(i) {
  
  purrr::map_dfr(1:itr, function(j) {
    
    message("\r> Variability: ", i, "/", nrow(variability_input), 
            " --- Iteration: ", j, "/", itr, "\t\t", 
            appendLF = FALSE)
    
    meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                              variability = as.numeric(variability_input[i, 1:2]),
                              input_mn = input_mn, freq_mn = freq_mn) %>% 
      meta.arrR::get_input_df() %>% 
      dplyr::select(Timestep, Gamma) %>% 
      dplyr::bind_cols(itr = j, .)}) %>% 
    dplyr::group_by(Timestep) %>% 
    dplyr::summarise(Mean = mean(Gamma), SD = sd(Gamma))}) %>% 
  purrr::set_names(variability_input$combined_label) %>% 
  dplyr::bind_rows(.id = "Input") %>% 
  dplyr::mutate(Input = factor(Input, levels = variability_input$combined_label, 
                               labels = paste0(variability_input$amplitude_label, " Amplitude // ", 
                                               variability_input$phase_label, " Phase")))

#### Save data ####

overwrite <- FALSE

suppoRt::save_rds(object = alpha_df, filename = "example_alpha.rds", 
                  path = "02_Data/", overwrite = overwrite)

suppoRt::save_rds(object = gamma_df, filename = "example_gamma.rds", 
                  path = "02_Data/", overwrite = overwrite)

#### Load data ####

alpha_df <- readRDS("02_Data/example_alpha.rds")

gamma_df <- readRDS("02_Data/example_gamma.rds")

#### Create ggplot ####

gg_input_local <- ggplot(data = alpha_df) +
  geom_line(aes(x = Timestep, y = Value, col = Meta)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  facet_wrap(. ~ Input) +
  scale_color_viridis_d(name = "Meta-ecosystems") +
  labs(x = "Time step", y = "Nutrient input") +
  guides(color = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")

gg_input_regional <- ggplot(data = gamma_df) +
  geom_ribbon(aes(x = Timestep, ymin = Mean - SD, ymax = Mean + SD, fill  = "Gamma"), 
              alpha = 0.25) + 
  geom_line(aes(x = Timestep, y = Mean, col = "Gamma")) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  facet_wrap(. ~ Input) +
  scale_fill_manual(name = "", values = "black") +
  scale_color_manual(name = "Meta-ecosystems (mean±sd)", values = "black") +
  guides(fill = "none") +
  labs(x = "Time step", y = "Nutrient input") +
  guides(color = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_input_local, filename = "gg_example-input_local.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_input_regional, filename = "gg_example-input_regional.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
