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

# run the model for n years
years <- 5
max_i <- (60 * 24 * 365 * years) / min_per_i

# set frequency of input peaks
freq_mn <- years * 1/2

#### Stable values ####

stable_values <- arrR::get_stable_values(bg_biomass = starting_list$bg_biomass,
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

# create variability data.frame with all combinations 
variability_input <- expand.grid(amplitude = c(0.05, 0.5, 1.0), 
                                 phase = c(0.05, 0.5, 1)) %>% 
  dplyr::mutate(amplitude_label = dplyr::case_when(amplitude == 0.05 ~ "Low", 
                                                   amplitude == 0.5 ~ "Medium", 
                                                   TRUE ~ "High"), 
                phase_label = dplyr::case_when(phase == 0.05 ~ "Low", 
                                               phase == 0.5 ~ "Medium", 
                                               TRUE ~ "High")) %>% 
  tidyr::unite("combined_label", amplitude_label, phase_label, sep = "_", remove = FALSE)

#### Simulate input ####

# simulate all treatment levels
alpha_df <- purrr::map(1:nrow(variability_input), function(i) {
  
  message("\r> Variability: ", i, "/", nrow(variability_input), "\t\t", 
          appendLF = FALSE)
  
  meta.arrR::sim_nutr_input(n = n, max_i = max_i, 
                            variability = as.numeric(variability_input[i, 1:2]),
                            input_mn = stable_values$nutrients_input, freq_mn = freq_mn) %>% 
    meta.arrR::get_input_df(gamma = FALSE, long = TRUE)}) %>% 
  purrr::set_names(variability_input$combined_label) %>% 
  dplyr::bind_rows(.id = "input") %>% 
  dplyr::mutate(input = factor(input, levels = variability_input$combined_label, 
                               labels = paste0(variability_input$amplitude_label, " Amplitude // ", 
                                               variability_input$phase_label, " Phase")))

# simulate all treatment levels
gamma_df <- purrr::map(1:nrow(variability_input), function(i) {
  
  purrr::map_dfr(1:iterations, function(j) {
    
    message("\r> Variability: ", i, "/", nrow(variability_input), 
            " --- Iteration: ", j, "/", iterations, "\t\t", appendLF = FALSE)
    
    meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                              variability = as.numeric(variability_input[i, 1:2]),
                              input_mn = stable_values$nutrients_input, freq_mn = freq_mn) %>% 
      meta.arrR::get_input_df() %>% 
      dplyr::select(timestep, gamma) %>% 
      dplyr::bind_cols(iterations = j, .)}) %>% 
    dplyr::group_by(timestep) %>% 
    dplyr::summarise(mean = mean(gamma), sd = sd(gamma))}) %>% 
  purrr::set_names(variability_input$combined_label) %>% 
  dplyr::bind_rows(.id = "input") %>% 
  dplyr::mutate(input = factor(input, levels = variability_input$combined_label, 
                               labels = paste0(variability_input$amplitude_label, " Amplitude // ", 
                                               variability_input$phase_label, " Phase")))

#### Create ggplot ####

gg_input_local <- ggplot(data = alpha_df) +
  geom_line(aes(x = timestep, y = value, col = meta)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  facet_wrap(. ~ input) +
  scale_color_viridis_d(name = "Meta-ecosystems") +
  labs(x = "Time step", y = "Nutrient input") +
  guides(color = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")

gg_input_regional <- ggplot(data = gamma_df) +
  geom_ribbon(aes(x = timestep, ymin = mean - sd, ymax = mean + sd, fill  = "Gamma"), 
              alpha = 0.25) + 
  geom_line(aes(x = timestep, y = mean, col = "Gamma")) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  facet_wrap(. ~ input) +
  scale_fill_manual(name = "", values = "black") +
  scale_color_manual(name = "Meta-ecosystems (mean±sd)", values = "black") +
  guides(fill = "none") +
  labs(x = "Time step", y = "Nutrient input") +
  guides(color = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")
