##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

source("05_Various/setup.R")
source("01_Functions/fit_nls.R")
source("01_Functions/predict_nls.R")

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

# set number of repetitions
itr <- 50

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

#### Variability for increasing scale ####

# sample variability for different treatment lvls and increasing scale
variability_input <- purrr::map_dfr(1:nrow(variability_experiment), function(i) {
  
    # print progress
    message("\r> Progress: ", i, " / ", nrow(variability_experiment), "\t\t", appendLF = FALSE)
    
    # simulate input
    input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                            variability = as.numeric(variability_experiment[i, 1:2]),
                                            input_mn = input_mn, freq_mn = freq_mn)
    
    # sample cv
    cv_temp <- meta.arrR::sample_variability(input_temp)
    
    # create data.frame
    dplyr::bind_cols(input = variability_experiment[i, 3], cv_temp)}) %>%
  dplyr::mutate(input = factor(input, levels = unique(variability_experiment$combined_label), 
                               labels = unique(paste0(variability_experiment$amplitude_label, " Amplitude // ",
                                                      variability_experiment$phase_label, " Phase"))), 
                n = factor(n, ordered = TRUE)) %>% 
  tidyr::pivot_longer(-c(input, n), names_to = "stat") %>% 
  dplyr::mutate(stat = factor(stat, levels = c("alpha", "beta", "gamma", "synchrony")))

#### Fit NLS decay model ####

# create vector with names
variability_input_names <- paste(rep(x = levels(variability_input$input), each = 2), 
                              c("gamma", "synchrony")) %>% 
  stringr::str_remove_all(pattern = "// ") %>%
  stringr::str_replace_all(c("Amplitude" = "-", "Phase " = "", " - " = "-", " " = "_")) %>% 
  stringr::str_to_lower()

# fit NLS model
variability_input_fit <- dplyr::filter(variability_input, stat %in% c("gamma", "synchrony")) %>% 
  dplyr::mutate(n = as.numeric(n)) %>% 
  dplyr::group_by(input, stat) %>% 
  dplyr::group_split() %>% 
  purrr::set_names(variability_input_names) %>% 
  purrr::map(fit_nls)

# predict data
variability_input_pred <- purrr::map_dfr(variability_input_fit, predict_nls, x = 1:9, .id = "input") %>% 
  tidyr::separate(col = input, sep = "_", into = c("input", "stat")) %>% 
  dplyr::mutate(input = factor(input, levels = c("low-low", "medium-low", "high-low", 
                                                 "low-medium", "medium-medium", "high-medium",
                                                 "low-high", "medium-high", "high-high"), 
                               labels = levels(variability_input$input))) 

#### Table with coefficients ####

variability_input_coef <- dplyr::bind_rows(variability_input_fit, .id = "input") %>% 
  tidyr::separate(col = input, sep = "_", into = c("input", "stat")) %>% 
  dplyr::mutate(input = factor(input, levels = c("low-low", "medium-low", "high-low", 
                                                 "low-medium", "medium-medium", "high-medium",
                                                 "low-high", "medium-high", "high-high"), 
                               labels = levels(variability_input$input)),
                model = factor(dplyr::case_when(is.na(p.value) ~ "const", 
                                                !is.na(p.value) ~ "nls"), 
                               levels = c("const", "nls")))

dplyr::filter(variability_input_coef, p.value < 0.05)

#### Create ggplot ####

size_text <- 2.75

x_text <- 8
y_text_a <- 0.9
y_text_b <- 0.8

digits_text <- 5

gg_variability_input <- dplyr::filter(variability_input, stat %in% c("gamma", "synchrony")) %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2, col = "grey") +
  geom_boxplot(aes(x = n, y = value, fill = stat), alpha = 0.25) +
  geom_line(data = variability_input_pred, aes(x = x, y = value, col = stat)) +
  geom_label(data = dplyr::filter(variability_input_coef, 
                                  stat == "gamma", term %in% c("n", "beta")), 
            aes(x = x_text, y = y_text_a, col = "gamma", 
                label = paste(term, ":", round(estimate, digits = digits_text))), 
            size = size_text, show.legend = FALSE) + 
  geom_label(data = dplyr::filter(variability_input_coef, 
                                  stat == "synchrony", term %in% c("n", "beta")), 
            aes(x = x_text, y = y_text_b, col = "synchrony", 
                label = paste(term, ":", round(estimate, digits = digits_text))), 
            size = size_text, show.legend = FALSE) + 
  facet_wrap(. ~ input, nrow = 3, ncol = 3) +
  scale_fill_manual(name = "", values = c("#ED5B66", "#60BAE4")) +
  scale_color_manual(name = "", values = c("#ED5B66", "#60BAE4")) +
  labs(x = "# local ecosystems", y = "Variability value") +
  theme_classic() +
  theme(legend.position = "bottom")

#### Save ggplot ####

overwrite <- FALSE

suppoRt::save_rds(object = variability_input, filename = "variability_example.rds", 
                  path = "02_Data/", overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_variability_input, filename = "gg_example-input_variability.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
