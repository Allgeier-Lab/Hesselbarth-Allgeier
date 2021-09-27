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

#### Load data ####

# variability_experiment <- readRDS("02_Data/variability_experiment.rds")

variability_nofish <- readRDS("02_Data/variability_nofish.rds")

#### Setup globals ####

size_text <- 2.75

x_text <- 8
y_text_a <- 0.9
y_text_b <- 0.8

digits_text <- 5

#### Preprocess data ####

parts <- c("bg_biomass", "ag_biomass", "bg_production", "ag_production")

# create vector with names
variability_names_input <- paste0(rep(x = unique(variability_nofish$combined_label), 
                                      each = 2), c("-gamma", "-synchrony")) %>% 
  stringr::str_to_lower()

variability_names_output <- paste0(rep(x = variability_names_input, each = 4), "-", 
                                   parts) %>% 
  stringr::str_to_lower()

# create labels for nice figure
variability_labels <- paste0(variability_nofish$amplitude_label, " Amplitude // ", 
                             variability_nofish$phase_label, " Phase") %>% 
  unique()

# convert combined as factor
variability_nofish <- dplyr::mutate(variability_nofish, 
                                    part = factor(part, levels = parts), 
                                    combined_label = factor(combined_label,
                                                            levels = unique(variability_nofish$combined_label), 
                                                            labels = variability_labels), 
                                    n = factor(n))

#### Input variability ####

# filter only input data, get required columns to reshape to long
variability_input <- dplyr::filter(variability_nofish, what == "input") %>% 
  dplyr::select(-c(input, what, part, alpha, beta, amplitude, phase)) %>% 
  tidyr::pivot_longer(cols = c(gamma, synchrony), names_to = "stat") %>% 
  dplyr::mutate(stat = factor(stat, levels = c("gamma", "synchrony")))

#### Fit decay model #####

# fit NLS model
variability_input_fit <- dplyr::mutate(variability_input, n = as.numeric(n)) %>% 
  dplyr::group_by(combined_label, stat) %>% 
  dplyr::group_split() %>% 
  purrr::set_names(variability_names_input) %>%
  purrr::map(fit_nls)

# predict data
variability_input_pred <- purrr::map_dfr(variability_input_fit, predict_nls, x = 1:9, 
                                         .id = "combined_label") %>% 
  tidyr::separate(col = combined_label, sep = "-", into = c("combined_label", "stat")) %>% 
  dplyr::mutate(combined_label = factor(combined_label, levels = c("low_low", "medium_low", "high_low", 
                                                                   "low_medium", "medium_medium", "high_medium",
                                                                   "low_high", "medium_high", "high_high"), 
                                        labels = variability_labels)) 

# get coefficients of model
variability_input_coef <- dplyr::bind_rows(variability_input_fit, .id = "combined_label") %>% 
  tidyr::separate(col = combined_label, sep = "-", into = c("combined_label", "stat")) %>% 
  dplyr::mutate(combined_label = factor(combined_label, levels = c("low_low", "medium_low", "high_low", 
                                                                   "low_medium", "medium_medium", "high_medium",
                                                                   "low_high", "medium_high", "high_high"), 
                                        labels = variability_labels))

#### Create ggplot ####

gg_variability_input <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, col = "grey") +
  geom_boxplot(data = variability_input, aes(x = n, y = value, fill = stat), alpha = 0.25) +
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
  facet_wrap(. ~ combined_label, nrow = 3, ncol = 3) +
  scale_fill_manual(name = "", values = c("#ED5B66", "#60BAE4")) +
  scale_color_manual(name = "", values = c("#ED5B66", "#60BAE4")) +
  labs(x = "# local ecosystems", y = "Variability value", title = "Nutrient input") +
  theme_classic() +
  theme(legend.position = "bottom")

#### Output variability ####

# filter only input data, get required columns to reshape to long
variability_output <- dplyr::filter(variability_nofish, what == "output") %>% 
  dplyr::select(-c(input, what, alpha, beta, amplitude, phase)) %>% 
  tidyr::pivot_longer(cols = c(gamma, synchrony), names_to = "stat") %>% 
  dplyr::mutate(stat = factor(stat, levels = c("gamma", "synchrony")))

#### Fit decay model #####

# fit NLS model
variability_output_fit <- dplyr::mutate(variability_output, n = as.numeric(n)) %>% 
  dplyr::group_by(combined_label, stat, part,) %>% 
  dplyr::group_split() %>%
  purrr::set_names(variability_names_output) %>%
  purrr::map(fit_nls)

# predict data
variability_output_pred <- purrr::map_dfr(variability_output_fit, predict_nls, x = 1:9, 
                                          .id = "combined_label") %>% 
  tidyr::separate(col = combined_label, sep = "-", into = c("combined_label", "stat", "part")) %>% 
  dplyr::mutate(combined_label = factor(combined_label, levels = c("low_low", "medium_low", "high_low", 
                                                                   "low_medium", "medium_medium", "high_medium",
                                                                   "low_high", "medium_high", "high_high"), 
                                        labels = variability_labels))

# get coefficients of model
variability_output_coef <- dplyr::bind_rows(variability_output_fit, .id = "combined_label") %>% 
  tidyr::separate(col = combined_label, sep = "-", into = c("combined_label", "stat", "part")) %>% 
  dplyr::mutate(combined_label = factor(combined_label, levels = c("low_low", "medium_low", "high_low", 
                                                                   "low_medium", "medium_medium", "high_medium",
                                                                   "low_high", "medium_high", "high_high"), 
                                        labels = variability_labels))

#### Create ggplot ####

gg_variability_output <- purrr::map(parts, function(i) {
  
  # filter data by part
  output_temp <- dplyr::filter(variability_output, part == i)
  
  pred_temp <- dplyr::filter(variability_output_pred, part == i)
  
  coef_temp <- dplyr::filter(variability_output_coef, part == i)
  
  ggplot() +
    geom_hline(yintercept = 0, linetype = 2, col = "grey") +
    geom_boxplot(data = output_temp, aes(x = n, y = value, fill = stat), alpha = 0.25) +
    geom_line(data = pred_temp, aes(x = x, y = value, col = stat)) +
    geom_label(data = dplyr::filter(coef_temp, stat == "gamma", term %in% c("n", "beta")),
               aes(x = x_text, y = y_text_a, col = "gamma",
                   label = paste(term, ":", round(estimate, digits = digits_text))),
               size = size_text, show.legend = FALSE) +
    geom_label(data = dplyr::filter(coef_temp, stat == "synchrony", term %in% c("n", "beta")),
               aes(x = x_text, y = y_text_b, col = "synchrony",
                   label = paste(term, ":", round(estimate, digits = digits_text))),
               size = size_text, show.legend = FALSE) +
    facet_wrap(. ~ combined_label, nrow = 3, ncol = 3) +
    scale_fill_manual(name = "", values = c("#ED5B66", "#60BAE4")) +
    scale_color_manual(name = "", values = c("#ED5B66", "#60BAE4")) +
    labs(x = "# local ecosystems", y = "Variability value", title = i) +
    theme_classic() +
    theme(legend.position = "bottom")
})
  
#### Save ggplot ####

overwrite <- FALSE

suppoRt::save_ggplot(plot = gg_variability_input, filename = "gg_input_variability.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

# save all output figures
purrr::walk(seq_along(parts), function(i) {
  
  part_temp <- stringr::str_replace(string = parts[i], pattern = "_", replacement = "-")
  
  name_temp <- paste0("gg_output_variability_", part_temp, ".png")
  
  suppoRt::save_ggplot(plot = gg_variability_output[[i]], filename = name_temp, 
                       path = "04_Figures/", width = height, height = width, dpi = dpi, 
                       units = units, overwrite = overwrite)
})
