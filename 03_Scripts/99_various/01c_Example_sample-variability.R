##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##


#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/fit_nls.R")
source("01_Functions/predict_nls.R")

#### Setup experiment ####

input_mn <- get_stable_values(starting_values = default_starting, 
                              parameters = default_parameters) %>% 
  magrittr::extract2("nutr_input")

# set number of repetitions
itr <- 50

# create variability data.frame with all combinations 
variability_experiment <- expand.grid(amplitude = c(0, 0.5, 1), 
                                phase = c(0, 0.5, 1)) %>% 
  dplyr::slice(rep(1:n(), each = itr))

#### Create HPC function ####

foo <- function(amplitude, phase) {
  
  # simulate input
  input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                          variability = c(amplitude, phase),
                                          input_mn = input_mn, freq_mn = freq_mn)
  
  # sample cv
  cv_temp <- meta.arrR::sample_variability(input_temp, itr = itr, verbose = FALSE)
  
  # combine to vector
  cbind(amplitude = amplitude, phase = phase, cv_temp)

}

#### Submit to HPC #### 

variability_sbatch <- rslurm::slurm_apply(f = foo, params = variability_experiment, 
                                          global_objects = c("n", "max_i", "input_mn", "freq_mn", "itr"),
                                          jobname = "example_sample_vari",
                                          nodes = nrow(variability_experiment), cpus_per_node = 1, 
                                          slurm_options = list("account" = account, 
                                                               "partition" = "standard",
                                                               "time" = "00:10:00"), ## hh:mm::ss
                                          pkgs = "meta.arrR",
                                          rscript_path = rscript_path, sh_template = sh_template, 
                                          submit = FALSE)

#### Collect results ####

variability_result <- rslurm::get_slurm_out(variability_sbatch, outtype = "table")

rslurm::cleanup_files(variability_sbatch)

#### Save data ####

suppoRt::save_rds(object = variability_result, filename = "01_example_sample-variability.rds", 
                  path = "02_Data/", overwrite = overwrite)

#### Load data ####

variability_result <- readRDS("02_Data/01_example_sample-variability.rds")

#### Pre-process data ####

variability_result <- dplyr::mutate(variability_result, 
                                    amplitude_label = dplyr::case_when(amplitude == 0 ~ "Low",
                                                                       amplitude == 0.5 ~ "Medium",
                                                                       TRUE ~ "High"),
                                    phase_label = dplyr::case_when(phase == 0 ~ "Low",
                                                                   phase == 0.5 ~ "Medium",
                                                                   TRUE ~ "High")) %>%
  tidyr::unite("input", amplitude_label, phase_label, sep = "_", remove = FALSE) %>%
  # tidyr::pivot_longer(c(alpha, beta, gamma, synchrony), names_to = "stat") %>% 
  dplyr::mutate(part = factor(part), 
                n = factor(n, ordered = TRUE),
                input = factor(input, levels = input, 
                               labels = paste0(amplitude_label, " Amplitude // ", 
                                               phase_label, " Phase")), 
                stat = factor(stat, levels = c("alpha", "beta", "gamma", "synchrony")))

#### Fit NLS decay model ####

# fit NLS model
variability_fit <- dplyr::filter(variability_result, stat %in% c("gamma", "synchrony")) %>% 
  dplyr::mutate(n = as.numeric(n)) %>% 
  dplyr::group_by(input, stat) %>% 
  dplyr::group_split() 

# create names
variability_names <- purrr::map_chr(variability_fit, function(i) {
  
  paste(unique(i$input), unique(i$stat), sep = "-") 
  
})

# set names and fit model
variability_fit <- purrr::set_names(variability_fit, variability_names) %>%
  purrr::map(fit_nls)

# predict data
variability_pred <- purrr::map_dfr(variability_fit, predict_nls, x = 1:9, .id = "input") %>% 
  tidyr::separate(col = input, sep = "-", into = c("input", "stat")) %>% 
  dplyr::mutate(input = factor(input, levels = unique(input))) 

#### Table with coefficients ####

variability_coef <- dplyr::bind_rows(variability_fit, .id = "input") %>% 
  tidyr::separate(col = input, sep = "-", into = c("input", "stat")) %>% 
  dplyr::mutate(input = factor(input, levels = unique(input)),
                model = factor(dplyr::case_when(is.na(p.value) ~ "const", 
                                                !is.na(p.value) ~ "nls"), 
                               levels = c("const", "nls")))

dplyr::filter(variability_coef, p.value < 0.05)

#### Create ggplot ####

size_text <- 2.75

x_text <- 8
y_text_a <- 0.9
y_text_b <- 0.8

digits_text <- 5

gg_variability <- dplyr::filter(variability_result, stat %in% c("gamma", "synchrony")) %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2, col = "grey") +
  geom_boxplot(aes(x = n, y = mean, fill = stat), alpha = 0.25) +
  geom_line(data = variability_pred, aes(x = x, y = value, col = stat)) +
  geom_label(data = dplyr::filter(variability_coef, 
                                  stat == "gamma", term %in% c("n", "beta")), 
            aes(x = x_text, y = y_text_a, col = "gamma", 
                label = paste(term, ":", round(estimate, digits = digits_text))), 
            size = size_text, show.legend = FALSE) + 
  geom_label(data = dplyr::filter(variability_coef, 
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

suppoRt::save_ggplot(plot = gg_variability, filename = "01_gg_example_sample-variability.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
