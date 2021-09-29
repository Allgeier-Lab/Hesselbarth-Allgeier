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
  dplyr::slice(rep(1:n(), each = itr))

#### HPC ####

#### Create function ####

foo <- function(amplitude, phase) {
  
  # simulate input
  input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                          variability = c(amplitude, phase),
                                          input_mn = input_mn, freq_mn = freq_mn)
  
  # sample cv
  cv_temp <- meta.arrR::sample_variability(input_temp)
  
  # combine to vector
  cbind(amplitude = amplitude, phase = phase, cv_temp)

}

#### Submit to HPC #### 

variability_sbatch <- rslurm::slurm_apply(f = foo, params = variability_experiment, 
                                          global_objects = c("n", "max_i", "input_mn", "freq_mn"),
                                          jobname = "exampl_vari",
                                          nodes = nrow(variability_experiment), cpus_per_node = 1, 
                                          slurm_options = list("account" = "jeallg1", 
                                                               "partition" = "standard",
                                                               "time" = "00:10:00", ## hh:mm::ss
                                                               "error" = "rslurm.log"),
                                          pkgs = "meta.arrR",
                                          rscript_path = rscript_path, sh_template = sh_template, 
                                          submit = FALSE)

#### Collect results ####

variability_result <- rslurm::get_slurm_out(variability_sbatch, outtype = "table")

rslurm::cleanup_files(variability_sbatch)

#### Save data ####

suppoRt::save_rds(object = variability_result, filename = "variability_example.rds", 
                  path = "02_Data/", overwrite = FALSE)

#### Load data ####

variability_result <- readRDS("02_Data/variability_example.rds")

variability_result <- dplyr::mutate(variability_result, 
                                    amplitude_label = dplyr::case_when(amplitude == 0 ~ "Low",
                                                                       amplitude == 0.5 ~ "Medium",
                                                                       TRUE ~ "High"),
                                    phase_label = dplyr::case_when(phase == 0 ~ "Low",
                                                                   phase == 0.5 ~ "Medium",
                                                                   TRUE ~ "High")) %>%
  tidyr::unite("input", amplitude_label, phase_label, sep = "_", remove = FALSE) %>% 
  dplyr::mutate(input = factor(input, levels = input, 
                               labels = paste0(amplitude_label, " Amplitude // ", 
                                               phase_label, " Phase")), 
                n = factor(n, ordered = TRUE)) %>% 
  tidyr::pivot_longer(c(alpha, beta, gamma, synchrony), names_to = "stat") %>% 
  dplyr::mutate(stat = factor(stat, levels = c("alpha", "beta", "gamma", "synchrony")))

#### Fit NLS decay model ####

# create vector with names
variability_input_names <- paste(rep(x = levels(variability_result$input), each = 2), 
                              c("gamma", "synchrony")) %>% 
  stringr::str_remove_all(pattern = "// ") %>%
  stringr::str_replace_all(c("Amplitude" = "-", "Phase " = "", " - " = "-", " " = "_")) %>% 
  stringr::str_to_lower()

# fit NLS model
variability_input_fit <- dplyr::filter(variability_result, stat %in% c("gamma", "synchrony")) %>% 
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
                               labels = levels(variability_result$input))) 

#### Table with coefficients ####

variability_input_coef <- dplyr::bind_rows(variability_input_fit, .id = "input") %>% 
  tidyr::separate(col = input, sep = "_", into = c("input", "stat")) %>% 
  dplyr::mutate(input = factor(input, levels = c("low-low", "medium-low", "high-low", 
                                                 "low-medium", "medium-medium", "high-medium",
                                                 "low-high", "medium-high", "high-high"), 
                               labels = levels(variability_result$input)),
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

gg_variability_input <- dplyr::filter(variability_result, stat %in% c("gamma", "synchrony")) %>%
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

suppoRt::save_ggplot(plot = gg_variability_input, filename = "gg_example-input_variability.png", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = FALSE)
