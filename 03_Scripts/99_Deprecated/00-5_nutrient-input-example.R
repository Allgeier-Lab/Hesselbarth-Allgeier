##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

years <- 50

max_i <- (60 * 24 * 365 * years) / min_per_i

##### Plot different input variability  #####

# gg_list <- purrr::map(seq(from = 0, to = 1.0, length.out = 4), function(i ) {
#   
#   meta.arrR::simulate_nutr_input(n = n, max_i = max_i, frequency = freq_mn, 
#                                  input_mn = nutrient_input, 
#                                  amplitude_mn = 0.5, amplitude_sd = i) %>% 
#     plot() + 
#     ggplot2::geom_hline(yintercept = 2 * nutrient_input, linetype = 2, color = "grey") +
#     ggplot2::geom_hline(yintercept = nutrient_input, linetype = 2, color = "grey") +
#     ggplot2::scale_y_continuous(limits = c(0, 2 * nutrient_input)) +
#     ggplot2::labs(subtitle = paste0("amplitude_sd=", i)) + 
#     ggplot2::guides(color = "none")
#   
# })
# 
# cowplot::plot_grid(plotlist = gg_list, nrow = 2, ncol = 2)

#### CV ####

amplitude_lvl <- c(0.05, 0.5, 0.95)

cv_input <- purrr::map_dfr(1:length(amplitude_lvl), function(i) {
  
  purrr::map_dfr(1:iterations, function(j) {
      
    message("\r> Progress: i=", i, "/", length(amplitude_lvl), "; j=", j , "/", 
            iterations, "\t\t", appendLF = FALSE)
      
    sd_temp <- runif(n = 1, min = 0.0, max = 1.0)
    
    # sd_temp <- sample(x = seq(from = 0, to = 1, length.out = 5), size = 1)
    # sd_temp <- sample(x = c(0.05, 0.5, 0.95), size = 1)
    
    df_amp <- meta.arrR::simulate_nutr_input(n = n, max_i = max_i, frequency = freq_mn, 
                                             input_mn = nutrient_input, 
                                             amplitude_mn = amplitude_lvl[i], amplitude_sd = sd_temp, 
                                             phase_mn = 0.0, phase_sd = 0.0) %>%
      meta.arrR::calc_variability() %>% 
      dplyr::mutate(amplitude_mn = amplitude_lvl[i], sd = sd_temp, stochastic = "amplitude")
    
    df_pha <- meta.arrR::simulate_nutr_input(n = n, max_i = max_i, frequency = freq_mn, 
                                             input_mn = nutrient_input, 
                                             amplitude_mn = amplitude_lvl[i], amplitude_sd = 0.0, 
                                             phase_mn = 0.0, phase_sd = sd_temp) %>%
      meta.arrR::calc_variability() %>% 
      dplyr::mutate(amplitude_mn = amplitude_lvl[i], sd = sd_temp, stochastic = "phase")
    
    df_sim <- meta.arrR::simulate_nutr_input(n = n, max_i = max_i, frequency = freq_mn, 
                                             input_mn = nutrient_input, 
                                             amplitude_mn = amplitude_lvl[i], amplitude_sd = sd_temp, 
                                             phase_mn = 0.0, phase_sd = sd_temp) %>%
      meta.arrR::calc_variability() %>% 
      dplyr::mutate(amplitude_mn = amplitude_lvl[i], sd = sd_temp, stochastic = "both")
    
    dplyr::bind_rows(df_amp, df_pha, df_sim)}, .id = "itr")}) %>%
  dplyr::select(-part) %>% 
  dplyr::mutate(itr = as.numeric(itr), measure = as.factor(measure), amplitude_mn = factor(amplitude_mn),
                stochastic = factor(stochastic, levels = c("amplitude", "phase", "both"))) %>%
  tibble::as_tibble()

dplyr::filter(cv_input, measure %in% c("alpha", "gamma", "synchrony")) %>% 
ggplot(aes(x = sd, y = value, color = measure)) + 
  geom_hline(yintercept = 0.0, linetype = 2, col = "grey") +
  geom_point(alpha = 0.25) + 
  # geom_smooth(method = "lm", se = FALSE, size = 0.65) +
  geom_smooth(method = "loess", se = FALSE, size = 0.65) +
  facet_grid(rows = vars(stochastic), cols = vars(amplitude_mn), labeller = label_both) +
  scale_x_continuous(breaks = seq(from = 0, to = 1.0, length.out = 5), limits = c(0, 1)) +
  scale_color_manual(name = "", values = c("#007AA1", "#DF4E25", "#41B282")) +
  labs(x = "Variability nutrients input (sd)", y = "CV nutrients input") + 
  theme_classic() + theme(legend.position = "top", 
                          panel.border = element_rect(fill = NA, size = 0.5),
                          strip.background = element_rect(fill = NA, size = 0.5))

# dplyr::filter(cv_input, measure == "beta") %>%
#   ggplot(aes(x = sd, y = log(value))) +
#   geom_hline(yintercept = 0.0, linetype = 2, col = "grey") +
#   geom_point(alpha = 0.5) +
#   # geom_smooth(method = "lm", se = FALSE, size = 0.5, col = "black") +
#   geom_smooth(method = "loess", se = FALSE, size = 0.5, col = "black") +
#   facet_grid(rows = vars(stochastic), cols = vars(amplitude_mn), labeller = label_both) +
#   scale_x_continuous(breaks = seq(from = 0, to = 1.0, length.out = 5), limits = c(0, 1)) +
#   labs(x = "Variability nutrients input (sd)", y = "PE nutrients input") +
#   theme_classic() + theme(legend.position = "bottom",
#                           panel.border = element_rect(fill = NA, size = 0.5),
#                           strip.background = element_rect(fill = NA, size = 0.5))
