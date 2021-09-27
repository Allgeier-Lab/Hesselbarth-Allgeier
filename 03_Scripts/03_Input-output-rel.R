##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

source("05_Various/setup.R")

#### Load data ####

sampled_cv_input <- readRDS("02_Data/sampled_cv_input.rds")

sampled_cv_nofish <- readRDS("02_Data/sampled_cv_nofish.rds")

#### Pre-process data ####

# calculate mean and sd for each input and local ecosystem
sampled_cv_input_sum <- dplyr::group_by(sampled_cv_input, label, n, stat) %>% 
  dplyr::summarise(value_mn = mean(value), 
                   value_sd = sd(value)) %>% 
  dplyr::ungroup()

# calculate mean and sd for each input and local ecosystem
sampled_cv_nofish_sum <- dplyr::group_by(sampled_cv_nofish, label, part, n, stat) %>% 
  dplyr::summarise(value_mn = mean(value), 
                   value_sd = sd(value)) %>%
  dplyr::ungroup()

# left join for full data.frame
sampled_cv_full <- dplyr::left_join(x = sampled_cv_nofish_sum, y = sampled_cv_input_sum, 
                                    by = c("label", "n", "stat"), 
                                    suffix = c(".out", ".in")) %>% 
  dplyr::mutate(mean_ratio = value_mn.out / value_mn.in, 
                sd_ratio = value_sd.out / value_sd.in) %>% 
  dplyr::select(label, part, n, stat, mean_ratio, sd_ratio) %>% 
  tidyr::separate(part, into = c("part", "measure"))

#### Combine data ####

dplyr::filter(sampled_cv_full, 
              part == "bg", stat %in% c("gamma", "synchrony"), n == max(n)) %>% 
  tidyr::replace_na(list(sd_ratio = 0)) %>% 
  dplyr::select(-sd_ratio, -n) %>% 
  tidyr::pivot_wider(names_from = stat, values_from = mean_ratio)


dplyr::filter(sampled_cv_full, 
              part == "ag", stat %in% c("gamma", "synchrony"), n == max(n)) %>% 
  tidyr::replace_na(list(sd_ratio = 0)) %>% 
  dplyr::select(-sd_ratio, -n) %>% 
  tidyr::pivot_wider(names_from = stat, values_from = mean_ratio)

dplyr::filter(sampled_cv_full, part == "ag", stat %in% c("gamma", "synchrony")) %>% 
  ggplot() + 
  geom_point(aes(x = n , y = mean_ratio, col = measure)) +
  geom_linerange(aes(x = n , ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio, 
                     col = measure), alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = 2, col = "grey") + 
  # geom_hline(yintercept = 0, linetype = 2, col = "grey") + 
  scale_color_manual(name = "", values = c("#EB2C88", "#3D98D3")) +
  facet_wrap(. ~ stat + label, ncol = 4, nrow = 2, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "bottom")
