##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

#### Load data #### 

cv_both <- readr::read_rds("02_Data/01-all-variability.rds") %>% 
  purrr::map_dfr(function(i) i$cv) %>% 
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass", 
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")), 
                amplitude_mean = factor(amplitude_mean, labels = c("Low mean", "Medium mean", "High mean"), 
                                        ordered = TRUE),
                stochastic = factor(stochastic, levels = c("amplitude_move", "phase_move", "all"))) %>% 
  tibble::tibble()

#### Variance partition #### 

variance_partion <- dplyr::filter(cv_both, part %in% c("ag_production", "bg_production", "ttl_production"),
                                  measure %in% c("alpha", "gamma", "beta"), stochastic == "all") %>% 
  dplyr::group_by(part, measure, amplitude_mean) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(i) {
    
    part_temp <- vegan::varpart(i$value, ~ amplitude_sd, ~ phase_sd, ~ move_meta_sd, data = i)
    
    tibble::tibble(part = unique(i$part), measure = unique(i$measure), amplitude_mean = unique(i$amplitude_mean),
                   explanatory = c("amplidue", "phase", "move", "amplitude_phase", 
                                   "phase_move", "amplitude_move", "amplitude_phase_move", 
                                   "residuals"), 
                   value = part_temp$part$indfract$Adj.R.square,
                   x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0, 2),
                   y = c(1.2, -0.6, 0.5, -0.6, 0.5, -1, 0, -2.2)) %>% 
      tidyr::separate(explanatory, sep = "_", into = c("first", "second", "third"), 
                      fill = "right")

  }) %>% 
  dplyr::mutate(value_label = dplyr::case_when(value > 0 ~ paste0(round(value * 100, 2), "%"), 
                                               TRUE ~ "---"), 
                y = dplyr::case_when(part == "ag_production" ~ y + 0.25, 
                                     part == "ttl_production" ~ y - 0.25, 
                                     TRUE ~ y))

#### Setup ggplot ####

df_circle <- tibble::tibble(x = rep(x = c(0, 0.866, -0.866), times = 3),
                            y = rep(x = c(1, -0.5, -0.5), times = 3), 
                            explanatory = rep(x = c("Amplitude", "Phase", "Move"), times = 3))

# create facet labels
label_measure <- c(alpha = "Local scale", gamma = "Meta-ecosystem scale", beta = "Portfolio effect")

label_amplitude <- c("Low mean" = "Low amplitude mean", "Medium mean" = "Medium amplitude mean",
                     "High mean" = "High amplitude mean")

#### Create ggplot ####

gg_var_partion <- ggplot() + 
  geom_circle(data = df_circle, aes(x0 = x, y0 = y, r = 1.5, fill = explanatory),
              alpha = 0.1, color = "black", size = 0.25) +
  geom_text(data = variance_partion, aes(x = x, y = y, label = value_label, color = part),
            size = 2.75) +
  annotate(geom = "text", x = 0.75, y = -2.25, label = "Residuals") +
  annotate(geom = "text", x = 0, y = 2.75, angle = 0, color = "#a71c00", label = "Amplitude") +
  annotate(geom = "text", x = 2.75, y = 0, angle = 270, color = "#138cc1", label = "Phase") +
  annotate(geom = "text", x = -2.75, y = 0, angle = 90, color = "#f2af35", label = "Move meta") +
  scale_fill_manual(name = "", values = c("#a71c00", "#f2af35", "#138cc1")) +
  scale_color_manual(name = "", values = c("#d05c4a", "#62782c", "#003666")) +
  guides(fill = "none", color = guide_legend(nrow = 3)) +
  facet_grid(cols  = vars(amplitude_mean), rows = vars(measure), 
             labeller = labeller(amplitude_mean = label_amplitude, measure = label_measure)) + 
  coord_equal() +
  theme_bw() + 
  theme(legend.position = "left", axis.title = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(), 
        panel.background = element_blank())
