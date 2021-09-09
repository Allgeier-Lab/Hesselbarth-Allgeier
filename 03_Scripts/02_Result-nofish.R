##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

source("05_Various/setup.R")

#### Load data ####

sampled_cv_nofish <- readRDS("02_Data/sampled_cv_nofish.rds")

#### Pre-process data ####

output_sampled <- tidyr::separate(sampled_cv_nofish, part, into = c("part", "measure"))

#### Create ggplot ####

gg_output_sampled_bg <- dplyr::filter(output_sampled, stat %in% c("gamma", "synchrony"), 
                                      part == "bg") %>%
  ggplot() +
  geom_boxplot(aes(x = n, y = value, col = measure), fill = "grey") +
  geom_hline(yintercept = 0, linetype = 2, col = "grey") +
  facet_wrap(. ~ stat + label,  ncol = 4, nrow = 2, scales = "fixed") +
  scale_color_manual(name = "", values = c("#EB2C88", "#3D98D3")) +
  labs(x = "# local ecosystems", y = "Variability value") +
  theme_classic() +
  theme(legend.position = "bottom")

gg_output_sampled_ag <- dplyr::filter(output_sampled, stat %in% c("gamma", "synchrony"), 
                                   part == "ag") %>%
  ggplot() +
  geom_boxplot(aes(x = n, y = value, col = measure), fill = "grey") +
  geom_hline(yintercept = 0, linetype = 2, col = "grey") +
  facet_wrap(. ~ stat + label,  ncol = 4, nrow = 2, scales = "fixed") +
  scale_color_manual(name = "", values = c("#EB2C88", "#3D98D3")) +
  labs(x = "# local ecosystems", y = "Variability value") +
  theme_classic() +
  theme(legend.position = "bottom")

#### Save ggplot ####

overwrite <- FALSE

suppoRt::save_ggplot(plot = gg_output_sampled_bg, filename = "gg_output_sampled_bg.pdf", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_output_sampled_ag, filename = "gg_output_sampled_ag.pdf", 
                     path = "04_Figures/", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
