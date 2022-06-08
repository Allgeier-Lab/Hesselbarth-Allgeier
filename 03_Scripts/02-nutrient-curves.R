##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Figure of nutrient input curves

#### Load setup ####

source("05_Various/setup.R")

extension <- ".pdf"

# simulate nutrient input
df_input_treatments <- purrr::map_dfr(c(low = 0.5, high = 0.95), function(i) {
  
  meta.arrR::simulate_nutr_input(n = 1, max_i = max_i, frequency = years,
                                 input_mn = nutrient_input, amplitude_mn = i,
                                 verbose = FALSE) %>% 
    meta.arrR::filter_meta(filter = c(0, (max_i / years) * 5)) %>% 
    meta.arrR::get_input_df(gamma = FALSE)}, .id = "amplitude") %>% 
  dplyr::mutate(amplitude = factor(amplitude, levels = c("base", "low", "high"), 
                                   labels = c("Average", "Low", "High")))

gg_nutrient_intput <- ggplot(data = df_input_treatments, 
                             aes(x = timestep, y = meta_1, color = amplitude)) +
  geom_hline(yintercept = nutrient_input, color = "grey", linetype = 2) +
  geom_line() +
  annotate(geom = "text", x = (max_i / years) * 4.65, y = nutrient_input * 1.05,  color = "grey", 
           label = "Input equals excretion") +
  scale_color_manual(name = "Amplitude treatment", values = c("#e0d145", "#a62b29" ,"#d29040")) +
  scale_x_continuous(breaks = seq(from = 0, to = (max_i / years) * 5, by = (max_i / years)), 
                     labels = 0:5) +
  scale_y_continuous(breaks = c(nutrient_input * 0.05, nutrient_input * 0.95,
                                nutrient_input * 1.05, nutrient_input * 1.95), 
                     labels = c("5%", "95%%", "105%", "195%"), limits = c(0.0, nutrient_input * 2)) + 
  labs(x = "Years", y = "Nutrient input relative to average fish excretion") + 
  theme_classic(base_size = 10) + 
  theme(legend.position = c(0.9, 0.9), legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

suppoRt::save_ggplot(plot = gg_nutrient_intput, filename = paste0("02-years-input", extension),
                     path = "04_Figures/", width = width, height = height * 0.40,
                     units = units, dpi = dpi, overwrite = overwrite)
