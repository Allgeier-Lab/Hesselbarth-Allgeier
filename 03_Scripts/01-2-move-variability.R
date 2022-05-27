##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/convert_label.R")

extension <- ".pdf"

#### Load/wrangle simulated data ####

file_names <- list.files(path = "02_Data/", pattern = "01-move-variability-*", 
                         full.names = TRUE) %>% 
  stringr::str_sort(numeric = TRUE) %>% 
  purrr::set_names(c("8", "16", "32", "64"))

df_cv_prod <- purrr::map_dfr(file_names, function(i) {
  readr::read_rds(i) %>% purrr::map_dfr(function(j) {
    dplyr::left_join(x = j$cv, y = j$production, by = c("part", "move_meta_mean", 
                                                        "move_meta_sd", "amplitude_mn"), 
                     suffix = c(".cv", ".prod")) %>% 
      dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))
    })
  }, .id = "pop_n") %>% 
  dplyr::mutate(pop_n = factor(as.numeric(pop_n), ordered = TRUE), 
                part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass", 
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")), 
                amplitude_mn = factor(amplitude_mn, ordered = TRUE, labels = c("Low", "Medium", "High")), 
                value.cv.log = log(value.cv), value.prod.log = log(value.prod), 
                value.move.log = log(value.move)) %>% 
  tibble::tibble()

df_dist_prod <- purrr::map_dfr(file_names, function(i) {
    readr::read_rds(i) %>% purrr::map_dfr(function(j) {
      dplyr::rename(j$dist, value.prod = value) %>% 
        dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))
    })
  }, .id = "pop_n") %>% 
  dplyr::mutate(pop_n = factor(as.numeric(pop_n), ordered = TRUE), 
                dist_class = factor(dist_class, ordered = TRUE),
                part = factor(part, levels = c("ag_production", "bg_production", "ttl_production",
                                               "nutrients_pool")),
                amplitude_mn = factor(amplitude_mn, ordered = TRUE, labels = c("Low", "Medium", "High")), 
                value.prod.log = log(value.prod), value.move.log = log(value.move)) %>%
  tibble::tibble()

#### Setup globals ####

parts <- c("ag_production", "bg_production", "ttl_production")

measures <- c("alpha", "gamma")

colors_ampl <- c("Low" = "#622567", "Medium" = "#007054", "High" = "#ec7122") 

colors_pop <- c("8" = "#e19ed5", "16" = "#7a92dc", "32" = "#ebcc60", "64" = "#cb3a54")

#### Total production ####

gg_pp_agbg <- dplyr::filter(df_cv_prod, part %in% c("ag_production", "bg_production"), 
                       measure == "alpha", value.move != 0) %>% 
  ggplot(aes(x = amplitude_mn, y = value.prod)) + 
  geom_jitter(aes(color = pop_n), alpha = 0.1, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5),
             shape = 19) +
  geom_boxplot(aes(fill = pop_n), position = position_dodge(width = 0.5), outlier.color = NA) +
  scale_color_manual(name = "Population size", values = colors_pop) +
  scale_fill_manual(name = "Population size", values = colors_pop) +
  facet_wrap(. ~ part, scales = "free_y", nrow = 2, ncol = 1,
             labeller = labeller(part = c(ag_production = "Aboveground", 
                                          bg_production = "Belowground ", 
                                          ttl_production = "Total"))) + 
  labs(x = "Amplitude treatment", y = "Total primary production") +
  theme_classic(base_size = 10) + 
  theme(legend.position = "none", strip.text = element_text(hjust = 0), 
        strip.background = element_blank(), panel.border = element_rect(size = 0.5, 
                                                                        fill = NA))

gg_pp_tt <- dplyr::filter(df_cv_prod, part == "ttl_production", 
                          measure == "alpha", value.move != 0) %>% 
  ggplot(aes(x = amplitude_mn, y = value.prod)) + 
  geom_jitter(aes(color = pop_n), alpha = 0.1, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5),
              shape = 19) +
  geom_boxplot(aes(fill = pop_n), position = position_dodge(width = 0.5), outlier.color = NA) +
  scale_color_manual(name = "Population size", values = c("8" = "#e19ed5", "16" = "#374a98", 
                                                          "32" = "#ebcc60", "64" = "#cb3a54")) +
  scale_fill_manual(name = "Population size", values = c("8" = "#e19ed5", "16" = "#374a98", 
                                                         "32" = "#ebcc60", "64" = "#cb3a54")) +
  facet_wrap(. ~ part, scales = "free_y", nrow = 1, ncol = 1,
             labeller = labeller(part = c(ag_production = "Aboveground", 
                                          bg_production = "Belowground ", 
                                          ttl_production = "Total"))) + 
  labs(x = "Amplitude treatment", y = "") +
  theme_classic(base_size = 10) + 
  theme(legend.position = "none", strip.text = element_text(hjust = 0), 
        strip.background = element_blank(), panel.border = element_rect(size = 0.5, 
                                                                        fill = NA))

gg_pp <- cowplot::plot_grid(gg_pp_agbg, gg_pp_tt, ncol = 2)

suppoRt::save_ggplot(plot = gg_pp, filename = paste0("00-amplitude-prod", extension), 
                     path = "04_Figures/", width = width, height = height / 2, 
                     units = units, dpi = dpi, overwrite = overwrite)

#### Mobility vs. CV (log-linear) ####

for (pop_i in c(8, 16, 32, 64)) {
  
 gg_parts <- purrr::map(parts, function(part_i) {
   
   gg_measures <- purrr::map(measures, function(measure_i) {
     
     gg_tile <- purrr::map(c("Low", "Medium", "High"), function(amp_i) {
       
       df_tile <- dplyr::filter(df_cv_prod, pop_n == pop_i, part == part_i, measure == measure_i, 
                                amplitude_mn == amp_i, value.move != 0)
       
       dplyr::group_by(df_tile, amplitude_mn, move_meta_mean, move_meta_sd) %>%
         dplyr::summarise(value.cv = mean(value.cv), .groups = "drop") %>%
         tibble::add_row(amplitude_mn = amp_i, move_meta_mean = 0.0, move_meta_sd = 0.0, value.cv = NA) %>%
         ggplot(aes(x = move_meta_mean, y = move_meta_sd, fill = value.cv)) +
         geom_tile() +
         scale_fill_viridis_c(option = "A", na.value = "grey") +
         scale_color_manual(values = c("black", "grey", "white")) +
         scale_x_continuous(breaks = c(0.0, 0.5, 1.0)) +
         scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
         labs(x = ifelse(amp_i == "Medium",
                         yes = expression(italic(move_meta_mean)), no = ""),
              y = ifelse(amp_i == "Low" && measure_i == "alpha",
                         yes = expression(italic(move_meta_sd)), no = "")) +
         labs(title = paste(amp_i, "amplitude")) +
         coord_equal() +
         theme_classic(base_size = 8) +
         theme(legend.position = "none", plot.title = element_text(size = 5.0),
               plot.margin = unit(c(-1.5, 2.5, -3.5, 2.5), "mm"), # t, r, b, l
               axis.line = element_blank(), panel.border = element_rect(size = 1, fill = NA), 
               strip.background = element_blank(), strip.text = element_text(hjust = 0))
       
     })
     
     gg_tile <- cowplot::plot_grid(plotlist = gg_tile, nrow = 1, ncol = 3)
    
     df_scatter <- dplyr::filter(df_cv_prod, pop_n == pop_i, part == part_i, 
                                 measure == measure_i, value.move != 0)
     
     df_regression <- dplyr::group_by(df_scatter, amplitude_mn) %>%
       dplyr::group_split() %>%
       purrr::map_dfr(function(df_temp) {

         lm_temp <- lm(value.cv.log ~ value.move, data = df_temp)
         
         pred <- exp(coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * df_temp$value.move)
         
         r.squared <- cor(df_temp$value.cv, pred)
         
         tibble::tibble(amplitude_mn = unique(df_temp$amplitude_mn), 
                        value.move = df_temp$value.move, value.cv.pred = pred, 
                        int = exp(lm_temp$coefficients[[1]]), coef = lm_temp$coefficients[[2]],
                        r.squared = r.squared, x = mean(df_temp$value.move), y = mean(df_temp$value.cv))
         
     })
     
     df_label <- dplyr::select(df_regression, -value.move, -value.cv.pred) %>% 
       dplyr::group_by(amplitude_mn) %>% 
       dplyr::distinct() %>% 
       dplyr::mutate(dir = dplyr::case_when(coef < 0.0 ~ "-", 
                                            coef > 0.0 ~ "+"),
                     int = convert_label(int), coef = convert_label(abs(coef)), 
                     r.squared = round(x = r.squared, digits = 2))

     gg_scatter <- ggplot(data = df_scatter, aes(x = value.move, y = value.cv, color = amplitude_mn)) +
       geom_hline(yintercept = 0, linetype = 2, color = "grey") +
       geom_point(shape = 19, alpha = 0.1) +
       geom_line(data = df_regression, aes(x = value.move, y = value.cv.pred, color = amplitude_mn)) +
       geom_label(data = df_label, parse = TRUE, size = 1.5, 
                  aes(x = x, y = y, color = amplitude_mn, 
                      label = paste0("italic(y)==", int, dir, coef, "*italic(x)*';'~R^2==", r.squared))) +
       scale_color_manual(name = "Amplitude", values = colors_ampl) +
       scale_x_continuous(breaks = seq(from = min(df_scatter$value.move), to = max(df_scatter$value.move),
                                       length.out = 5),
                          limits = c(min(df_scatter$value.move), max(df_scatter$value.move)),
                          labels = function(x) round(x, 1)) +
       scale_y_continuous(breaks = seq(from = 0.0, to = max(df_scatter$value.cv), length.out = 5),
                          limits = c(0.0, max(df_scatter$value.cv)),
                          labels = function(x) round(x, 3)) +
       labs(x = "Mean cross-ecosystem movement", 
            y = ifelse(measure_i == "alpha", yes = "Coefficient of variation", no = "")) +
       theme_classic(base_size = 8) + 
       theme(legend.position = "none", plot.margin = unit(c(-3.5, 2.5, 3.5, 2.5), "mm"), # t r b l
             strip.text = element_text(hjust = 0.0), strip.background = element_blank())
     
     cowplot::plot_grid(gg_tile, gg_scatter, nrow = 2, ncol = 1)
    
    })
   
   if (part_i == "ag_production") {
     
     cowplot::plot_grid(plotlist = gg_measures, nrow = 1, ncol = 2, 
                        labels = c("Local scale", "Meta scale"), hjust = c(-3, -3), 
                        label_fontface = "italic", label_size = 10)
     
   } else {
     
    cowplot::plot_grid(plotlist = gg_measures, nrow = 1, ncol = 2)
     
    }
  })
 
 gg_dummy <- ggplot(data = dplyr::filter(df_cv_prod, pop_n == 8),
                    aes(x = value.move, y = value.cv, color = amplitude_mn)) +
   geom_point(shape = 19) +
   scale_color_manual(name = "Amplitude", values = colors_ampl) +
   theme_classic(base_size = 8) + 
   theme(legend.position = "bottom")
 
 gg_legend <- cowplot::get_legend(gg_dummy)
 
 gg_mobility_cv <- cowplot::plot_grid(plotlist = gg_parts, nrow = 3, ncol = 1, 
                                      labels = c("a)", "b)", "c)"), 
                                      label_fontface = "italic", label_size = 10)
 
 gg_mobility_cv <- cowplot::plot_grid(gg_mobility_cv, gg_legend, ncol = 1, 
                                      rel_heights = c(0.95, 0.05))
 
 suppoRt::save_ggplot(plot = gg_mobility_cv, filename = paste0("01-mobility-cv-", pop_i, extension), 
                      path = "04_Figures/", height = height, width = width, 
                      units = units,  dpi = dpi, overwrite = overwrite)
  
}

#### CV vs. PP (log-log) ####

for (pop_i in c(8, 16, 32, 64)) {
  
  gg_parts <- purrr::map(parts, function(part_i) {
    
    gg_measures <- purrr::map(measures, function(measure_i) {
      
      df_scatter <- dplyr::filter(df_cv_prod, pop_n == pop_i, measure == measure_i, 
                                  part %in% part_i, value.move != 0)
      
      df_regression <- dplyr::group_by(df_scatter, amplitude_mn) %>%
        dplyr::group_split() %>%
        purrr::map_dfr(function(df_temp) {

          lm_temp <- lm(value.prod.log ~ value.cv.log, data = df_temp)
          
          pred <- exp(coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * df_temp$value.cv.log)
          
          r.squared <- cor(df_temp$value.prod, pred)
          
          tibble::tibble(amplitude_mn = unique(df_temp$amplitude_mn), 
                         value.cv = df_temp$value.cv, value.prod.pred = pred,
                         int = exp(lm_temp$coefficients[[1]]), coef = exp(lm_temp$coefficients[[2]]),
                         r.squared = r.squared)
          
      })
      
      if (part_i == "ag_production") {
        
        y_placement <- c(0.1, 0.2, 0.3)
        
      } else {
        
        y_placement <- c(0.7, 0.8, 0.9)
        
      }
      
      df_label <- dplyr::select(df_regression, -value.cv, -value.prod.pred) %>% 
        dplyr::group_by(amplitude_mn) %>% 
        dplyr::distinct() %>% 
        dplyr::mutate(dir = dplyr::case_when(coef < 0.0 ~ "-", coef > 0.0 ~ "+"),
                      int = convert_label(int), coef = convert_label(abs(coef)), 
                      r.squared = round(x = r.squared, digits = 2))
      
      df_label$x <- min(df_scatter$value.cv) + 
        (max(df_scatter$value.cv) - min(df_scatter$value.cv)) * 0.8
      
      df_label$y <- min(df_scatter$value.prod) + 
        (max(df_scatter$value.prod) - min(df_scatter$value.prod)) * y_placement
      
      if (measure_i == "alpha") {
        
        y_axis <- NULL
        
      } else {
        
        y_axis <- element_blank()
        
      }
      
      ggplot(data = df_scatter, aes(x = value.cv, y = value.prod, color = amplitude_mn)) +
        geom_point(shape = 19, alpha = 0.05) +
        geom_line(data = df_regression, aes(x = value.cv, y = value.prod.pred, color = amplitude_mn)) +
        geom_label(data = df_label, parse = TRUE, size = 2.0, 
                   aes(x = x, y = y, color = amplitude_mn, 
                       label = paste0("italic(y)==", int, dir, coef, "*italic(x)*';'~R^2==", r.squared))) +
        scale_color_manual(name = "Amplitude treatment", values = colors_ampl) +
        scale_x_continuous(breaks = seq(from = min(df_scatter$value.cv), to = max(df_scatter$value.cv),
                                        length.out = 5),
                           limits = c(min(df_scatter$value.cv), max(df_scatter$value.cv)),
                           labels = function(x) round(x, 3)) +
        scale_y_continuous(breaks = seq(from = min(df_scatter$value.prod), to = max(df_scatter$value.prod),
                                        length.out = 5),
                           limits = c(min(df_scatter$value.prod), max(df_scatter$value.prod)),
                           labels = function(x) round(x, 1)) +
        labs(x = "", y = "") +
        theme_classic(base_size = 10) + 
        theme(legend.position = "none", strip.text = element_text(hjust = 0.0), 
              strip.background = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
              axis.text.y = y_axis, plot.margin = unit(c(5.0, 2.5, 1.0, 2.5), "mm")) # t r b l
  
    })
    
    if (part_i == "ag_production") {
    
    cowplot::plot_grid(plotlist = gg_measures, ncol = 2, labels = c("Local scale", "Meta scale"),  
                       hjust = c(-4.0, -4.0), label_fontface = "italic", label_size = 10) 
    
    } else {
      
      cowplot::plot_grid(plotlist = gg_measures, ncol = 2)
      
    }
  })
  
  gg_dummy <- ggplot(data = df_scatter, aes(x = value.cv, y = value.prod, color = amplitude_mn)) +
    geom_point(shape = 19) + geom_line(size = 0.5) +
    scale_color_manual(name = "Amplitude treatment", values = colors_ampl) + 
    theme_classic(base_size = 10) + 
    theme(legend.position = "bottom")
  
  gg_legend <- cowplot::get_legend(gg_dummy)
  
  gg_scatter <- cowplot::plot_grid(plotlist = gg_parts, nrow = 3, 
                                   labels = c("a)", "b)", "c)"), label_fontface = "italic", 
                                   label_size = 10) + 
    draw_label("Coeffiecent of variation", x = 0.5, y = 0, vjust = -0.5, angle = 0, 
               size = 10) +
    draw_label("Total primary production", x = 0, y = 0.5, vjust = 1.25, angle = 90, 
               size = 10)
  
  gg_scatter <- plot_grid(gg_scatter, gg_legend, ncol = 1, rel_heights = c(0.95, 0.05))
  
  suppoRt::save_ggplot(plot = gg_scatter, filename = paste0("01-cv-pp-", pop_i, extension), 
                       path = "04_Figures/", width = height, height = width, units = units,  dpi = dpi, 
                       overwrite = overwrite)
  
}

#### PE vs. PP (log-log) ####

for (pop_i in c(8, 16, 32, 64)) {
  
  gg_parts <- purrr::map(parts, function(part_i) {
    
    gg_amplitude <- purrr::map(c("Low", "Medium", "High"), function(amplitude_i) {
      
      df_scatter <- dplyr::filter(df_cv_prod, pop_n == pop_i, measure == "beta", 
                                  amplitude_mn == amplitude_i, part %in% part_i, 
                                  value.move != 0)
      
      df_regression <- dplyr::group_by(df_scatter, amplitude_mn) %>%
        dplyr::group_split() %>%
        purrr::map_dfr(function(df_temp) {
          
          lm_temp <- lm(value.prod.log ~ value.cv.log, data = df_temp)
          
          pred <- exp(coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * df_temp$value.cv.log)
          
          r.squared <- cor(df_temp$value.prod, pred)
          
          tibble::tibble(amplitude_mn = unique(df_temp$amplitude_mn), 
                         value.cv = df_temp$value.cv, value.prod.pred = pred,
                         int = exp(lm_temp$coefficients[[1]]), coef = exp(lm_temp$coefficients[[2]]),
                         r.squared = r.squared)
          
      })
      
      y_placement <- ifelse(test = part_i == "ag_production", yes = 0.1, no = 0.9)
      
      df_label <- dplyr::select(df_regression, -value.cv, -value.prod.pred) %>% 
        dplyr::group_by(amplitude_mn) %>% 
        dplyr::distinct() %>% 
        dplyr::mutate(x = min(df_scatter$value.cv) + 
                        (max(df_scatter$value.cv) - min(df_scatter$value.cv)) * 0.8,
                      y = min(df_scatter$value.prod) + 
                        (max(df_scatter$value.prod) - min(df_scatter$value.prod)) * y_placement,
                      dir = dplyr::case_when(coef < 0.0 ~ "-", coef > 0.0 ~ "+"),
                      int = convert_label(int), coef = convert_label(abs(coef)), 
                      r.squared = round(x = r.squared, digits = 2))
      
      if (amplitude_i == "Low") {
        
        y_axis <- NULL
        
      } else {
        
        y_axis <- element_blank()
        
      }
      
      ggplot(data = df_scatter, aes(x = value.cv, y = value.prod, color = amplitude_mn)) +
        geom_point(shape = 19, alpha = 0.1) +
        geom_line(data = df_regression, aes(x = value.cv, y = value.prod.pred, color = amplitude_mn)) +
        geom_label(data = df_label, parse = TRUE, size = 1.75,
                   aes(x = x, y = y, color = amplitude_mn,
                       label = paste0("italic(y)==", int, dir, coef, "*italic(x)*';'~R^2==", r.squared))) +
        scale_color_manual(name = "Amplitude treatment", values = colors_ampl) +
        scale_x_continuous(breaks = seq(from = min(df_scatter$value.cv), to = max(df_scatter$value.cv),
                                        length.out = 5),
                           limits = c(min(df_scatter$value.cv), max(df_scatter$value.cv)),
                           labels = function(x) round(x, 2)) +
        scale_y_continuous(breaks = seq(from = min(df_scatter$value.prod), to = max(df_scatter$value.prod),
                                        length.out = 5),
                           limits = c(min(df_scatter$value.prod), max(df_scatter$value.prod)),
                           labels = function(x) round(x, 1)) +
        labs(x = "", y = "") +
        theme_classic(base_size = 10) + 
        theme(legend.position = "none", strip.text = element_text(hjust = 0.0), 
              strip.background = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
              strip.text.x = element_blank(), axis.text.y = y_axis, 
              plot.margin = unit(c(5.0, 2.5, 1.0, 2.5), "mm"))
  
    })
    
    cowplot::plot_grid(plotlist = gg_amplitude, ncol = 3)
    
  })
  
  gg_dummy <- ggplot(data = df_scatter, aes(x = value.cv, y = value.prod, color = amplitude_mn)) +
    geom_point(shape = 19) + geom_line(size = 0.5) +
    scale_color_manual(name = "Amplitude treatment", values = colors_ampl) + 
    theme_classic(base_size = 10) + 
    theme(legend.position = "bottom")
  
  gg_legend <- cowplot::get_legend(gg_dummy)
  
  gg_scatter <- cowplot::plot_grid(plotlist = gg_parts, nrow = 3, labels = c("a)", "b)", "c)"), 
                                   label_fontface = "italic", label_size = 10) + 
    draw_label("Portfolio effect", x = 0.5, y = 0, vjust = -1.0, angle = 0, 
               size = 10) +
    draw_label("Total primary production", x = 0, y = 0.5, vjust = 1.25, angle = 90, 
               size = 10)
  
  gg_scatter <- plot_grid(gg_scatter, gg_legend, nrow = 2, rel_heights = c(0.95, 0.05))
  
  suppoRt::save_ggplot(plot = gg_scatter, filename = paste0("01-pe-pp-", pop_i, extension), 
                       path = "04_Figures/", width = height, height = width, units = units,  dpi = dpi, 
                       overwrite = overwrite)
  
}

#### Dist vs. PP ####

gg_dist <- dplyr::group_by(df_dist_prod, pop_n, dist_class, part, amplitude_mn) %>% 
  dplyr::summarise(value.prod = mean(value.prod), .groups = "drop") %>%
  dplyr::filter(part != "nutrients_pool") %>%
  ggplot(aes(x = as.numeric(dist_class), y = value.prod, color = amplitude_mn, group = amplitude_mn)) + 
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  facet_grid(rows = dplyr::vars(part), cols = dplyr::vars(pop_n), scales = "free_y", 
             labeller = labeller(part = c(ag_production = "Aboveground production", 
                                          bg_production = "Belowground  production", 
                                          ttl_production = "Total production", 
                                          nutrients_pool = "Nutrients"), 
                                 pop_n = c("8" = "Local population n=8", 
                                           "16" = "Local population n=16", 
                                           "32" = "Local population n=32", 
                                           "64" = "Local population n=64"))) + 
  scale_color_manual(name = "Amplitude treatment", values = colors_ampl) +
  scale_x_continuous(breaks = seq(from = 0, to = length(levels(df_dist_prod$dist_class)), 
                                  by = 5)) +
  labs(x = "Distance to AR [m]", y = "Mean value per cell") +
  theme_classic(base_size = 10) + 
  theme(legend.position = "bottom", panel.border = element_rect(size = 0.5, fill = NA))

suppoRt::save_ggplot(plot = gg_dist, filename = paste0("00-dist-pp", extension), 
                     path = "04_Figures/", width = height, height = width, units = units,  dpi = dpi, 
                     overwrite = overwrite)
