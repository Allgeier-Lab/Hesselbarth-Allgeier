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
  purrr::set_names(c("8", "16", "32", "64", "128"))

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
 
# df_dist_prod <- purrr::map_dfr(file_names, function(i) {
#     readr::read_rds(i) %>% purrr::map_dfr(function(j) {
#       dplyr::rename(j$dist, value.prod = value) %>%
#         dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))
#     })
#   }, .id = "pop_n") %>%
#   dplyr::mutate(pop_n = factor(as.numeric(pop_n), ordered = TRUE),
#                 dist_class = factor(dist_class, ordered = TRUE),
#                 part = factor(part, levels = c("ag_production", "bg_production", "ttl_production",
#                                                "nutrients_pool")),
#                 amplitude_mn = factor(amplitude_mn, ordered = TRUE, labels = c("Low", "Medium", "High")),
#                 value.prod.log = log(value.prod), value.move.log = log(value.move)) %>%
#   tibble::tibble()

#### Setup globals ####

parts <- c("ag_production", "bg_production", "ttl_production")

measures <- c("alpha", "gamma")

colors_ampl <- c("Low" = "#622567", "Medium" = "#007054", "High" = "#ec7122") 

colors_pop <- c("8" = "#0c7156", "16" = "#e2998a", "32" = "#ea7428", "64" = "#cf3a36", "128" = "#663171")

# #### Total production ####

gg_pp_agbg <- dplyr::filter(df_cv_prod, part %in% c("ag_production", "bg_production"),
                            measure == "alpha") %>%
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

gg_pp_tt <- dplyr::filter(df_cv_prod, part == "ttl_production", measure == "alpha") %>%
  ggplot(aes(x = amplitude_mn, y = value.prod)) +
  geom_jitter(aes(color = pop_n), alpha = 0.1, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5),
              shape = 19) +
  geom_boxplot(aes(fill = pop_n), position = position_dodge(width = 0.5), outlier.color = NA) +
  scale_color_manual(name = "Population size", values = colors_pop) +
  scale_fill_manual(name = "Population size", values = colors_pop) +
  facet_wrap(. ~ part, scales = "free_y", nrow = 1, ncol = 1,
             labeller = labeller(part = c(ag_production = "Aboveground",
                                          bg_production = "Belowground ",
                                          ttl_production = "Total"))) +
  labs(x = "Amplitude treatment", y = "") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(hjust = 0),
        strip.background = element_blank(), panel.border = element_rect(size = 0.5,
                                                                        fill = NA))

gg_pp <- cowplot::plot_grid(gg_pp_agbg, gg_pp_tt, ncol = 2)

suppoRt::save_ggplot(plot = gg_pp, filename = paste0("00-amplitude-prod", extension),
                     path = "04_Figures/", width = width, height = height / 2,
                     units = units, dpi = dpi, overwrite = overwrite)

#### Mobility vs. CV (log-linear) ####

for (pop_i in c(8, 16, 32, 64, 128)) {
  
 gg_parts <- purrr::map(parts, function(part_i) {
   
   gg_measures <- purrr::map(measures, function(measure_i) {
     
     gg_tile <- purrr::map(c("Low", "Medium", "High"), function(amp_i) {
       
       df_tile <- dplyr::filter(df_cv_prod, pop_n == pop_i, part == part_i, 
                                measure == measure_i, amplitude_mn == amp_i)
       
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
                                 measure == measure_i)
     
     # log-linear
     df_regression <- dplyr::group_by(df_scatter, amplitude_mn) %>%
       dplyr::group_split() %>%
       purrr::map_dfr(function(df_temp) {

         lm_temp <- lm(value.cv.log ~ value.move, data = df_temp)
         
         pred <- coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * df_temp$value.move
         
         r.squared <- cor(df_temp$value.cv.log, pred)
         
         p.value <- summary(lm_temp)[["coefficients"]][2, 4]
         
         tibble::tibble(amplitude_mn = unique(df_temp$amplitude_mn), 
                        value.move = df_temp$value.move, value.cv.log.pred = pred, 
                        int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                        r.squared = r.squared, p.value = p.value,
                        x = mean(df_temp$value.move), y = mean(df_temp$value.cv.log))
         
     })
     
     df_label <- dplyr::select(df_regression, -value.move, -value.cv.log.pred) %>% 
       dplyr::group_by(amplitude_mn) %>% 
       dplyr::distinct() %>% 
       dplyr::mutate(dir = dplyr::case_when(coef < 0.0 ~ "-", 
                                            coef > 0.0 ~ "+"),
                     int = round(int, digits = 2), coef = convert_label(abs(coef), digits = 2), 
                     r.squared = round(x = r.squared, digits = 2), 
                     p.value = dplyr::case_when(p.value <  0.001 ~  "'***'", 
                                                p.value <  0.01 ~  "'**'",
                                                p.value <  0.05 ~  "'*'", 
                                                p.value >  0.05 ~  "n.s."), 
                     label = paste0("italic(y)==", int, dir, coef, "*italic(x)*';'~R^2==", r.squared, "*';'~", p.value),
                     label = dplyr::case_when(p.value == "n.s." ~ "n.s.", 
                                              TRUE ~ label))

     gg_scatter <- ggplot(data = df_scatter, aes(x = value.move, y = value.cv.log, 
                                                 color = amplitude_mn)) +
       geom_point(shape = 19, alpha = 0.25) +
       geom_line(data = df_regression, aes(x = value.move, y = value.cv.log.pred, color = amplitude_mn)) +
       geom_label(data = df_label, parse = TRUE, size = 1.75, aes(x = x, y = y, color = amplitude_mn, label = label)) +
       scale_color_manual(name = "Amplitude", values = colors_ampl) +
       scale_x_continuous(breaks = seq(from = min(df_scatter$value.move), to = max(df_scatter$value.move),
                                       length.out = 5),
                          limits = c(min(df_scatter$value.move), max(df_scatter$value.move)),
                          labels = function(x) round(x, 1)) +
       scale_y_continuous(breaks = seq(from = min(df_scatter$value.cv.log),
                                       to = max(df_scatter$value.cv.log), length.out = 5),
                          limits = c(min(df_scatter$value.cv.log), max(df_scatter$value.cv.log)),
                          labels = function(x) round(x, 3)) +
       labs(x = "Mean cross-ecosystem movement", 
            y = ifelse(measure_i == "alpha", yes = "log(Coefficient of variation)", no = "")) +
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
 
 gg_dummy <- dplyr::filter(df_cv_prod, part == "ag_production", measure == "alpha", 
                           amplitude_mn == "Low") %>% 
   ggplot(aes(x = value.move, y = value.cv, color = amplitude_mn)) +
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

for (pop_i in c(8, 16, 32, 64, 128)) {
  
  gg_parts <- purrr::map(parts, function(part_i) {
    
    gg_measures <- purrr::map(measures, function(measure_i) {
      
      df_scatter <- dplyr::filter(df_cv_prod, pop_n == pop_i, 
                                  measure == measure_i, part %in% part_i)
      
      df_regression <- dplyr::group_by(df_scatter, amplitude_mn) %>%
        dplyr::group_split() %>%
        purrr::map_dfr(function(df_temp) {

          lm_temp <- lm(value.prod.log ~ value.cv.log, data = df_temp)
          
          pred <- coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * df_temp$value.cv.log
          
          r.squared <- cor(df_temp$value.prod.log, pred)
          
          p.value <- summary(lm_temp)[["coefficients"]][2, 4]
          
          tibble::tibble(amplitude_mn = unique(df_temp$amplitude_mn), 
                         value.cv.log = df_temp$value.cv.log, value.prod.log.pred = pred,
                         int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                         r.squared = r.squared, p.value = p.value)
          
      })
      
      if (part_i != "ag_production") {
        
        y_placement <- c(0.05, 0.2, 0.35)
        
      } else {
        
        y_placement <- c(0.65, 0.8, 0.95)
        
      }
      
      df_label <- dplyr::select(df_regression, -value.cv.log, -value.prod.log.pred) %>% 
        dplyr::group_by(amplitude_mn) %>% 
        dplyr::distinct() %>% 
        dplyr::mutate(dir = dplyr::case_when(coef < 0.0 ~ "-", coef > 0.0 ~ "+"),
                      int = round(int, digits = 2), coef = convert_label(abs(coef), digits = 2), 
                      r.squared = round(x = r.squared, digits = 2), 
                      p.value = dplyr::case_when(p.value <  0.001 ~  "'***'", 
                                                 p.value <  0.01 ~  "'**'",
                                                 p.value <  0.05 ~  "'*'", 
                                                 p.value >  0.05 ~  "n.s."), 
                      label = paste0("italic(y)==", int, dir, coef, "*italic(x)*';'~R^2==", r.squared, "*';'~", p.value),
                      label = dplyr::case_when(p.value == "n.s." ~ "n.s.", 
                                               TRUE ~ label))
      
      df_label$x <- min(df_scatter$value.cv.log) + 
        (max(df_scatter$value.cv.log) - min(df_scatter$value.cv.log)) * 0.15
      
      df_label$y <- min(df_scatter$value.prod.log) + 
        (max(df_scatter$value.prod.log) - min(df_scatter$value.prod.log)) * y_placement
      
      if (measure_i == "alpha") {
        
        y_axis <- NULL
        
      } else {
        
        y_axis <- element_blank()
        
      }
      
      ggplot(data = df_scatter, aes(x = value.cv.log, y = value.prod.log, color = amplitude_mn)) +
        geom_point(shape = 19, alpha = 0.25) +
        geom_line(data = df_regression, aes(x = value.cv.log, y = value.prod.log.pred, color = amplitude_mn)) +
        geom_label(data = df_label, parse = TRUE, size = 2.0,
                   aes(x = x, y = y, color = amplitude_mn, label = label)) +
        scale_color_manual(name = "Amplitude treatment", values = colors_ampl) +
        scale_x_continuous(breaks = seq(from = min(df_scatter$value.cv.log), 
                                        to = max(df_scatter$value.cv.log), length.out = 5),
                           limits = c(min(df_scatter$value.cv.log), max(df_scatter$value.cv.log)),
                           labels = function(x) round(x, 2)) +
        scale_y_continuous(breaks = seq(from = min(df_scatter$value.prod.log), 
                                        to = max(df_scatter$value.prod.log), length.out = 5),
                           limits = c(min(df_scatter$value.prod.log), max(df_scatter$value.prod.log)),
                           labels = function(x) round(x, 2)) +
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
  
  gg_dummy <- dplyr::filter(df_cv_prod, part == "ag_production", measure == "alpha",
                            amplitude_mn == "Low") %>% 
    ggplot(aes(x = value.cv, y = value.prod, color = amplitude_mn)) +
    geom_point(shape = 19) + geom_line(size = 0.5) +
    scale_color_manual(name = "Amplitude treatment", values = colors_ampl) + 
    theme_classic(base_size = 10) + 
    theme(legend.position = "bottom")
  
  gg_legend <- cowplot::get_legend(gg_dummy)
  
  gg_scatter <- cowplot::plot_grid(plotlist = gg_parts, nrow = 3, 
                                   labels = c("a)", "b)", "c)"), label_fontface = "italic", 
                                   label_size = 10) + 
    draw_label("log(Coeffiecent of variation)", x = 0.5, y = 0, vjust = -0.5, angle = 0, 
               size = 10) +
    draw_label("log(Total primary production)", x = 0, y = 0.5, vjust = 1.25, angle = 90, 
               size = 10)
  
  gg_scatter <- plot_grid(gg_scatter, gg_legend, ncol = 1, rel_heights = c(0.95, 0.05))
  
  suppoRt::save_ggplot(plot = gg_scatter, filename = paste0("01-cv-pp-", pop_i, extension), 
                       path = "04_Figures/", width = height, height = width, units = units,  dpi = dpi, 
                       overwrite = overwrite)
  
}

#### PE vs. PP (log-log) ####

for (pop_i in c(8, 16, 32, 64, 128)) {
  
  gg_parts <- purrr::map(parts, function(part_i) {
    
    gg_amplitude <- purrr::map(c("Low", "Medium", "High"), function(amplitude_i) {
      
      df_scatter <- dplyr::filter(df_cv_prod, pop_n == pop_i, measure == "beta", 
                                  amplitude_mn == amplitude_i, part %in% part_i)
      
      df_regression <- dplyr::group_by(df_scatter, amplitude_mn) %>%
        dplyr::group_split() %>%
        purrr::map_dfr(function(df_temp) {
          
          lm_temp <- lm(value.prod.log ~ value.cv.log, data = df_temp)
          
          pred <- coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * df_temp$value.cv.log
          
          r.squared <- cor(df_temp$value.prod.log, pred)
          
          p.value <- summary(lm_temp)[["coefficients"]][2, 4]
          
          tibble::tibble(amplitude_mn = unique(df_temp$amplitude_mn), 
                         value.cv.log = df_temp$value.cv.log, value.prod.log.pred = pred,
                         int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                         r.squared = r.squared, p.value = p.value)
      })
      
      y_placement <- ifelse(test = part_i == "ag_production", yes = 0.1, no = 0.9)
      
      df_label <- dplyr::select(df_regression, -value.cv.log, -value.prod.log.pred) %>% 
        dplyr::group_by(amplitude_mn) %>% 
        dplyr::distinct() %>% 
        dplyr::mutate(x = min(df_scatter$value.cv.log) + 
                        (max(df_scatter$value.cv.log) - min(df_scatter$value.cv.log)) * 0.8,
                      y = min(df_scatter$value.prod.log) + 
                        (max(df_scatter$value.prod.log) - min(df_scatter$value.prod.log)) * y_placement,
                      dir = dplyr::case_when(coef < 0.0 ~ "-", coef > 0.0 ~ "+"),
                      int = round(int, digits = 2), coef = convert_label(abs(coef), digits = 2), 
                      r.squared = round(x = r.squared, digits = 2), 
                      p.value = dplyr::case_when(p.value <  0.001 ~  "'***'", 
                                                 p.value <  0.01 ~  "'**'",
                                                 p.value <  0.05 ~  "'*'", 
                                                 p.value >  0.05 ~  "n.s."), 
                      label = paste0("italic(y)==", int, dir, coef, "*italic(x)*';'~R^2==", r.squared, "*';'~", p.value),
                      label = dplyr::case_when(p.value == "n.s." ~ "n.s.", 
                                               TRUE ~ label))
      
      if (amplitude_i == "Low") {
        
        y_axis <- NULL
        
      } else {
        
        y_axis <- element_blank()
        
      }
      
      ggplot(data = df_scatter, aes(x = value.cv.log, y = value.prod.log, color = amplitude_mn)) +
        geom_point(shape = 19, alpha = 0.25) +
        geom_line(data = df_regression, aes(x = value.cv.log, y = value.prod.log.pred, color = amplitude_mn)) +
        geom_label(data = df_label, parse = TRUE, size = 1.75, aes(x = x, y = y, color = amplitude_mn, label = label)) +
        scale_color_manual(name = "Amplitude treatment", values = colors_ampl) +
        scale_x_continuous(breaks = seq(from = min(df_scatter$value.cv.log), 
                                        to = max(df_scatter$value.cv.log), length.out = 5),
                           limits = c(min(df_scatter$value.cv.log), max(df_scatter$value.cv.log)),
                           labels = function(x) round(x, 2)) +
        scale_y_continuous(breaks = seq(from = min(df_scatter$value.prod.log),
                                        to = max(df_scatter$value.prod.log), length.out = 5),
                           limits = c(min(df_scatter$value.prod.log), max(df_scatter$value.prod.log)),
                           labels = function(x) round(x, 2)) +
        labs(x = "", y = "") +
        theme_classic(base_size = 10) + 
        theme(legend.position = "none", strip.text = element_text(hjust = 0.0), 
              strip.background = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
              strip.text.x = element_blank(), axis.text.y = y_axis, 
              plot.margin = unit(c(5.0, 2.5, 1.0, 2.5), "mm"))
  
    })
    
    cowplot::plot_grid(plotlist = gg_amplitude, ncol = 3)
    
  })
  
  gg_dummy <- dplyr::filter(df_cv_prod, part == "ag_production", measure == "alpha", 
                            amplitude_mn == "Low") %>% 
    ggplot(aes(x = value.cv, y = value.prod, color = amplitude_mn)) +
    geom_point(shape = 19) + geom_line(size = 0.5) +
    scale_color_manual(name = "Amplitude treatment", values = colors_ampl) + 
    theme_classic(base_size = 10) + 
    theme(legend.position = "bottom")
  
  gg_legend <- cowplot::get_legend(gg_dummy)
  
  gg_scatter <- cowplot::plot_grid(plotlist = gg_parts, nrow = 3, labels = c("a)", "b)", "c)"), 
                                   label_fontface = "italic", label_size = 10) + 
    draw_label("log(Portfolio effect)", x = 0.5, y = 0, vjust = -1.0, angle = 0, 
               size = 10) +
    draw_label("log(Total primary production)", x = 0, y = 0.5, vjust = 1.25, angle = 90, 
               size = 10)
  
  gg_scatter <- plot_grid(gg_scatter, gg_legend, nrow = 2, rel_heights = c(0.95, 0.05))
  
  suppoRt::save_ggplot(plot = gg_scatter, filename = paste0("01-pe-pp-", pop_i, extension), 
                       path = "04_Figures/", width = height, height = width, units = units,  dpi = dpi, 
                       overwrite = overwrite)
  
}

#### Heatmap ####

df_regression <- dplyr::filter(df_cv_prod, part %in% parts, measure %in% measures, pop_n != 128) %>% 
  dplyr::group_by(pop_n, part, measure, amplitude_mn) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(df_temp) {
    
    lm_temp <- lm(value.prod.log ~ value.cv.log, data = df_temp)
  
    pred <- coef(lm_temp)[[1]] + coef(lm_temp)[[2]] * df_temp$value.cv.log
    
    r.squared <- cor(df_temp$value.prod.log, pred)
    
    p.value <- summary(lm_temp)[["coefficients"]][2, 4]
    
    tibble::tibble(
      pop_n = unique(df_temp$pop_n), part = unique(df_temp$part), 
      measure = unique(df_temp$measure), amplitude_mn = unique(df_temp$amplitude_mn), 
      int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]], 
      r.squared = r.squared, p.value = p.value)}) %>% 
  dplyr::mutate(dir = dplyr::case_when(coef < 0.0 ~ "-", coef > 0.0 ~ "+"),
                # int = round(int, digits = 2), # coef = convert_label(abs(coef), digits = 2), 
                # r.squared = round(x = r.squared, digits = 2), 
                p.value = dplyr::case_when(p.value <  0.001 ~  "'***'", p.value <  0.01 ~  "'**'",
                                           p.value <  0.05 ~  "'*'", p.value >  0.05 ~  "n.s."), 
                # label = paste0("int==", int, "*'\n'~R^2==", r.squared, "*';'~", p.value),
                label_coef = paste0("coef=", round(coef,2), "\n R2=",  round(r.squared,2), "\n", p.value),
                label_int = paste0("int=", round(int,2), "\n R2=",  round(r.squared,2), "\n", p.value),
                label_coef = dplyr::case_when(p.value == "n.s." ~ "", TRUE ~ label_coef), 
                label_int = dplyr::case_when(p.value == "n.s." ~ "", TRUE ~ label_int))

# brks <- c(-max(c(abs(min(x$coef)), abs(max(x$coef)))), 0.0, max(c(abs(min(x$coef)), abs(max(x$coef)))))

gg_int <- ggplot(data = x, aes(x = pop_n, y = amplitude_mn)) + 
  geom_tile(aes(fill = int)) + 
  geom_text(aes(label = label_coef), size = 2.0, color = "white") +
  facet_grid(rows = vars(part), cols = vars(measure)) +
  # scale_fill_gradientn(breaks = brks, limits = c(brks[[1]], brks[3]), 
  #                      colors = c("darkorange1", "white", "navyblue")) +
  scale_fill_viridis_c(name = "Intercept", option = "C") + 
  labs(x = "Population size", y = "Amplitude treatment") +
  coord_equal() + 
  theme_classic() +
  theme(legend.position = "bottom", 
        axis.line = element_blank(), panel.border = element_rect(size = 1, fill = NA), 
        strip.background = element_blank(), strip.text.y = element_blank(), 
        strip.text.x = element_text(hjust = 0), legend.key.width = unit(10.0,"mm"))

gg_coef <- ggplot(data = x, aes(x = pop_n, y = amplitude_mn)) + 
  geom_tile(aes(fill = coef)) + 
  geom_text(aes(label = label_int), size = 2.0, color = "white") +
  facet_grid(rows = vars(part), cols = vars(measure)) +
  # scale_fill_gradientn(breaks = brks, limits = c(brks[[1]], brks[3]), 
  #                      colors = c("darkorange1", "white", "navyblue")) +
  scale_fill_viridis_c(name = "Coeffiecent", option = "C") + 
  labs(x = "Population size", y = "") +
  coord_equal() + 
  theme_classic() +
  theme(legend.position = "bottom", 
        axis.line = element_blank(), axis.text.y = element_blank(),
        panel.border = element_rect(size = 1, fill = NA), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0), 
        plot.margin = unit(c(5.0, 5.0, 5.0, 5.0), "mm"), legend.key.width = unit(10.0,"mm"))

gg_regression <- cowplot::plot_grid(gg_int, gg_coef, ncol = 2)

suppoRt::save_ggplot(plot = gg_regression, filename = paste0("01-regression", extension), 
                     path = "04_Figures/", width = height, height = width, units = units,  dpi = dpi, 
                     overwrite = T)

# #### Dist vs. PP ####
# 
# gg_dist <- dplyr::group_by(df_dist_prod, pop_n, dist_class, part, amplitude_mn) %>% 
#   dplyr::summarise(value.prod = mean(value.prod), .groups = "drop") %>%
#   dplyr::filter(part != "nutrients_pool") %>%
#   ggplot(aes(x = as.numeric(dist_class), y = value.prod, color = amplitude_mn, group = amplitude_mn)) + 
#   geom_line() +
#   geom_hline(yintercept = 0, linetype = 2, color = "grey") +
#   facet_grid(rows = dplyr::vars(part), cols = dplyr::vars(pop_n), scales = "free_y", 
#              labeller = labeller(part = c(ag_production = "Aboveground production", 
#                                           bg_production = "Belowground  production", 
#                                           ttl_production = "Total production", 
#                                           nutrients_pool = "Nutrients"), 
#                                  pop_n = c("8" = "Local population n=8", 
#                                            "16" = "Local population n=16", 
#                                            "32" = "Local population n=32", 
#                                            "64" = "Local population n=64"))) + 
#   scale_color_manual(name = "Amplitude treatment", values = colors_ampl) +
#   scale_x_continuous(breaks = seq(from = 0, to = length(levels(df_dist_prod$dist_class)), 
#                                   by = 5)) +
#   labs(x = "Distance to AR [m]", y = "Mean value per cell") +
#   theme_classic(base_size = 10) + 
#   theme(legend.position = "bottom", panel.border = element_rect(size = 0.5, fill = NA))
# 
# suppoRt::save_ggplot(plot = gg_dist, filename = paste0("00-dist-pp", extension), 
#                      path = "04_Figures/", width = height, height = width, units = units,  dpi = dpi, 
#                      overwrite = overwrite)
