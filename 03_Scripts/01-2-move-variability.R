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
                amplitude_mn = factor(amplitude_mn, ordered = TRUE, labels = c("Low", "Medium", "High"))) %>% 
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
                amplitude_mn = factor(amplitude_mn, ordered = TRUE, labels = c("Low", "Medium", "High"))) %>%
  tibble::tibble()

#### Setup globals ####

parts <- c("ag_production", "bg_production", "ttl_production")

measures <- c("alpha", "gamma")

colors_ampl <- c("Low" = "#622567", "Medium" = "#007054", "High" = "#ec7122") 

colors_pop <- c("8" = "#e19ed5", "16" = "#7a92dc", "32" = "#ebcc60", "64" = "#cb3a54")

#### Total production

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
  theme(legend.position = "bottom", strip.text = element_text(hjust = 0))

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
  theme(legend.position = "none", strip.text = element_text(hjust = 0))

gg_pp <- cowplot::plot_grid(gg_pp_agbg, gg_pp_tt, ncol = 2)

suppoRt::save_ggplot(plot = gg_pp, filename = paste0("00-amplitude-prod.pdf"), 
                     path = "04_Figures/", width = width, height = height / 2, 
                     units = units,  dpi = dpi,overwrite = overwrite)

#### Mobility vs. CV ####

for (pop_i in c(8, 16, 32, 64)) {
  
 df_pop_n <- dplyr::filter(df_cv_prod, pop_n == pop_i, value.move != 0)
 
 gg_parts <- purrr::map(parts, function(part_i) {
   
   gg_measures <- purrr::map(measures, function(measure_i) {
     
     gg_tile <- purrr::map(c("Low", "Medium", "High"), function(amp_i) {
       
       df_tile <- dplyr::filter(df_pop_n, part == part_i, measure == measure_i, 
                                amplitude_mn == amp_i)
       
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
    
     df_scatter <- dplyr::filter(df_pop_n, part == part_i, measure == measure_i)
     
     df_scatter_label <- dplyr::group_by(df_scatter, amplitude_mn) %>% 
       dplyr::summarise(label = paste(unique(amplitude_mn), "amplitude"), 
                        x = min(value.move) + (max(value.move) - min(value.move)) * 0.15,
                        y = mean(value.cv), .groups = "drop")
     
     gg_scatter <- ggplot(data = df_scatter, aes(x = value.move, y = value.cv, 
                                                 color = amplitude_mn)) +
       geom_hline(yintercept = 0, linetype = 2, color = "grey") +
       geom_point(shape = 19, alpha = 0.1) +
       geom_smooth(formula = "y ~ x", method = "lm", se = FALSE, size = 0.5) +
       geom_label(data = df_scatter_label, aes(x = x, y = y, label = label), size = 1.75) +
       scale_color_manual(name = "Amplitude", values = colors_ampl) +
       scale_x_continuous(breaks = seq(from = min(df_scatter$value.move), to = max(df_scatter$value.move), 
                                       length.out = 5), 
                          limits = c(min(df_scatter$value.move), max(df_scatter$value.move)), 
                          labels = function(x) round(x, 2)) +
       scale_y_continuous(breaks = seq(from = 0.0, to = max(df_scatter$value.cv), length.out = 5), 
                          limits = c(0.0, max(df_scatter$value.cv)), labels = function(x) round(x, 2)) +
       labs(x = "Mean cross-ecosystem movement", 
            y = ifelse(measure_i == "alpha", yes = "Coefficient of variation", no = "")) +
       theme_classic(base_size = 8) + 
       theme(legend.position = "none", plot.margin = unit(c(-3.5, 2.5, 5.0, 2.5), "mm"), # t r b l
             strip.text = element_text(hjust = 0.0), strip.background = element_blank())
     
     cowplot::plot_grid(gg_tile, gg_scatter, nrow = 2, ncol = 1)
    
    })
   
   if (part_i == "ag_production") {
     
     cowplot::plot_grid(plotlist = gg_measures, nrow = 1, ncol = 2, labels = c("Local scale", "Meta scale"),  
                        hjust = c(-3, -3), label_fontface = "italic", label_size = 10)
     
   } else {
     
    cowplot::plot_grid(plotlist = gg_measures, nrow = 1, ncol = 2)
     
   }
  })
 
  gg_mobility_cv <- cowplot::plot_grid(plotlist = gg_parts, nrow = 3, ncol = 1, 
                                       labels = c("ag)", "bg)", "ttl)"), label_fontface = "italic", 
                                       label_size = 10)
  
  suppoRt::save_ggplot(plot = gg_mobility_cv, filename = paste0("01-mobility-cv-", pop_i, ".pdf"), 
                       path = "04_Figures/", height = height, width = width, 
                       units = units,  dpi = dpi, overwrite = overwrite)
  
}

df_regression <- dplyr::filter(df_cv_prod, part %in% c("ag_production", "bg_production", "ttl_production"), 
              measure %in% c("alpha", "gamma")) %>% 
  dplyr::group_by(pop_n, part, measure, amplitude_mn) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(df_temp) {
    
    lm_temp <- lm(value.cv ~ value.move, data = df_temp)
  
    tibble::tibble(pop_n = unique(df_temp$pop_n), part = unique(df_temp$part), 
                   measure = unique(df_temp$measure), amplitude_mn = unique(df_temp$amplitude_mn),
                   int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                   rr = summary(lm_temp)[["adj.r.squared"]], 
                   p = summary(lm_temp)[["coefficients"]][2, 4]) 
  })

#### CV vs. PP ####

for (pop_i in c(8, 16, 32, 64)) {
  
  gg_parts <- purrr::map(parts, function(part_i) {
    
    gg_measures <- purrr::map(measures, function(measure_i) {
      
      df_scatter <- dplyr::filter(df_cv_prod, pop_n == pop_i, measure == measure_i, 
                                  part %in% part_i, value.move != 0)
      
      df_regression <- dplyr::group_by(df_scatter, amplitude_mn) %>%
        dplyr::group_split() %>%
        purrr::map_dfr(function(lm_i) {
          
          lm_temp <- lm(value.prod ~ value.cv, data = lm_i)
          
          tibble::tibble(int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                         rr = summary(lm_temp)[["adj.r.squared"]], p =  summary(lm_temp)[["coefficients"]][[2, 4]]) %>%
            dplyr::mutate(rr = dplyr::case_when(rr <= 0.0 ~ 0.0, 
                                                rr > 0.0 ~ rr),
                          dir = dplyr::case_when(coef < 0 ~ "-", TRUE ~ "+"),
                          signif = dplyr::case_when(p < 0.05 & p >= 0.01 ~ "'*'",
                                                    p < 0.01 & p >= 0.001 ~ "'**'",
                                                    p < 0.001 ~ "'***'", TRUE ~ "' '"),
                          int = convert_label(int), coef = convert_label(abs(coef)), rr = signif(rr, 1),
                          part = unique(lm_i$part), measure = unique(lm_i$measure), 
                          amplitude_mn = unique(lm_i$amplitude_mn))
        }) %>% 
        dplyr::mutate(x = min(df_scatter$value.cv) + 
                        (max(df_scatter$value.cv) - min(df_scatter$value.cv)) * 0.175)
      
      if (dplyr::filter(df_scatter, value.cv == min(value.cv))[["value.prod"]] >
          mean(df_scatter$value.prod)) {
        
        df_regression$y <- min(df_scatter$value.prod) + 
          (max(df_scatter$value.prod) - min(df_scatter$value.prod)) * c(0.25, 0.15, 0.05)
        
      } else {
        
        df_regression$y <- min(df_scatter$value.prod) + 
          (max(df_scatter$value.prod) - min(df_scatter$value.prod)) * c(0.75, 0.85, 0.95)
        
      }
      
      if (measure_i == "alpha") {
        y_axis <- NULL
      } else {
        y_axis <- element_blank()
      }
      
      ggplot(data = df_scatter, aes(x = value.cv, y = value.prod, color = amplitude_mn)) +
        geom_point(shape = 19, alpha = 0.05) +
        geom_smooth(formula = "y ~ x", method = "lm", se = FALSE, size = 0.5) +
        geom_text(data = df_regression, parse = TRUE, size = 2.75,
                  aes(x = x, y = y, label = paste0("italic(y)==", int, dir, coef, "*italic(x)*';'~R[adj]^2==", rr),
                      color = amplitude_mn)) +
        scale_color_manual(name = "Amplitude treatment", values = colors_ampl) +
        scale_x_continuous(breaks = seq(from = min(df_scatter$value.cv), to = max(df_scatter$value.cv), 
                                        length.out = 5), 
                           limits = c(min(df_scatter$value.cv), max(df_scatter$value.cv)), 
                           labels = function(x) round(x, 2)) +
        scale_y_continuous(breaks = seq(from = min(df_scatter$value.prod), to = max(df_scatter$value.prod), 
                                        length.out = 5), 
                           limits = c(min(df_scatter$value.prod), max(df_scatter$value.prod)), 
                           labels = function(x) round(x, 2)) +
        labs(x = "", y = "") +
        theme_classic(base_size = 10) + 
        theme(legend.position = "none", strip.text = element_text(hjust = 0.0), 
              strip.background = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
              axis.text.y = y_axis)
  
    })
    
    cowplot::plot_grid(plotlist = gg_measures, ncol = 2)
    
  })
  
  gg_dummy <- ggplot(data = df_scatter, aes(x = value.cv, y = value.prod, color = amplitude_mn)) +
    geom_point(shape = 19) + geom_line(size = 0.5) +
    scale_color_manual(name = "Amplitude treatment", values = colors_ampl) + 
    theme_classic(base_size = 10) + 
    theme(legend.position = "bottom")
  
  gg_legend <- cowplot::get_legend(gg_dummy)
  
  gg_scatter <- cowplot::plot_grid(plotlist = gg_parts, nrow = 3) + 
    draw_label("Coeffiecent of variation", x = 0.5, y = 0, vjust = -1.25, angle = 0, 
               size = 10) +
    draw_label("Total primary production", x = 0, y = 0.5, vjust = 1.25, angle = 90, 
               size = 10)
  
  gg_scatter <- plot_grid(gg_scatter, gg_legend, ncol = 1, rel_heights = c(0.95, 0.05))
  
  suppoRt::save_ggplot(plot = gg_scatter, filename = paste0("01-cv-pp-", pop_i, ".pdf"), 
                       path = "04_Figures/", width = height, height = width, units = units,  dpi = dpi, 
                       overwrite = overwrite)
  
}


#### PE vs. PP ####

for (pop_i in c(8, 16, 32, 64)) {
  
  gg_parts <- purrr::map(parts, function(part_i) {
    
    gg_amplitude <- purrr::map(c("Low", "Medium", "High"), function(amplitude_i) {
      
      df_scatter <- dplyr::filter(df_cv_prod, pop_n == pop_i, measure == "beta", 
                                  amplitude_mn == amplitude_i, part %in% part_i, value.move != 0)
      
      df_regression <- dplyr::group_by(df_scatter, amplitude_mn) %>%
        dplyr::group_split() %>%
        purrr::map_dfr(function(lm_i) {
          
          lm_temp <- lm(value.prod ~ value.cv, data = lm_i)
          
          tibble::tibble(int = lm_temp$coefficients[[1]], coef = lm_temp$coefficients[[2]],
                         rr = summary(lm_temp)[["adj.r.squared"]], p =  summary(lm_temp)[["coefficients"]][[2, 4]]) %>%
            dplyr::mutate(rr = dplyr::case_when(rr <= 0.0 ~ 0.0,
                                                rr > 0.0 ~ rr),
                          dir = dplyr::case_when(coef < 0 ~ "-", TRUE ~ "+"),
                          signif = dplyr::case_when(p < 0.05 & p >= 0.01 ~ "'*'",
                                                    p < 0.01 & p >= 0.001 ~ "'**'",
                                                    p < 0.001 ~ "'***'", TRUE ~ "' '"),
                          int = convert_label(int), coef = convert_label(abs(coef)), rr = signif(rr, 1),
                          part = unique(lm_i$part), measure = unique(lm_i$measure),
                          amplitude_mn = unique(lm_i$amplitude_mn))
        }) %>%
        dplyr::mutate(x = min(df_scatter$value.cv) +
                        (max(df_scatter$value.cv) - min(df_scatter$value.cv)) * 0.75)
      
      if (dplyr::filter(df_scatter, value.cv == min(value.cv))[["value.prod"]] <
          mean(df_scatter$value.prod)) {
        
        df_regression$y <- min(df_scatter$value.prod) +
          (max(df_scatter$value.prod) - min(df_scatter$value.prod)) * 0.1
        
      } else {
        
        df_regression$y <- min(df_scatter$value.prod) +
          (max(df_scatter$value.prod) - min(df_scatter$value.prod)) * 0.9
        
      }
      
      if (amplitude_i == "Low") {
        y_axis <- NULL
      } else {
        y_axis <- element_blank()
      }
      
      ggplot(data = df_scatter, aes(x = value.cv, y = value.prod, color = amplitude_mn)) +
        geom_point(shape = 19, alpha = 0.1) +
        geom_smooth(formula = "y ~ x", method = "lm", se = FALSE, size = 0.5) +
        geom_text(data = df_regression, parse = TRUE, size = 2.75,
                  aes(x = x, y = y, label = paste0("italic(y)==", int, dir, coef, "*italic(x)*';'~R[adj]^2==", rr),
                      color = amplitude_mn), color = "black") +
        scale_color_manual(name = "Amplitude treatment", values = colors_ampl) +
        scale_x_continuous(breaks = seq(from = min(df_scatter$value.cv), to = max(df_scatter$value.cv),
                                        length.out = 5),
                           limits = c(min(df_scatter$value.cv), max(df_scatter$value.cv)),
                           labels = function(x) round(x, 2)) +
        scale_y_continuous(breaks = seq(from = min(df_scatter$value.prod), to = max(df_scatter$value.prod),
                                        length.out = 5),
                           limits = c(min(df_scatter$value.prod), max(df_scatter$value.prod)),
                           labels = function(x) round(x, 2)) +
        labs(x = "", y = "") +
        theme_classic(base_size = 10) + 
        theme(legend.position = "none", strip.text = element_text(hjust = 0.0), 
              strip.background = element_blank(), panel.border = element_rect(size = 0.5, fill = NA), 
              strip.text.x = element_blank(),
              axis.text.y = y_axis)
  
    })
    
    cowplot::plot_grid(plotlist = gg_amplitude, ncol = 3)
    
  })
  
  gg_dummy <- ggplot(data = df_scatter, aes(x = value.cv, y = value.prod, color = amplitude_mn)) +
    geom_point(shape = 19) + geom_line(size = 0.5) +
    scale_color_manual(name = "Amplitude treatment", values = colors_ampl) + 
    theme_classic(base_size = 10) + 
    theme(legend.position = "bottom")
  
  gg_legend <- cowplot::get_legend(gg_dummy)
  
  gg_scatter <- cowplot::plot_grid(plotlist = gg_parts, nrow = 3) + 
    draw_label("Portfolio effect", x = 0.5, y = 0, vjust = -1.25, angle = 0, 
               size = 10) +
    draw_label("Total primary production", x = 0, y = 0.5, vjust = 1.25, angle = 90, 
               size = 10)
  
  gg_scatter <- plot_grid(gg_scatter, gg_legend, nrow = 2, rel_heights = c(0.95, 0.05))
  
  suppoRt::save_ggplot(plot = gg_scatter, filename = paste0("01-pe-pp-", pop_i, ".pdf"), 
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
             labeller = labeller(part = c(ag_production = "Aboveground PP", 
                                          bg_production = "Belowground  PP", 
                                          ttl_production = "Total PP", 
                                          nutrients_pool = "Nutrients"), 
                                 pop_n = c("8" = "Local population n: 8", 
                                           "16" = "Local population n: 16", 
                                           "32" = "Local population n: 32", 
                                           "64" = "Local population n: 64"))) + 
  scale_color_manual(name = "Amplitude_treatment", values = colors_ampl) +
  scale_x_continuous(breaks = seq(from = 0, to = length(levels(df_dist_prod$dist_class)), 
                                  by = 5)) +
  labs(x = "Distance to AR [m]", y = "Mean value per cell") +
  theme_classic(base_size = 10) + 
  theme(legend.position = "bottom", panel.border = element_rect(size = 0.5, fill = NA))

suppoRt::save_ggplot(plot = gg_dist, filename = "01-dist-pp-.pdf", 
                     path = "04_Figures/", width = height, height = width, units = units,  dpi = dpi, 
                     overwrite = overwrite)

