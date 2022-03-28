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

paths <- list(nofish = "02_Data/01_VarAmp-CV-nofish.rds", local = "02_Data/01_VarAmp-CV-local.rds")

#### Combine to one data.frame ####

df_results <- purrr::map(paths, function(i) {
  
  readr::read_rds(i) %>% 
    dplyr::group_by(enrichment_lvl, amplitude_lvl, n_diff, type, part, measure) %>%
    dplyr::summarise(value_mn = mean(value), value_sd = sd(value), .groups = "drop") %>% 
    dplyr::mutate(enrichment_lvl = factor(enrichment_lvl, levels = c(0.75, 1.0, 1.25), labels = c("low", "medium", "high")), 
                  amplitude_lvl = factor(amplitude_lvl, levels = c(0.05, 0.5, 1.0), labels = c("low", "medium", "high")),
                  n_diff = as.integer(n_diff),
                  type = factor(type, levels = c("cv", "absolute")),
                  part = factor(part, levels = c("ag_biomass", "ag_production", "ttl_biomass",
                                                 "bg_biomass", "bg_production", "ttl_production")),
                  measure = factor(measure, levels = c("alpha", "beta", "gamma", "synchrony")))}) %>% 
  purrr::reduce(dplyr::left_join, by = c("enrichment_lvl", "amplitude_lvl", "n_diff", 
                                         "type", "part", "measure"), suffix = c(".nf", ".im")) %>% 
  dplyr::mutate(mn_diff = (value_mn.im - value_mn.nf) / value_mn.nf * 100,
                sd_diff = (value_sd.im - value_sd.nf) / value_sd.nf * 100)

#### Create ggplot ####

col_palette <- c("#5ABCD6", "#FAD510", "#F22301")

legend_position <- c("none", "none", "bottom")

x_axis <- c(" ", " ", "Variability of amplitude")

x_labels <- list(element_blank(), element_blank(), NULL)

labels_facet <- list(c(low = "Low enrichment", medium = "Medium enrichment", high = "High enrichment"), 
                     c(low = "", medium = "", high = ""), 
                     c(low = "", medium = "", high = ""))

rel_heights <- c(1, 1, 1.25)

# loop through production and biomass output
gg_results <- purrr::map(c("production", "biomass"), function(i) {
  
  # create parts to loop through
  parts <- paste0(c("ag_", "bg_", "ttl_"), i)
  
  # create names for plot labelling
  names(parts) <- paste(c("Aboveground", "Belowground", "Total"), i)
  
  gg_scale <- purrr::map(c("gamma", "beta"), function(j) {
    
    # y_axis <- c("", expression(paste("Coefficient of variation ", gamma)), "")
    y_axis <- c(" ", paste0("Difference CV [%] (", j, ")"), " ")
    
    gg_temp <- purrr::map(seq_along(parts), function(k) {
    
      df_temp <- dplyr::filter(df_results, type == "cv", measure == j, part == parts[[k]], 
                               is.finite(mn_diff), is.finite(sd_diff))
    
        ggplot(data = df_temp) +
        geom_point(aes(x = n_diff, y = mn_diff, col = amplitude_lvl)) +
        geom_line(aes(x = n_diff, y = mn_diff, col = amplitude_lvl)) +
        # geom_linerange(aes(x = n_diff, ymin = mn_diff - sd_diff,
        #                    ymax = mn_diff + sd_diff, col = amplitude_lvl, group = part)) +
        geom_hline(yintercept = 0, linetype = 2, col = "grey") +
        facet_wrap(. ~ enrichment_lvl, ncol = 3, nrow = 1, 
                   labeller = labeller(enrichment_lvl = labels_facet[[k]])) +
        scale_x_continuous(breaks = seq(from = 0, to = 9, by = 1)) +
        scale_y_continuous(limits = c(-max(abs(df_temp$mn_diff)), max(abs(df_temp$mn_diff)))) +
        scale_color_manual(name = "Amplitude treatment", values = col_palette) +
        labs(y = y_axis[k], x = x_axis[k], subtitle = names(parts)[k]) +
        theme_classic(base_size = base_size) + 
        theme(legend.position = legend_position[k], strip.background = element_blank(), 
              strip.text = element_text(hjust = 0), 
              axis.text.x = x_labels[[k]])
      
    })
    
    # combine to one ggplot
    gg_temp <- cowplot::plot_grid(plotlist = gg_temp, ncol = 1, nrow = 3, 
                                  rel_heights = rel_heights)
    
  })
  
  # set names of scale level
  names(gg_scale) <- c("gamma", "beta")
  
  return(gg_scale)
  
})

names(gg_results) <- c("production", "biomass")

# loop through output level
purrr::walk(seq_along(gg_results), function(i) {
  
  # loop through scale level
  purrr::walk(seq_along(gg_results[[i]]), function(j) {
    
    # create file name
    filename_temp <- paste0("01_VarAmp_CV_", names(gg_results)[[i]], "_", 
                            names(gg_results[[i]])[[j]], "_diff.png")
    
    # save ggplot
    suppoRt::save_ggplot(plot = gg_results[[i]][[j]], filename = filename_temp,
                         path = "04_Figures", width = height, height = width, dpi = dpi, 
                         units = units, overwrite = overwrite)
    
  })
})
 
