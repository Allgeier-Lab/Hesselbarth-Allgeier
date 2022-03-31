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

paths <- list(nofish =  "02_Data/01_VarAmp-CV-nofish.rds", 
              local = "02_Data/01_VarAmp-CV-local.rds", 
              mobile = "02_Data/01_VarAmp-CV-mobile.rds")

#### Combine to one data.frame ####

result_rmse <- purrr::map_dfr(paths, readr::read_rds, .id = "run") %>%
  tidyr::pivot_wider(names_from = run, values_from = value, values_fn = list) %>% 
  tidyr::unnest(cols = c(nofish, local, mobile)) %>% 
  dplyr::mutate(nofish_local = (local - nofish) / nofish * 100, 
                nofish_mobile = (mobile - nofish) / nofish * 100, 
                local_mobile = (mobile - local) / local * 100) %>% 
  dplyr::group_by(enrichment_lvl, amplitude_lvl, n_diff, type, part, measure) %>%
  dplyr::summarise(nofish_local_mn = mean(nofish_local), nofish_local_sd = sd(nofish_local), 
                   nofish_mobile_mn = mean(nofish_mobile), nofish_mobile_sd = sd(nofish_mobile), 
                   local_mobile_mn = mean(local_mobile), local_mobile_sd = sd(local_mobile), 
                   .groups = "drop") %>% 
  dplyr::mutate(enrichment_lvl = factor(enrichment_lvl, levels = c(0.75, 1.0, 1.25), labels = c("low", "medium", "high")), 
                amplitude_lvl = factor(amplitude_lvl, levels = c(0.05, 0.5, 1.0), labels = c("low", "medium", "high")),
                n_diff = as.integer(n_diff),
                type = factor(type, levels = c("cv", "absolute")),
                part = factor(part, levels = c("ag_biomass", "ag_production", "ttl_biomass",
                                               "bg_biomass", "bg_production", "ttl_production")),
                measure = factor(measure, levels = c("alpha", "beta", "gamma", "synchrony")))

#### Create ggplot ####

# specify what to plot
comparisons <- c("nofish_local", "nofish_mobile", "local_mobile")

parts <- c(rep(x = "ag_production", times = 3), rep(x = "bg_production", times = 3))

enrichments <- rep(x = c("low", "medium", "high"), times = 2)

# setup some plot options
col_palette <- c("#5ABCD6", "#FAD510", "#F22301")

# setup labels
label_parts <- c("Aboveground production", "" , "", 
                 "Belowground production", "", "")

label_enrich <- c("Low enrichment", "Medium enrichment" , "High enrichment",
                  "", "", "")

# create plots in loop
gg_comparisons <- purrr::map(comparisons, function(i) {
  
  # get mean and sd of current comparison
  comparison_temp <- paste0(i, c("_mn", "_sd"))

  # loop through all parts and treatments
  gg_treatments <- purrr::map(seq_along(parts), function(j) {
    
    dplyr::filter(result_rmse, enrichment_lvl == enrichments[j], type == "cv", 
                  measure == "gamma", part == parts[j]) %>% 
      ggplot() +
      geom_hline(yintercept = 0, linetype = 2, col = "grey") +
      geom_point(aes(x = n_diff, y = get(comparison_temp[1]), col = amplitude_lvl)) +
      geom_line(aes(x = n_diff, y = get(comparison_temp[1]), col = amplitude_lvl), 
                alpha = 0.2) +
      geom_errorbar(aes(x = n_diff, ymin = get(comparison_temp[1]) - get(comparison_temp[2]),
                        ymax = get(comparison_temp[1]) + get(comparison_temp[2]), 
                        col = amplitude_lvl, group = part), width = 0.2) +
      scale_x_continuous(breaks = seq(from = 0, to = 9, by = 1)) +
      scale_color_manual(values = col_palette) +
      labs(y = "", x = "", subtitle = label_parts[j], title = label_enrich[j]) +
      theme_classic(base_size = base_size) +
      theme(legend.position = "none")
    
    })
  
  # create dummy plot to grab legend
  gg_dummy <- ggplot(data = result_rmse) + 
    geom_point(aes(x = 1, y = 1, col = amplitude_lvl)) +
    scale_color_manual(name = "Amplitude treatment", values = col_palette) + 
    theme(legend.position = "bottom")
  
  # grab legend
  legend <- cowplot::get_legend(gg_dummy)
  
  # create x/y labs
  x_lab <- "Variability of amplitude"
  
  y_lab <- stringr::str_replace(i, pattern = "_", replacement = " - ")
  y_lab <- expr(paste(Delta, "CV ", !!y_lab, " [%]"))
  
  # add axsis labels
  gg_temp <- cowplot::plot_grid(plotlist = gg_treatments, nrow = 2, ncol = 3) + 
    cowplot::draw_label(label = x_lab, x = 0.5, y = 0, vjust = -0.5, angle = 0) + 
    cowplot::draw_label(label = y_lab, x = 0, y = 0.5, vjust = 1.5, angle = 90)
  
  # add legend
  plot_grid(gg_temp, legend, rel_heights = c(0.9, 0.1), nrow = 2)
  
})
    
names(gg_comparisons) <- comparisons

