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

col_palette <- c("#5ABCD6", "#FAD510", "#F22301")

comparison <- c("nofish_local", "nofish_mobile", "local_mobile")

labels_col <- c(low = "Low enrichment", medium = "Medium enrichment", high = "High enrichment")

labels_row <- c(ag_production = "Aboveground production", bg_production = "Belowground production", 
                ttl_production = "Total production")

gg_diff <- purrr::map(comparison, function(i) {
  
  comparison_temp <- paste0(i, c("_mn", "_sd"))
  
  y_lab <- paste("Diff CV", stringr::str_replace(i, pattern = "_", replacement = "-"), "[%]")
  
  dplyr::filter(result_rmse, type == "cv", measure == "gamma", 
                part %in% c("ag_production", "bg_production", "ttl_production")) %>% 
    dplyr::select(enrichment_lvl, amplitude_lvl, n_diff, part, comparison_temp) %>% 
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2, col = "grey") +
    geom_point(aes(x = n_diff, y = get(comparison_temp[1]), col = amplitude_lvl)) +
    geom_line(aes(x = n_diff, y = get(comparison_temp[1]), col = amplitude_lvl), 
              alpha = 0.2) +
    geom_errorbar(aes(x = n_diff, ymin = get(comparison_temp[1]) - get(comparison_temp[2]),
                      ymax = get(comparison_temp[1]) + get(comparison_temp[2]), 
                      col = amplitude_lvl, group = part), width = 0.2) +
    facet_grid(rows = vars(part), cols = vars(enrichment_lvl), scales = "free_y", 
               labeller = labeller(enrichment_lvl = labels_col, part = labels_row)) +
    scale_x_continuous(breaks = seq(from = 0, to = 9, by = 1)) +
    scale_color_manual(name = "Amplitude treatment", values = col_palette) +
    labs(y = y_lab, x = "Variability of amplitude") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "bottom")
})

names(gg_diff) <- comparison
gg_diff$local_mobile

