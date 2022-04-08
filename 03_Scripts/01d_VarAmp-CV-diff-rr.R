##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

source("01_Functions/log_response.R")

#### Load data ####

paths <- c(nofish = "02_Data/01_VarAmp-CV-nofish.rds", 
           local = "02_Data/01_VarAmp-CV-local.rds", 
           mobile = "02_Data/01_VarAmp-CV-mobile.rds") %>% 
  suppoRt::expand_grid_unique(x = ., y = .) %>%
  tibble::as_tibble() %>% 
  magrittr::set_names(c("control", "treatment")) 

comparisons <- c("nofish_local", "nofish_mobile", "local_mobile")

#### Combine to one data.frame ####

response_rations <- purrr::map(1:nrow(paths), function(i) {
  
  control <- readr::read_rds(paths[[i, 1]]) %>% 
    tibble::add_column(id_row = 1:nrow(.), .before = "enrichment_lvl")
  
  treatment <- readr::read_rds(paths[[i, 2]]) %>% 
    tibble::add_column(id_row = 1:nrow(.), .before = "enrichment_lvl")
  
  treatments_list <- dplyr::inner_join(x = control, y = treatment, 
                                by = c("enrichment_lvl", "amplitude_lvl", "n_diff", 
                                       "type", "part", "measure", "id_row"), 
                                suffix = c(".ctrl", ".trtm")) %>% 
    dplyr::filter(type == "cv", measure == "gamma") %>% 
    dplyr::group_by(enrichment_lvl, amplitude_lvl, n_diff, part) %>% 
    dplyr::group_split()
  
  purrr::map_dfr(seq_along(treatments_list), function(j) {
      
    df_temp <- treatments_list[[j]]
      
    message("\r> Progres: i = ", i , "/", nrow(paths), " || j = ", j , "/",  
            length(treatments_list), "\t\t", appendLF = FALSE)
      
    if (any(df_temp$value.ctrl == 0) || any(df_temp$value.trtm == 0)) {
        
      tibble(enrichment_lvl = unique(df_temp$enrichment_lvl), amplitude_lvl = unique(df_temp$amplitude_lvl), 
             n_diff = unique(df_temp$n_diff), part = unique(df_temp$part), mean = as.numeric(NA), 
             lo = as.numeric(NA), hi = as.numeric(NA))
        
        
    } else {
        
      bootstrap <- boot::boot(data = tibble::tibble(ctrl = df_temp$value.ctrl, 
                                                    trtm = df_temp$value.trtm),
                              statistic = log_response, R = 1000)
        
      bootstrap_ci <- boot::boot.ci(bootstrap, type = "norm", conf = 0.95)
        
      tibble(enrichment_lvl = unique(df_temp$enrichment_lvl), amplitude_lvl = unique(df_temp$amplitude_lvl), 
             n_diff = unique(df_temp$n_diff), part = unique(df_temp$part),
             mean = mean(bootstrap$t[, 1]), lo = bootstrap_ci$normal[2], hi = bootstrap_ci$normal[3])
        
    }
  })
}) 

names(response_rations) <- comparisons

response_rations <- dplyr::bind_rows(response_rations, .id = "id") %>% 
  dplyr::mutate(id = factor(id, levels = comparisons),
                enrichment_lvl = factor(enrichment_lvl, levels = c(0.75, 1.0, 1.25), 
                                        labels = c("low", "medium", "high")), 
                amplitude_lvl = factor(amplitude_lvl, levels = c(0.05, 0.5, 1.0), 
                                       labels = c("low", "medium", "high")),
                n_diff = as.integer(n_diff), 
                part = factor(part, levels = c("ag_biomass", "ag_production", "ttl_biomass",
                                               "bg_biomass", "bg_production", "ttl_production")))

#### Create ggplot ####

parts <- c(rep(x = "ag_biomass", times = 3), rep(x = "bg_biomass", times = 3), 
           rep(x = "ttl_biomass", times = 3))

# parts <- c(rep(x = "ag_production", times = 3), rep(x = "bg_production", times = 3), 
#            rep(x = "ttl_production", times = 3))

enrichments <- rep(x = c("low", "medium", "high"), times = 3)

# setup some plot options
col_palette <- c("#5ABCD6", "#FAD510", "#F22301")

# setup labels
label_parts <- c("Aboveground production", "" , "", 
                 "Belowground production", "", "", 
                 "Total production", "", "")

label_enrich <- c("Low enrichment", "Medium enrichment" , "High enrichment",
                  "", "", "", 
                  "", "", "")

label_y <- list(NULL, element_blank(), element_blank()) %>% 
  rep(times = 3)

label_x <- list(element_blank(), element_blank(), element_blank(),
                element_blank(), element_blank(), element_blank(), 
                NULL, NULL, NULL)

# create plots in loop
gg_comparisons <- purrr::map(comparisons, function(i) {
  
  # # get mean and sd of current comparison
  # comparison_temp <- paste0(i, c("_mn", "_sd"))
  
  # loop through all parts and treatments
  gg_treatments <- purrr::map(seq_along(parts), function(j) {
    
    df_temp <- dplyr::filter(response_rations, id == i, enrichment_lvl == enrichments[j], 
                             part == parts[j])
    
    y_range <- dplyr::filter(response_rations, part == parts[j]) %>%
      dplyr::select(lo, hi) %>% 
      range()
    
    y_range[1] <- ifelse(test = y_range[1] > 0, yes = 0, no = y_range[1])
    
    y_range[2] <- ifelse(test = y_range[2] < 0, yes = 0, no = y_range[2])
    
    ggplot(data = df_temp) +
      geom_hline(yintercept = 0, linetype = 2, col = "grey") +
      geom_point(aes(x = n_diff, y = mean, col = amplitude_lvl)) +
      geom_line(aes(x = n_diff, y = mean, col = amplitude_lvl), 
                alpha = 0.5) +
      geom_errorbar(aes(x = n_diff, ymin = lo, ymax = hi,
                        col = amplitude_lvl, group = part), width = 0.2) +
      scale_x_continuous(breaks = seq(from = 0, to = 9, by = 1)) +
      scale_y_continuous(limits = y_range) +
      scale_color_manual(values = col_palette) +
      labs(y = "", x = "", subtitle = label_parts[j], title = label_enrich[j]) +
      theme_classic(base_size = base_size) +
      theme(legend.position = "none", axis.text.y = label_y[[j]], axis.text.x = label_x[[j]])
    
  })
  
  # create dummy plot to grab legend
  gg_dummy <- ggplot(data = response_rations) + 
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
  gg_temp <- cowplot::plot_grid(plotlist = gg_treatments, nrow = 3, ncol = 3) + 
    cowplot::draw_label(label = x_lab, x = 0.5, y = 0, vjust = -0.5, angle = 0) + 
    cowplot::draw_label(label = y_lab, x = 0, y = 0.5, vjust = 1.5, angle = 90)
  
  # add legend
  plot_grid(gg_temp, legend, rel_heights = c(0.9, 0.1), nrow = 2)
  
})

names(gg_comparisons) <- comparisons

