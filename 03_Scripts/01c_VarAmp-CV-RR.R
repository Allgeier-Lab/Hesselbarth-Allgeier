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

paths <- list(local = "02_Data/01_VarAmp-CV-local.rds", 
              mobile = "02_Data/01_VarAmp-CV-mobile.rds")

#### Bootstrap RR ####

control <- readr::read_rds("02_Data/01_VarAmp-CV-local.rds") %>% 
  tibble::add_column(id_row = 1:nrow(.), .before = "enrichment_lvl")

treatment <- readr::read_rds("02_Data/01_VarAmp-CV-mobile.rds") %>% 
  tibble::add_column(id_row = 1:nrow(.), .before = "enrichment_lvl")
  
treatments_list <- dplyr::inner_join(x = control, y = treatment,
                                     by = c("id_row", "enrichment_lvl", "amplitude_lvl", "n_diff", 
                                            "type", "part", "measure"), 
                                     suffix = c(".ctrl", ".trtm")) %>% 
  dplyr::filter(type == "cv", measure == "gamma") %>% 
  dplyr::group_by(enrichment_lvl, amplitude_lvl, n_diff, part) %>% 
  dplyr::group_split()
  
response_ratios <- purrr::map_dfr(seq_along(treatments_list), function(i) {
  
  df_temp <- treatments_list[[i]]
  
  message("\r> Progres: i = ", i , "/", length(treatments_list), "\t\t", appendLF = FALSE)
      
  if (all(df_temp$value.ctrl == 0) || all(df_temp$value.trtm == 0)) {
    
    tibble(enrichment_lvl = unique(df_temp$enrichment_lvl), amplitude_lvl = unique(df_temp$amplitude_lvl), 
           n_diff = unique(df_temp$n_diff), part = unique(df_temp$part), mean = as.numeric(NA), 
           lo = as.numeric(NA), hi = as.numeric(NA))
    
  } else {
        
    bootstrap <- boot::boot(data = data.frame(ctrl = df_temp$value.ctrl, 
                                              trtm = df_temp$value.trtm),
                            statistic = log_response, R = 10000, relative = TRUE)
      
    bootstrap_ci <- boot::boot.ci(bootstrap, type = "norm", conf = 0.95)
      
    tibble(enrichment_lvl = unique(df_temp$enrichment_lvl), amplitude_lvl = unique(df_temp$amplitude_lvl), 
           n_diff = unique(df_temp$n_diff), part = unique(df_temp$part),
           mean = mean(bootstrap$t[, 1]), lo = bootstrap_ci$normal[2], hi = bootstrap_ci$normal[3])
      
  }}) %>% 
  dplyr::mutate(enrichment_lvl = factor(enrichment_lvl, levels = c(0.75, 1.0, 1.25), 
                                        labels = c("low", "medium", "high")), 
                amplitude_lvl = factor(amplitude_lvl, levels = c(0.05, 0.5, 1.0), 
                                       labels = c("low", "medium", "high")),
                n_diff = as.integer(n_diff), 
                part = factor(part, levels = c("ag_biomass", "ag_production", "ttl_biomass",
                                               "bg_biomass", "bg_production", "ttl_production")))

#### Create ggplot ####

parts <- rep(x = c("ag_production", "bg_production", "ttl_production"), each = 3)

enrichments <- rep(x = c("low", "medium", "high"), times = 3)

# setup labels
label_parts <- c("Aboveground production", "" , "", 
                 "Belowground production", "", "", 
                 "Total production", "", "")

label_enrich <- c("Low enrichment", "Medium enrichment" , "High enrichment",
                  "", "", "", 
                  "", "", "")

label_y <- list(NULL, element_blank(), element_blank(),
                NULL, element_blank(), element_blank(),
                NULL, element_blank(), element_blank())

label_x <- list(element_blank(), element_blank(), element_blank(),
                element_blank(), element_blank(), element_blank(), 
                NULL, NULL, NULL)

# setup some plot options
col_palette <- c("#5ABCD6", "#FAD510", "#F22301")

# create plots in loop

# loop through all parts and treatments
gg_treatments <- purrr::map(seq_along(parts), function(i) {
  
  df_temp <- dplyr::filter(response_ratios, enrichment_lvl == enrichments[i], part == parts[i])
  
  y_range <- dplyr::filter(response_ratios, part == parts[i]) %>%
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
    labs(y = "", x = "", subtitle = label_parts[i], title = label_enrich[i]) +
    theme_classic(base_size = base_size) +
    theme(legend.position = "none", axis.text.y = label_y[[i]], axis.text.x = label_x[[i]], 
          plot.title = element_text(size = base_size), plot.subtitle = element_text(size = base_size))
  
})
  
# create dummy plot to grab legend
gg_dummy <- ggplot(data = response_ratios) + 
  geom_point(aes(x = 1, y = 1, col = amplitude_lvl)) +
  scale_color_manual(name = "Amplitude treatment", values = col_palette) + 
  theme(legend.position = "bottom")
  
# grab legend
gg_legend <- cowplot::get_legend(gg_dummy)

# create x/y labs
x_lab <- "Variability of amplitude"

y_lab <- expr(paste("log(response ratios) CV", gamma, " local vs. mobile"))

# add axsis labels
gg_treatments <- cowplot::plot_grid(plotlist = gg_treatments, nrow = 3, ncol = 3) + 
  cowplot::draw_label(label = x_lab, x = 0.5, y = 0, vjust = -0.5, angle = 0) + 
  cowplot::draw_label(label = y_lab, x = 0, y = 0.5, vjust = 1.5, angle = 90)

# add legend
gg_treatments <- plot_grid(gg_treatments, gg_legend, rel_heights = c(0.9, 0.1), nrow = 2)
