##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-cv.R")
source("01_Functions/import-mortality.R")
source("01_Functions/plot-lm.R")
source("01_Functions/cut-borders.R")

#### Load/wrangle simulated data ####

results_noise_df <- import_cv(path = "02_Data/result-noise.rds")

#### Load mortality data ####

mortality_noise_df <- import_mortality(path = "02_Data/result-noise.rds") |> 
  dplyr::filter(pop_n != 0) |> 
  dplyr::group_by(row_id, pop_n, nutrient_input) |>
  dplyr::summarise(died_total = mean(died_total), .groups = "drop")

#### Filter abundances/mortality #### 

results_final_df <- dplyr::left_join(x = results_noise_df, y = mortality_noise_df, 
                                     by = c("row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(died_total > threshold_mort ~ "no", TRUE ~ "yes")) |> 
  dplyr::filter(include == "yes") |> 
  dplyr::mutate(treatment = dplyr::case_when(abiotic != 0.0 & biotic == 0.0 ~ "subsidies", 
                                             abiotic == 0.0 & biotic != 0.0 ~ "connectivity", 
                                             abiotic != 0.0 & biotic != 0.0 ~ "combined")) |>
  dplyr::mutate(part = factor(part, labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                               "ttl_productio" = "Total")),
                treatment = factor(treatment, levels = c("subsidies", "connectivity", "combined")))

#### Join PE and PP values ####

results_pe_pp_df <- dplyr::left_join(x = dplyr::filter(results_final_df, measure != "synchrony", treatment == "combined", 
                                                       part != "Belowground") |> 
                                       dplyr::select(measure, row_id, part, pop_n, nutrient_input, biotic, abiotic, value.cv),
                                     y = dplyr::filter(results_final_df, measure == "gamma", treatment == "combined", 
                                                       part != "Belowground") |> 
                                       dplyr::select(row_id, part, pop_n, nutrient_input, biotic, abiotic, value.prod), 
                                     by = c("row_id", "part", "pop_n", "nutrient_input", "biotic", "abiotic")) |> 
  dplyr::mutate(measure = factor(measure, levels = c("alpha", "gamma" , "beta"), 
                                 labels = c("alpha" = "Local scale", "gamma" = "Meta-ecosystem scale", 
                                            "beta" = "Portfolio effect")))

#### Split data #### 

results_final_list <- dplyr::group_by(results_pe_pp_df, part, measure) |>
  dplyr::group_split()

names_list <- purrr::map(results_final_list, function(i) as.character(c(unique(i$part), unique(i$measure))))

#### Fit regression model and dredge ####

best_lm_list <- vector(mode = "list", length = length(results_final_list))

lm_summary_list <- vector(mode = "list", length = length(results_final_list))

for(i in 1:length(results_final_list)) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  df_temp <- results_final_list[[i]]
  
  if (unique(df_temp$part) == "Aboveground") {
    
    df_temp$value.prod <- sqrt(df_temp$value.prod)
    
  } else if (unique(df_temp$part) == "Total") {
    
    df_temp$value.prod <- log(df_temp$value.prod)
    
  }
  
  df_temp <- dplyr::mutate(df_temp, value.cv = (value.cv - mean(value.cv)) / sd(value.cv),
                           biotic = (biotic - mean(biotic)) / sd(biotic),
                           abiotic = (abiotic - mean(abiotic)) / sd(abiotic))
  
  lm_temp <- lm(value.prod ~ value.cv + abiotic + biotic + nutrient_input + pop_n,
                data = df_temp, na.action = "na.fail")
  
  lm_dredge <- MuMIn::dredge(lm_temp, extra = c("R^2"))
  
  lm_summary_list[[i]] <- subset(lm_dredge, subset = 1:3) |>
    tibble::as_tibble() |>
    dplyr::mutate(part = names_list[[i]][[1]], measure = names_list[[i]][[2]], .before = "(Intercept)")
  
  best_lm_list[[i]] <- get.models(lm_dredge, subset = 1)
  
}

# convert to data.frame
lm_summary_df <- dplyr::bind_rows(lm_summary_list) |> 
  dplyr::rename("intercept" = "(Intercept)", "spatialVariation" = "abiotic", 
                "connectivity" = "biotic", "enrichment" = "nutrient_input", 
                "popSize" = "pop_n", "portfolioEffect" = "value.cv") |> 
  dplyr::mutate_if(is.numeric, round, digits = 5)

#### Rel imp ####

rel_importance_df <- purrr::map_dfr(seq_along(best_lm_list), function(i) {
  
  message("> Progress: ", i, "/", length(best_lm_list))
  
  rel_r2 <- relaimpo::boot.relimp(best_lm_list[[i]][[1]], type = "lmg", b = 500, level = 0.95, 
                                  fixed = FALSE) |>
    relaimpo::booteval.relimp(bty = "basic")
  
  tibble::tibble(
    part = names_list[[i]][[1]], measure = names_list[[i]][[2]],
    beta = c(rel_r2@namen[-1], "residual"), mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)),
    lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA)) |>
    dplyr::mutate(dplyr::across(mean:higher, ~ .x * 100), 
                  lower = dplyr::case_when(lower < 0 ~ 0.0, TRUE ~ lower),
                  beta = factor(beta, levels = c("abiotic", "biotic", "nutrient_input", "pop_n", "value.cv",
                                                 "residual"),
                                labels = c("abiotic" = "Enrich. var.", "biotic" = "Connectivity", 
                                           "nutrient_input" = "Enrichment", "pop_n" = "Pop. size",  
                                           "value.cv" = "Stability", "residual" = "Residuals"))) |> 
                    tidyr::complete(beta, fill = list(part = names_list[[i]][[1]], measure = names_list[[i]][[2]],
                                                      mean = 0, lower = 0, higher = 0))})

#### Setup ggplot ####

base_size <- 8.0

#### Create ggplot ####

gg_part <- purrr::map(c(Aboveground = "Aboveground", Total = "Total"), function(part_i) {

  gg_scale <- purrr::map(c(alpha = "Local scale", gamma = "Meta-ecosystem scale", beta = "Portfolio effect"), function(scale_i) {
    
    df_sum <- dplyr::filter(results_pe_pp_df, measure == scale_i, part == part_i) |> 
      dplyr::mutate(pe_class = cut(value.cv, breaks = seq(from = min(value.cv), 
                                                          to = max(value.cv), length.out = 20), 
                                   include.lowest = TRUE)) |> 
      dplyr::group_by(pe_class) |> 
      dplyr::summarise(value.prod = mean(value.prod), n = dplyr::n(), .groups = "drop")
    
    # label_temp <- paste0(
    # ifelse(test = part_i == "Aboveground", yes = "A", no = "B"),
    # ifelse(test = scale_i == "Local scale", yes = ".1)", no = ifelse(test = scale_i == "Meta-ecosystem scale", 
    #                                                                 yes = ".2)", no = ".3)"))
    # )
    
    if (scale_i == "Local scale") label_temp <- "i)"
    if (scale_i == "Meta-ecosystem scale") label_temp <- "ii)"
    if (scale_i == "Portfolio effect") label_temp <- "iii)"
  
    breaks_temp <- df_sum$pe_class[seq(2, length(df_sum$pe_class), 4)]
    
    breaks_label <- cut_borders(breaks_temp) |> 
      apply(MARGIN = 1, FUN = function(i) sprintf("%.2f", median(i)))
      # apply(MARGIN = 2, FUN = function(i) sprintf("%.3f", round(i, digits = 3)))
    
    # breaks_label <- paste0("(", breaks_label[, 1], "-", breaks_label[, 2], "]")
    
    rel_imp_temp <- dplyr::filter(rel_importance_df, part == part_i, measure == scale_i, 
                                  beta != "Residuals")
    
    breaks_col <- seq(min(df_sum$value.prod) * 1.1, max(df_sum$value.prod) * 0.95, length.out = 2)
    
    gg_pp <- ggplot(data = df_sum) +
      
      # adding tiles with PP
      geom_col(aes(x = pe_class, y = n, fill = value.prod)) +
      
      # # adding facet labels
      annotate(geom = "text", x = 1, y = max(df_sum$n) * 1.03, label = label_temp, size = 2.5) +
      
      # adding box
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set scales 
      # scale_fill_gradient2(name = expression(paste("P.P.")),
      #                      low = "white", mid = "#fadea3", high = "#224b20") +
      scale_fill_gradientn(name = "", # name = expression(paste("PP [", gDW~d^-1~m^-2, "]")),
                           colours = MetBrewer::met.brewer(name = "VanGogh3", n = 255), 
                           breaks = breaks_col, labels = round(breaks_col, 1)) +
      scale_x_discrete(breaks = breaks_temp, labels = breaks_label) +
      
      # theming
      labs(x = "", y = "") +
      guides(fill = guide_colorbar(direction = "horizontal",  title.position = "top")) +
      theme_classic(base_size = base_size) + 
      theme(axis.line = element_blank(),
            legend.position = c(0.8, 0.4), legend.title = element_text(size = base_size * 0.5),
            legend.text = element_text(size = base_size * 0.65), legend.key.height = unit(3.5, "mm"),
            legend.key.size = unit(0.25, "cm"),
            legend.background = element_rect(colour = NA, fill = NA))
    
    gg_inset <- ggplot(data = rel_imp_temp, aes(x = beta, y = mean, group = beta)) + 
      
      geom_errorbar(aes(ymin = lower, ymax = higher, group = beta), 
                    width = 0, linewidth = 0.5, color = "grey25") +
      
      # add columns with rel imp
      geom_col(width = 0.75, color = "grey25") +
      
      # scale fill, color, and y-axis
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 100)) +
      
      # general theming
      coord_flip() +
      labs(y = "Rel. importance [%]") +
      theme_classic(base_size = base_size * 0.65) +
      theme(legend.position = "none", axis.title.y = element_blank(), 
            axis.text.y = element_text(size = base_size * 0.65))
    
    ggdraw(gg_pp) + 
      draw_plot(gg_inset, x = 0.495, y = 0.61, width = 0.475, height = 0.3525)
  
  })
  
 cowplot::plot_grid(plotlist = gg_scale, ncol = 3)
  
})

gg_pe_pp <- cowplot::plot_grid(gg_part$Aboveground, gg_part$Total, nrow = 2, 
                               labels = c("A)", "B)"), label_fontface = "plain", 
                               label_size = 12.0, hjust = 0.5) |> 
  cowplot::ggdraw(xlim = c(-0.035, 1.0), ylim = c(0.0, 1.0)) +
  cowplot::draw_label("AG Number of observations", x = 0.015, y = 0.275, angle = 90, size = base_size) + 
  cowplot::draw_label("TTL Number of observations", x = 0.015, y = 0.78, angle = 90, size = base_size) + 
  # cowplot::draw_label("AG", x = -0.02, y = 0.8, angle = 90, size = base_size) + 
  # cowplot::draw_label("TTL", x = -0.02, y = 0.275, angle = 90, size = base_size) + 
  cowplot::draw_label(expression(paste("Local scale (", italic(CV), italic(alpha), ")")), 
                      x = 0.175, y = 0.02, size = base_size) +
  cowplot::draw_label(expression(paste("Meta-ecosystem scale (", italic(CV), italic(gamma), ")")), 
                      x = 0.525, y = 0.02, size = base_size) +
  cowplot::draw_label("Portfolio effect (PE)", x = 0.875, y = 0.02, size = base_size)

#### Save ggplots ####

suppoRt::save_ggplot(plot = gg_pe_pp, filename = "Figure-4.png",
                     path = "04_Figures/", width = width, height = height * 1/3,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_pe_pp, filename = "Figure-4.pdf",
                     path = "04_Figures/", width = width, height = height * 1/3,
                     units = units, dpi = dpi, overwrite = FALSE)
