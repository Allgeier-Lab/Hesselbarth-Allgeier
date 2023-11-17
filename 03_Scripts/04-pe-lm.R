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
  dplyr::mutate(part = factor(part, levels = c("ag_production", "bg_production", "ttl_production"),
                              labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                         "ttl_production" = "Total")),
                treatment = factor(treatment, levels = c("subsidies", "connectivity", "combined")))

#### Means of raw ####

cv_means <- dplyr::filter(results_final_df, measure != "synchrony", part != "Belowground", 
                          treatment == "combined") |> 
  dplyr::mutate(facet = dplyr::case_when(measure == "beta" ~ "PE", TRUE ~ "cv")) |> 
  dplyr::group_by(facet, part, measure, nutrient_input) |> 
  dplyr::summarise(mn = mean(value.cv), sd = sd(value.cv), .groups = "drop") |> 
  dplyr::mutate(measure = factor(measure, levels = c("alpha", "gamma", "beta"), 
                                 labels = c("Local scale", "Meta-ecosystem scale", "Portfolio effect")))

cv_means_nf <- dplyr::filter(results_final_df, measure != "synchrony", part != "Belowground", 
              treatment == "subsidies") |> 
  dplyr::mutate(facet = dplyr::case_when(measure == "beta" ~ "PE", TRUE ~ "cv")) |> 
  dplyr::group_by(facet, part, measure, nutrient_input) |> 
  dplyr::summarise(mn = mean(value.cv), sd = sd(value.cv), .groups = "drop") |> 
  dplyr::mutate(measure = factor(measure, levels = c("alpha", "gamma", "beta"), 
                                 labels = c("Local scale", "Meta-ecosystem scale", "Portfolio effect")))

#### Split data #### 

results_final_list <- dplyr::filter(results_final_df, measure != "synchrony", 
                                    part %in% c("Aboveground", "Total"), treatment == "combined") |>  
  dplyr::group_by(part, measure) |>
  dplyr::group_split()

names_list <- purrr::map(results_final_list, function(i) as.character(c(unique(i$part), unique(i$measure))))

#### Fit regression model and dredge ####

best_lm_list <- vector(mode = "list", length = length(results_final_list))

for (i in 1:length(results_final_list)) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  df_temp <- results_final_list[[i]] |> 
    dplyr::mutate(value.cv = log(value.cv), abiotic = log(abiotic), biotic = log(biotic)) |>
    dplyr::mutate(value.cv = (value.cv - mean(value.cv)) / sd(value.cv),
                  biotic = (biotic - mean(biotic)) / sd(biotic),
                  abiotic = (abiotic - mean(abiotic)) / sd(abiotic))
  
  lm_temp <- lm(value.cv ~ nutrient_input + abiotic + pop_n + biotic + 
                  nutrient_input*pop_n + nutrient_input*biotic + abiotic*pop_n + abiotic*biotic,
                data = df_temp, na.action = "na.fail")
  
  # nutrient_input*pop_n THIS IS THE ISSUE

  lm_dredge <- MuMIn::dredge(lm_temp, extra = c("R^2"))

  best_lm_list[[i]] <- get.models(lm_dredge, subset = 1)
  
}

#### Relative importance ####

# rel_importance_df <- purrr::map_dfr(seq_along(results_final_list), function(i) {
# 
#   message("> Progress: ", i, "/", length(results_final_list))
# 
#   rel_r2 <- relaimpo::boot.relimp(best_lm_list[[i]][[1]], type = "lmg", b = 500, level = 0.95,
#                                   fixed = FALSE) |>
#     relaimpo::booteval.relimp(bty = "basic")
# 
#   tibble::tibble(
#     part = names_list[[i]][[1]], measure = names_list[[i]][[2]],
#     beta = c(rel_r2@namen[-1], "residual"), mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)),
#     lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA)) |>
#     dplyr::mutate(dplyr::across(mean:higher, ~ .x * 100),
#                   lower = dplyr::case_when(lower < 0 ~ 0.0, TRUE ~ lower))}) |>
#   dplyr::mutate(measure = factor(measure, levels = c("alpha", "gamma" , "beta"),
#                                  labels = c("alpha" = "Local scale", "gamma" = "Meta-ecosystem scale",
#                                             "beta" = "Portfolio effect")),
#                 beta = factor(beta, levels = c("abiotic", "biotic", "nutrient_input", "pop_n",
#                                                "nutrient_input:pop_n", "biotic:nutrient_input",
#                                                "abiotic:pop_n", "abiotic:biotic", "residual"),
#                               labels = c("abiotic" = "Enrich. var.", "biotic" = "Connectivity",
#                                          "nutrient_input" = "Enrichment", "pop_n" = "Pop. size",
#                                          "nutrient_input:pop_n" = "Enrich:Pop", "biotic:nutrient_input" = "Enrich:Con",
#                                          "abiotic:pop_n" = "Var:Pop", "abiotic:biotic" = "Var:Con",
#                                          "residual" = "Residuals")))
# 
# suppoRt::save_rds(object = rel_importance_df, filename = "rel-importance-df.rds",
#                   path = "02_Data/", overwrite = FALSE)

rel_importance_df <- readr::read_rds("02_Data/rel-importance-df.rds")

#### Marginal means ####

marginal_means <- purrr::map_dfr(seq_along(results_final_list), function(i) {
  
  message("> Progress: ", i, "/", length(results_final_list))
  
  modelbased::estimate_means(best_lm_list[[i]][[1]], at = c("pop_n", "nutrient_input"), 
                             ci = 0.95) |> 
    tibble::as_tibble() |> 
    dplyr::mutate(part = names_list[[i]][[1]], measure = names_list[[i]][[2]], .before = "pop_n")}) |> 
  dplyr::mutate(measure = factor(measure, levels = c("alpha", "gamma" , "beta"), 
                                 labels = c("alpha" = "Local scale", "gamma" = "Meta-ecosystem scale", 
                                            "beta" = "Portfolio effect")))

#### Create dummy plot #### 

base_size <- 13.0

w <- 1.0

text_size <- 5

color_enrich <- c("low" = "#009E73", "medium" = "#F0E442", "high" = "#CC79A7")

#### Create ggplot ####

gg_importance_part <- dplyr::filter(rel_importance_df, beta != "Residuals") |> 
  dplyr::group_by(part) |> 
  dplyr::group_split() |> 
  purrr::map(function(i) {
    
    label_temp <- data.frame(label = "i)", x = 1, y = 60)
    lgd_pos <- c(0.8, 0.9)
    y_just <- 0.1
    y_push <- 1.5
    y_temp <- "Relative importance [%]" # ifelse(test = i == "Aboveground", yes = "Stability value", no = "")
    
    if(unique(i$part) == "Total") {
      
      label_temp <- data.frame(label = "ii)", x = 1, y = 60)
      lgd_pos <- "none"
      y_just <- 0.35
      y_temp <- ""
      y_push <- 4.25
      
    }
    
    ggplot(data = i) +
      
      # add error bars
      geom_errorbar(aes(x = beta, ymin = lower, ymax = higher, group = measure),
                    position = position_dodge(w * 0.75), width = 0, linewidth = 0.5) +
      
      # relative importance bars
      geom_col_pattern(aes(x = beta, y = mean, pattern = measure), pattern_density = 0.01, pattern_spacing = 0.02,
                       pattern_color = "grey", position = position_dodge(w * 0.75), width = w * 0.7) +
      
      geom_text(data = label_temp, aes(x = x, y = y, label = label), 
                size = text_size, nudge_x = -0.25) +
      
      # adding box
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set scales
      scale_pattern_manual(name = "", values = c("none", "stripe", "circle"),
                           labels = c("Local scale", "Meta-Ecosys.", "Portfolio effect")) + 
      scale_y_continuous(limits = c(0, max(rel_importance_df$higher))) +
      
      # labels and themes
      guides(pattern = guide_legend(nrow = 3, override.aes = list(color = "black"))) +
      labs(y = y_temp, x = "") +
      theme_classic(base_size = base_size) +
      theme(legend.position = lgd_pos, legend.background = element_blank(), legend.key.size = unit(0.5, "cm"),
            axis.line = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y = element_text(hjust = y_just, margin = unit(c(0, y_push, 0, 0), "mm")), 
            strip.text = element_blank())

  })

gg_importance <- cowplot::plot_grid(plotlist = gg_importance_part, nrow = 1)

# figure raw values
gg_raw_part <- purrr::map(c("Aboveground", "Total"), function(i) {

  df_temp <- dplyr::filter(cv_means, part == i)

  y_temp <- ifelse(test = i == "Aboveground", yes = "Stability value", no = "")

  label_greek <- data.frame(facet = c("cv", "PE"), label = c("i)", "ii)"),
                            x = c(1, 1), y = c(0.325, 7.15))

  lgd_pos <- c(0.325, 0.9)
  
  nudge_temp <- c(-0.4, -0.45)

  if (i != "Aboveground") {

    label_greek <- data.frame(facet = c("cv", "PE"), label = c("iii)", "iv)"),
                              x = c(1, 1), y = c(0.04, 9.2))

    lgd_pos <- "none"
    
    nudge_temp <- c(-0.3, -0.45)

  }

  # filter part
  ggplot(data = df_temp) +

    # mean stability values bars
    geom_col(aes(x = measure, y = mn, fill = nutrient_input), position = position_dodge(w * 0.75),
           width = w * 0.70, color = NA) +

    geom_errorbar(aes(x = measure, ymax = mn + sd, ymin = mn * 0.9, color = nutrient_input),
                position = position_dodge(w * 0.75), width = 0.0, linewidth = 0.5) +

    # mean and error of nutrients only model
    geom_point(data = dplyr::filter(cv_means_nf, part == i),
               aes(x = measure, y = mn, group = nutrient_input), position = position_dodge(w * 0.75),
               color = "black", shape = 8, size = 2.5) +

    # adding facet labels
    geom_text(data = label_greek, aes(x = x , y = y, label = label), size = text_size,
              nudge_x = nudge_temp) +

    # facet wrap by alpha/gamma and PE
    facet_wrap(. ~ facet, scales = "free", ncol = 2) +

    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +

    # set scales
    scale_x_discrete(labels = c("Local scale" = "Local", "Meta-ecosystem scale" = "Meta-Ecosys.",
                                "Portfolio effect" = "Portfolio effect")) +

    scale_fill_manual(name = "", values = color_enrich) +
    scale_color_manual(name = "", values = color_enrich) +

    # change theme
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "", y = y_temp) +
    theme_classic(base_size = base_size) +
    theme(legend.position = lgd_pos, legend.background = element_blank(),
          legend.text = element_text(size = base_size * 0.85), legend.key.size = unit(0.35, "cm"),
          strip.text = element_blank(), axis.line = element_blank())
})

gg_raw <- cowplot::plot_grid(plotlist = gg_raw_part, nrow = 1)

# figure marginal means
gg_means_part <- dplyr::filter(marginal_means, measure == "Portfolio effect") |> 
  dplyr::group_by(part) |> 
  dplyr::group_split() |> 
  purrr::map(function(i) {
    
    y_max <- exp(max(marginal_means$CI_high, na.rm = TRUE))
    
    label_temp <- data.frame(label = "i)", x = 1, y = y_max * 1)
    lgd_pos <- c(0.125, 0.8)
    y_temp <- "Portfolio effect (PE)"
    y_push <- 2.5
    
    if(unique(i$part) == "Total") {
      
      label_temp <- data.frame(label = "ii)", x = 1, y = y_max * 1)
      lgd_pos <- "none"
      y_temp <- ""
      y_push <- 5.15
      
    }
    
    dplyr::mutate(i, nutrient_input = factor(nutrient_input, ordered = TRUE, 
                                             labels = c("low", "medium", "high"))) |> 
      ggplot() +
     
      # adding mean and error bars
      geom_point(aes(x = pop_n, y = exp(Mean), color = nutrient_input), size = 3.5) + 
      geom_errorbar(aes(x = pop_n, ymin = exp(CI_low), ymax = exp(CI_high), color = nutrient_input), 
                    width = 0, linewidth = 1.5) +
      
      geom_line(aes(x = pop_n, y = exp(Mean), color = nutrient_input, group = nutrient_input), 
                alpha = 0.25, linewidth = 1.5) + 
      
      # adding facet labels
      
      geom_text(data = label_temp, aes(x = x, y = y, label = label),
                size = text_size, nudge_x = -0.35) +
      
      # adding box
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
      annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
      
      # set colors
      scale_color_manual(name = "", values = color_enrich) +
      scale_y_continuous(limits = c(0.0, y_max), breaks = seq(0, ceiling(y_max), 1)) +
      guides(color = guide_legend(nrow = 3)) +
      
      # change theme
      labs(x = "", y = y_temp) +
      theme_classic(base_size = base_size) +
      theme(legend.position = lgd_pos, legend.background = element_blank(),
            legend.text = element_text(size = base_size * 0.85), legend.key.size = unit(0.35, "cm"),
            axis.line = element_blank(), 
            axis.title.y = element_text(hjust = 0.0, margin = unit(c(0, y_push, 0, 0), "mm")))
    
  })

gg_means <- cowplot::plot_grid(plotlist = gg_means_part, ncol = 2)

gg_means <- cowplot::ggdraw(gg_means) +
  cowplot::draw_label("Population size", x = 0.3, y = 0.05, size = base_size) + 
  cowplot::draw_label("Population size", x = 0.8, y = 0.05, size = base_size)

gg_final <- cowplot::plot_grid(gg_importance, gg_raw, gg_means, nrow = 3, rel_heights = c(0.4, 0.3, 0.3), 
                               labels = c("A)", "B)", "C)"), label_fontface = "plain", 
                               vjust = c(1.5, 1.5, 1.0)) |> 
  cowplot::ggdraw(ylim = c(0.0, 1.025)) +
  cowplot::draw_line(x = c(0.52, 0.52), y = c(0.05, 1.025 - 0.05), color = "grey") +
  cowplot::draw_label("Aboveground", x = 0.275, y = 1.01, size = base_size) + 
  cowplot::draw_label("Total", x = 0.8, y = 1.01, size = base_size)

#### Save table and ggplot #### 

suppoRt::save_ggplot(plot = gg_final, filename = "Figure-3.png",
                     path = "04_Figures/", width = width, height = height * 0.75,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_final, filename = "Figure-3.pdf",
                     path = "04_Figures/", width = width, height = height * 0.75,
                     units = units, dpi = dpi, overwrite = FALSE)
