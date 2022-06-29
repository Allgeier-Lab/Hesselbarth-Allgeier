library(ecofolio)

# source setup file
source("05_Various/setup.R")

# create own function based on ecofolio plot_mv
plot_port <- function(x, title = "") {

  # calculate alpha scale mean
  mean_alpha <- apply(X = x, FUN = function(i) mean(i), MARGIN = 2)

  # calculate alpha scale var
  var_alpha <- apply(X = x, FUN = function(i) var(i), MARGIN = 2)

  # combine to one data.frame
  df_alpha <- data.frame(mn = mean_alpha, vr = var_alpha)

  # calculate gamma scale sum of each timestep
  sum_gamma <- apply(X = x, FUN = sum, MARGIN = 1)

  # calculate mean and variance of gamma scale
  df_gamma <- data.frame(mn = mean(sum_gamma), vr = var(sum_gamma))

  # linear regression of variance and mean (log-log)
  lm_gamma <- lm(log(vr) ~ log(mn), data = df_alpha)
  # lm_gamma <- robustbase::lmrob(log(vr) ~ log(mn), data = df_alpha)
  
  # predict variance
  df_pred <- data.frame(mn = c(mean_alpha, gamma = df_gamma$mn),
                        vr = exp(predict(lm_gamma, newdata = data.frame(mn = c(mean_alpha, sum = df_gamma$mn)))), 
                        row.names = 1:(ncol(x) + 1))

  # get id of sum
  id_gamma <- ncol(x) + 1
  
  # calculate portfolio effect
  pe <- (sqrt(df_pred$vr[id_gamma]) / df_pred$mn[id_gamma]) / (sqrt(df_gamma$vr) / df_gamma$mn)

  # combine to one data.frame
  df_gamma_ttl <- cbind(type = rbind("observed", "predicted"),
                        rbind(df_gamma, df_pred[id_gamma, ]))
  
  y_pred <- ifelse(test = df_gamma_ttl[df_gamma_ttl$type == "predicted", "vr"] < max(df_alpha$vr),
                   yes = min(df_pred$vr[-id_gamma]), no = max(df_pred$vr[-id_gamma]))
  
  # MH: Issue in the geom_segment blue if slope is negative!
  ggplot() +
    geom_line(data = df_pred[-ncol(x) - 1, ], aes(x = mn, y = vr), col = "blue") +
    geom_segment(aes(x = max(df_pred$mn[-id_gamma]), xend = df_gamma_ttl$mn[2],
                     y = y_pred, yend = df_gamma_ttl$vr[2]),
                 linetype = 2, col = "blue") +
    geom_segment(aes(x = df_gamma_ttl$mn[1], xend = df_gamma_ttl$mn[2],
                     y = df_gamma_ttl$vr[[1]], yend = df_gamma_ttl$vr[[2]]),
                 linetype = 2) +
    geom_point(data = df_alpha, aes(x = mn, y = vr)) +
    geom_point(data = df_gamma_ttl, aes(x = mn, y = vr, col = type), pch = 4, size = 2.5) +
    geom_text(aes(x = min(df_alpha$mn) * 1.25, y = max(df_gamma_ttl$vr),
                  label = paste0("pe=",round(pe, 5)))) +
    scale_color_manual(name = "", values = c("black", "red")) +
    scale_y_log10() + scale_x_log10() +
    labs(x = "log(mean)", y = "log(variance)", title = title) +
    theme_classic() +
    theme(legend.position = "bottom")

}

# load parameters
list_starting <- readRDS("02_Data/list_starting.rds")

list_parameters <- readRDS("02_Data/list_parameters.rds")

#### Adapt parameters ####

list_starting$pop_n <- 0

list_parameters$nutrients_diffusion <- 0.0
list_parameters$detritus_diffusion <- 0.0
list_parameters$detritus_fish_diffusion <- 0.0

#### Stable values ####

stable_values <- arrR::get_stable_values(bg_biomass = list_starting$bg_biomass,
                                         ag_biomass = list_starting$ag_biomass,
                                         parameters = list_parameters)

list_starting$nutrients_pool <- stable_values$nutrients_pool

list_starting$detritus_pool <- stable_values$detritus_pool

#### Setup input and metaecosyst ####

# simulate input
input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                        input_mn = stable_values$nutr_input, freq_mn = freq_mn,
                                        amplitude_mod = seq(from = 0.01, to = 1.0, length.out = n),
                                        phase_mod = runif(n = n, min = 0.0, max = 0.0))

plot(input_temp, gamma = FALSE) +
  geom_hline(yintercept = stable_values$nutr_input, linetype = 2)

# setup metaecosystems
metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i,
                                       starting_values = list_starting,
                                       parameters = list_parameters,
                                       dimensions = dimensions, grain = grain,
                                       reef = NULL)

# run model
result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutr_input = input_temp,
                                   parameters = list_parameters,
                                   max_i = max_i, min_per_i = min_per_i,
                                   seagrass_each = seagrass_each, save_each = save_each)

# result_temp <- readRDS("02_Data/result_temp.rds")
# input_temp <- result_temp$nutr_input

# result_temp <- meta.arrR::filter_meta(x = result_temp, filter = c(max_i / 2, max_i))

plot(result_temp, summarize = TRUE)
plot_meta_production(result_temp, lag = TRUE)

result_sum <- purrr::map_dfr(result_temp$seafloor, function(i) {
  
  dplyr::group_by(i, timestep) %>% 
    dplyr::summarise(bg_biomass = sum(bg_biomass), ag_biomass = sum(ag_biomass),
                     bg_production = sum(bg_production), ag_production = sum(ag_production), 
                     .groups = "drop") %>% 
    dplyr::mutate(ttl_biomass = bg_biomass + ag_biomass, 
                  ttl_production = bg_production + ag_production)}, .id = "meta")

# ggplot(data = results_sum) + 
#   geom_line(aes(x = timestep, y = ttl_biomass, col = meta)) + 
#   scale_color_viridis_d(name = "", option = "A") + 
#   theme_classic() + theme(legend.position = "bottom")

# result_prod <- meta.arrR::get_meta_production(result_temp, lag = TRUE)

#### PE ####

col_name <- "ag_biomass"

df_input <- dplyr::select(result_sum, meta, timestep, all_of(col_name)) %>%
  tidyr::pivot_wider(values_from = all_of(col_name), names_from = meta) %>% 
  as.data.frame()

# plot(log(apply(X = df_input[, -1], FUN = mean, MARGIN = 2)), 
#      log(apply(X = df_input[, -1], FUN = var, MARGIN = 2)))

# df_input <- dplyr::filter(result_prod, part == "bg_production") %>% 
#   dplyr::select(-part) %>% 
#   dplyr::mutate(meta = paste0("Meta_", meta)) %>% 
#   tidyr::pivot_wider(values_from = value, names_from = meta) %>% 
#   dplyr::slice(-1) %>% 
#   as.data.frame()

# df_input <- meta.arrR::get_input_df(x = input_temp, gamma = FALSE)

# mean(apply(X = df_input[, -1], FUN = terra::cv, MARGIN = 2)) /
#   terra::cv(apply(X = df_input[, -1], FUN = sum, MARGIN = 1))
# 
# calc_variability(result_temp, lag = FALSE)

ecofolio::pe_avg_cv(x = df_input[, -1])
ecofolio::pe_mv(x = df_input[, -1])

plot_port(x = df_input[, -1], title = col_name)
# ecofolio::plot_mv(df_input[, -1], show = c("linear", "robust"))

