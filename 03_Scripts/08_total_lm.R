##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose:

#### Load setup ####

source("05_Various/setup.R")

amplitude <- "095"

file_path <- c(phase = paste0("02_Data/05-variability-phase-", amplitude, ".rds"), 
               noise = paste0("02_Data/05-variability-noise-", amplitude, ".rds"))

df_cv_prod <- purrr::map_dfr(file_path, function(path_i) {
  readr::read_rds(path_i) %>% 
    purrr::map_dfr(function(result_j){
      
      abiotic_sd <- names(result_j$cv)[[6]]
      
      df_temp <- dplyr::left_join(x = result_j$cv, y = result_j$production, 
                       by = c("part", "pop_n", "move_meta_sd", abiotic_sd), 
                       suffix = c(".cv", ".prod")) %>% 
        dplyr::mutate(value.move = mean(result_j$moved$moved, na.rm = TRUE))
      
      names(df_temp)[[6]] <- "abiotic_sd"
      
      return(df_temp)})}, .id = "abiotic_trtmnt") %>%
  dplyr::mutate(abiotic_trtmnt = factor(abiotic_trtmnt), 
                part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass",
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                pop_n = factor(as.numeric(pop_n), ordered = TRUE)) %>% 
  dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), 
                measure %in% c("alpha", "gamma")) %>% 
  tibble::tibble()

summary_lm_ <- dplyr::group_by(df_cv_prod, part, measure) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(df_temp) {
    
    lm_temp <- dplyr::mutate(df_temp, pop_n = as.numeric(pop_n)) %>%
      lm(formula = log(value.cv) ~ pop_n + move_meta_sd + abiotic_trtmnt + abiotic_sd, data = .)
    
    broom::tidy(lm_temp) %>% 
      dplyr::mutate(part = unique(df_temp$part), measure = unique(df_temp$measure), 
                    .before = term)}) %>% 
  dplyr::mutate(term = dplyr::case_when(term == "(Intercept)" ~ "Intercept", 
                                        TRUE ~ term), 
                p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                                 p.value < 0.05 ~ "*", p.value >= 0.05 ~ "n.s."))

summary_lm <- dplyr::group_by(df_cv_prod, abiotic_trtmnt, part, measure) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(df_temp) {
    
    lm_temp <- dplyr::mutate(df_temp, pop_n = as.numeric(pop_n)) %>%
      lm(formula = value.cv ~ pop_n + move_meta_sd + abiotic_sd, data = .)
    
    broom::tidy(lm_temp) %>% 
      dplyr::mutate(abiotic_trtmnt = unique(df_temp$abiotic_trtmnt), part = unique(df_temp$part), 
                    measure = unique(df_temp$measure), 
                    .before = term)}) %>% 
  dplyr::mutate(term = dplyr::case_when(term == "(Intercept)" ~ "Intercept", 
                                        TRUE ~ term), 
                p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                                 p.value < 0.05 ~ "*", p.value >= 0.05 ~ "n.s."))


