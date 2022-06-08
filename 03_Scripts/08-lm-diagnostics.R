##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Regression model diagnostics

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/model_assumptions.R")

extension <- ".pdf"

colors_pop <- c("8" = "#ebcc60", "16" = "#459b75", "32" = "#374a98", "64" = "#913f98")

#### Load/wrangle simulated data ####

file_names <- list.files(path = "02_Data/", pattern = "05-move-variability-*", 
                         full.names = TRUE)

names(file_names) <- stringr::str_sub(file_names, start = -10, end = -5)

df_cv_prod <- purrr::map_dfr(file_names, function(i) {
  readr::read_rds(i) %>% 
    purrr::map_dfr(function(j) {
      dplyr::left_join(x = j$cv, y = j$production, by = c("part", "pop_n", "move_meta_mean", 
                                                          "move_meta_sd"), 
                       suffix = c(".cv", ".prod")) %>% 
        dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))
    })
}, .id = "experiment") %>%
  tidyr::separate(col = experiment, sep = "-", into = c("amplitude", "cycles")) %>% 
  dplyr::mutate(amplitude = factor(amplitude, levels = c("005", "095"), labels = c("5%", "95%")), 
                cycles = factor(cycles, levels = c("50", "12"), labels = c("50", "12")),
                part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass",
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                pop_n = factor(as.numeric(pop_n), ordered = TRUE),
                value.cv.log = log(value.cv), value.prod.log = log(value.prod), 
                value.move.log = log(value.move)) %>% 
  dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), 
                measure %in% c("alpha", "gamma", "beta")) %>% 
  tibble::tibble()

list_assumptions <- dplyr::group_by(df_cv_prod, amplitude, cycles, part, measure, pop_n) %>% 
  dplyr::group_split() %>% 
  purrr::map(function(i) {
    
    title_temp <- paste0("Amp=", unique(i$amplitude), "; Cycles=", unique(i$cycles), 
                         "; Part=", unique(i$part), "; Scale=", unique(i$measure), 
                         "; Pop n=", unique(i$pop_n))
    
    lm_temp <- lm(value.prod.log ~ value.cv.log, data = i)
    
    model_assumptions(lm_temp, title = title_temp)
  })
