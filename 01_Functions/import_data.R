import_data <- function(path) {
  
  readr::read_rds(path) %>% 
    purrr::map_dfr(function(j) {
      dplyr::left_join(x = j$cv, y = j$prod, by = c("part", "measure", "pop_n", "biotic", "abiotic"), 
                       suffix = c(".cv", ".prod")) %>% 
        dplyr::rename("value.cv.sd" = "sd", "value.cv.mn" = "mean") %>% 
        dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))
    }) %>% 
    dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass",
                                                 "ag_production", "bg_production", "ttl_production")), 
                  measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                  pop_n = factor(as.numeric(pop_n), ordered = TRUE))
  
}
