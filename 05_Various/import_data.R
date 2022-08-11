import_data <- function(path) {
  
  abiotic <- paste0(strsplit(x = path, split = "-")[[1]][3], "_sd")
  
  readr::read_rds(path) %>% 
    purrr::map_dfr(function(j) {
      dplyr::left_join(x = j$cv, y = j$prod, by = c("part", "measure", "pop_n", "move_meta_sd", abiotic), 
                       suffix = c(".cv", ".prod")) %>% 
        dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))
    }) %>% 
    dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass",
                                                 "ag_production", "bg_production", "ttl_production")), 
                  measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                  pop_n = factor(as.numeric(pop_n), ordered = TRUE)) %>% 
    dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), 
                  measure %in% c("alpha", "gamma", "beta")) %>% 
    tibble::tibble()
  
}
