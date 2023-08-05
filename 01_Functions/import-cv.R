import_cv <- function(path, near = FALSE) {
  
  readr::read_rds(path) |> 
    purrr::map_dfr(function(j) {
      
      if (near) {
        
        x <- dplyr::left_join(x = j$cv_near, y = j$prod_near, suffix = c(".cv", ".prod"),
                              by = c("part", "measure", "pop_n", "nutrient_input",
                                     "biotic", "abiotic"))
        
      } else {
        
        x <- dplyr::left_join(x = j$cv, y = j$prod, suffix = c(".cv", ".prod"),
                              by = c("part", "measure", "pop_n", "nutrient_input",
                                     "biotic", "abiotic"))
        
      }
    
      dplyr::rename(x, "value.sd" = "sd", "value.mn" = "mean")
      
    }, .id = "row_id") |> 
    dplyr::mutate(row_id = as.numeric(row_id), 
                  part = factor(part, levels = c("ag_production", "bg_production", "ttl_production")),
                  measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                  pop_n = factor(as.numeric(pop_n), levels = c(0,8, 16, 32, 64, 128)), 
                  nutrient_input = factor(nutrient_input, levels = c(6.18915377290985e-05, 0.000121005, 0.000191781), 
                                          labels = c("low", "medium", "high"))) |>
    dplyr::select("row_id", "part", "measure", "pop_n", "nutrient_input",
                  "biotic", "abiotic","value.cv", "value.mn", "value.sd", "value.prod") |> 
    tibble::tibble()
  
}
