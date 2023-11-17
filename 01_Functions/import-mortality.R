import_mortality <- function(path) {
  
  readr::read_rds(path) |> 
    purrr::map_dfr(function(j) {
      
      dplyr::left_join(x = j$connectivity, y = j$mortality, 
                       by = c("id", "pop_n", "nutrient_input", "biotic", "abiotic"),
                       suffix = c(".connectivity", ".mortality"), multiple = "all")}, .id = "row_id") |> 
    dplyr::mutate(row_id = as.numeric(row_id), pop_n = factor(as.numeric(pop_n), levels = c(0,8, 16, 32, 64, 128)), 
                  nutrient_input = factor(nutrient_input, levels = c(6.18915377290985e-05, 0.000121005, 0.000191781), 
                                          labels = c("low", "medium", "high"))) |>
    tibble::tibble()
  
}

import_mortality_null <- function(path) {
  
  readr::read_rds(path) |> 
    purrr::map_dfr(function(j) {
      dplyr::mutate(j$mortality, pop_n = factor(as.numeric(pop_n), levels = c(0,8, 16, 32, 64, 128)), 
                    nutrient_input = factor(nutrient_input, levels = c(6.18915377290985e-05, 0.000121005, 0.000191781), 
                                            labels = c("low", "medium", "high")))
      }, .id = "row_id") |> 
    dplyr:::mutate(row_id = as.numeric(row_id))
  
}
