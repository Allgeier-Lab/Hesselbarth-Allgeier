import_mortality <- function(path) {
  
  readr::read_rds(path) |> 
    purrr::map_dfr(function(j) {
      
      dplyr::left_join(x = j$connectivity, y = j$mortality, by = c("id", "pop_n", "nutrient_input", 
                                                                   "biotic", "abiotic"),
                       suffix = c(".connectivity", ".mortality"), multiple = "all")}, .id = "row_id") |> 
    dplyr::mutate(row_id = as.numeric(row_id), pop_n = factor(as.numeric(pop_n), ordered = TRUE),
                  nutrient_input = factor(nutrient_input, ordered = TRUE, labels = c("low", "medium", "high"))) |>
    tibble::tibble()
  
}
