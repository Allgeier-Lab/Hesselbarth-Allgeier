import_abundance <- function(path) {
  
  readr::read_rds(path) |> 
    purrr::map_dfr(function(j) j$abundance, .id = "row_id") |> 
    dplyr::mutate(row_id = as.numeric(row_id), pop_n = factor(as.numeric(pop_n), ordered = TRUE),
                  nutrient_input = factor(nutrient_input, ordered = TRUE, labels = c("low", "medium", "high"))) |>
    tibble::tibble()
  
}
