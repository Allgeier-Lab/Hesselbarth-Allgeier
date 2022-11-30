import_cv <- function(path) {
  
  readr::read_rds(path) |> 
    purrr::map_dfr(function(j) {
      dplyr::left_join(x = j$cv, y = j$prod_cumulative, by = c("part", "measure", "pop_n", "nutrient_input",
                                                               "biotic", "abiotic"), 
                       suffix = c(".cv", ".prod")) |> 
        dplyr::rename("value.sd" = "sd", "value.mn" = "mean")
    }, .id = "row_id") |> 
    dplyr::mutate(row_id = as.numeric(row_id),part = factor(part, levels = c("ag_production", "bg_production", "ttl_production")),
                  measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                  pop_n = factor(as.numeric(pop_n), ordered = TRUE),
                  nutrient_input = factor(nutrient_input, ordered = TRUE, labels = c("low", "medium", "high"))) |>
    dplyr::select("row_id", "part", "measure", "pop_n", "nutrient_input",
                  "biotic", "abiotic","value.cv", "value.mn", "value.sd", "value.prod") |> 
    tibble::tibble()
  
}
