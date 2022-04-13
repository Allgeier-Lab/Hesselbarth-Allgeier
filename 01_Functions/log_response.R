log_response <- function(data, indices, relative = FALSE) {
  
  if (!all(names(data) == c("ctrl", "trtm"))) stop("Columns must be named 'ctrl' and 'trtm'.")
  
  df_boot <- data[indices, ] 
  
  rr <- log(mean(df_boot$trtm) / mean(df_boot$ctrl))
  
  if (relative) {
    
    rr <- 100 * (exp(rr) - 1)
    
  }
  
  return(rr)
}
