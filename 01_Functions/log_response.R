log_response <- function(data, indices) {
  
  x <- data[indices, ] 
  
  f <- mean(log(x$trtm)) - mean(log(x$ctrl))
  
  # f <- (exp(f) - 1) * 100
  
  return(f)
} 
