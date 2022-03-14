# n: Number of modifiers that are returned
# local: Current local modifier level
# modifier: vector with amplitude levels
# df: Return dataframe or list
# mehod: Which method to use
get_modifier <- function(n, local, modifier, df = FALSE, method = "runif") {
  
  result <- purrr::map(seq(from = n, to = 0, by = -1), function(i) {
    
    # repeat modifier of local metaecosyst n to 0 times
    local_mod <- rep(x = local, times = i)
    
    if (method == "runif") {
    
      # sample others
      ampl_mod <- runif(n = n - i)
      
      phase_mod <- runif(n = n - i)
      
    } else if (method == "sample") {
      
      ampl_mod <- sample(x = modifier, size = n - i, replace = TRUE)
      
      phase_mod <- sample(x = modifier, size = n - i, replace = TRUE)
      
      # ampl_mod <- sample(x = modifier[modifier != local], size = n - i, replace = TRUE) 
      
      # phase_mod <- sample(x = modifier[modifier != local], size = n - i, replace = TRUE) 
      
    } else{
      
      stop("Wrong method.", call. = FALSE)
      
    }
  
    # combine to one data.frame
    result_temp <- data.frame(n_diff = n - i, 
                              amplitude = c(local_mod, ampl_mod), 
                              phase = c(local_mod, phase_mod))
    
    if (nrow(result_temp) != n) {
      
      stop("Each local metaecosystem needs a modifier value.", call. = FALSE)
      
    }
    
    return(result_temp)
    
  })
  
  if (df) {
    
    result <- do.call(what = "rbind", args = result)
    
  }
  
  return(result)

}
