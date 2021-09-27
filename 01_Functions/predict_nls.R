# alpha: coefficients$estimate[1]
# beta: coefficients$estimate[]

predict_nls <- function(x, coefficients) {
  
  # predict values using lm formula
  if (all(coefficients$model == "lm")) {
    
    # predict value ~ n
    predictions <- coefficients$estimate[1] + (coefficients$estimate[2] * x)
    
  # predict values for NLS formula
  } else if (all(coefficients$model == "nls")) {
    
    # predict value ~ (alpha * exp(beta * n)) + theta
    predictions <- (coefficients$estimate[1] * exp(coefficients$estimate[2] * x)) +
      coefficients$estimate[3]
  }
  
  # this shouldnt happen
  else {
    
    stop("This should not happen...", call. = FALSE)
    
  }
  
  # return value
  return(data.frame(x = x, value = predictions))
  
}
