# @references
# https://rpubs.com/mengxu/exponential-model

# alpha:  starting value
# beta:   rate of decay
# theta:  asymptotic (?)

fit_nls <- function(x) {
  
  # fit linear model for starting value estimate
  linear_model <- lm(mean ~ n, data = x)
  
  # estimate starting value for theta
  theta_start <- min(x$mean) * 0.5
  
  alpha_start <- exp(coef(linear_model)[1])
  
  beta_start <- coef(linear_model)[2]
  
  # create list with starting values
  start_values <- list(alpha = alpha_start, beta = beta_start, theta = theta_start)

  # fit NLS model and get coefficients
  coefficients <- tryCatch(expr = dplyr::bind_cols(broom::tidy(stats::nls(mean ~ (alpha * exp(beta * n)) + theta,
                                                         start = start_values, data = x)),
                                                   model = "nls"),
                           error = function(e) dplyr::bind_cols(broom::tidy(linear_model), 
                                                                model = "lm"))

  # return result
  return(coefficients)

}
