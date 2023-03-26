plot_lm <- function(x) {
  
  gg_fit_res <- ggplot2::ggplot(data.frame(fitted = fitted(x), resid = resid(x))) + 
    ggplot2::geom_point(ggplot2::aes(x = fitted, y = resid), shape = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "#e76254") +
    ggplot2::labs(x = "Fitted values", y = "Residuals") +
    ggplot2::theme_classic()
  
  gg_qq <- ggplot2::ggplot(data = x, aes(sample = rstandard(x))) + 
    ggplot2::stat_qq(shape = 1) +     
    ggplot2::stat_qq_line(linetype = 2, color = "#e76254") +   
    ggplot2::labs(x = "Theoretical Quantiles", y = "Standardized Residuals") +
    ggplot2::theme_classic()
  
  gg_res_hist <- ggplot2::ggplot(data.frame(resid = resid(x)), aes(x = resid)) + 
    ggplot2::geom_density() +
    ggplot2::geom_histogram(aes(y = ..density..), fill = NA, color = "black") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "#e76254") +
    ggplot2::labs(x = "Residuals", y = "Density") +
    ggplot2::theme_classic()
  
  cowplot::plot_grid(gg_fit_res, gg_qq, gg_res_hist, ncol = 3)
   
}
