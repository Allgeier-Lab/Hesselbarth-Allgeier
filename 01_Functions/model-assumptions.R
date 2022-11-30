model_assumptions <- function(x, title = "") {
  
  # df_residuals <- data.frame(residuals = resid(x), fitted = fitted(x), stand = rstandard(x), 
  #                            leverage = hatvalues(x))
  # 
  # gg_topleft <- ggplot(data = df_residuals, aes(x = fitted, y = residuals)) + 
  #   geom_point(shape = 1) + 
  #   geom_smooth(se = FALSE, linetype = 2, color = "red", size = 0.5) + 
  #   labs(x = "Fitted values", y = "Residuals", title = "Residuals vs. Fitted") + 
  #   theme_classic(base_size = 10)
  # 
  # gg_topright <- ggplot(data = df_residuals) 
  # 
  # gg_bottomleft <- ggplot(data = df_residuals, aes(x = fitted, y = sqrt(stand))) + 
  #   geom_point(shape = 1) + 
  #   geom_smooth(se = FALSE, linetype = 2, color = "red", size = 0.5) + 
  #   labs(x = "Fitted values", y = expression(paste(sqrt("Standardized residuals"))), title = "Scale-Location") + 
  #   theme_classic()
  # 
  # gg_bottomright <- ggplot(data = df_residuals, aes(x = leverage, y = stand)) + 
  #   geom_point(shape = 1) + 
  #   geom_smooth(se = FALSE, linetype = 2, color = "red", size = 0.5) + 
  #   labs(x = "Leverage", y = "Standardized residuals", title = "Residuals vs. Leverage") + 
  #   theme_classic(base_size = 10)
  # 
  # gg_freq <- ggplot(data = df_residuals, aes(residuals)) + 
  #   geom_histogram(color = "black", fill = NA) + 
  #   labs(x = "Residuals", y = "Frequency") + 
  #   theme_classic(base_size = 10)
  # 
  # gg_test <- ggplot() +
  #   annotate(geom = "text", x = 0.0, y = 1.0, 
  #            label = paste("Shapiro = ", round(shapiro.test(df_residuals$residuals)[[2]], 3))) + 
  #   annotate(geom = "text", x = 0.0, y = 0.0, 
  #            label = paste("Lillie = ", round(nortest::lillie.test(df_residuals$residuals)[[2]], 3))) + 
  #   scale_y_continuous(limits = c(-5, 5)) +
  #   theme_void(base_size = 10)
  # 
  # gg_titel <- ggplot() + labs(title = title) + theme_void(base_size = 10)
  # 
  # gg_model <- cowplot::plot_grid(gg_topleft, gg_topright, 
  #                                gg_bottomleft, gg_bottomright, 
  #                                gg_freq, gg_test, nrow = 3, ncol = 2)
  # 
  # cowplot::plot_grid(gg_titel, gg_model, ncol = 1, rel_heights = c(0.1, 1))
  
  #checks model assumptions
  #quartz() #opens new window
  #par(mfrow=c(3,2))
  rm=resid(x)
  fm=fitted(x)
  model2=lm(rm~fm)
  par(mfrow=c(3,2))
  plot(model2)
  hist(rm)
  plot(c(0,5),c(0,5))
  text(3,3.5,paste("shapiro = ",round(shapiro.test(rm)[[2]],3)),cex=2)
  #library(nortest)
  text(3,2,paste("lillie = ",round(nortest::lillie.test(rm)[[2]],3)),cex=2)
  
}
