convert_label <- function(x, exp = TRUE, digits = 3) {
  
  if (exp) {
    
    sprintf(paste0("%.", digits, "f*'×'*10^%d"), x/10 ^ floor(log10(abs(x))), 
            floor(log10(abs(x))))
    
  } else {
    
    as.character(round(x, digits = digits + 1))
    
  }
}
