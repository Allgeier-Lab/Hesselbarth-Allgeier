convert_label <- function(x, digits = 2) {
  
  if (x < 0.0001 || x > 10000) {
  
  sprintf(paste0("%.", digits, "f*'×'*10^%d"), x/10 ^ floor(log10(abs(x))), 
          floor(log10(abs(x))))
    
  } else {
    
    as.character(round(x, digits = digits + 1))
    
  }
}
