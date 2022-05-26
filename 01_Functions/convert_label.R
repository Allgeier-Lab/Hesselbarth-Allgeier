convert_label <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f*'×'*10^%d"), x/10 ^ floor(log10(abs(x))), floor(log10(abs(x))))
}
