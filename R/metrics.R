metrics <- function(y, y_pred, quantiles, tau=0.5) {
  mse <- mean((quantiles - y_pred)^2, na.rm = TRUE)
  mad <- 0.5*mean(abs(quantiles - y_pred), na.rm = TRUE)
  quantile_loss <- mean((y - y_pred)*(tau - 1*(y - y_pred < 0)), na.rm = TRUE)
  return(c(
    mse = mse,
    mad = mad,
    quantile_loss = quantile_loss
    ))
}
