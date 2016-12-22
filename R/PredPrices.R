PredPrices <- function(HF_pred, LF_pred, Trend_pred, TestPrices = NULL, plot = TRUE){
  # Arguments:
  # HF_pred: vector of predictions of the high frequency time series.
  # LF_pred: vector of predictions of the low frequency time series.
  # Trend_pred: vector of predictions of the trend time series.
  # TestPrices: time series of prices to test the prediction.
  # plot: flag indicating whether or not plot original and fitted data.
  #
  # Values:
  # PredPrices: the vector of predictions of the price time series.
  # PRMSE: percentage root mean square errors at different steps in the future.
  # RMSE: root mean square errors at different steps in the future.

  if(!is.list(HF_pred) | !is.list(LF_pred) | !is.list(Trend_pred))
    stop("HF_pred, LF_pred e Trend_pred devono essere delle liste")

  sh <- min(length(HF_pred[[1]]$mean),
            length(LF_pred[[1]]$predictions),
            length(Trend_pred[[1]]$predictions))

  PredPrices <- c()
  for(i in 1:length(HF_pred)){
    tmp <- HF_pred[[i]]$mean[1:sh]
    PredPrices <- rbind(PredPrices, tmp)
  }
  for(j in 1:length(LF_pred)){
    tmp <- LF_pred[[j]]$predictions[1:sh]
    PredPrices <- rbind(PredPrices, tmp)
  }
  PredPrices <- colSums(rbind(PredPrices, Trend_pred[[1]]$predictions[1:sh]))

  if(!is.null(TestPrices)){
    PRMSE <- c()
    RMSE <- c()
    for (h in 1:sh){
      RMSE <- c(RMSE,rmse(TestPrices[1:h],PredPrices[1:h] ))
      PRMSE <- c(PRMSE, RMSE[h]/sqrt(mean(TestPrices[1:h]^2)))
    }
    if(plot){
      plot(TestPrices[1:sh], type = "l",
           ylim = c(min(min(PredPrices),min(TestPrices[1:sh])),
                    max(max(PredPrices),max(TestPrices[1:sh]))),
           ylab = "TestPrices", xlab = "Steps-ahead")
      lines(PredPrices, col = "red")
      plot(PRMSE, ylab = "Percentage Root Mean Squared Error", xlab = "Number of steps ahead")
      plot(RMSE, ylab = "Root Mean Squared Error", xlab = "Number of steps ahead")
    }
    return(list(PredPrices = PredPrices, PRMSE = PRMSE, RMSE = RMSE))
  }
  return(PredPrices)
}
