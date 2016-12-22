HighFrequencyTester <- function(HF_pred, lags, type = c("Box-Pierce", "Ljung-Box")){
  # Arguments:
  # HF_pred: list, output from HighFrequencyTester
  # lags: number of lags to be used in the chosen test
  # type: type of test to be used

  # Values:
  # The function will print the outputs of the chosen test
  # and all the plot returned by tsdiag, for all the elements of the input list.

  if(!is.list(HF_pred))
    stop("HF_pred must be a list")
  if(!("method" %in% names(HF_pred)) |
     !("residuals" %in% names(HF_pred)) |
     !("coef" %in% names(HF_pred)) |
     !("model" %in% names(HF_pred)))
    stop("HF_pred must contain 'method','residuals','coef' and 'model'")

  for (j in 1:length(HF_pred)){
    print(paste(HF_pred[[j]]$method))
    print(Box.test(HF_pred[[j]]$residuals, lag = lags, type = type,
                   fitdf = length(HF_pred[[j]]$model$coef)))
    tsdiag(HF_pred[[j]]$model)
  }
}
