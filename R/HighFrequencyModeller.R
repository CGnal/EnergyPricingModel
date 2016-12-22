HighFrequencyModeller <- function(
  TS, sh = 10, level = c(80,95), bootstrap = FALSE, npaths = 5000,
  d = NA, D = NA, max.p = 10, max.q = 10, max.d = 2, max.P = 2, max.Q = 2, max.D = 2,
  max.order = 20, start.p = 2, start.q = 2, start.P = 1, start.Q = 1,
  stationary = FALSE, seasonal = TRUE, ic = c("aicc", "aic", "bic"),
  stepwise = TRUE, approximation = FALSE, truncate = NULL, xreg = NULL,
  test = c("kpss","adf","pp"), seasonal.test = c("ocsb","ch"),
  allowdrift = TRUE, allowmean = TRUE, lambda = NULL, biasadj = FALSE,
  parallel = FALSE, num.cores = NULL, plot = FALSE){

  # Arguments:
  # TS: the (matrix) time series to be analyzed
  #
  # See forecast.Arima for the description of:
  # sh, level, bootstrap, npaths
  #
  # See auto.arima for the description of:
  # d, D, max.p, max.d, max.q, max.P, max.Q,
  # max.D, max.order, start.p, start.q, start.P, start.Q, stationary, seasonal,
  # ic, test, biasadj, stepwise, allowdrift, allowmean, lambda, seasonal.test, approximation
  #
  # parallel: if TRUE, the search for the best arima model is done in parallel
  # num.cores: number of cores to be used in parallel computing. If NULL
  #            all the available cores will be used.
  # plot: flag indicating whether or not plot original and fitted data.

  # Values:
  # The function returns a list with length equal to the smallest between the number
  # of rows and columns of the input TS matrix.
  # Each element of the list contains:
  # method: a string describing the used model
  # mode: a list containing the details of the used model (as returned by auto.arima)
  # level: the confidence level to be used for prediction
  # mean: the estimate of the forecast
  # lower: the lower extreme of the confidence interval of the prediction (at the given level)
  # upper: the upper extreme of the confidence interval of the prediction (at the given level)
  # x: the time series in input
  # xname: the names associated to the elements of the time series in input
  # fitted: the fitted values of the time series
  # residuals: the series of the residuals
  # fit_rmse: root mean squared error of the fit

  require(forecast)
  require(Metrics)

  if(!is.numeric(TS))
    stop("TS must be numeric")

  if(is.matrix(TS)){
    if(nrow(TS) < ncol(TS)){
      TS <- t(TS)
    }
    J <- ncol(TS)
  } else if(is.ts(TS)){
    TS <- matrix(TS, nrow = length(TS), ncol = 1)
    J <- 1
  }else if(is.vector(TS)){
    TS <- matrix(TS, nrow = length(TS), ncol = 1)
    J <- 1
  }else{
    stop("TS must be a matrix, a vector or a time series")
  }

  if(parallel){
    require(parallel)
    require(doParallel)

    if(is.null(num.cores)){
      registerDoParallel(makeCluster(detectCores()))
    }else{
      registerDoParallel(makeCluster(num.cores))
    }

    TS_pred <-
      foreach(j = 1:J, .packages = c("forecast","Metrics")) %dopar%{
        tmp <- forecast::forecast.Arima(
          forecast::auto.arima(TS[,j], d = d, D = D, max.p = max.p, max.d = max.d,
                               max.q = max.q, max.P = max.P,max.Q = max.Q, max.D = max.D, max.order = max.order,
                               start.p = start.p, start.q = start.q, start.P = start.P,
                               start.Q = start.Q, stationary = stationary, seasonal = seasonal,
                               ic = ic, test = test, biasadj = biasadj, stepwise = stepwise,
                               allowdrift = allowdrift, allowmean = allowmean, lambda = lambda,
                               seasonal.test = seasonal.test, approximation = approximation),
          h = sh, level = level, biasadj = biasadj, bootstrap = bootstrap, npaths = npaths,
          xreg = xreg, lambda = lambda)

        append(tmp, list(fit_rmse = Metrics::rmse(TS[,j], tmp$fitted)))
      }
    if(plot){
      for(j in 1:J){
        plot(TS[,j], type = "l", ylab = paste0("TS[,",j,"]"), xlab = "Time",
             main = paste0("HF series: ",j))
        lines(TS_pred[[j]]$fitted, col = "red")
        legend("topleft", legend = c("True","Fitted"), col = c("black", "red"),
               lty = 1, lwd = 1)
      }
    }
  }else{
    TS_pred <- list()
    for(j in 1:J){
      TS_pred[[j]] <-
        forecast::forecast.Arima(
          forecast::auto.arima(TS[,j], d = d, D = D, max.p = max.p, max.d = max.d,
                               max.q = max.q, max.P = max.P,max.Q = max.Q, max.D = max.D, max.order = max.order,
                               start.p = start.p, start.q = start.q, start.P = start.P,
                               start.Q = start.Q, stationary = stationary, seasonal = seasonal,
                               ic = ic, test = test, biasadj = biasadj, stepwise = stepwise,
                               allowdrift = allowdrift, allowmean = allowmean, lambda = lambda,
                               seasonal.test = seasonal.test, approximation = approximation),
          h = sh, level = level, biasadj = biasadj, bootstrap = bootstrap, npaths = npaths,
          xreg = xreg, lambda = lambda)
      TS_pred[[j]]$fit_rmse <- Metrics::rmse(TS[,j], TS_pred[[j]]$fitted)
    }

    if(plot){
      for(j in 1:J){
        plot(TS[,j], type = "l", ylab = paste0("TS[,",j,"]"), xlab = "Time",
             main = paste0("HF series: ",j))
        lines(TS_pred[[j]]$fitted, col = "red")
        legend("topleft", legend = c("True","Fitted"), col = c("black", "red"),
               lty = 1, lwd = 1)
      }
    }
  }
  return(TS_pred)
}
