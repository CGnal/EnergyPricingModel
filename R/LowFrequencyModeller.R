LowFrequencyModeller <- function(
  TS, range = 20, by = 1, metric = c("rmse"), sh = 10, parallel = FALSE, plot = TRUE){
  # Arguments:
  # TS: a (matrix) time series
  # range: the range of the interval around the minimal number of lags required to
  #        annihilate the acf to be spanned to look for the optimal number of
  #        regression components
  # by: the length of the steps used to span the space of possible numbers of
  #     regression components
  # metric: the metric to be used to define the best number of regression components
  # sh: the number of steps ahead to forecast and to be used to evaluate the
  #     goodness of a regression model.
  # parallel: a flag indicating whether or not to  use parallel computing.
  #
  # Values:

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

  svar <- function(TS, p, sh){
    # Support vector autoregression
    #
    # Arguments:
    # TS: the input time series (matrix or vector)
    # p: the number of lags to be considered
    # sh: the number of steps ahead to be considered for the prediction
    #
    # Values:
    # The function returns a list of four elements
    # mod: a list containing the details of the svm model used
    # fit_rmse: the root means squared error of the fitted values
    # predictions: the predicted values
    # num_reg: the number p of regressors considered

    require(Metrics)
    require(e1071)

    lagger <- function(x, p){
      # Arguments:
      # x: a vector
      # p: the considered number of lags
      #
      # Values:
      # A matrix containing with nrow = length(x), ncol = p and LagMat[,j] = Lag(x, shift = j)

      require(Hmisc)
      if(!is.vector(x))
        x <- as.vector(x)

      LagMat <- matrix(NA, nrow = length(x), ncol = p)
      for (j in 1:p)
        LagMat[,j] <- Hmisc::Lag(x, shift = j)
      return(LagMat)
    }

    X <- data.frame(cbind(TS, lagger(TS, p)))
    colnames(X) <- c("y",paste0(rep("x",p),1:p))
    X <- X[(p+1):nrow(X),]
    mod <- e1071::svm(as.formula(paste0("y~",paste0(colnames(X)[-1], collapse = "+"))),
               data = X, kernel = "radial")
    OSApred <- cbind()
    for (h in 1:sh){
      if(ncol(X) - h > 0){
        if(h == 1){
          OSAreg <- X[nrow(X), 1:(ncol(X) - h)]
        }else{
          OSAreg <- cbind(OSApred, X[nrow(X),1:(ncol(X)-h)])
        }
        colnames(OSAreg) <- paste0(rep("x",p),1:p)
        OSApred <- cbind(OSApred, predict(mod,OSAreg))
      }else{
        OSAreg <- matrix(OSApred[,(h-ncol(X) + 1) : ncol(OSApred)],
                         nrow = 1, ncol = ncol(X) - 1)
        colnames(OSAreg) <- paste0(rep("x",p),1:p)
        OSApred <- cbind(OSApred, predict(mod,OSAreg))
      }
    }
    return(list(mod = mod, fit_rmse = Metrics::rmse(TS[(p+1):length(TS)],mod$fitted), predictions = as.vector(OSApred), num_reg = p))
  }

  if(metric == "rmse"){
    val_metric <- function(actual, predicted) {
      require(Metrics)
      Metrics::rmse(actual, predicted)
      }
  }else
    stop("metric must be one of the predefined")

  train <- 1:(nrow(TS) - sh)
  test <- (nrow(TS) - sh + 1):nrow(TS)

  if (!parallel){
    Out <- list()
    for(j in 1:J){
      X <- TS[,j]
      bp <- which(acf(X, length(X), plot = FALSE)$acf < 0.025)[1]
      best_eval <- Inf
      for(n in seq(from = max((bp-range),4), to = min((bp+range),length(X)), by = by)){
        temp_pred <- svar(X[train], p = n, sh = sh)
        eval_pred <- val_metric(X[test][1:sh], temp_pred$predictions)
        if(eval_pred < best_eval){
          best_eval <- eval_pred
          best_p <- n
        }
      }
      Out[[j]] <- svar(X, p = best_p, sh = sh)
      Out[[j]]$pred_eval <- best_eval
    }
  }else{
    require(foreach)
    require(parallel)
    require(doParallel)
    registerDoParallel(makeCluster(detectCores()))
    Out <-
      foreach::foreach::foreach(j = 1:J, .export = c("svar")) %dopar% {
        X <- TS[,j]
        bp <- which(acf(X, length(X), plot = FALSE)$acf < 0.025)[1]
        best_eval <- Inf
        for(n in seq(from = max((bp-range),4), to = min((bp+range),length(X)), by = by)){
          temp_pred <- svar(X[train], p = n, sh = sh)
          eval_pred <- val_metric(X[test][1:sh], temp_pred$predictions)
          if(eval_pred < best_eval){
            best_eval <- eval_pred
            best_p <- n
          }
        }
        append(svar(X, p = best_p, sh = sh), list(pred_eval = best_eval))
      }
  }
  if(plot){
    for(j in 1:J){
      plot(TS[Out[[j]]$num_reg:nrow(TS),j], type = "l", main = paste0("LF-IMF ",j))
      lines(Out[[j]]$mod$fitted, col = "red")
    }
  }
  return(Out)
}
