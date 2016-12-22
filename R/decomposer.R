# Decompose a complex time series in simpler parts
#
# Author: Nicola Donelli
# Date: 21/12/2016

decomposer <- function(TS, method = c("emd","eemd","ceemdan"), plot = FALSE, ...){
  # Arguments:
  # TS: Vector of length N. The input signal to decompose.
  # method: the method of decomposition.
  # plot: flag indicating whether or not to plot the results of the decomposition
  # ...: arguments to be supplied to the chosen method
  #      (see the respective descriptions on kernlab package documentation).
  #
  # Value:
  # A list containing:
  #   - Time series object of class "mts" where series corresponds to IMFs of the input
  #     signal, with the last series being the final residual.
  #   - the maximum error associated with the chosen decomposition
  #
  #
  # Details:
  # Decompose input data to Intrinsic Mode Functions (IMFs) with the
  # Complete Ensemble Empirical Mode Decomposition with Adaptive Noise (CEEMDAN) algorithm,
  # an adaptive variant of EEMD.


  require(Rlibeemd)
  varargs = list(...)

  if(!is.vector(TS))
    stop("TS must be a vector")

  if(method == "ceemdan"){
    IMFs <- Rlibeemd::ceemdan(TS,
                              num_imfs = ifelse("num_imfs" %in% names(varargs),varargs[["num_imfs"]],0),
                              ensemble_size = ifelse("ensemble_size" %in% names(varargs),varargs[["ensemble_size"]],250),
                              noise_strength = ifelse("noise_strength" %in% names(varargs),varargs[["noise_strength"]],0.2),
                              S_number = ifelse("S_number" %in% names(varargs),varargs[["S_number"]],4),
                              num_siftings = ifelse("num_siftings" %in% names(varargs),varargs[["num_siftings"]],50),
                              rng_seed = ifelse("rng_seed" %in% names(varargs),varargs[["rng_seed"]],0),
                              threads = ifelse("threads" %in% names(varargs),varargs[["threads"]],0))
  }else if(method == "eemd"){
    IMFs <- Rlibeemd::eemd(TS,
                 num_imfs = ifelse("num_imfs" %in% names(varargs),varargs[["num_imfs"]],0),
                 ensemble_size = ifelse("ensemble_size" %in% names(varargs),varargs[["ensemble_size"]],250),
                 noise_strength = ifelse("noise_strength" %in% names(varargs),varargs[["noise_strength"]],0.2),
                 S_number = ifelse("S_number" %in% names(varargs),varargs[["S_number"]],4),
                 num_siftings = ifelse("num_siftings" %in% names(varargs),varargs[["num_siftings"]],50),
                 rng_seed = ifelse("rng_seed" %in% names(varargs),varargs[["rng_seed"]],0),
                 threads = ifelse("threads" %in% names(varargs),varargs[["threads"]],0))

  }else if(method == "emd"){
    IMFs <- Rlibeemd::eemd(TS,
                 num_imfs = ifelse("num_imfs" %in% names(varargs),varargs[["num_imfs"]],0),
                 S_number = ifelse("S_number" %in% names(varargs),varargs[["S_number"]],4),
                 num_siftings = ifelse("num_siftings" %in% names(varargs),varargs[["num_siftings"]],50))
  }else{
    stop("'method' must be one of 'emd','ceemdan','eemd'")
  }
  if(plot){
    for(j in 1:ncol(IMFs)){
      if(j == ncol(IMFs)){
        plot(IMFs[,j], type = "l", main = "Residuals")
      }else{
        plot(IMFs[,j], type = "l", main = paste0("IMF ",j))
      }
    }
    plot(TS, type = "l", ylab="", xlab = "time")
    lines(IMFs[,ncol(IMFs)], col = "red")
    legend("topleft", legend = c("Prices","Residuals"), col = c("black", "red"),
           lty = 1, lwd = 1)
  }

  return(list(IMFs = IMFs, MaxErr = max(abs(TS - rowSums(IMFs)))))
}
