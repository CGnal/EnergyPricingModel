SeriesClassificator <- function(IMFs, trend = NULL, pval = 0.05){
  # Arguments:
  # IMFs: matrix containing the time series to be classified
  # trend: number indicating the column/row of IMFs in which is stored the trend component
  # pval: the level of the pvalue of the t-test used for classification

  # Values:
  # A list containing three elements:
  #   -HF: a matrix containing the time series classified as "high frequency" (i.e. with mean == 0)
  #   -LF: a matrix containing the time series classified as "low frequency" (i.e. with mean != 0)
  #   -Trend: a vector containing the trend series

  if(!is.matrix(IMFs))
    stop("IMFs must be a matrix")

  if(nrow(IMFs) < ncol(IMFs))
    IMFs <- t(IMFs)
  if(is.null(trend))
    trend <- ncol(IMFs)

  p_val_t <- c()
  for (j in setdiff(1:ncol(IMFs), trend)){
    if (j==1){
      p_val_t[j] <- t.test(IMFs[,1:j], alternative = "two.sided")$p.value
    }else{
      p_val_t[j] <- t.test(rowSums(IMFs[,1:j]), alternative = "two.sided")$p.value
    }
  }
  HF <- IMFs[,-trend][,p_val_t > pval]
  LF <- IMFs[,-trend][,p_val_t <= pval]
  Trend <- IMFs[,trend]

  return(list(HF = HF, LF = LF, Trend = Trend, pvals = p_val_t))
}
