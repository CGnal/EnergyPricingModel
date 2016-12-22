\name{HighFrequencyModeller}
\alias{HighFrequencyModeller}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Automatic model choice and forecasting high frequency time series using ARIMA models
}
\description{
This function performs model choice and forecasting using ARIMA models.
}
\usage{
HighFrequencyModeller(TS, sh = 10, level = c(80,95), bootstrap = FALSE, npaths = 5000,
  d = NA, D = NA, max.p = 10, max.q = 10, max.d = 2, max.P = 2, max.Q = 2, max.D = 2,
  max.order = 20, start.p = 2, start.q = 2, start.P = 1, start.Q = 1,
  stationary = FALSE, seasonal = TRUE, ic = c("aicc", "aic", "bic"),
  stepwise = TRUE, approximation = FALSE, truncate = NULL, xreg = NULL,
  test = c("kpss","adf","pp"), seasonal.test = c("ocsb","ch"),
  allowdrift = TRUE, allowmean = TRUE, lambda = NULL, biasadj = FALSE,
  parallel = FALSE, num.cores = NULL, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TS}{the (matrix) time series to be analyzed.}
  \item{sh}{number of periods for forecasting. If xreg is used, sh is ignored and the number of forecast periods is set to the number of rows of xreg.}
  \item{level}{confidence level for prediction interval.}
  \item{bootstrap}{if TRUE, then prediction intervals computed using simulation with resampled errors.}
  \item{npaths}{number of sample paths used in computing simulated prediction intervals when bootstrap=TRUE.}
  \item{d}{order of first-differencing. If missing, will choose a value based on KPSS test.}
  \item{D}{order of seasonal-differencing. If missing, will choose a value based on OCSB test.}
  \item{max.p}{maximum value of p.}
  \item{max.q}{maximum value of q.}
  \item{max.d}{maximum value of d.}
  \item{max.P}{maximum value of P.}
  \item{max.Q}{maximum value of Q.}
  \item{max.D}{maximum value of D.}
  \item{max.order}{maximum value of p+q+P+Q if model selection is not stepwise.}
  \item{start.p}{starting value of p in stepwise procedure.}
  \item{start.q}{starting value of q in stepwise procedure.}
  \item{start.P}{starting value of P in stepwise procedure.}
  \item{start.Q}{starting value of Q in stepwise procedure.}
  \item{stationary}{if TRUE, restricts search to stationary models.}
  \item{seasonal}{if FALSE, restricts search to non-seasonal models.}
  \item{ic}{information criterion to be used in model selection.}
  \item{stepwise}{if TRUE, will do stepwise selection (faster). Otherwise, it searches over all models. Non-stepwise selection can be very slow, especially for seasonal models.}
  \item{approximation}{if TRUE, estimation is via conditional sums of squares andthe information criteria used for model selection are approximated. The final model is still computed using maximum likelihood estimation. Approximation should be used for long time series or a high seasonal period to avoid excessive computation times.}
  \item{truncate}{an integer value indicating how many observations to use in model selection. The last truncate values of the series are used to select a model when truncate is not NULL and approximation=TRUE. All observations are used if either truncate=NULL or approximation=FALSE.}
  \item{xreg}{optionally, a vector or matrix of external regressors, which must have the same number of rows as y.}
  \item{test}{type of unit root test to use. See ndiffs for details.}
  \item{seasonal.test}{this determines which seasonal unit root test is used. See nsdiffs for details.}
  \item{allowdrift}{if TRUE, models with drift terms are considered.}
  \item{allowmean}{If TRUE, models with a non-zero mean are considered.}
  \item{lambda}{Box-Cox transformation parameter. Ignored if NULL. Otherwise, data transformed before model is estimated.}
  \item{biasadj}{use adjusted back-transformed mean for Box-Cox transformations. If TRUE, point forecasts and fitted values are mean forecast. Otherwise, these points can be considered the median of the forecast densities.}
  \item{parallel}{if TRUE, the search for the best arima model is done in parallel.}
  \item{num.cores}{number of cores to be used in parallel computing. If NULL all the available cores will be used.}
  \item{plot}{if TRUE, plot original and fitted data.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{method}{a string describing the used model}
\item{mod}{a list containing the details of the used model (as returned by auto.arima)}
\item{level}{the confidence level to be used for prediction}
\item{mean}{the estimate of the forecast}
\item{lower}{the lower extreme of the confidence interval of the prediction (at the given level)}
\item{upper}{the upper extreme of the confidence interval of the prediction (at the given level)}
\item{x}{the time series in input}
\item{xname}{the names associated to the elements of the time series in input}
\item{fitted}{the fitted values of the time series}
\item{residuals}{the series of the residuals}
\item{fit_rmes}{root mean squared error of the fit}
}
\references{

}
\author{
Nicola Donelli
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
auto.arima, forecast
}
\examples{
load(data/HenryHubTrainPrices.RData)

Dec <- decomposer(TS = TrainPrices, method = "ceemdan", plot = T, num_imfs = 0,
                  ensemble_size = 250L, noise_strength = 0.2, S_number = 4L,
                  num_siftings = 50L, rng_seed = 0L, threads = 0L)

#### Find High and Low frequency time series
classifiedSeries <- SeriesClassificator(Dec$IMFs)

#### Model, asses goodness of fit and predict High Frequency IMFs
system.time(
  HF_pred <- HighFrequencyModeller(
    classifiedSeries$HF, sh = 10, level = 95,
    max.p = 10, max.q = 10, max.d = 2, max.P = 2, max.Q = 2, max.D = 2, max.order = 15,
    start.p = 2, start.q = 2, start.P = 1, start.Q = 1, seasonal = T, ic = "aic",
    seasonal.test = "ocsb", stepwise = F, approximation = F, test = "adf", parallel = T,
    biasadj = T, plot = T)
)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
