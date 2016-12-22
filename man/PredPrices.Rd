\name{PredPrices}
\alias{PredPrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictions of the original time series
}
\description{
This function assembles the outputs of HighFrequencyModeller and LowFrequencyModeller.
}
\usage{
PredPrices(HF_pred, LF_pred, Trend_pred, TestPrices = NULL, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HF_pred}{vector of predictions of the high frequency time series.}
  \item{LF_pred}{vector of predictions of the low frequency time series.}
  \item{Trend_pred}{vector of predictions of the trend time series.}
  \item{TestPrices}{a time series of real values used to evaluate the quality of the prediction.}
  \item{plot}{if TRUE and TestPrices != NULL, plot original and fitted data.}
}
\details{
This function assembles the predictions obtained for high-frequency, low-frequency and trend time series and eventually check their conformance with true test data.
}
\value{
\item{PredPrices}{the vector of predictions of the price time series.}
\item{PRMSE}{if TestPrices != NULL, percentage root mean square errors at different steps in the future.}
\item{RMSE}{ if TestPrices != NULL, root mean square errors at different steps in the future.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicola Donelli
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
load(data/HenryHubTrainPrices.RData)
load(data/HenryHubTestPrices.RData)

Dec <- decomposer(TS = TrainPrices, method = "ceemdan", plot = T, num_imfs = 0,
                  ensemble_size = 250L, noise_strength = 0.2, S_number = 4L,
                  num_siftings = 50L, rng_seed = 0L, threads = 0L)

#### Find High and Low frequency time series
classifiedSeries <- SeriesClassificator(Dec$IMFs)

#### Model, asses goodness of fit and predict High Frequency IMFs
HF_pred <- HighFrequencyModeller(
  classifiedSeries$HF, sh = 10, level = 95,
  max.p = 10, max.q = 10, max.d = 2, max.P = 2, max.Q = 2, max.D = 2, max.order = 15,
  start.p = 2, start.q = 2, start.P = 1, start.Q = 1, seasonal = T, ic = "aic",
  seasonal.test = "ocsb", stepwise = F, approximation = F, test = "adf", parallel = T,
  biasadj = T, plot = T)

HighFrequencyTester(HF_pred, lags = 30, type = "Box-Pierce")

#### Model, asses goodness of fit and predict Low Frequency IMFs
LF_pred <-
  LowFrequencyModeller(classifiedSeries$LF, range = 200, by = 4, metric = "rmse",
                       sh = 10, parallel = T, plot = T)
#### Model, asses goodness of fit and predict Trend
Trend_pred <-
  LowFrequencyModeller(classifiedSeries$Trend, range = 200, by = 4, metric = "rmse",
                       sh = 10, parallel = F, plot = T)

#### Predictions
PredPrices <- PredPrices(HF_pred, LF_pred, Trend_pred, TestPrices, plot = TRUE)
}
