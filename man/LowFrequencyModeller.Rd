\name{LowFrequencyModeller}
\alias{LowFrequencyModeller}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model low frequency time series
}
\description{
This function implements a support vector auto regression to fit and forecast low frequency time series.
}
\usage{
LowFrequencyModeller(TS, range = 20, by = 1, metric = c("rmse"), sh = 10, parallel = FALSE, plot = TRUE)
}

\arguments{
  \item{TS}{a (matrix) time series.}
  \item{range}{the range of the interval around the minimal number of lags required to annihilate the acf to be spanned to look for the optimal number of regression components.}
  \item{by}{the length of the steps used to span the space of possible numbers of regression components.}
  \item{metric}{the metric to be used to define the best number of regression components.}
  \item{sh}{the number of steps ahead to forecast and to be used to evaluate the goodness of a regression model.}
  \item{parallel}{if TRUE, the search for the best svar model is done in parallel.}
  \item{plot}{if TRUE, plot original and fitted data.}
}
\details{
The function runs an automatic selection of the best support vector autoregression model. It performs a grid search of the best number of autoregressive components: starting from the first value of lags for which the autocorrelation function becomes zero, the function spans an interval of range 2*range around the selected number of lags with the given length of steps. The evaluation of the models is based on the root mean square error of the prediction of sh steps ahead.
}
\value{
The function returns a list containing, for each considered time series:
\item{mod}{a list containing the details of the used model (as returned by svr).}
\item{fit_rmse}{the root mean squared error of the fit.}
\item{predictions}{the sh predictions.}
\item{num_reg}{the number of autoregressive components chosen for the best model.}
\item{pred_eval}{the root mean square error of the predictions obtained with the chosen model.}
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
svr
}
\examples{
load(data/HenryHubTrainPrices.RData)

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


}
