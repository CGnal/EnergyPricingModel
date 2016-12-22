\name{HighFrequencyTester}
\alias{HighFrequencyTester}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test the goodness of fit of the ARIMA model
}
\description{
This function runs diagnostic tests for ARIMA models
}
\usage{
HighFrequencyTester(HF_pred, lags, type = c("Box-Pierce", "Ljung-Box"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HF_pred}{list, output from HighFrequencyTester}
  \item{lags}{number of lags to be used in the chosen test}
  \item{type}{type of test to be used}
}
\details{
The function will print the outputs of the chosen test and all the plot returned by tsdiag, for all the elements of the input list.
}
\value{
Ouputs of the Box.test and ts.diag functions.
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
}
