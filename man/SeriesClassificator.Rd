\name{SeriesClassificator}
\alias{SeriesClassificator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classifies the time series in High or Low frequency
}
\description{
We test the zero mean hypothesis to identify zero-mean (high frequecy) and non-zero mean (low frequency) time series
}
\usage{
SeriesClassificator(IMFs, trend = NULL, pval = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IMFs}{matrix of time series to be tested}
  \item{trend}{number indicating the column/row of IMFs in which is stored the trend component}
  \item{pval}{the level of the pvalue of the t-test used for classification}
}
\details{
For each series a twp-sided t-test is run to test the hypothesis of zero mean.
}
\value{
\item{HF}{a matrix containing the time series classified as "high frequency" (i.e. with mean == 0)}
\item{LF}{a matrix containing the time series classified as "low frequency" (i.e. with mean != 0)}
\item{Trend}{a vector containing the trend series}
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
\code{t.test}
}
\examples{
load("data/HenryHubTrainPrices.RData")
Dec <- decomposer(TS = TrainPrices, method = "ceemdan", plot = T, num_imfs = 0,
                  ensemble_size = 250L, noise_strength = 0.2, S_number = 4L,
                  num_siftings = 50L, rng_seed = 0L, threads = 0L)

classifiedSeries <- SeriesClassificator(Dec$IMFs)
}
