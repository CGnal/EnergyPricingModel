\name{decomposer}
\alias{decomposer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to decompose complex price series
}
\description{
This function implements different possible decomposition methods
}
\usage{
decomposer(TS, method = c("emd","eemd","ceemdan"), plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TS}{Vector of length N. The input signal to decompose.}
  \item{method}{Decomposition method to be used.}
  \item{plot}{Logical. Indicating whether or not to plot the results of the decomposition}
  \item{...}{arguments to be supplied to the chosen method (see the respective descriptions on kernlab package documentation).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{IMFs}{Time series object of class "mts" where series corresponds to IMFs of the input signal, with the last series being the final residual.}
\item{MaxErr}{The maximum absolute error associated with the chosen decomposition}
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
\code{\link{ceemdan}, \link{eemd}, \link{emd}}
}
\examples{
attach(HenryHubDailyPrices)
#### Decomposition in IMFs
Dec <- decomposer(TS = TrainPrices, method = "ceemdan", plot = T, num_imfs = 0,
                  ensemble_size = 250L, noise_strength = 0.2, S_number = 4L,
                  num_siftings = 50L, rng_seed = 0L, threads = 0L)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
