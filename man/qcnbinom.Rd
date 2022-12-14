% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cN.R
\name{qcnbinom}
\alias{qcnbinom}
\title{Quantile function for the continuous Negative Binomial distribution with parameters r and beta}
\usage{
qcnbinom(p, r, beta)
}
\arguments{
\item{p}{vector of probabilities}

\item{r}{parameter from the hierarchical gamma distribution}

\item{beta}{parameter from the hierarchical gamma distribution}
}
\value{
vector of quantiles
}
\description{
Quantile function for the continuous Negative Binomial distribution with parameters r and beta
}
\examples{
qcnbinom(p = 0.3, r = 2, beta = 5)
}
