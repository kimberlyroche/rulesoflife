% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_model.R
\name{get_Gamma}
\alias{get_Gamma}
\title{Define a kernel (function) over samples}
\usage{
get_Gamma(
  kernel_scale,
  diet_weight,
  min_correlation = 0.1,
  days_to_min_autocorrelation = 90
)
}
\arguments{
\item{kernel_scale}{total variance for the composite kernel}

\item{diet_weight}{relative contribution of first three diet PCs to
covariance across samples}

\item{rho}{bandwidth for SE kernel}
}
\value{
list containing kernel function and bandwidth parameter
}
\description{
Define a kernel (function) over samples
}
\details{
Composite kernel is built from (1) squared exponential kernel (base
autocorrelation component) and (2) seasonal kernel (periodic)
}
