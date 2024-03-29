% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{pull_series}
\alias{pull_series}
\title{Pull the predicted host x taxon series from a data.frame of predictions}
\usage{
pull_series(prediction_obj, sname, idx)
}
\arguments{
\item{prediction_obj}{a long tibble containing aligned centered trajectories}

\item{sname}{host short name}

\item{idx}{taxon index}
}
\value{
filtered tibble
}
\description{
Pull the predicted host x taxon series from a data.frame of predictions
}
