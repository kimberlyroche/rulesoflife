% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_data.R
\name{load_predictions}
\alias{load_predictions}
\title{Pull a prediction object from the specified directory}
\usage{
load_predictions(sname, output_dir, generate = FALSE)
}
\arguments{
\item{sname}{host short name indicating which baboon's series to fit}

\item{output_dir}{specifies the subdirectory of the model fits directory in
which to save the predictions}

\item{generate}{flag indicating whether or not to generate prediction data
if not found}
}
\value{
a named list of predicted values or NULL if no such object exists
}
\description{
Pull a prediction object from the specified directory
}
