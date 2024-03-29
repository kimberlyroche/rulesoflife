% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_data.R
\name{load_fit}
\alias{load_fit}
\title{Pull a fitted model object from the specified directory}
\usage{
load_fit(sname, MAP, output_dir)
}
\arguments{
\item{sname}{host short name indicating which baboon's series to fit}

\item{MAP}{flag indicating whether or not to return (single) MAP estimate}

\item{output_dir}{specifies the subdirectory of the model fits directory in
which to save the predictions}
}
\value{
a bassetfit object or NULL if no such object exists
}
\description{
Pull a fitted model object from the specified directory
}
