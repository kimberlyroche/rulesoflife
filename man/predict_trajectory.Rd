% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visuals.R
\name{predict_trajectory}
\alias{predict_trajectory}
\title{Draw predictions from a fido::basset model in the CLR coordinate system.
This function interpolates 5% of missing observations, which gives good quick
visuals in the ABRP data.}
\usage{
predict_trajectory(sname, output_dir, resolution = 5)
}
\arguments{
\item{sname}{host short name indicating which baboon's series to fit}

\item{output_dir}{specifies the subdirectory of the model fits directory in
which to save the predictions}

\item{resolution}{percent of missing days to interpolate; default is 5% which
is plenty good enough for plotting}
}
\description{
Draw predictions from a fido::basset model in the CLR coordinate system.
This function interpolates 5% of missing observations, which gives good quick
visuals in the ABRP data.
}
