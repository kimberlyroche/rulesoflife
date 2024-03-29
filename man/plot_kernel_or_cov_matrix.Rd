% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visuals.R
\name{plot_kernel_or_cov_matrix}
\alias{plot_kernel_or_cov_matrix}
\title{Renders a kernel or covariance matrix as a square heatmap.}
\usage{
plot_kernel_or_cov_matrix(K, save_name = NULL)
}
\arguments{
\item{K}{symmetric matrix object}

\item{save_name}{if not NULL, a filename under which to save the plot}
}
\value{
a ggplot object
}
\description{
Renders a kernel or covariance matrix as a square heatmap.
}
