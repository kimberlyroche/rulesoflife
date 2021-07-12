#' Generate baseline abundances and cooperative/competitive parameters for use a generalized
#' stochastic Lotka-Volterra
#' This code is adapted directly from Tredennick et al. (2017)
#'
#' @param S number of species to simulate
#' @param start_time time point at which to begin "recording" the series (previous time points will be
#' discarded as burn-in)
#' @param end_time full length of series to simulate (including burn-in)
#' @param PhiE value from -1 to 1 defining average "synchrony" of species in the system;
#' -1 indicates strong asynchrony, 1 indicates strong synchrony, and 0, neutrality
#' @param almean value from -1 to 1 defining average competition or cooperation of species in
#' the system; I believe -1.0 indicates strong cooperation, 1.0 indicates strong competition
#' @return named list of parameters
#' @import LaplacesDemon
#' @export
generate_innate_params <- function(S, start_time, end_time, PhiE, almean) {
  if(PhiE != 1/S) {
    beta = (PhiE - 1 + sqrt(PhiE*(1 - PhiE)*(S - 1)))/(S*PhiE - 1)
  } else {
    beta = 0.5
  }

  # K <- rnorm(S, mean = 5000, sd = 0.7)
  K <- round(rdirichlet(1, rep(1, S))*10000)[1,] + runif(S, min = 1, max = 20) # prevent zeros here

  alsi <- matrix(rnorm(S*S, mean = almean, sd = runif(1, min = 0.001, max = 0.05)), nrow=S, ncol=S)
  al <- alsi*matrix(K, S, S)/t(matrix(K, S, S))
  diag(al) <- 1
  NiTh <- solve(al, K)

  StbEq <- F
  NtryStbEq <- 0
  while((StbEq == F) & (NtryStbEq < 10000)) {
    rmbound <- c(runif(1, min = 0, max = 1), runif(1, min = 0, max = 3))
    rm <- runif(S, min = min(rmbound), max = max(rmbound))
    A <- -rm*NiTh/K*al
    diag(A) <- 1 + diag(A)
    domEig <- eigen(A)$values[1]
    if(abs(domEig) < 1) {
      StbEq <- T
    }
    NtryStbEq <- NtryStbEq + 1
  }

  # sigO <- runif(S, min = 0, max = 0.2)
  sigO <- 0

  return(list(PhiE = PhiE,
              almean = almean,
              K = K,
              NiTh = NiTh,
              beta = beta,
              rm = rm,
              al = al,
              sigO = sigO))
}

#' Generate matrix that encodes per-species response to environmental perturbation
#' This code is adapted directly from Tredennick et al. (2017)
#'
#' @param S number of species to simulate
#' @param beta parameter returned by generate_innate_params
#' @return matrix
#' @export
generate_env_link <- function(S, beta) {
  EsigE <- runif(S, min = 0, max = 1)
  AsigE <- matrix(-1/S, nrow = S, ncol = S)
  diag(AsigE) <- beta - 1/S
  AsigE <- AsigE * EsigE/sqrt(beta^2 - 2*beta/S + 1/S)
  return(AsigE)
}

#' Use a generalized stochastic Lotka-Volterra model to (mechanistically) simulate count data
#' This code is adapted directly from Tredennick et al., 2017
#'
#' @param S number of species to simulate
#' @param NiTh parameter returned by generate_innate_params (near-steady state abundances?)
#' @param AsigE parameter returned by generate_env_link (per-species response to environmental perturbation)
#' @param rm parameter returned by generate_env_link (per-species growth rates)
#' @param al parameter returned by generate_env_link (cooperative/competitive dynamics)
#' @param K per-species carrying capacities
#' @param sigO magnitude of observational noise
#' @param start_time time point at which to begin "recording" the series (previous time points will be
#' discarded as burn-in)
#' @param end_time full length of series to simulate (including burn-in)
#' @param perturbations_global global environmental perturbation matrix
#' @param perturbations_sd standard deviation of zero-mean perturbations
#' @param alpha proportion of perturbations_global in incorporate into the local perturbations matrix
#' @return named list of simulated species series and generated parameters
#' @examples
#' @export
generate_series <- function(S, NiTh, AsigE, rm, al, K, sigO, start_time, end_time, perturbations_global, perturbation_sd, alpha) {
  perturbations_local <- matrix(rnorm(S*(end_time + 10 + 1), 0, perturbation_sd), S, end_time + 10 + 1)
  perturbations_local <- alpha*perturbations_global + (1 - alpha)*perturbations_local

  TimeSeries <- array(0, dim = c(S, end_time + 10 + 1))
  TimeSeries[,1] <- NiTh
  t <- 1
  while(t <= end_time + 10) {
    # epsilon <- AsigE %*% rnorm(S)
    epsilon <- AsigE %*% perturbations_local[,t]
    Nt <- TimeSeries[,t]
    TimeSeries[,t+1] <- Nt + rm*(1 - (al%*%Nt)/K) + epsilon
    if(any(TimeSeries[,t+1] <= 0)) {
      TimeSeries[which(TimeSeries[,t+1] <= 0),t+1] <- 1
    }
    t <- t + 1
  }

  # Adding observation error
  TimeSeriesObs <- array(TimeSeries, dim = c(S, end_time + 10 + 1)) +
    array(rep(sigO), dim = c(S, end_time + 10 + 1))*array(rnorm(S*(end_time + 10 + 1)), dim = c(S, end_time + 10 + 1))

  TS <- TimeSeriesObs[,start_time:end_time]
  return(TS)
}

#' CLR-convert a set of species abundance series, calculate covariance across species, and return
#' the upper triangular portion of this covariance matrix
#'
#' @param TS multi-species time series returned by generate_series
#' @return atomic vector
#' @export
convert_series_clr <- function(TS) {
  clr_TS <- clr_array(TS + 0.5, parts = 1)
  temp <- cov2cor(cov(t(clr_TS)))
  return(temp[upper.tri(temp, diag = FALSE)]) # row of heatmap
}

#' @param S number of species to simulate
#' @param H number of hosts to simulate
#' @param shared_param_level see details
#' @param shared_noise_proportion proportion of environmental perturbations "in-common" to the population
#' @param start_time time point at which to begin "recording" the series (previous time points will be
#' discarded as burn-in)
#' @param end_time full length of series to simulate (including burn-in)
#' @param PhiE value from -1 to 1 defining average "synchrony" of species in the system;
#' -1 indicates strong asynchrony, 1 indicates strong synchrony, and 0, neutrality
#' @param almean value from -1 to 1 defining average competition or cooperation of species in
#' the system; I believe -1.0 indicates strong cooperation, 1.0 indicates strong competition
#' @return named list of parameters
#' @details
#' The shared parameter flag can take values
#' 0: Nothing in common between hosts
#' 1: "Innate" parameters (baseline abundances and cooperative/competitive dynamics) in common
#' 2: "Innate" parameters plus response to environment in common
#' @export
simulate_scenario <- function(S, H, shared_param_level, shared_noise_proportion, start_time = 1000, end_time = 2000,
                              PhiE = 0.5, almean = 0.1) {
  # Original simulation used PhiE = 0.5 and swept through almean in c(-0.1, 0.1, 0.5, 0.9)
  # Generate global perturbations
  perturbation_sd <- 100
  perturbations_global <- matrix(rnorm(S*(end_time + 10 + 1), 0, perturbation_sd), S, end_time + 10 + 1)

  series_list <- list()
  heatmap <- matrix(NA, H, (S^2 - S)/2)

  global_innate_params <- NULL
  global_env_link <- NULL

  for(h in 1:H) {
    if(shared_param_level == 0 | is.null(global_innate_params)) {
      innate_params <- generate_innate_params(S, start_time, end_time, PhiE, almean)
      if(is.null(global_innate_params)) {
        global_innate_params <- innate_params
      }
    }
    if(shared_param_level == 1) {
      innate_params <- global_innate_params
    }
    if(shared_param_level < 2 | is.null(global_env_link)) {
      env_link <- generate_env_link(S, innate_params$beta)
      if(is.null(global_env_link)) {
        global_env_link <- env_link
      }
    }
    if(shared_param_level == 2) {
      env_link <- global_env_link
    }
    series <- generate_series(S, innate_params$NiTh, env_link, innate_params$rm, innate_params$al, innate_params$K, innate_params$sigO,
                              start_time, end_time, perturbations_global, perturbation_sd, shared_noise_proportion)
    series_list[[h]] <- series
    heatmap[h,] <- convert_series_clr(series)
  }
  return(list(series = series_list, heatmap = heatmap))
}

#' Generate a random high-contrast palette (for bar plotting)
#'
#' @param S desired palette size
#' @return vector of hex color strings
#' @import RColorBrewer
#' @export
generate_highcontrast_palette <- function(S) {
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  sample(getPalette(S))
}

#' Plot a random 50-day segment from a gLV simulation
#'
#' @param series matrix of observations from a single host (rows are species, columns are time points)
#' @param palette optional color palette (as a vector of hex color strings)
#' @return NULL
#' @import driver
#' @import ggplot2
#' @export
plot_gLV_bars <- function(series, palette = NULL) {
  idx_start <- sample(1:(ncol(series)-51))[1]
  idx_end <- idx_start + 50
  idx <- idx_start:idx_end
  df <- gather_array(series[,idx], "value", "taxon", "time")
  df$taxon <- factor(df$taxon)
  if(is.null(palette)) {
    palette <- generate_highcontrast_palette(length(levels(df$taxon)))
  }
  ggplot(df) +
    geom_bar(aes(x = time, y = value, fill = taxon), position = "fill", stat = "identity") +
    scale_fill_manual(values = palette) +
    theme(legend.position="right")
}

#' Plot a gLV simulation as lines
#'
#' @param series matrix of observations from a single host (rows are species, columns are time points)
#' @param palette optional color palette (as a vector of hex color strings)
#' @return NULL
#' @import driver
#' @import ggplot2
#' @export
plot_gLV_lines <- function(series, palette = NULL) {
  df <- gather_array(series, "value", "taxon", "time")
  df$taxon <- factor(df$taxon)
  if(is.null(palette)) {
    palette <- generate_highcontrast_palette(length(levels(df$taxon)))
  }
  ggplot(df) +
    geom_path(aes(x = time, y = value, color = taxon)) +
    scale_color_manual(values = palette) +
    theme(legend.position="right")
}

#' Plot a gLV simulation as a heatmap
#'
#' @param heatmap matrix of correlations (rows are hosts, columns are taxon-taxon pairs)
#' @return NULL
#' @import driver
#' @import ggplot2
#' @export
plot_gLV_heatmap <- function(heatmap) {
  d <- dist(heatmap)
  clustering.hosts <- hclust(d)
  d <- dist(t(heatmap))
  clustering.interactions <- hclust(d)
  interactions.reordered <- heatmap[clustering.hosts$order,]
  interactions.reordered <- interactions.reordered[,clustering.interactions$order]
  df <- gather_array(interactions.reordered, "value", "host", "pair")
  ggplot(df) +
    geom_tile(aes(x = pair, y = host, fill = value)) +
    scale_fill_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    theme(legend.position="right")
}
