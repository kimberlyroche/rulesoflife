#' Generate a color palette suitable for stacked bar plots given a feature number (S)
#'
#' @param S number of features (colors to generate)
#' @return list of hex colors
#' @import RColorBrewer
#' @export
generate_highcontrast_palette <- function(S) {
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  sample(getPalette(S))
}

#' Generate a color palette suitable for stacked bar plots given a feature number (S)
#'
#' @param hex_min hex code for min color in ramp
#' @param hex_mid hex code for median color in ramp
#' @param hex_max hex code for max color in ramp
#' @param S number of colors to generate
#' @return list of hex colors
#' @import RColorBrewer
#' @export
generate_palette <- function(hex_min, hex_mid, hex_max, S) {
  getPalette <- colorRampPalette(c(hex_min, hex_mid, hex_max))
  getPalette(S)
}

#' Renders a kernel or covariance matrix as a square heatmap.
#'
#' @param K symmetric matrix object
#' @param save_name if not NULL, a filename under which to save the plot
#' @return NULL
#' @import ggplot2
#' @import tidyr
#' @export
plot_kernel_or_cov_matrix <- function(K, save_name = NULL) {
  K <- cbind(1:nrow(K), K)
  colnames(K) <- c("sample1", 1:nrow(K))
  K <- pivot_longer(as.data.frame(K), !sample1, names_to = "sample2", values_to = "covariance")
  K$sample2 <- as.numeric(K$sample2)
  p <- ggplot(K, aes(x = sample1, y = sample2)) +
    geom_raster(aes(fill = covariance)) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                         midpoint = 0)
  if(is.null(save_name)) {
    show(p)
  } else {
    output_dir <- check_dir(c("output", "figures"))
    ggsave(file.path(output_dir, paste0(save_name, ".png")),
           p,
           units = "in",
           height = 3,
           width = 4.25)
  }
}

#' Draw predictions from a fido::basset model in the CLR coordinate system.
#' This function interpolates 5% of missing observations, which gives good quick
#' visuals in the ABRP data.
#'
#' @param sname host short name indicating which baboon's series to fit
#' @param output_dir specifies the subdirectory of the model fits directory in
#' which to save the predictions
#' @return NULL
#' @import fido
#' @export
predict_trajectory <- function(sname, output_dir) {
  fit <- load_fit(sname = sname, MAP = FALSE, output_dir = output_dir)
  if(fit$coord_system != "clr") {
    fit <- to_clr(fit)
  }
  first_day <- min(fit$X[1,])
  last_day <- max(fit$X[1,])
  n_days <- last_day
  resolution <- 5
  span <- round(seq(from = first_day, to = last_day, length.out = round(n_days * resolution/100)))
  span <- sort(c(span, c(fit$X[1,]))) # include observations

  newdata <- matrix(NA, nrow(fit$X), length(span))
  newdata[1,] <- span
  # Linearly interpolate the diet PCs for lack of a better option
  x <- fit$X[1,]
  for(cov_idx in 2:4) {
    y <- fit$X[cov_idx,]
    newdata[cov_idx,] <- approx(x = x, y = y, xout = span, ties = "ordered")$y
  }

  pred <- predict(fit, newdata = newdata, response = "Eta") # or LambdaX or Y
  predictions <- list(Eta = pred, # taxa x time points x n_iter
                      span = c(span),
                      resolution = resolution)

  output_dir <- check_dir(c("output", "model_fits", output_dir, "predictions"))
  filename <- paste0(fit$sname, ".rds")
  saveRDS(predictions, file = file.path(output_dir, filename))
  return(list(fit = fit, predictions = predictions))
}

#' Visualize inferred logratio trajectories (CLR by default)
#'
#' @param sname host short name indicating which baboon's series to fit
#' @param output_dir specifies the subdirectory of the model fits directory in
#' which to save the predictions
#' @param coord index of logratio to plot
#' @param coord_label optional coord label (e.g. "family Lachnospiraceae")
#' @param show_observations flag indicating whether or not to logratio transform
#' and render the observed data; this can be useful for judging sensitivity of
#' fit
#' @param save_file flag indicating whether or not to save rendered plot to file
#' @return NULL
#' @import driver
#' @import fido
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @export
plot_trajectory <- function(sname, output_dir, coord, coord_label = NULL,
                            show_observations = FALSE, save_file = FALSE) {
  predictions <- load_predictions(sname, output_dir)
  if(is.null(predictions)) {
    fit_and_predictions <- predict_trajectory(fit)
    fit <- fit_and_predictions$fit
    predictions <- fit_and_predictions$predictions
  }
  # Gather predictions into intervals
  days <- predictions$span
  eta_df <- gather_array(predictions$Eta, val, coord, sample, iteration)
  map <- data.frame(sample = 1:max(eta_df$sample), day = days)
  eta_df <- left_join(eta_df, map, by = "sample")
  eta_df <- eta_df %>%
    select(coord, day, iteration, val) %>%
    group_by(coord, day) %>%
    summarize(p2.5 = quantile(val, probs = c(0.025)),
              p25 = quantile(val, probs = c(0.25)),
              mean = mean(val),
              p75 = quantile(val, probs = c(0.75)),
              p97.5 = quantile(val, probs = c(0.975)), .groups = "keep")
  plot_df <- eta_df[eta_df$coord == coord,]

  if(show_observations) {
    # Logratio transform real data
    observations <- fit$Y
    observations <- t(driver::clr((t(observations) + 0.5)))
    observations <- data.frame(day = c(fit$X), logratio = observations[coord,])

    plot_df <- left_join(plot_df, observations, by = "day")
  }
  if(save_file) {
    filename <- paste0("GP_", fit$sname, "_", "coord-", coord, ".png")
  }

  p <- ggplot(plot_df) +
    geom_ribbon(aes(x = day, ymin = p2.5, ymax = p97.5), fill = "grey80") +
    geom_ribbon(aes(x = day, ymin = p25, ymax = p75), fill = "grey60") +
    geom_line(aes(x = day, y = mean))
  if(!is.null(coord_label)) {
    p <- p +
      ylab(coord_label)
  } else {
    p <- p +
      ylab(paste0("logratio index ", coord))
  }
  if(show_observations) {
    p <- p +
      geom_point(aes(x = day, y = logratio), size = 1, na.rm = TRUE)
  }
  if(save_file) {
    output_dir <- check_dir(c("output", "figures"))
    ggsave(file.path(output_dir, filename), p, units = "in", height = 3, width = 10)
  } else {
    show(p)
  }
}

#' Summarize the pairwise correlations (optionally via proportionality) from all
#' fitted models of a given type. Uses MAP estimates if available, otherwise,
#' summarizes full posterior model fits to give a "rug" matrix where rows are
#' hosts and columns are unique pairwise interactions between taxa.
#'
#' @param output_dir specifies the subdirectory of the model fits directory in
#' which to look for fitted model output
#' @param use_proportionality flag indicating whether to compute proportionality
#' between features; if FALSE, CLR correlations are computed
#' @return named list of "rug" and associated labels
#' @import fido
#' @export
summarize_Sigmas <- function(output_dir, use_proportionality = FALSE) {
  # Get all fitted model objects
  output_dir_full <- check_dir(c("output", "model_fits", output_dir, "MAP"))
  file_list <- list.files(path = output_dir_full, pattern = "*.rds")
  if(length(file_list) == 0) {
    output_dir_full <- check_dir(c("output", "model_fits", output_dir, "full_posterior"))
    file_list <- list.files(path = output_dir_full, pattern = "*.rds")
  }
  if(length(file_list) == 0) {
    stop("No model output found!")
  }
  # Get taxa number and posterior sample number
  fit <- readRDS(file.path(output_dir_full, file_list[1]))
  D <- fit$D
  iter <- fit$iter
  # Initialize the stuff we'll return
  pairs <- combn(1:(fit$D), m = 2)
  # Exclude the "other" category by default (the reference taxon, D)
  include_pairs <- which(pairs[1,] != fit$D & pairs[2,] != fit$D)
  pair1 <- pairs[1,include_pairs]
  pair2 <- pairs[2,include_pairs]
  rug <- matrix(NA, length(file_list), (D-1)*(D - 2)/2) # minus reference
  hosts <- character(length(file_list))
  if(iter == 1) {
    # MAP estimates
    for(f in 1:length(file_list)) {
      file <- file_list[f]
      fit <- readRDS(file.path(output_dir_full, file))
      # Convert to CLR
      if(fit$coord_system != "clr") {
        fit <- to_clr(fit)
      }
      if(use_proportionality) {
        Sigma <- fit$Sigma[,,1]
        Sigma <- Sigma[1:(D-1),1:(D-1)]
        for(m in 1:length(pair1)) {
          j <- pair1[m]
          k <- pair2[m]
          var_j <- Sigma[j,j] # diagonal element
          var_k <- Sigma[k,k] # diagonal element
          var_j_minus_k <- var_j + var_k - 2*Sigma[j,k]
          rho_jk <- 1 - var_j_minus_k / (var_j + var_k)
          rug[f,m] <- rho_jk
        }
      } else {
        Sigma_correlation <- cov2cor(fit$Sigma[,,1])
        Sigma_correlation <- Sigma_correlation[1:(D-1),1:(D-1)]
        vector_Sigma <- Sigma_correlation[lower.tri(Sigma_correlation,
                                                    diag = FALSE)]
        rug[f,] <- vector_Sigma
      }
      hosts[f] <- fit$sname
    }
  } else {
    # Full posteriors
    for(f in 1:length(file_list)) {
      file <- file_list[f]
      fit <- readRDS(file.path(output_dir_full, file))
      # Convert to CLR
      if(fit$coord_system != "clr") {
        fit <- to_clr(fit)
      }
      if(use_proportionality) {
        Sigma <- apply(fit$Sigma, c(1,2), mean)
        Sigma <- Sigma[1:(D-1),1:(D-1)]
        for(m in 1:length(pair1)) {
          j <- pair1[m]
          k <- pair2[m]
          var_j <- Sigma[j,j] # diagonal element
          var_k <- Sigma[k,k] # diagonal element
          var_j_minus_k <- var_j + var_k - 2*Sigma[j,k]
          rho_jk <- 1 - var_j_minus_k / (var_j + var_k)
          rug[f,m] <- rho_jk
        }
      } else {
        Sigma_correlation <- fit$Sigma
        for(i in 1:fit$iter) {
          Sigma_correlation[,,i] <- cov2cor(Sigma_correlation[,,i])
        }
        Sigma_summary <- apply(Sigma_correlation, c(1,2), mean)
        Sigma_summary <- Sigma_summary[1:(D-1),1:(D-1)]
        vector_Sigma <- Sigma_summary[lower.tri(Sigma_summary, diag = FALSE)]
        rug[f,] <- vector_Sigma
      }
      hosts[f] <- fit$sname
    }
  }
  return(list(hosts = hosts, rug = rug, tax_idx1 = pair1, tax_idx2 = pair2))
}

#' Renders the "rug" (host x pairwise association matrix) as a heatmap.
#'
#' @param rug heatmap matrix output from summarize_Sigmas()
#' @param canonical_col_order if not NULL, order of pairs (columns) to use
#' @param canonical_row_order if not NULL, order of rows (hosts) to use
#' @param row_labels optional parameter with row labels (ordered host short
#' names)
#' @param save_name name with which to save heatmap file
#' @return named list with column and row ordering
#' @import ggplot2
#' @import tidyr
#' @export
plot_rug <- function(rug, canonical_col_order = NULL,
                     canonical_row_order = NULL, row_labels = NULL,
                     save_name = NULL) {
  # Cluster
  if(is.null(canonical_row_order)) {
    d <- dist(rug)
    canonical_row_order <- hclust(d)$order
  }
  if(is.null(canonical_col_order)) {
    d <- dist(t(rug))
    canonical_col_order <- hclust(d)$order
  }
  rug <- rug[canonical_row_order,canonical_col_order]

  rug <- cbind(1:nrow(rug), rug)
  colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
  rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
  rug$pair <- as.numeric(rug$pair)

  p <- ggplot(rug, aes(x = pair, y = host)) +
    geom_raster(aes(fill = correlation)) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                         midpoint = 0)
  if(!is.null(row_labels)) {
    p <- p +
      scale_y_continuous(breaks = 1:length(row_labels),
                         labels = row_labels) +
      theme(axis.text = element_text(size = 4))
  }
  if(is.null(save_name)) {
    show(p)
  } else {
    output_dir <- check_dir(c("output", "figures"))
    ggsave(file.path(output_dir, paste0(save_name, ".png")),
           p,
           units = "in",
           height = 3,
           width = 8)
  }
  return(list(row_order = canonical_row_order,
              col_order = canonical_col_order))
}

#' Renders the host x pairwise associations as a histogram
#'
#' @param rug heatmap matrix output from summarize_Sigmas()
#' @param save_name name with which to save heatmap file
#' @return named list with column and row ordering
#' @import ggplot2
#' @export
plot_correlation_histogram <- function(rug, save_name = NULL) {
  plot_df <- data.frame(x = c(rug))
  p <- ggplot(plot_df, aes(x = x)) +
    geom_histogram(bins = 30, color = "white") +
    xlim(c(-1,1)) +
    xlab("CLR correlation")
  if(is.null(save_name)) {
    show(p)
  } else {
    output_dir <- check_dir(c("output", "figures"))
    ggsave(file.path(output_dir, paste0(save_name, ".png")),
           p,
           units = "in",
           height = 3,
           width = 4)
  }
}

#' Utility function for pulling collection date-aligned trajectories from a
#' given host's prediction set.
#'
#' @import fido
get_trajectory_df <- function(dates, coord, common_baseline, predictions,
                              center = NULL) {
  # Get days from common baseline
  days <- round(unname(sapply(dates, function(x) {
      difftime(x, common_baseline, units = "days")
    }))) + 1
  # Recalculate the "days" index for this host given the new baseline
  offset <- days[1] - 1

  days <- predictions$span + offset
  eta_df <- gather_array(predictions$Eta, val, coord, sample, iteration)
  map <- data.frame(sample = 1:max(eta_df$sample), day = days)
  eta_df <- left_join(eta_df, map, by = "sample")
  eta_df <- eta_df %>%
    select(coord, day, iteration, val) %>%
    group_by(coord, day) %>%
    summarize(p2.5 = quantile(val, probs = c(0.025)),
              p25 = quantile(val, probs = c(0.25)),
              mean = mean(val),
              p75 = quantile(val, probs = c(0.75)),
              p97.5 = quantile(val, probs = c(0.975)), .groups = "keep")
  plot_df <- eta_df[eta_df$coord == coord,]
  if(!is.null(center)) {
    new_center <- mean(plot_df$mean) + center
    plot_df <- plot_df %>%
      mutate(p2.5 = p2.5 - new_center,
             p25 = p25 - new_center,
             p75 = p75 - new_center,
             p97.5 = p97.5 - new_center,
             mean = mean - new_center)
  }
  return(plot_df)
}

#' Utility function for pulling collection date-aligned trajectories for a pair
#' of taxa from a given host's prediction set.
get_paired_trajectories <- function(dates, coord1, coord2, host, common_baseline,
                                    predictions, center = NULL) {
  pair_df_tax1 <- get_trajectory_df(dates,
                                    coord1,
                                    common_baseline,
                                    predictions,
                                    center = center)
  pair_df_tax1 <- cbind(host = host, pair_df_tax1)
  pair_df_tax2 <- get_trajectory_df(dates,
                                    coord2,
                                    common_baseline,
                                    predictions,
                                    center = center)
  pair_df_tax2 <- cbind(host = host, pair_df_tax2)
  pair_df <- rbind(pair_df_tax1, pair_df_tax2)
  return(pair_df)
}

#' Utility function for pulling collection date-aligned trajectories for a pair
#' of taxa from a set of hosts.
get_predictions_host_list <- function(ref_hosts, output_dir, metadata) {
  prediction_list <- list()
  date_list <- list()
  for(host in ref_hosts) {
    predictions <- load_predictions(host, output_dir, generate = TRUE)
    dates <- metadata %>%
      filter(sname == host) %>%
      pull(collection_date)
    prediction_list[[host]] <- predictions
    date_list[[host]] <- dates
  }
  # Find common baseline date
  common_baseline <- min(unname(unlist(date_list)))
  return(list(common_baseline = common_baseline,
              predictions = prediction_list,
              dates = date_list))
}

#' Render date-aligned trajectories for a pair of taxa and a given set of hosts.
#' This allows us to eyeball whether taxa strongly correlated within hosts show
#' similar dynamics across hosts.
#'
#' @param output_dir specifies the subdirectory of the model fits directory in
#' which to save the predictions
#' @param tax_idx1 index of first taxon
#' @param tax_idx2 index of second taxon
#' @param tax_label1 readable label for first taxon
#' @param tax_label2 readable label for second taxon
#' @param metadata metadata data.frame (with sname and collection_date column)
#' @param bind_taxa plot all trajectories for taxon 1 together; ditto taxon 2
#' @param save_file flag indicating whether or not to save rendered plot to file
#' @return NULL
#' @import dplyr
#' @export
plot_aligned_trajectories <- function(output_dir, tax_idx1, tax_idx2,
                                      tax_label1, tax_label2, metadata,
                                      bind_taxa = FALSE, save_file = FALSE) {
  # Pull "reference" hosts
  host_list <- get_reference_hosts()$sname
  # Pull models/predictions
  pred_obj <- get_predictions_host_list(host_list, output_dir, metadata)
  if(bind_taxa) {
    host_centers <- rep(0, length(host_list))
  } else {
    host_centers <- seq(from = 0, by = 10, length.out = length(host_list))
  }
  names(host_centers) <- sort(host_list, decreasing = TRUE)
  pair_df <- NULL
  for(host in host_list) {
    pair_df_host <- get_paired_trajectories(pred_obj$dates[[host]],
                                            tax_idx1,
                                            tax_idx2,
                                            host,
                                            pred_obj$common_baseline,
                                            pred_obj$predictions[[host]],
                                            center = host_centers[[host]])
    if(is.null(pair_df)) {
      pair_df <- pair_df_host
    } else {
      pair_df <- rbind(pair_df, pair_df_host)
    }
  }

  # The pair_df data.frame contains the interpolated trajectories (50% CI and 95%
  # CI) for a pair of taxa in two hosts.

  p <- ggplot()
  if(bind_taxa) {
    for(ref_host in host_list) {
      p <- p +
        geom_ribbon(data = pair_df[pair_df$host == ref_host & pair_df$coord == tax_idx1, ],
                    aes(x = day, ymin = p25, ymax = p75, fill = factor(host)),
                    # fill = "red",
                    alpha = 0.33) +
        geom_ribbon(data = pair_df[pair_df$host == ref_host & pair_df$coord == tax_idx2, ],
                    aes(x = day, ymin = p25 - 10, ymax = p75 - 10, fill = factor(host)),
                    # fill = "blue",
                    alpha = 0.33) +
        labs(fill = "Host")
    }
    p <- p + scale_y_continuous(breaks = c(0, -10),
                                labels = c("taxon 1", "taxon 2")) +
      labs(title = paste0("Correlation: ", tax_label1, " x ", tax_label2))
    if(save_file) {
      filename <- paste0("aligned_series_coords_", tax_idx1, "-", tax_idx2, "_",
                         paste0(host_list, collapse = "-"), "_boundtaxa.png")
      save_dir <- check_dir(c("output", "figures"))
      ggsave(file.path(save_dir, filename),
             p,
             units = "in",
             height = 3,
             width = 10)
    } else {
      show(p)
    }
  } else {
    for(ref_host in host_list) {
      p <- p +
        geom_ribbon(data = pair_df[pair_df$host == ref_host & pair_df$coord == tax_idx1, ],
                    aes(x = day, ymin = p25, ymax = p75),
                    fill = "red",
                    alpha = 0.33) +
        geom_ribbon(data = pair_df[pair_df$host == ref_host & pair_df$coord == tax_idx2, ],
                    aes(x = day, ymin = p25, ymax = p75),
                    fill = "blue",
                    alpha = 0.33)
    }
    p <- p + scale_y_continuous(breaks = -unname(host_centers),
                                labels = sort(host_list, decreasing = TRUE)) +
      labs(title = paste0("Correlation: ", tax_label1, " x ", tax_label2))
    if(save_file) {
      filename <- paste0("aligned_series_coords_", tax_idx1, "-", tax_idx2, "_",
                         paste0(host_list, collapse = "-"), ".png")
      save_dir <- check_dir(c("output", "figures"))
      ggsave(file.path(save_dir, filename),
             p,
             units = "in",
             height = length(host_list),
             width = 10)
    } else {
      show(p)
    }
  }
}

#' Utility function to plot aligned series for the pair of taxa with the largest
#' universality score from among a subset of association indices
render_universal_pairs <- function(output_dir, select_idx, scores,
                                   tax_idx1, tax_idx2, taxonomy) {
  select_scores <- scores[select_idx]
  select_assoc_tax1 <- tax_idx1[select_idx]
  select_assoc_tax2 <- tax_idx2[select_idx]

  max_idx <- which(select_scores == max(select_scores))
  pair <- c(select_assoc_tax1[max_idx], select_assoc_tax2[max_idx])

  for(bind_taxa in c(FALSE, TRUE)) {
    plot_aligned_trajectories(output_dir = output_dir,
                              tax_idx1 = pair[1],
                              tax_idx2 = pair[2],
                              tax_label1 = get_tax_label(taxonomy,
                                                         pair[1],
                                                         "CLR"),
                              tax_label2 = get_tax_label(taxonomy,
                                                         pair[2],
                                                         "CLR"),
                              metadata = metadata,
                              bind_taxa = bind_taxa,
                              save_file = TRUE)
  }
}

#' Render date-aligned trajectories for the most universal positive and negative
#' pairs of taxa in the top-sampled hosts within each social group
#'
#' @param output_dir specifies the subdirectory of the model fits directory in
#' which to save the predictions
#' @param metadata metadata data.frame (with sname and collection_date column)
#' @param taxonomy appropriate taxonomy object for this level
#' @return NULL
#' @export
plot_trajectories_top_pairs <- function(output_dir, metadata, taxonomy) {
  rug_obj <- summarize_Sigmas(output_dir = output_dir)
  rug <- rug_obj$rug
  scores <- apply(rug, 2, calc_universality_score)
  # Get prevailing sign of association
  agreement_signs <- apply(rug, 2, function(x) {
    sum(sign(x))
  })
  consensus_signs <- sign(agreement_signs)
  # Most universal positive associations
  render_universal_pairs(select_idx = which(consensus_signs > 0),
               scores = scores,
               tax_idx1 = rug_obj$tax_idx1,
               tax_idx2 = rug_obj$tax_idx2,
               taxonomy = taxonomy,
               output_dir = output_dir)

  # Most universal negative associations
  render_universal_pairs(select_idx = which(consensus_signs < 0),
               scores = scores,
               tax_idx1 = rug_obj$tax_idx1,
               tax_idx2 = rug_obj$tax_idx2,
               taxonomy = taxonomy,
               output_dir = output_dir)
}
