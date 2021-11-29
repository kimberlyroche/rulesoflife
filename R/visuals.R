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
#' @return a ggplot object
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
    return(p)
  } else {
    output_dir <- check_dir(c("output", "figures"))
    ggsave(file.path(output_dir, paste0(save_name, ".png")),
           plot = p,
           units = "in",
           dpi = 100,
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
#' @param resolution percent of missing days to interpolate; default is 5% which
#' is plenty good enough for plotting
#' @return NULL
#' @import fido
#' @export
predict_trajectory <- function(sname, output_dir, resolution = 5) {
  fit <- load_fit(sname = sname, MAP = FALSE, output_dir = output_dir)
  if(fit$coord_system != "clr") {
    fit <- to_clr(fit)
  }
  first_day <- min(fit$X[1,])
  last_day <- max(fit$X[1,])
  n_days <- last_day
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
#' @return NULL
#' @import driver
#' @import fido
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @export
plot_trajectory <- function(sname, output_dir, coord, coord_label = NULL,
                            show_observations = FALSE) {
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
  filename <- paste0("GP_", fit$sname, "_", "coord-", coord, ".png")

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
  output_dir <- check_dir(c("output", "figures"))
  ggsave(file.path(output_dir, filename), p, units = "in", height = 3, width = 10)
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
#' @param cluster_obj an optional pre-computed hierarchical clustering; if
#' included, a dendrogram will be rendered from this
#' @param save_name name with which to save heatmap file
#' @return named list with column and row ordering
#' @import ggplot2
#' @import tidyr
#' @import grid
#' @import ggdendro
#' @import cowplot
#' @export
plot_rug <- function(rug,
                     canonical_col_order = NULL,
                     canonical_row_order = NULL,
                     row_labels = NULL,
                     cluster_obj = NULL,
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
  row_labels <- row_labels[canonical_row_order]

  rug <- cbind(1:nrow(rug), rug)
  colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
  rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
  rug$pair <- as.numeric(rug$pair)

  p <- ggplot(rug, aes(x = pair, y = host)) +
    geom_raster(aes(fill = correlation)) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                         midpoint = 0) +
    labs(y = "")

  if(is.null(row_labels)) {
    row_labels <- 1:length(rug$host)
  }
  p <- p +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:length(row_labels),
                       labels = row_labels,
                       expand = c(0, 0)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(fill = "correlation\nmetric")

  if(!is.null(cluster_obj)) {
    p <- p +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 14, face = "plain"),
            legend.title = element_text(size = 14))
    dhc <- as.dendrogram(cluster_obj)
    ddata <- dendro_data(dhc, type = "rectangle")
    pd <- ggplot(segment(ddata)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      coord_flip() +
      scale_y_reverse() +
      theme_nothing()
    if(!is.null(save_name)) {
      png(save_name, height = 600, width = 1200)
    }
    grid.newpage()
    print(pd, vp = viewport(x = 0.10, y = 0.515, width = 0.15, height = 0.995))
    print(p, vp = viewport(x = 0.57, y = 0.5, width = 0.8, height = 0.97))
    if(!is.null(save_name)) {
      dev.off()
    }
  } else {
    p <- p +
      theme(axis.text = element_text(size = 7),
            axis.title = element_text(size = 12, face = "plain"),
            legend.title = element_text(size = 12))
    if(is.null(save_name)) {
      show(p)
    } else {
      ggsave(file.path(output_dir, save_name),
             plot = p,
             dpi = 100,
             units = "in",
             height = 6,
             width = 12)
    }
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

#' Pull collection date-aligned trajectories from a given host's prediction set
#'
#' @param dates named list of dates in YYYY-MM-DD format indexed by host short
#' name
#' @param coord index of taxon
#' @param common_baseline earliest date to use for all hosts (YYYY-MM-DD format)
#' @param predictions named list of predicted (interpolated) Eta matrices
#' indexed by host short name
#' @param center center to use for taxa (e.g. 0)
#' @param mean_only return mean rather than quantiles
#' @return NULL
#' @import driver
#' @import fido
#' @export
get_trajectory_df <- function(dates, coord, common_baseline, predictions,
                              center = NULL, mean_only = FALSE) {
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
    group_by(coord, day)
  if(mean_only) {
    eta_df <- eta_df %>%
      summarize(mean = mean(val), .groups = "keep")
  } else {
    eta_df <- eta_df %>%
      summarize(p2.5 = quantile(val, probs = c(0.025)),
                p25 = quantile(val, probs = c(0.25)),
                mean = mean(val),
                p75 = quantile(val, probs = c(0.75)),
                p97.5 = quantile(val, probs = c(0.975)), .groups = "keep")
  }
  plot_df <- eta_df[eta_df$coord == coord,]
  if(!is.null(center)) {
    new_center <- mean(plot_df$mean) + center
    if(mean_only) {
      plot_df <- plot_df %>%
        mutate(mean = mean - new_center)
    } else {
      plot_df <- plot_df %>%
        mutate(p2.5 = p2.5 - new_center,
               p25 = p25 - new_center,
               p75 = p75 - new_center,
               p97.5 = p97.5 - new_center,
               mean = mean - new_center)
    }
  }
  return(plot_df)
}

#' Pull a collection date-aligned trajectories for a pair of taxa from a given
#' host's prediction set (list returned by get_predictions_host_list())
#'
#' @param dates named list of dates in YYYY-MM-DD format indexed by host short
#' name
#' @param coord1 index of taxon 1
#' @param coord2 index of taxon 2
#' @param host host short name
#' @param common_baseline earliest date to use for all hosts (YYYY-MM-DD format)
#' @param predictions named list of predicted (interpolated) Eta matrices
#' indexed by host short name
#' @param center center to use for taxa (e.g. 0)
#' @param mean_only return mean rather than quantiles
#' @return NULL
#' @import dplyr
#' @export
get_paired_trajectories <- function(dates, coord1, coord2, host, common_baseline,
                                    predictions, center = NULL, mean_only = FALSE) {
  pair_df_tax1 <- get_trajectory_df(dates,
                                    coord1,
                                    common_baseline,
                                    predictions,
                                    center = center,
                                    mean_only = mean_only)
  pair_df_tax1 <- cbind(host = host, pair_df_tax1)
  pair_df_tax2 <- get_trajectory_df(dates,
                                    coord2,
                                    common_baseline,
                                    predictions,
                                    center = center,
                                    mean_only = mean_only)
  pair_df_tax2 <- cbind(host = host, pair_df_tax2)
  pair_df <- rbind(pair_df_tax1, pair_df_tax2)
  return(pair_df)
}

#' Pulling a collection date-aligned trajectories for a pair of taxa from a set
#' of hosts
#'
#' @param host_list list of host short names
#' @param output_dir specifies the subdirectory of the model fits directory in
#' which to save the predictions
#' @param metadata metadata data.frame (with sname and collection_date column)
#' @return NULL
#' @import dplyr
#' @export
get_predictions_host_list <- function(host_list, output_dir, metadata) {
  prediction_list <- list()
  date_list <- list()
  for(host in host_list) {
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
#' @param return_plot if TRUE, returns the ggplot2 object instead of saving it
#' @param file_tag if not NULL, a tag to append to the saved filename
#' @return NULL
#' @import dplyr
#' @export
plot_aligned_trajectories <- function(output_dir, tax_idx1, tax_idx2,
                                      tax_label1, tax_label2, metadata,
                                      bind_taxa = FALSE, return_plot = FALSE,
                                      file_tag = NULL) {
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
    cat("Getting paired trajectory for",host,"\n")
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
    if(return_plot) {
      return(p)
    } else {
      filename <- paste0("aligned_series_coords_", tax_idx1, "-", tax_idx2, "_",
                         paste0(host_list, collapse = "-"), "_boundtaxa")
      if(!is.null(file_tag)) {
        filename <- paste0(filename, "_", file_tag)
      }
      filename <- paste0(filename, ".png")
      save_dir <- check_dir(c("output", "figures"))
      ggsave(file.path(save_dir, filename),
             p,
             units = "in",
             height = 3,
             width = 10)
    }
  } else {
    p <- ggplot()
    for(ref_host in host_list) {
      p <- p +
        geom_ribbon(data = pair_df[pair_df$host == ref_host & pair_df$coord == tax_idx1, ],
                    aes(x = day, ymin = p25, ymax = p75),
                    fill = "#2a9d8f",
                    alpha = 0.8) +
        geom_ribbon(data = pair_df[pair_df$host == ref_host & pair_df$coord == tax_idx2, ],
                    aes(x = day, ymin = p25, ymax = p75),
                    fill = "#fb8500",
                    alpha = 0.8)
    }
    p <- p + scale_y_continuous(breaks = -unname(host_centers),
                                labels = sort(host_list, decreasing = TRUE)) +
      labs(title = paste0("Correlation: ", tax_label1, " x ", tax_label2),
           x = "day from common baseline")

    if(return_plot) {
      return(p)
    } else {
      filename <- paste0("aligned_series_coords_", tax_idx1, "-", tax_idx2, "_",
                         paste0(host_list, collapse = "-"))
      if(!is.null(file_tag)) {
        filename <- paste0(filename, "_", file_tag)
      }
      filename <- paste0(filename, ".png")
      save_dir <- check_dir(c("output", "figures"))
      ggsave(file.path(save_dir, filename),
             p,
             units = "in",
             height = length(host_list),
             width = 10)
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
                              bind_taxa = bind_taxa)
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
  consensus_signs <- apply(rug, 2, calc_consensus_sign)
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

#' Render line plot of a given CLR taxon against diet PC1 for a given host.
#' This plot is meant to demonstrate the largely uncorrelated nature of the diet
#' data with the observed logratio abundances.
#'
#' @param sname host short name indicating which baboon's series to fit
#' @param tax_idx symmetric matrix object
#' @param counts filtered 16S count table (taxa x samples)
#' @param metadata annotations data.frame
#' @return NULL
#' @import driver
#' @import ggplot2
#' @export
plot_clr_vs_diet <- function(sname, tax_idx, counts, metadata) {
  sname_idx <- which(metadata$sname == sname)
  sub_md <- metadata[sname_idx,]
  sub_counts <- counts[,sname_idx]
  clr_counts <- clr_array(sub_counts + 0.5, parts = 1)
  plot_df <- data.frame(time = rep(1:nrow(sub_md), 2),
                        value = c(scale(sub_md$diet_PC1, center = TRUE, scale = FALSE),
                                  scale(clr_counts[tax_idx,], center = TRUE, scale = FALSE)),
                        type = c(rep("diet", nrow(sub_md)), rep("taxon", nrow(sub_md))))
  plot_df$type <- factor(plot_df$type)
  p <- ggplot(plot_df, aes(x = time, y = value, color = type)) +
    geom_line(size = 0.5) +
    geom_point() +
    xlab("sample index")
  filename <- paste0("dietPC1_vs_tax", tax_idx, "_", sname, ".png")
  output_dir <- check_dir(c("output", "figures"))
  ggsave(file.path(output_dir, filename),
         p,
         units = "in",
         height = 3,
         width = 6)
}

#' Plot relative representation of ASV pairs (using family labels)
#'
#' @param frequencies_subset table of ASV pairs (labeled as "Family - Family")
#' and their frequency in some subset of case of interest
#' @param frequencies table of ASV pairs (labeled as "Family - Family")
#' and their pverall frequency
#' @param plot_height plot height in inches
#' @param plot_width plot width in inches
#' @param legend_topmargin top margin of legend in points
#' @param legend_leftmargin left margin of legend in inches
#' @param use_pairs if TRUE, renders relative abundances of pairs of taxa,
#' otherwise, renders families
#' @param save_name if NULL, ggplot object is returned
#' @return NULL
#' @import ggplot2
#' @import cowplot
#' @export
plot_enrichment <- function(frequencies_subset, frequencies, plot_height, plot_width,
                            legend_topmargin, legend_leftmargin, use_pairs = TRUE, save_name = NULL) {
  # Define a huge color palette over all observed family-family pairs
  if(use_pairs) {
    palette_fn <- file.path("output", "family-family_palette.rds")
    if(file.exists(palette_fn)) {
      fam_palette <- readRDS(palette_fn)
    } else {
      fam_palette <- generate_highcontrast_palette(length(frequencies))
      names(fam_palette) <- names(frequencies)
      saveRDS(fam_palette, palette_fn)
    }
  } else {
    palette_fn <- file.path("output", "family_palette.rds")
    if(file.exists(palette_fn)) {
      fam_palette <- readRDS(palette_fn)
    } else {
      fam_palette <- generate_highcontrast_palette(length(frequencies))
      names(fam_palette) <- names(frequencies)
      saveRDS(fam_palette, palette_fn)
    }
  }

  observed_relative <- data.frame(pair = names(frequencies_subset),
                                  count = unname(c(frequencies_subset)))
  observed_relative$prop <- observed_relative$count / sum(observed_relative$count)

  observed_relative$pair <- factor(observed_relative$pair,
                                   levels = c(sort(observed_relative$pair[observed_relative$pair != "Other"]), "Other"))

  fam_palette_sub <- fam_palette[names(fam_palette) %in% observed_relative$pair]

  p2 <- ggplot(observed_relative, aes(x = 1, y = prop, fill = pair)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = fam_palette_sub) +
    theme_bw() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.justification = c("top"),
          legend.margin = margin(t = legend_topmargin, r = 0, b = 0, l = 0, unit = "pt"))
  if(use_pairs) {
    p2 <- p2 +
    labs(x = "\nlow phylogenetic distance,\nhigh median assoc. strength",
         y = "proportion shared family ASV pair",
         fill = "Shared family ASV pair")
  } else {
    p2 <- p2 +
      labs(x = "\nlow phylogenetic distance,\nhigh median assoc. strength",
           y = "proportion shared ASVs from family",
           fill = "ASV family")
  }

  # Stacked bar: render these pairs in the baseline representation of all
  # family-family pairs; all other family-family pairs will be colored gray

  baseline_relative <- observed_relative %>%
    full_join(data.frame(pair = names(frequencies),
                         count = unname(c(frequencies))), by = "pair") %>%
    select(pair, count.y) %>%
    arrange(pair)
  colnames(baseline_relative) <- c("pair", "count")
  baseline_relative$prop <- baseline_relative$count / sum(baseline_relative$count)

  fam_palette_gray <- fam_palette
  repl_names <- baseline_relative$pair[which(!(baseline_relative$pair %in% observed_relative$pair))]
  repl_idx <- which(names(fam_palette_gray) %in% repl_names)
  grays <- sample(c("#dddddd", "#d5d5d5", "#cccccc", "#c5c5c5"), replace = TRUE, size = length(repl_idx))
  fam_palette_gray[repl_idx] <- grays

  p1 <- ggplot(baseline_relative, aes(x = 1, y = prop, pair, fill = pair)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = fam_palette_gray) +
    theme_bw() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  if(use_pairs) {
    p1 <- p1 +
      labs(x = "\noverall\n",
           y = "proportion shared family ASV pair")
  } else {
    p1 <- p1 +
      labs(x = "\noverall\n",
           y = "proportion shared ASVs from family")
  }

  legend <- get_legend(p2)
  p2 <- p2 +
    theme(legend.position = "none")

  p <- plot_grid(p1, NULL, p2, NULL, legend, ncol = 5, rel_widths = c(1, 0.35, 1, legend_leftmargin, 4.5), scale = 1)

  if(!is.null(save_name)) {
    ggsave(file.path("output", "figures", save_name),
           p,
           dpi = 100,
           units = "in",
           height = plot_height,
           width = plot_width)
  } else {
    p
  }
}


