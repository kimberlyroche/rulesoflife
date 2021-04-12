#' Draw predictions from a fido::basset model in the CLR coordinate system.
#' This function interpolates 5% of missing observations, which gives good quick
#' visuals in the ABRP data.
#'
#' @param fit bassetfit object
#' @return NULL
#' @import fido
#' @export
predict_trajectory <- function(fit) {
  if("bassetfit" %in% class(fit)) {
    if(fit$coord_system != "clr") {
      fit <- to_clr(fit)
    }
    first_day <- min(fit$X)
    last_day <- max(fit$X)
    n_days <- last_day
    resolution <- 5
    span <- round(seq(from = first_day, to = last_day, length.out = round(n_days * resolution/100)))
    span <- sort(c(span, c(fit$X))) # include observations
    dim(span) <- c(1, length(span))
    pred <- predict(fit, newdata = span, response = "Eta") # or LambdaX or Y
    predictions <- list(Eta = pred, # taxa x time points x n_iter
                        span = c(span),
                        resolution = resolution)
    output_dir <- check_dir(c("output", "GP_fits", "predictions"))
    filename <- paste0(fit$sname, ".rds")
    saveRDS(predictions, file = file.path(output_dir, filename))
    return(predictions)
  } else {
    stop("Fit object is not of bassetfit class!")
  }
}

#' Visualize inferred logratio trajectories. GP objects will be visualized in
#' the CLR coordinate system; DLM objects will be visualized in the ALR (for
#' now).
#'
#' @param fit bassetfit or labraduckfit object (from fido)
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
plot_trajectory <- function(fit, coord, coord_label = NULL,
                            show_observations = FALSE, save_file = FALSE) {
  if("bassetfit" %in% class(fit)) {
    output_dir <- check_dir(c("output", "GP_fits", "predictions"))
    filename <- paste0(fit$sname, ".rds")
    if(file.exists(file.path(output_dir, filename))) {
      predictions <- readRDS(file.path(output_dir, filename))
    } else {
      predictions <- predict_trajectory(fit)
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
  } else if("labraduckfit" %in% class(fit)) {
    output_dir <- check_dir(c("output", "DLM_fits", "predictions"))
    filename <- paste0(fit$sname, "_", "coord-", coord, ".rds")
    if(file.exists(file.path(output_dir, filename))) {
      predictions <- readRDS(file.path(output_dir, filename))
    } else {
      output_dir <- check_dir(c("output", "DLM_fits", "predictions"))
      filename <- paste0(fit$sname, ".rds")
      if(file.exists(file.path(output_dir, filename))) {
        predictions <- readRDS(file.path(output_dir, filename))
      } else {
        if(max(fit$Eta_DLM) == 0) {
          stop("Visualization of filtered models not implemented yet!")
        }
        # Otherwise the "prediction" has already been done by the smoother
        # We just need to extract the appropriate quantities for the coordinate
        # of interest
        F <- fit$F
        Y <- fit$Y
        T <- max(fit$observations)
        N <- ncol(Y)
        D <- nrow(Y)
        n_samples <- dim(fit$Eta)[3]
        Ft <- t(F)
        Q <- nrow(F)
        Eta_samples <- matrix(NA, T, n_samples)
        for(k in 1:n_samples) {
          EtasS_1T <- fit$Eta_DLM[,k]
          dim(EtasS_1T) <- c(Q, D-1, T)
          Eta_samples[,k] <- EtasS_1T[,coord,]
        }
        # Thin samples for easier plotting
        resolution <- 20
        span <- round(seq(from = min(fit$observations),
                          to = T,
                          length.out = round(T * resolution/100)))
        span <- sort(c(span, c(fit$observations))) # include observations
        predictions <- list(Eta = Eta_samples[span,], # time points x n_iter
                            span = span,
                            resolution = resolution)
        saveRDS(predictions, file = file.path(output_dir, filename))
      }
    }

    # Gather predictions into intervals
    days <- predictions$span
    eta_df <- gather_array(predictions$Eta, val, sample, iteration)
    map <- data.frame(sample = 1:max(eta_df$sample), day = days)
    eta_df <- left_join(eta_df, map, by = "sample")
    eta_df <- eta_df %>%
      select(day, iteration, val) %>%
      group_by(day) %>%
      summarize(p2.5 = quantile(val, probs = c(0.025)),
                p25 = quantile(val, probs = c(0.25)),
                mean = mean(val),
                p75 = quantile(val, probs = c(0.75)),
                p97.5 = quantile(val, probs = c(0.975)), .groups = "keep")
    plot_df <- eta_df

    if(show_observations) {
      # Logratio transform real data
      observations <- fit$Y
      observations <- t(driver::alr((t(observations) + 0.5), d = fit$D))
      observations <- data.frame(day = c(fit$observations), logratio = observations[coord,])

      plot_df <- left_join(plot_df, observations, by = "day")
    }
    filename <- paste0("DLM_", fit$sname, "_", "coord-", coord, ".png")
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
    output_dir <- check_dir(c("output", "images"))
    ggsave(file.path(output_dir, filename), p, units = "in", height = 3, width = 10)
  } else {
    show(p)
  }
}

