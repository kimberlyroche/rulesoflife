#' Draw predictions from a bassetfit or labraduckfit object over a given day
#' range
#'
#' @param fit bassetfit or labraduckfit object (from fido)
#' @param span optional day range over which to predict
#' @return NULL
#' @import fido
#' @export
predict_span <- function(fit, span = NULL) {
  if("bassetfit" %in% class(fit)) {
    if(is.null(span)) {
      span <- min(fit$X):max(fit$X)
    }
    # Check to see if this prediction has already been rendered
    output_dir <- check_dir(c("output", "GP_fit", "predictions"))
    filename <- paste0(fit$sname, "_", min(span), "-", max(span), ".rds")
    dim(span) <- c(1, length(span))
    pred <- predict(fit, newdata = span, response = "Eta") # or LambdaX or Y
    predictions <- list(Eta = pred,
                        span = c(span),
                        coord_system = fit$coord_system)
    saveRDS(predictions, file = file.path(output_dir, filename))
    return(predictions)
  } else {
    stop("predict_span not implemented for types other than bassetfit yet!")
  }
}

#' Visualize inferred logratio trajectories
#'
#' @param fit bassetfit or labraduckfit object (from fido)
#' @param predictions named list of predictions from a bassetfit or
#' labraduckfit object
#' @param taxon index of taxon to plot
#' @param taxon_label optional taxon label (e.g. "family Lachnospiraceae")
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
plot_trajectory <- function(fit, predictions, taxon, taxon_label = NULL,
                            show_observations = FALSE, save_file = FALSE) {
  if(dim(predictions$Eta)[3] < 100) {
    stop("Fit object doesn't have at least 100 posterior samples!\n")
  }
  if("bassetfit" %in% class(fit)) {
    # Gather predictions into intervals
    days <- predictions$span
    eta_df <- gather_array(predictions$Eta, val, taxon, sample, iteration)
    map <- data.frame(sample = 1:max(eta_df$sample), day = days)
    eta_df <- left_join(eta_df, map, by = "sample")
    eta_df <- eta_df %>%
      select(taxon, day, iteration, val) %>%
      group_by(taxon, day) %>%
      summarize(p2.5 = quantile(val, probs = c(0.025)),
                p25 = quantile(val, probs = c(0.25)),
                mean = mean(val),
                p75 = quantile(val, probs = c(0.75)),
                p97.5 = quantile(val, probs = c(0.975)), .groups = "keep")
    plot_df <- eta_df[eta_df$taxon == taxon,]

    if(show_observations) {
      # Logratio transform real data
      observations <- fit$Y
      if(predictions$coord_system == "clr") {
        observations <- t(driver::clr((t(observations) + 0.5)))
      } else if(predictions$coord_system == "alr") {
        observations <- t(driver::alr((t(observations) + 0.5)))
      } else {
        stop("Predictions are in unrecognized coordinate system!")
      }
      observations <- data.frame(day = c(fit$X), logratio = observations[taxon,])

      plot_df <- left_join(plot_df, observations, by = "day")
    }
    p <- ggplot(plot_df) +
      geom_ribbon(aes(x = day, ymin = p2.5, ymax = p97.5), fill = "grey80") +
      geom_ribbon(aes(x = day, ymin = p25, ymax = p75), fill = "grey60") +
      geom_line(aes(x = day, y = mean))
    if(!is.null(taxon_label)) {
      p <- p +
        ylab(paste0("CLR(", taxon_label, ")"))
    } else {
      p <- p +
        ylab(paste0("CLR(taxon ", taxon, ")"))
    }
    if(show_observations) {
      p <- p +
        geom_point(aes(x = day, y = logratio), size = 1, na.rm = TRUE)
    }
    if(save_file) {
      filename <- paste0(fit$sname, "_", min(span), "-", max(span), "_",
                         "taxon", taxon, ".png")
      output_dir <- check_dir(c("output", "images"))
      width <- 10 * (max(predictions$span) / 4500)
      if(width < 4) {
        width <- 4
      }
      ggsave(file.path(output_dir, filename), p, units = "in", height = 3, width = width)
    } else {
      show(p)
    }
  }
}


