source("path_fix.R")

library(fido)
library(ggplot2)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure SX - consistency of parameterizations
#
# ------------------------------------------------------------------------------

alt_sigmas <- function(path, fits) {
  Sigmas <- NULL
  for(i in 1:length(fits)) {
    fit <- fits[i]
    cat(paste0("Loading ", fit, "...\n"))
    fit <- readRDS(file.path(path, fit))
    mean_clr_sigma <- apply(to_clr(fit)$Sigma, c(1,2), mean)
    mean_clr_cor <- cov2cor(mean_clr_sigma)
    mean_clr_cor <- mean_clr_cor[1:(nrow(mean_clr_cor)-1),1:(nrow(mean_clr_cor)-1)]
    upper_only <- mean_clr_cor[upper.tri(mean_clr_cor)]
    if(is.null(Sigmas)) {
      Sigmas <- matrix(NA, length(fits), length(upper_only))
    }
    Sigmas[i,] <- upper_only
  }
  Sigmas
}

parameterizations <- data.frame(days = c(30, 90, 90, 90),
                                diet = c(0, 0, 0, 50),
                                scale = c(1, 1, 2, 1),
                                plot = c(F, F, F, T))

# Canonical parameterization
folder <- paste0("asv_days90_diet25_scale1")
path <- file.path("output", "model_fits", folder, "full_posterior")
fits <- list.files(path)
canonical_Sigmas <- alt_sigmas(path, fits)

r_squared_vals <- c()
for(p in 1:length(parameterizations)) {
  folder <- paste0("asv_days", parameterizations$days[p], "_diet", parameterizations$diet[p], "_scale", parameterizations$scale[p])
  path <- file.path("output", "model_fits", folder, "full_posterior")
  fits <- list.files(path)
  Sigmas <- alt_sigmas(path, fits)
  plot_df <- data.frame(x = c(canonical_Sigmas), y = c(Sigmas))
  if(parameterizations$plot[p]) {
    plot <- ggplot(plot_df, aes(x = x, y = y)) +
      geom_point(size = 2, shape = 21, fill = "#bbbbbb") +
      theme_bw() +
      labs(x = "canonical model CLR correlation estimates",
           y = "alternative model CLR correlation estimates")
    ggsave(file.path("output", "figures", "SX.png"),
           plot,
           units = "in",
           dpi = 100,
           height = 8,
           width = 8)
  }
  cat(paste0("R^2 ", folder, ": ", round(cor(plot_df$x, plot_df$y)**2, 8), "\n"))
}
