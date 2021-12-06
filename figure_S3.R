source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S3 - scrambled/permuted "rugs"
#
# ------------------------------------------------------------------------------

D_combos <- NULL
i <- 1 # permuted sample to use
pdir <- paste0("asv_days90_diet25_scale1_scramble-sample-", sprintf("%02d", i))
fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
for(j in 1:length(fits)) {
  Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
  if(is.null(D_combos)) {
    D <- nrow(Sigma)
    D_combos <- (D^2 - D)/2
    rug <- matrix(NA, length(fits), D_combos)
  }
  rug[j,] <- Sigma[upper.tri(Sigma)]
}

rug <- cbind(1:nrow(rug), rug)
colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
rug$pair <- as.numeric(rug$pair)

p <- ggplot(rug, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0,
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(fill = "Correlation",
       x = "ASV-ASV pairs",
       y = "hosts") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14, face = "plain"),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.title = element_text(size = 14),
        # legend.position = "bottom",
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())

ggsave(file.path("output", "figures", "S3_ASV.png"),
       p,
       dpi = 100,
       units = "in",
       height = 6,
       width = 8)
