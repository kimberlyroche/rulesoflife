library(rulesoflife)
library(fido)
library(ggplot2)
library(dplyr)

output_dir_full <- check_dir(c("output", "model_fits", "asv_days90_diet25_scale1", "MAP"))

hosts <- c("DUI", "DUX", "LIW", "PEB", "VET")

# Get taxa number and posterior sample number
for(host in hosts) {
  fit <- readRDS(file.path(output_dir_full, paste0(host, ".rds")))
  fit.clr <- to_clr(fit)
  p <- plot_kernel_or_cov_matrix(cov2cor(fit.clr$Sigma[,,1]))
  p <- p +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  ggsave(paste0("corrmat_", host, ".png"),
         p,
         units = "in",
         dpi = 100,
         height = 6,
         width = 6)
}

data <- load_data("ASV")
data$metadata %>%
  filter(sname %in% hosts) %>%
  select(sname, sex) %>%
  distinct()
