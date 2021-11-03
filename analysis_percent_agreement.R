source("path_fix.R")

library(rulesoflife)
library(tidyverse)

source("ggplot_fix.R")

map <- data.frame(level = c("phylum", "family", "ASV"),
                  short_name = c("phy", "fam", "asv"))

for(level in map$level) {
  cat("Evaluating", level, "\n")
  short_name <- map[map$level == level,]$short_name
  rug_fn <- file.path("output", paste0("rug_", short_name, ".rds"))
  if(file.exists(rug_fn)) {
    rug_obj <- readRDS(rug_fn)
  } else {
    rug_obj <- summarize_Sigmas(output_dir = paste0(short_name, "_days90_diet25_scale1"))
    saveRDS(rug_obj, rug_fn)
  }
  rug <- rug_obj$rug
  weak_idx <- which(rug < 0.25 & rug > -0.25, arr.ind = TRUE)
  rug[weak_idx] <- NA

  agree_count <- apply(rug, 2, function(x) {
    y <- sign(x[!is.na(x)])
    canonical_sign <- sign(mean(y))
    sum(y == canonical_sign)
  })

  cat(paste0("Percent agreement, thresholded correlations: ",
             round(sum(agree_count) / sum(!is.na(rug)), 3)*100, "\n"))

  cat(paste0("Percent agreement, all correlations: ",
             round(sum(agree_count) / length(rug), 3)*100, "\n"))
}

