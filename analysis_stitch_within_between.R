source("path_fix.R")

library(tidyverse)

files_to_stitch <- list.files(path = "output", pattern = "within_between_distros_.*?\\.rds", full.names = TRUE)
df <- readRDS(files_to_stitch[1])
for(f in files_to_stitch[2:length(files_to_stitch)]) {
	df <- bind_rows(df, readRDS(f))
}
saveRDS(df, file.path("output", "within_between_distros.rds"))

cat(paste0("Found ", length(unique(df$pair)), " pairs\n"))
