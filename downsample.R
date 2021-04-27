source("path_fix.R")

library(rulesoflife)
library(uuid)

results <- downsample_counts()
for(iter in 1:4) {
  results <- rbind(results, downsample_counts())
}

saveRDS(results, paste0("output/simdata_", UUIDgenerate(), ".rds"))
