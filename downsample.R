if(R.version$major == 4 & Sys.info()[[1]] != "Windows") {
  cat("Updating lib paths...\n")
  .libPaths(c("/gpfs/fs1/data/mukherjeelab/roche/Rlibs", .libPaths()[2]))
}

library(rulesoflife)
library(uuid)

results <- downsample_counts()
for(iter in 1:4) {
  results <- rbind(results, downsample_counts())
}

saveRDS(results, paste0("output/simdata_", UUIDgenerate(), ".rds"))
