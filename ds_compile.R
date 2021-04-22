if(R.version$major == 4) {
  cat("Updating lib paths...\n")
  .libPaths(c("/gpfs/fs1/data/mukherjeelab/roche/Rlibs", .libPaths()[2]))
}

library(stringr)

# Generalize to take a `p` argument
files <- list.files(path = "output", pattern = "simdata_.*")
results_all <- NULL
for(file in files) {
  results <- readRDS(file.path("output", file))
  if(is.null(results_all)) {
    results_all <- results
  } else {
    results_all <- rbind(results_all, results)
  }
}
saveRDS(results_all, file = file.path("output", "simdata_results.rds"))
