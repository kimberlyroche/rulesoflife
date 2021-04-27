if(R.version$major == 4 & Sys.info()[[1]] != "Windows") {
  cat("Updating lib paths...\n")
  .libPaths(c("/gpfs/fs1/data/mukherjeelab/roche/Rlibs", .libPaths()[2]))
}

