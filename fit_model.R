if(R.version$major == 4 & Sys.info()[[1]] != "Windows") {
  cat("Updating lib paths...\n")
  .libPaths(c("/gpfs/fs1/data/mukherjeelab/roche/Rlibs", .libPaths()[2]))
}

library(optparse)
library(rulesoflife)

option_list = list(
  make_option(c("--sname"), type = "character", default = NULL,
              help = "short name of host series to fit", metavar = "character"),
  make_option(c("--method"), type = "character", default = "GP",
              help = "model: GP or DLM", metavar = "character"),
  make_option(c("--MAP"), type = "logical", default = FALSE,
              help = "MAP estimation flag", metavar = "logical"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "output subdirectory of model fit directory", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$sname)) {
  stop("No short name provided!")
}

# Load
data <- load_data(tax_level = "family")

point_est <- opt$MAP
if(opt$method == "DLM") {
  fit <- fit_DLM(sname = opt$sname,
                 counts = data$counts,
                 metadata = data$metadata,
                 point_est = point_est,
                 smoothed = TRUE)
} else if(opt$method == "GP") {
  fit <- fit_GP(sname = opt$sname,
                counts = data$counts,
                metadata = data$metadata,
                point_est = point_est)
}
