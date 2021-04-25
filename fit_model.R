if(R.version$major == 4 & Sys.info()[[1]] != "Windows") {
  cat("Updating lib paths...\n")
  .libPaths(c("/gpfs/fs1/data/mukherjeelab/roche/Rlibs", .libPaths()[2]))
}

library(optparse)
library(rulesoflife)

option_list = list(
  make_option(c("--sname"), type = "character", default = NULL,
              help = "short name of host series to fit", metavar = "character"),
  make_option(c("--tax_level"), type = "character", default = NULL,
              help = "taxonomic level: phylum, family, ASV", metavar = "character"),
  make_option(c("--MAP"), type = "logical", default = FALSE,
              help = "MAP estimation flag", metavar = "logical"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "output subdirectory of model fit directory", metavar = "character"),
  make_option(c("--days_min_cor"), type = "numeric", default = 0,
              help = "days to minimum autocorrelation", metavar = "numeric"),
  make_option(c("--diet_weight"), type = "numeric", default = 0,
              help = "proportion of variance contributed to sample kernel by diet PCs", metavar = "numeric"),
  make_option(c("--var_scale_taxa"), type = "numeric", default = 1,
              help = "scale associated with taxonomic covariance matrix", metavar = "numeric"),
  make_option(c("--var_scale_samples"), type = "numeric", default = 1,
              help = "scale associated with sample covariance matrix", metavar = "numeric"),
  make_option(c("--scramble"), type = "logical", default = FALSE,
              help = "use permuted data", metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$sname)) {
  stop("No short name provided!")
}

if(!(opt$tax_level %in% c("phylum", "family", "ASV"))) {
  stop("Invalid taxonomic level!")
}

# Load
if(scramble) {
  data <- load_scrambled_data(tax_level = opt$tax_level)
} else {
  data <- load_data(tax_level = opt$tax_level)
}

MAP <- opt$MAP
fit <- fit_GP(sname = opt$sname,
              counts = data$counts,
              metadata = data$metadata,
              output_dir = opt$output_dir,
              MAP = MAP,
              diet_weight = opt$diet_weight,
              days_to_min_autocorrelation = opt$days_min_cor)
