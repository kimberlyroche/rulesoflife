library(optparse)
library(rulesoflife)

option_list = list(
  make_option(c("--sname"), type = "character", default = NULL,
              help = "short name of host series to fit", metavar = "character"),
  make_option(c("--method"), type = "character", default = "GP",
              help = "model: GP or DLM", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$sname)) {
  stop("No short name provided!")
}

# Load
data <- load_data(tax_level = "family")

if(opt$method == "DLM") {
  fit <- fit_DLM(sname = opt$sname,
                 counts = data$counts,
                 metadata = data$metadata,
                 point_est = FALSE,
                 smoothed = TRUE)
} else if(opt$method == "GP") {
  fit <- fit_GP(sname = opt$sname,
                counts = data$counts,
                metadata = data$metadata,
                point_est = FALSE)
}
