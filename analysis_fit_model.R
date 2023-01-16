source("path_fix.R")

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
  make_option(c("--days_min_cor"), type = "numeric", default = 90,
              help = "days to minimum autocorrelation", metavar = "numeric"),
  make_option(c("--diet_weight"), type = "numeric", default = 0.25,
              help = "proportion of variance contributed to sample kernel by diet PCs", metavar = "numeric"),
  make_option(c("--var_scale_taxa"), type = "numeric", default = 1,
              help = "scale associated with taxonomic covariance matrix", metavar = "numeric"),
  make_option(c("--var_scale_samples"), type = "numeric", default = 1,
              help = "scale associated with sample covariance matrix", metavar = "numeric"),
  make_option(c("--concentration"), type = "numeric", default = 10,
              help = "concentration parameter for IW prior over taxon-taxon covariance", metavar = "numeric"),
  make_option(c("--age_min"), type = "numeric", default = -Inf,
              help = "minimum baboon age", metavar = "numeric"),
  make_option(c("--age_max"), type = "numeric", default = Inf,
              help = "maximum baboon age", metavar = "numeric"),
  make_option(c("--scramble_sample"), type = "logical", default = FALSE,
              help = "scramble taxa within samples", metavar = "logical"),
  make_option(c("--scramble_spacing"), type = "logical", default = FALSE,
              help = "scrambled sampling frequency", metavar = "logical"),
  make_option(c("--scramble_order"), type = "logical", default = FALSE,
              help = "scrambled sampling order", metavar = "logical"),
  make_option(c("--use_adam"), type = "logical", default = FALSE,
              help = "use Adam for optimization", metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$sname)) {
  stop("No short name provided!")
}

if(!(opt$tax_level %in% c("phylum", "family", "ASV"))) {
  stop("Invalid taxonomic level!")
}

cat(paste0("Tax level is: ", opt$tax_level, "\n"))

data <- load_data(tax_level = opt$tax_level)

fit <- fit_GP(sname = opt$sname,
              counts = data$counts,
              metadata = data$metadata,
              output_dir = opt$output_dir,
              MAP = opt$MAP,
              days_to_min_autocorrelation = opt$days_min_cor,
              diet_weight = opt$diet_weight,
              concentration = opt$concentration,
              age_min = opt$age_min,
              age_max = opt$age_max,
              scramble_sample = opt$scramble_sample,
              scramble_spacing = opt$scramble_spacing,
              scramble_order = opt$scramble_order,
              use_adam = opt$use_adam)
