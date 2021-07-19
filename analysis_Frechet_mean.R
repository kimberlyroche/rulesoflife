source("path_fix.R")

library(rulesoflife)
library(fido)
library(shapes)
library(optparse)

option_list = list(
  make_option(c("--corr"), type = "logical", default = FALSE,
              help = "convert covariance to correlation", metavar = "logical"),
  make_option(c("--clr"), type = "logical", default = FALSE,
              help = "convert to CLR", metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

use_corr <- opt$corr
use_clr <- opt$clr

# ------------------------------------------------------------------------------
#   Get MAP covariance or correlation for all hosts
# ------------------------------------------------------------------------------

output_dir <- "asv_days90_diet25_scale1"
output_dir_full <- check_dir(c("output", "model_fits", output_dir, "MAP"))
file_list <- list.files(path = output_dir_full, pattern = "*.rds")

# Get taxa number and posterior sample number
fit <- readRDS(file.path(output_dir_full, file_list[1]))
D <- fit$D
iter <- fit$iter
limit <- length(file_list)

if(use_clr) {
  S <- array(NA, dim = c(D, D, length(file_list)))
} else {
  S <- array(NA, dim = c(D-1, D-1, length(file_list)))
}

for(i in 1:length(file_list)) {
  cat(paste0("Parsing file #", i, "\n"))
  fit <- readRDS(file.path(output_dir_full, file_list[i]))
  if(use_clr) {
    fit <- to_clr(fit)
  }
  Sigma <- fit$Sigma[,,1] + diag(nrow(fit$Sigma[,,1]))*1e-06
  if(use_corr) {
    Sigma <- cov2cor(Sigma)
  }
  S[,,i] <- Sigma
}

cat("Estimating Frechet mean...\n")
start <- Sys.time()
Frechet_mean <- estcov(S, method = "Riemannian", weights = 1)
diff <- Sys.time() - start
cat(paste0("Estimate took ",
           round(as.numeric(diff), 3),
           " ",
           attr(diff, "units"),
           "\n"))

saveRDS(Frechet_mean, file = file.path("output", paste0("Frechet_corr",
                                                        ifelse(use_corr, 1, 0),
                                                        "_",
                                                        ifelse(use_clr, "clr", "alr"),
                                                        ".rds")))
