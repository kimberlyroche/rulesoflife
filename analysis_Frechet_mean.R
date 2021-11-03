source("path_fix.R")

library(rulesoflife)
library(fido)
library(dplyr)
library(shapes)
library(optparse)

option_list = list(
  make_option(c("--corr"), type = "logical", default = FALSE,
              help = "convert covariance to correlation", metavar = "logical"),
  make_option(c("--clr"), type = "logical", default = FALSE,
              help = "convert to CLR", metavar = "logical"),
  make_option(c("--group"), type = "numeric", default = 0,
              help = "social group for which to calculate mean", metavar = "number"),
  make_option(c("--sex"), type = "character", default = NULL,
              help = "sex for which to calculate mean", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

use_corr <- opt$corr
use_clr <- opt$clr
use_group <- opt$group
use_sex <- opt$sex

if(!(use_group %in% c(0, 1.1, 1.21, 1.22, 2.1, 2.2))) {
  stop(paste0("Invalid social group: ", use_group, "!\n"))
}

if(!is.null(use_sex)) {
  if(!(use_sex %in% c("F", "M"))) {
    stop(paste0("Invalid sex: ", use_sex, "!\n"))
  }
}

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

host_list <- c()
for(i in 1:length(file_list)) {
  fn <- file_list[i]
  host <- substr(fn, 1, 3)
  host_list <- c(host_list, host)
  cat(paste0("Parsing file for host ", host, "\n"))
  fit <- readRDS(file.path(output_dir_full, fn))
  if(use_clr) {
    fit <- to_clr(fit)
  }
  Sigma <- fit$Sigma[,,1] + diag(nrow(fit$Sigma[,,1]))*1e-06
  if(use_corr) {
    Sigma <- cov2cor(Sigma)
  }
  S[,,i] <- Sigma
}

# Filter to social group if necessary
if(use_group != 0) {
  hosts_grp <- get_host_social_groups(host_list) %>%
    filter(grp == use_group) %>%
    pull(sname)
  S <- S[,,which(host_list %in% hosts_grp)]
}

# Filter to sex if necessary
if(!is.null(use_sex)) {
  metadata <- load_data()$metadata
  host_sex <- metadata %>%
    select(sname, sex) %>%
    distinct() %>%
    filter(sex == use_sex) %>%
    pull(sname)
  S <- S[,,which(host_list %in% host_sex)]
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

saveRDS(Frechet_mean,
        file = file.path("output", paste0("Frechet_corr",
                                          ifelse(use_corr, 1, 0),
                                          "_",
                                          ifelse(use_clr, "clr", "alr"),
                                          ifelse(use_group > 0, paste0("_group-", use_group), ""),
                                          ifelse(!is.null(use_sex), paste0("_sex-", use_sex), ""),
                                          ".rds")))

