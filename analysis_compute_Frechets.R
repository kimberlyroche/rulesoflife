source("path_fix.R")

library(rulesoflife)
library(driver)
library(fido)
library(shapes)
library(matrixsampling)

# 1: Compute Frechet mean of true (MAP) dynamics estimates
# 2: Permute true (MAP) dynamics estimates; compute Frechet mean
# 3: Draw random "dynamics"; compute Frechet mean
opt <- 3

output_dir <- "asv_days90_diet25_scale1"
output_dir_full <- check_dir(c("output", "model_fits", output_dir, "MAP"))

# Load a few posterior samples for a few hosts
data <- load_data(tax_level = "ASV")
md <- data$metadata
hosts <- sort(unique(md$sname))
N <- length(hosts)

# Pull feature number
fit <- readRDS(file.path(output_dir_full, paste0(hosts[1], ".rds")))
D <- fit$D

if(opt == 3) {
  Sigmas <- rinvwishart(N, D+2, diag(D))
  for(i in 1:N) {
    Sigmas[,,i] <- cov2cor(Sigmas[,,i])
  }
} else {
  Sigmas <- array(NA, dim = c(D, D, N))
  for(i in 1:length(hosts)) {
    host <- hosts[i]
    cat(paste0("Parsing host ", host, "...\n"))
    fit <- readRDS(file.path(output_dir_full, paste0(host, ".rds")))
    fit <- to_clr(fit)
    Sigma <- fit$Sigma[,,1] + diag(nrow(fit$Sigma[,,1]))*1e-06
    Sigma <- cov2cor(Sigma)
    if(opt == 2) {
      # Introduce permutations
      idx <- sample(1:D)
      Sigma <- Sigma[idx,idx]
    }
    Sigmas[,,i] <- Sigma
  }
}

cat(paste0("Estimating Frechet mean...\n"))
start <- Sys.time()
Frechet_mean <- estcov(Sigmas, method = "Riemannian", weights = 1)
run_diff <- Sys.time() - start
cat(paste0("\tTime elapsed: ", round(run_diff, 2), " ", attr(run_diff, "units"), "\n"))
saveRDS(list(Sigmas = Sigmas,
             mean = Frechet_mean$mean),
        file.path("output", paste0("Frechet_", opt, ".rds")))

