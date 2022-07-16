source("path_fix.R")

library(tidyverse)
library(codaDE)

# ------------------------------------------------------------------------------
#   ASV pairs sharing same family (in top 2.5%)
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")
tax <- data$taxonomy
md <- data$metadata
data <- data$data

tax_families <- tax$family[1:134]
tax_families[is.na(tax_families)] <- "unknown"
tax_families <- data.frame(index = 1:134, family = tax_families)

combos <- combn(1:134, m = 2)
temp <- data.frame(idx1 = combos[1,],
                   idx2 = combos[2,])
temp <- temp %>%
  left_join(tax_families, by = c("idx1" = "index")) %>%
  left_join(tax_families, by = c("idx2" = "index"))

rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
rug_scores <- data.frame(idx1 = rug_obj$tax_idx1,
                         idx2 = rug_obj$tax_idx2,
                         score = apply(rug_obj$rug, 2, calc_universality_score))

temp <- temp %>%
  left_join(rug_scores, by = c("idx1", "idx2")) %>%
  arrange(desc(score))

top_n <- round(nrow(temp)*0.025)
x <- sum(temp$family.x[1:top_n] == temp$family.y[1:top_n] & (temp$family.x[1:top_n] != "unknown"))
p <- sum(temp$family.x == temp$family.y & (temp$family.x != "unknown"))/nrow(temp)
binom.test(x, top_n, p)

# Results:
#   Null % same-family pairs 13.7%
#   Observed % same-family pairs 35.9% (80 / 223)
#   Reject the null, p-value = 0

# ------------------------------------------------------------------------------
#   ...and Lachno-Lachno in particular
# ------------------------------------------------------------------------------

x <- sum(temp$family.x[1:top_n] == "Lachnospiraceae" & temp$family.y[1:top_n] == "Lachnospiraceae")
p <- sum(temp$family.x == "Lachnospiraceae" & temp$family.y == "Lachnospiraceae")/nrow(temp)
binom.test(x, top_n, p)

# Results:
#   Null % LL-family pairs 5.6%
#   Observed % same-family pairs 20.2% (45 / 223)
#   Reject the null, p-value = 0

# ------------------------------------------------------------------------------
#   Median gap in sampling within a host
# ------------------------------------------------------------------------------

hosts <- sort(unique(md$sname))
all_gaps <- c()
for(host in hosts) {
  cat(paste0("Host: ", host, "\n"))
  samples <- md %>%
    filter(sname == host) %>%
    arrange(collection_date) %>%
    pull(collection_date)
  gaps <- numeric(length(samples)-1)
  for(i in 2:length(samples)) {
    gaps[i] <- difftime(samples[i], samples[i-1], units = "days")
  }
  all_gaps <- c(all_gaps, gaps)
}
median(all_gaps)
min(all_gaps)
max(all_gaps)

# ------------------------------------------------------------------------------
#   Median universality
# ------------------------------------------------------------------------------

# Assume ASV loaded
median(rug_obj$rug)

rug_obj <- summarize_Sigmas(output_dir = "fam_days90_diet25_scale1")
median(rug_obj$rug)

rug_obj <- summarize_Sigmas(output_dir = "phy_days90_diet25_scale1")
median(rug_obj$rug)
