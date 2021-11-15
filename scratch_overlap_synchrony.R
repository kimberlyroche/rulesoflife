# How much overlap exists between host sampling at the weekly level?
# Monthy level?

library(rulesoflife)
library(tidyverse)
library(foreach)
library(fido)

data <- load_data(tax_level = "ASV")

# Get sampling dates for all hosts
md <- data$metadata
hosts <- sort(unique(md$sname))
host_dates <- list()
for(host in hosts) {
  host_dates[[host]] <- md %>%
    filter(sname == host) %>%
    select(collection_date, sample_id)
}

pairs <- combn(1:length(hosts), 2)

overlap <- "daily"
save_fn <- file.path("output", paste0(overlap, "_host-host_overlap.rds"))

# ------------------------------------------------------------------------------
#   Calculate (or load) all overlaps of a given frequency
# ------------------------------------------------------------------------------

if(!file.exists(save_fn)) {
  overlap_obj <- data.frame(host1 = c(), host2 = c(),
                            overlap_date1 = c(), overlap_date2 = c(),
                            sample_id1 = c(), sample_id2 = c())
  start <- Sys.time()
  for(i in 1:ncol(pairs)) {
    h1 <- hosts[pairs[1,i]]
    h2 <- hosts[pairs[2,i]]

    cat(paste0("Host pair: ", h1, " x ", h2, " (", i, " / ", ncol(pairs), ")\n"))

    d1 <- host_dates[[h1]]$collection_date
    d2 <- host_dates[[h2]]$collection_date
    start <- min(c(d1, d2))
    if(d1[1] == start) {
      # Start with d1
      ref_h <- h1
      alt_h <- h2
    } else {
      # Start with d2
      ref_h <- h2
      alt_h <- h1
    }

    ref_d <- host_dates[[ref_h]]$collection_date
    alt_d <- host_dates[[alt_h]]$collection_date
    ref_s <- host_dates[[ref_h]]$sample_id
    alt_s <- host_dates[[alt_h]]$sample_id

    for(j in 1:length(ref_d)) {
      delta <- abs(sapply(alt_d, function(x) difftime(x, ref_d[j], units = "days")))
      if(overlap == "daily") {
        hits <- unname(which(delta <= 1))
      }
      if(overlap == "weekly") {
        hits <- unname(which(delta <= 7))
      }
      if(overlap == "monthly") {
        hits <- unname(which(delta <= 30))
      }
      for(hit in hits) {
        date <- names(date)
        overlap_obj <- rbind(overlap_obj,
                             data.frame(host1 = h1, host2 = h2,
                                        overlap_date1 = ref_d[j], overlap_date2 = alt_d[hit],
                                        sample_id1 = ref_s[j], sample_id2 = alt_s[hit]))
      }
    }
  }
  rt <- Sys.time() - start
  cat("Full run time:", rt, attr(rt, "units"), "\n")

  saveRDS(overlap_obj, paste0(overlap, "_host-host_overlap.rds"))
} else {
  overlap_obj <- readRDS(save_fn)
}

# ------------------------------------------------------------------------------
#   Plot the alignment of these overlaps; how clustered are they?
# ------------------------------------------------------------------------------

all_data <- md %>%
  select(collection_date, sname)

# For a random pair
h <- sample(hosts, size = 2, replace = FALSE)
h <- c("OPH", "ORI")
sampled <- overlap_obj %>%
  filter(host1 == h[1], host2 == h[2])

ggplot() +
  geom_point(data = all_data,
             mapping = aes(x = collection_date, y = sname),
             size = 2,
             shape = 21,
             fill = "#dddddd") +
  geom_point(data = sampled,
             mapping = aes(x = overlap_date1, y = host1),
             size = 3,
             shape = 21,
             fill = "#ff0000") +
  geom_point(data = sampled,
             mapping = aes(x = overlap_date2, y = host2),
             size = 3,
             shape = 21,
             fill = "#0000ff") +
  theme(axis.text.x = element_blank()) +
  labs(x = "sample date",
       y = "host")

# ------------------------------------------------------------------------------
#   For each host-pair (1540), pull an aligned date. For a given taxon, sample
#   a 1540 length vector of these the estimated Etas pairs at this date for
#   these hosts. Calculate the correlation of these to get an idea of
#   "synchrony."
# ------------------------------------------------------------------------------

# Pull host Etas
Etas <- list()
for(h in 1:length(hosts)) {
  host <- hosts[h]
  cat(paste0("Loading Etas for host: ", host, "\n"))
  fit <- readRDS(file.path("output",
                           "model_fits",
                           "asv_days90_diet25_scale1",
                           "full_posterior",
                           paste0(hosts[h], ".rds")))
  fit.clr <- to_clr(fit)
  Eta <- apply(fit.clr$Eta, c(1,2), mean)
  Etas[[host]] <- Eta
}

n_tax <- nrow(Etas[[1]])

sampled_overlap <- NULL
for(i in 1:ncol(pairs)) {
  h1 <- hosts[pairs[1,i]]
  h2 <- hosts[pairs[2,i]]

  cat(paste0("Host pair: ", h1, " x ", h2, " (", i, " / ", ncol(pairs), ")\n"))

  temp <- overlap_obj %>%
    filter((host1 == h1 & host2 == h2) | (host1 == h2 & host2 == h1))
  if(nrow(temp) > 0) {
    temp <- temp %>%
      arrange(sample(1:nrow(temp))) %>%
      slice(1)
    if(is.null(sampled_overlap)) {
      sampled_overlap <- temp
    } else {
      sampled_overlap <- rbind(sampled_overlap, temp)
    }
  }
}

# Scrambled version!!!
for(i in 1:nrow(sampled_overlap)) {
  h1 <- sampled_overlap[i,]$host1
  h2 <- sampled_overlap[i,]$host2
  s1 <- host_dates[[h1]]
  s1 <- s1[sample(1:nrow(s1), size = 1),]
  s2 <- host_dates[[h2]]
  s2 <- s2[sample(1:nrow(s2), size = 1),]
  sampled_overlap[i,]$overlap_date1 <- s1$collection_date
  sampled_overlap[i,]$sample_id1 <- s1$sample_id
  sampled_overlap[i,]$overlap_date2 <- s2$collection_date
  sampled_overlap[i,]$sample_id2 <- s2$sample_id
}

head(sampled_overlap)

# Looks like the host #1, host #2 assignment is wrong; fix that!
# For now, we can ballpark it, because we know each host will have a sample on
#   approximately the same day.

starts <- seq(from = 1, to = nrow(sampled_overlap), by = 100)
sampled_list <- foreach(k = 1:length(starts), .combine = rbind) %dopar% {
  start <- starts[k]
  end <- min(c(nrow(sampled_overlap), start + 99))
  subset_data <- sampled_overlap[start:end,]
  row_combos <- nrow(subset_data)*(n_tax-1)*2
  sampled_Etas <- data.frame(Eta = numeric(row_combos),
                             tax_idx = numeric(row_combos),
                             partner = numeric(row_combos))
  row_counter <- 1
  for(i in 1:nrow(subset_data)) {
    # cat(paste0("Processing row ", i, " / ", nrow(subset_data), "\n"))
    s_row <- subset_data[i,]
    h1 <- s_row$host1
    h2 <- s_row$host2
    Eta1 <- Etas[[h1]]
    Eta2 <- Etas[[h2]]

    s1_idx <- which(host_dates[[h1]]$collection_date %in% c(s_row$overlap_date1,
                                                            s_row$overlap_date2))
    s2_idx <- which(host_dates[[h2]]$collection_date %in% c(s_row$overlap_date1,
                                                            s_row$overlap_date2))
    s1_idx <- s1_idx[1]
    s2_idx <- s2_idx[1]

    for(j in 1:(n_tax-1)) {
      sampled_Etas[row_counter,] <- data.frame(Eta = Eta1[j,s1_idx],
                                               tax_idx = j,
                                               partner = 1)
      row_counter <- row_counter + 1
      sampled_Etas[row_counter,] <- data.frame(Eta = Eta2[j,s2_idx],
                                               tax_idx = j,
                                               partner = 2)
      row_counter <- row_counter + 1
    }
  }
  sampled_Etas
}

# ------------------------------------------------------------------------------
#   Plot: Mean correlations for each taxon (averaged over host pairs)
# ------------------------------------------------------------------------------

correlations <- c()
for(i in 1:(n_tax-1)) {
  x <- sampled_list %>%
    filter(partner == 1 & tax_idx == i) %>%
    pull(Eta)
  y <- sampled_list %>%
    filter(partner == 2 & tax_idx == i) %>%
    pull(Eta)
  correlations <- c(correlations, cor(x,y))
}

# plot_df <- data.frame(tax_idx = 1:(n_tax-1), mean_cor = correlations)
# ggplot(plot_df, aes(x = mean_cor)) +
#   geom_histogram(color = "white") +
#   theme_bw() +
#   labs(x = "CLR(taxon) correlation across hosts")

# ------------------------------------------------------------------------------
#   Plot: Mean correlations for each taxon pair x universality scores
# ------------------------------------------------------------------------------

rug_obj <- readRDS(file.path("output", "rug_asv.rds"))
scores <- apply(rug_obj$rug, 2, calc_universality_score)

plot_df <- data.frame(synchrony = c(),
                      universality = c())
for(i in 1:length(rug_obj$tax_idx1)) {
  t1 <- rug_obj$tax_idx1[i]
  t2 <- rug_obj$tax_idx2[i]
  plot_df <- rbind(plot_df,
                   data.frame(synchrony = mean(c(correlations[t1], correlations[t2])),
                              universality = scores[i]))
}
ggplot(plot_df, aes(x = synchrony, y = universality, fill = universality)) +
  geom_point(size = 2, shape = 21) +
  scale_fill_gradientn(colours = rainbow(4)) +
  xlim(c(min(correlations), 0.46)) +
  theme_bw() +
  theme(legend.position = "none")

cat(paste0("R^2: ", round(cor(plot_df$synchrony, plot_df$universality)^2, 3), "\n"))

