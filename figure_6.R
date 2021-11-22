source("path_fix.R")

library(driver)
library(tidyverse)
library(rulesoflife)
library(cowplot)
library(gridGraphics)
library(doParallel)
library(foreach)
library(fido)

registerDoParallel(6)

# ------------------------------------------------------------------------------
#
#   Figure 6 - paired synchrony vs. universality plots, plus enrichment barplots
#              for relative abundances of family-family pairs in subregions
#
# ------------------------------------------------------------------------------

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

# Get all pairs of hosts
pairs <- combn(1:length(hosts), 2)

overlap <- "daily"

# ------------------------------------------------------------------------------
#   Calculate (or load) all overlaps of a given frequency
# ------------------------------------------------------------------------------

save_fn <- file.path("output", paste0(overlap, "_host-host_overlap.rds"))
if(!file.exists(save_fn)) {
  # This takes almost 90 min. to run from scratch
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
        overlap_obj <- rbind(overlap_obj,
                             data.frame(host1 = ref_h, host2 = alt_h,
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
#   For each host-pair (1540), pull an aligned date. For a given taxon, sample
#   a 1540 length vector of these the estimated Etas pairs at this date for
#   these hosts. Calculate the correlation of these to get an idea of
#   "synchrony."
#
#   This uses ~27-28 samples per host on average.
# ------------------------------------------------------------------------------

# Pull host Etas
# This takes < 30 sec.
Etas <- list()
for(h in 1:length(hosts)) {
  host <- hosts[h]
  fit <- readRDS(file.path("output",
                           "model_fits",
                           "asv_days90_diet25_scale1",
                           "MAP",
                           paste0(hosts[h], ".rds")))
  fit.clr <- to_clr(fit)
  Etas[[host]] <- fit.clr$Eta[,,1]
}

n_tax <- nrow(Etas[[1]])

# Build data.frame of overlaps
# This takes < 30 sec.
sampled_overlap <- NULL
for(i in 1:ncol(pairs)) {
  h1 <- hosts[pairs[1,i]]
  h2 <- hosts[pairs[2,i]]

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

if(FALSE) {
  # Scrambled/permuted version
  # Haven't run this in a long time; need to check
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
}

# Parallelized over 6 cores this takes < 1 min.
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
#   Plot: Mean correlations for each taxon pair x universality scores
# ------------------------------------------------------------------------------

rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
scores <- apply(rug_asv$rug, 2, calc_universality_score)
consensus_signs <- apply(rug_asv$rug, 2, calc_consensus_sign)

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

# Print out the most synchronous pairs
cat(paste0("Most synchronous pairs:\n"))
# top5 <- order(correlations, decreasing = TRUE)[1:5]
top_synchronous <- which(correlations > 0.4)
for(synchronous_idx in top_synchronous) {
  x <- data$taxonomy[synchronous_idx,2:7]
  max_level <- max(which(!is.na(x)))
  cat(paste0("\tASV ", synchronous_idx, ", in ", colnames(data$taxonomy)[max_level], " ", unname(x[max_level]), "\n"))
}
# Print out the least synchronous pairs
cat(paste0("Least synchronous pairs:\n"))
# top5 <- order(correlations, decreasing = FALSE)[1:5]
bottom_synchronous <- which(correlations < 0.05)
for(synchronous_idx in bottom_synchronous) {
  x <- data$taxonomy[synchronous_idx,2:7]
  max_level <- max(which(!is.na(x)))
  cat(paste0("\tASV ", synchronous_idx, ", in ", colnames(data$taxonomy)[max_level], " ", unname(x[max_level]), "\n"))
}

plot_df <- data.frame(synchrony = c(),
                      universality = c())
for(i in 1:length(rug_asv$tax_idx1)) {
  t1 <- rug_asv$tax_idx1[i]
  t2 <- rug_asv$tax_idx2[i]
  plot_df <- rbind(plot_df,
                   data.frame(synchrony = mean(c(correlations[t1], correlations[t2])),
                              universality = scores[i],
                              sign = consensus_signs[i]))
}
plot_df$sign <- factor(plot_df$sign, levels = c(1, -1))
levels(plot_df$sign) <- c("positive", "negative")

p <- ggplot(plot_df %>% filter(!is.na(sign)), aes(x = synchrony, y = universality, fill = sign)) + #, color = sign)) +
  geom_point(size = 2, shape = 21) +
  scale_fill_manual(values = c("#F25250", "#34CCDE")) +
  theme_bw() +
  xlim(c(min(correlations), 0.46)) +
  # xlim(c(min(correlations), 0.1)) + # for scrambled version
  theme_bw() +
  labs(fill = "Consensus\ncorrelation sign",
       x = "synchrony score",
       y = "universality score")
show(p)
ggsave(file.path(plot_dir, "F6.svg"),
       p,
       units = "in",
       dpi = 100,
       height = 5,
       width = 8)

cat(paste0("R^2: ", round(cor(plot_df$synchrony, plot_df$universality)^2, 3), "\n"))

# ------------------------------------------------------------------------------
#   Enrichment of top center and top right-hand parts
# ------------------------------------------------------------------------------

topcenter_pairs <- which(plot_df$synchrony < 0.3 & plot_df$universality > 0.5)
topright_pairs <- which(plot_df$synchrony > 0.3 & plot_df$universality > 0.4)

all_pairs <- data.frame(idx1 = rug_asv$tax_idx1,
                        idx2 = rug_asv$tax_idx2,
                        topcenter = FALSE,
                        topright = FALSE)
all_pairs$tax1 <- sapply(1:nrow(all_pairs), function(x) {
  paste0(data$taxonomy[all_pairs$idx1[x],6], collapse = "/")
})
all_pairs$tax2 <- sapply(1:nrow(all_pairs), function(x) {
  paste0(data$taxonomy[all_pairs$idx2[x],6], collapse = "/")
})

for(i in 1:length(topcenter_pairs)) {
  all_pairs$topcenter[all_pairs$idx1 == rug_asv$tax_idx1[topcenter_pairs[i]] &
                        all_pairs$idx2 == rug_asv$tax_idx2[topcenter_pairs[i]]] <- TRUE
}
for(i in 1:length(topright_pairs)) {
  all_pairs$topright[all_pairs$idx1 == rug_asv$tax_idx1[topright_pairs[i]] &
                       all_pairs$idx2 == rug_asv$tax_idx2[topright_pairs[i]]] <- TRUE
}

# 6105 of 8911 family-family pairs don't involve "NA" families
all_pairs_noNA <- all_pairs %>%
  filter(tax1 != "NA" & tax2 != "NA")

all_pairs_noNA <- all_pairs_noNA %>%
  mutate(taxpair = paste0(tax1, " - ", tax2))

frequencies <- table(all_pairs_noNA$taxpair)

# ------------------------------------------------------------------------------
#   Visualize differences in relative abundance of family pairs for top center
#   clump of pairs
# ------------------------------------------------------------------------------

# Stacked bar: plot observed relative representation of family-family pairs
frequencies_subset <- table(all_pairs_noNA$taxpair[all_pairs_noNA$topcenter == TRUE])

plot_enrichment(frequencies_subset,
                frequencies,
                plot_height = 6,
                plot_width = 8,
                legend_topmargin = 20,
                legend_leftmargin = -0.5,
                save_name = "F6-enrichment1.svg")

# ------------------------------------------------------------------------------
#   Visualize differences in relative abundance of family pairs for top center
#   clump of pairs
# ------------------------------------------------------------------------------

# Stacked bar: plot observed relative representation of family-family pairs
frequencies_subset <- table(all_pairs_noNA$taxpair[all_pairs_noNA$topright == TRUE])

plot_enrichment(frequencies_subset,
                            frequencies,
                            plot_height = 6,
                            plot_width = 9,
                            legend_topmargin = 60,
                            legend_leftmargin = 0.5,
                            save_name = "F6-enrichment2.svg")

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S7 - synchrony of most synchronous taxon across a sample
#                            of hosts
#
# ------------------------------------------------------------------------------

tax_idx <- 23
hosts <- c("DUI", "DUX", "LIW", "PEB", "VET")
host_y_offset <- 10

data <- load_data(tax_level = "ASV")

# Find common baseline
pred_objs <- NULL
min_date <- NULL
max_date <- NULL
for(host in hosts) {
  pred_obj <- get_predictions_host_list(host_list = host,
                                        output_dir = "asv_days90_diet25_scale1",
                                        metadata = data$metadata)
  if(is.null(min_date)) {
    min_date <- min(pred_obj$dates[[host]])
    max_date <- max(pred_obj$dates[[host]])
  } else {
    min_date <- min(min_date, min(pred_obj$dates[[host]]))
    max_date <- max(max_date, max(pred_obj$dates[[host]]))
  }
  pred_objs[[host]] <- pred_obj
}

plot_df <- NULL
y_offset_counter <- 0
for(host in hosts[c(2,4,3,5,1)]) {
  pred_obj <- pred_objs[[host]]
  pred_df <- gather_array(pred_obj$predictions[[host]]$Eta,
                          val,
                          coord,
                          sample,
                          iteration) %>%
    filter(coord == tax_idx) %>%
    group_by(coord, sample) %>%
    summarize(p25 = quantile(val, probs = c(0.25)),
              mean = mean(val),
              p75 = quantile(val, probs = c(0.75)))
  pred_df$host <- host

  # Calculate day from (shared) baseline
  addend <- as.numeric(difftime(min(pred_obj$dates[[host]]), min_date, units = "day"))

  day_span <- pred_obj$predictions[[host]]$span
  pred_df <- pred_df %>%
    left_join(data.frame(sample = 1:length(day_span), day = day_span), by = "sample")

  pred_df$day <- pred_df$day + addend

  if(y_offset_counter > 0) {
    pred_df$p25 <- pred_df$p25 + host_y_offset*y_offset_counter
    pred_df$mean <- pred_df$mean + host_y_offset*y_offset_counter
    pred_df$p75 <- pred_df$p75 + host_y_offset*y_offset_counter
  }
  y_offset_counter <- y_offset_counter + 1

  plot_df <- rbind(plot_df, pred_df)
}

plot(plot_df$day)

alpha <- 0.6

p <- ggplot()
for(this_host in hosts) {
  p <- p +
    geom_ribbon(data = plot_df %>% filter(host == this_host),
                mapping = aes(x = day, ymin = p25, ymax = p75),
                fill = "#fdbf6f",
                alpha = alpha) +
    geom_line(data = plot_df %>% filter(host == this_host),
              mapping = aes(x = day, y = mean),
              color = "#ff7f00",
              size = 1,
              alpha = alpha)
}
p <- p +
  theme_bw() +
  labs(x = paste0("days from first sample (", min_date, ")"),
       y = "CLR abundance") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path("output", "figures", "F6_most-synchronous.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 5,
       width = 7)

cat("Order of host series (from bottom to top):\n")
hosts[c(2,4,3,5,1)]
