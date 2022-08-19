source("path_fix.R")

library(driver)
library(tidyverse)
library(rulesoflife)
library(cowplot)
library(gridGraphics)
library(doParallel)
library(foreach)
library(fido)
library(magrittr)

palette_fn <- file.path("output", "family_palette.rds")
fpalette <- readRDS(palette_fn)

# registerDoParallel(detectCores())

null_case <- TRUE # needed if rendering Figure S9

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

pull_Etas <- function(sample_obj, n_tax, Etas, host_dates) {
  # Parallelized over 6 cores this takes < 1 min.
  starts <- seq(from = 1, to = nrow(sample_obj), by = 100)
  sampled_list <- foreach(k = 1:length(starts), .combine = rbind) %dopar% {
    start <- starts[k]
    end <- min(c(nrow(sample_obj), start + 99))
    subset_data <- sample_obj[start:end,]
    row_combos <- nrow(subset_data)*(n_tax-1)*2
    sampled_Etas <- data.frame(Eta = numeric(row_combos),
                               tax_idx = numeric(row_combos),
                               partner = numeric(row_combos))
    row_counter <- 1
    for(i in 1:nrow(subset_data)) {
      s_row <- subset_data[i,]
      h1 <- s_row$host1
      h2 <- s_row$host2
      Eta1 <- Etas[[h1]]
      Eta2 <- Etas[[h2]]

      s1_idx <- which(host_dates[[h1]]$collection_date == s_row$overlap_date1)
      s2_idx <- which(host_dates[[h2]]$collection_date == s_row$overlap_date2)
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
  return(sampled_list)
}

# ------------------------------------------------------------------------------
#   Perform the actual synchrony estimates
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")

tax <- data$taxonomy

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
#   This uses ~27-28 samples per host on average. We could further thin this.
# ------------------------------------------------------------------------------

save_fn <- file.path("output", "saved_synchrony_samples.rds")
if(!file.exists(save_fn)) {
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
  permuted_overlap <- NULL
  for(it in 1:10) {
    for(i in 1:ncol(pairs)) {
      h1 <- hosts[pairs[1,i]]
      h2 <- hosts[pairs[2,i]]

      if(it == 1) {
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

      # For "null" distribution -- take a random (likely non-overlapping) pair of
      # samples from these hosts
      h1_samples <- md %>%
        filter(sname == h1) %>%
        select(sample_id, collection_date)
      h1_sample <- h1_samples[sample(1:nrow(h1_samples), size = 1),]
      h2_samples <- md %>%
        filter(sname == h2) %>%
        select(sample_id, collection_date)
      h2_sample <- h2_samples[sample(1:nrow(h2_samples), size = 1),]
      permuted_overlap <- rbind(permuted_overlap,
                                data.frame(host1 = h1,
                                           host2 = h2,
                                           overlap_date1 = h1_sample$collection_date,
                                           overlap_date2 = h2_sample$collection_date,
                                           sample_id1 = h1_sample$sample_id,
                                           sample_id2 = h2_sample$sample_id,
                                           iteration = it))
    }
  }

  sampled_list <- pull_Etas(sampled_overlap, n_tax, Etas, host_dates)
  sampled_list_permuted <- NULL
  if(null_case) {
    for(it in 1:10) {
      cat(paste0("Iteration ", it, "\n"))
      sampled_list_permuted <- rbind(sampled_list_permuted,
                                     cbind(pull_Etas(permuted_overlap %>% filter(iteration == it), n_tax, Etas, host_dates),
                                           iteration = it))
    }
  }

  saveRDS(list(observed = sampled_list,
               permuted = sampled_list_permuted), save_fn)
} else {
  sampled_obj <- readRDS(save_fn)
  sampled_list <- sampled_obj$observed
  sampled_list_permuted <- sampled_obj$permuted
  n_tax <- length(unique(sampled_list$tax_idx)) + 1
}

rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
scores <- apply(rug_asv$rug, 2, calc_universality_score)
mcs <- apply(rug_asv$rug, 2, function(x) median(abs(x)))
consensus_signs <- apply(rug_asv$rug, 2, calc_consensus_sign)

save_fn <- file.path("output", "saved_synchrony_correlations.rds")
if(!file.exists(save_fn)) {
  paired_samples_x <- c()
  paired_samples_y <- c()
  correlations <- c()
  correlations_permuted <- c()
  for(i in 1:(n_tax-1)) {
    cat(paste0("Taxon iteration ", i, "\n"))
    x <- sampled_list %>%
      filter(partner == 1 & tax_idx == i) %>%
      pull(Eta)
    y <- sampled_list %>%
      filter(partner == 2 & tax_idx == i) %>%
      pull(Eta)
    if(i == 21) {
      paired_samples_x <- c(paired_samples_x, x)
      paired_samples_y <- c(paired_samples_y, y)
    }
    correlations <- c(correlations, cor(x,y))
    if(null_case) {
      for(it in 1:10) {
        cat(paste0("\tPermutation iteration ", it, "\n"))
        x <- sampled_list_permuted %>%
          filter(iteration == it & partner == 1 & tax_idx == i) %>%
          pull(Eta)
        y <- sampled_list_permuted %>%
          filter(iteration == it & partner == 2 & tax_idx == i) %>%
          pull(Eta)
        correlations_permuted <- c(correlations_permuted, cor(x,y))
      }
    }
  }

  saveRDS(list(correlations = correlations,
               correlations_permuted = correlations_permuted),
          save_fn)
} else {
  corr_obj <- readRDS(save_fn)
  correlations <- corr_obj$correlations
  correlations_permuted <- corr_obj$correlations_permuted

  paired_samples_x <- sampled_list %>%
    filter(partner == 1 & tax_idx == 23) %>%
    pull(Eta)
  paired_samples_y <- sampled_list %>%
    filter(partner == 2 & tax_idx == 23) %>%
    pull(Eta)
}

# ------------------------------------------------------------------------------
#   Write out table of per-taxon synchrony estimates
# ------------------------------------------------------------------------------

write_df <- data.frame(synchrony = correlations,
                       tax_idx = 1:length(correlations),
                       taxonomy = sapply(1:(nrow(data$taxonomy)-1), function(x) {
                         paste0(data$taxonomy[x,2:ncol(data$taxonomy)], collapse = " / ")
                       }))
write_df <- write_df %>%
  arrange(desc(synchrony))
write_df <- cbind(rank = 1:nrow(write_df), write_df)
write.table(write_df,
            file.path("output", "synchrony_table.tsv"),
            sep = "\t",
            quote = F,
            row.names = F)

# Print out the most synchronous pairs
cat(paste0("Most synchronous pairs:\n"))
top_synchronous <- which(correlations > 0.4)
for(synchronous_idx in top_synchronous) {
  x <- data$taxonomy[synchronous_idx,2:7]
  max_level <- max(which(!is.na(x)))
  cat(paste0("\tASV ", synchronous_idx, ", in ", colnames(data$taxonomy)[max_level], " ", unname(x[max_level]), "\n"))
}

# Print out the least synchronous pairs
cat(paste0("Least synchronous pairs:\n"))
bottom_synchronous <- which(correlations < 0.05)
for(synchronous_idx in bottom_synchronous) {
  x <- data$taxonomy[synchronous_idx,2:7]
  max_level <- max(which(!is.na(x)))
  cat(paste0("\tASV ",
             synchronous_idx,
             ", in ",
             colnames(data$taxonomy)[max_level],
             " ",
             unname(x[max_level]), "\n"))
}

# ------------------------------------------------------------------------------
#
#   Figures
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#   Synchrony versus universality
# ------------------------------------------------------------------------------

save_fn <- file.path("output", "synchrony_plot_obj.rds")
if(file.exists(save_fn)) {
  plot_obj <- readRDS(save_fn)
  plot_df <- plot_obj$plot_df
  plot_df_permuted <- plot_obj$plot_df_permuted
} else {
  plot_df <- NULL
  plot_df_permuted <- NULL
  for(i in 1:length(rug_asv$tax_idx1)) {
    if(i %% 1000 == 0) {
      cat(paste0("Taxon pair iteration ", i, " / ", length(rug_asv$tax_idx1), "\n"))
    }
    t1 <- rug_asv$tax_idx1[i]
    t2 <- rug_asv$tax_idx2[i]
    plot_df <- rbind(plot_df,
                     data.frame(synchrony = mean(c(correlations[t1], correlations[t2])),
                                universality = scores[i],
                                mcs = mcs[i],
                                sign = consensus_signs[i]))
    if(null_case) {
      for(it in 1:10) {
        t1_idx <- (it-1)*(n_tax-1) + t1
        t2_idx <- (it-1)*(n_tax-1) + t2
        plot_df_permuted <- rbind(plot_df_permuted,
                                  data.frame(synchrony = mean(c(correlations_permuted[t1_idx],
                                                                correlations_permuted[t2_idx])),
                                             universality = scores[i],
                                             sign = consensus_signs[i]))
      }
    }
  }
  plot_df$sign <- factor(plot_df$sign, levels = c(1, -1))
  levels(plot_df$sign) <- c("positive", "negative")
  plot_df_permuted$sign <- factor(plot_df_permuted$sign, levels = c(1, -1))
  levels(plot_df_permuted$sign) <- c("positive", "negative")
  saveRDS(list(plot_df = plot_df, plot_df_permuted = plot_df_permuted), save_fn)
}

p0 <- ggplot() +
  geom_point(data = plot_df %>% filter(!is.na(sign)),
             mapping = aes(x = synchrony, y = universality, fill = sign),
             size = 2, shape = 21) +
  geom_segment(data = data.frame(x = 0, xend = 0.5, y = 0.4, yend = 0.4),
               mapping = aes(x = x, y = y, xend = xend, yend = yend),
               color = "black",
               linetype = 2,
               size = 0.8) +
  geom_segment(data = data.frame(x = 0.3, xend = 0.3, y = 0.4, yend = 0.8),
               mapping = aes(x = x, y = y, xend = xend, yend = yend),
               color = "black",
               linetype = 2,
               size = 0.8) +
  geom_text(data = data.frame(x = 0, y = 0.78, label = "High universality\nLow synchrony"),
            mapping = aes(x = x, y = y, label = label),
            size = 5,
            hjust = 0,
            color = "black") +
  geom_text(data = data.frame(x = 0.50, y = 0.78, label = "High universality\nHigh synchrony"),
            mapping = aes(x = x, y = y, label = label),
            size = 5,
            hjust = 1,
            color = "black") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(fill = "Consensus\ncorrelation sign",
       x = "synchrony score",
       y = "universality score")

cat(paste0("R^2 (variance explained): ", round(cor(plot_df$synchrony, plot_df$universality)^2, 3), "\n"))
fit <- lm(scale(plot_df$synchrony) ~ scale(plot_df$universality))
cat(paste0("Beta: ", round(coef(summary(fit))[2,1], 3), "\n"))
cat(paste0("\tp-value: ", round(coef(summary(fit))[2,4], 3), "\n"))

median(plot_df$synchrony)
min(plot_df$synchrony)
max(plot_df$synchrony)

if(null_case) {
  p1p <- ggplot() +
    geom_point(data = plot_df_permuted %>% filter(!is.na(sign)),
               mapping = aes(x = synchrony, y = universality, fill = sign),
               size = 2, shape = 21) +
    scale_fill_manual(values = c("#F25250", "#34CCDE")) +
    theme_bw() +
    labs(fill = "Consensus\ncorrelation sign",
         x = "synchrony score",
         y = "universality score")

  cat(paste0("R^2: ", round(cor(plot_df_permuted$synchrony, plot_df_permuted$universality)^2, 3), "\n"))

  ggsave(file.path("output", "figures", paste0("synchrony", ifelse(null_case, "_null", ""), ".svg")),
         p1p,
         dpi = 100,
         units = "in",
         height = 6,
         width = 9)
}

# ------------------------------------------------------------------------------
#   Are Johannes' seasonal taxa more synchronous or universal than others?
#   (OMITTED)
# ------------------------------------------------------------------------------

seasonal_families <- list(wet = c("Helicobacteraceae",
                                  "Coriobacteriaceae",
                                  "Burkholderiaceae",
                                  "Eggerthellaceae",
                                  "Atopobiaceae",
                                  "Erysipelotrichaceae",
                                  "Lachnospiraceae"),
                          dry = c("Baceroidales RF16 group",
                                  "vadinBE97",
                                  "Spirochaetaceae",
                                  "Campylobacteraceae",
                                  "Christensenellaceae",
                                  "Syntrophomonadaceae"))

plot_df$idx1 <- rug_asv$tax_idx1
plot_df$idx2 <- rug_asv$tax_idx2
plot_df$fam1 <- tax[plot_df$idx1,6]
plot_df$fam2 <- tax[plot_df$idx2,6]

plot_df$fam1[!(plot_df$fam1 %in% seasonal_families$wet |
                 plot_df$fam1 %in% seasonal_families$dry)] <- NA
plot_df$fam2[!(plot_df$fam2 %in% seasonal_families$wet |
                 plot_df$fam2 %in% seasonal_families$dry)] <- NA

if(FALSE) {
  x <- plot_df %>%
    filter(fam1 == "Lachnospiraceae" & fam2 == "Lachnospiraceae") %>%
    dplyr::select(c(synchrony, universality))
  cat(paste0("R^2: ", round(cor(x$synchrony, x$universality)^2, 3), "\n"))

  x <- plot_df %>%
    filter(!is.na(fam1) & fam1 != "Lachnospiraceae" &
             !is.na(fam2) & fam2 != "Lachnospiraceae") %>%
    dplyr::select(c(synchrony, universality))
  cat(paste0("R^2: ", round(cor(x$synchrony, x$universality)^2, 3), "\n"))

  p1 <- ggplot() +
    geom_point(data = plot_df %>% filter(is.na(fam1) | fam1 == "Lachnospiraceae"),
               mapping = aes(x = synchrony, y = universality),
               color = "gray") +
    geom_point(data = plot_df %>% filter(!is.na(fam1) & fam1 != "Lachnospiraceae"),
               mapping = aes(x = synchrony, y = universality, fill = fam1),
               size = 2,
               shape = 21) +
    scale_fill_manual(values = fpalette) +
    theme_bw() +
    theme(legend.position = "none")
  p2 <- ggplot() +
    geom_point(data = plot_df %>% filter(is.na(fam1)),
               mapping = aes(x = synchrony, y = universality),
               color = "gray") +
    geom_point(data = plot_df %>% filter(fam1 == "Lachnospiraceae"),
               mapping = aes(x = synchrony, y = universality, fill = fam1),
               size = 2,
               shape = 21) +
    scale_fill_manual(values = fpalette) +
    theme_bw() +
    theme(legend.position = "none")
  p3 <- ggplot() +
    geom_point(data = plot_df %>% filter(is.na(fam2) | fam2 == "Lachnospiraceae"),
               mapping = aes(x = synchrony, y = universality),
               color = "gray") +
    geom_point(data = plot_df %>% filter(!is.na(fam2) & fam2 != "Lachnospiraceae"),
               mapping = aes(x = synchrony, y = universality, fill = fam2),
               size = 2,
               shape = 21) +
    scale_fill_manual(values = fpalette) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(fill = "family")
  legend <- get_legend(p3)
  p3 <- p3 +
    theme(legend.position = "none")
  p4 <- ggplot() +
    geom_point(data = plot_df %>% filter(is.na(fam2)),
               mapping = aes(x = synchrony, y = universality),
               color = "gray") +
    geom_point(data = plot_df %>% filter(fam2 == "Lachnospiraceae"),
               mapping = aes(x = synchrony, y = universality, fill = fam2),
               size = 2,
               shape = 21) +
    scale_fill_manual(values = fpalette) +
    theme_bw() +
    theme(legend.position = "none")
  r1 <- plot_grid(p1, p2, rel_widths = c(1, 1))
  r2 <- plot_grid(p3, p4, rel_widths = c(1, 1))
  p <- plot_grid(r1, r2, legend, ncol = 1, rel_heights = c(1, 1, 0.3))
  ggsave(file.path("output", "figures", "temp.svg"),
         p,
         dpi = 100,
         units = "in",
         height = 10,
         width = 10)

  # # Purposely overplot
  # # Labeled: any single partner is seasonal, non-Lachno
  p1 <- ggplot() +
    geom_point(data = plot_df %>% filter(is.na(fam1) | fam1 == "Lachnospiraceae" |
                                           is.na(fam2) | fam2 == "Lachnospiraceae"),
               mapping = aes(x = synchrony, y = universality),
               color = "gray") +
    geom_point(data = plot_df %>% filter(!is.na(fam1) & fam1 != "Lachnospiraceae"),
               mapping = aes(x = synchrony, y = universality, fill = fam1),
               size = 2,
               shape = 21) +
    geom_point(data = plot_df %>% filter(!is.na(fam2) & fam2 != "Lachnospiraceae"),
               mapping = aes(x = synchrony, y = universality, fill = fam2),
               size = 2,
               shape = 21) +
    scale_fill_manual(values = fpalette) +
    theme_bw() +
    theme(legend.position = "none")
  # Labeled: any single partner is Lachnospiraceae
  p2 <- ggplot() +
    geom_point(data = plot_df %>% filter(is.na(fam1) | fam1 != "Lachnospiraceae" |
                                           is.na(fam2) | fam2 != "Lachnospiraceae"),
               mapping = aes(x = synchrony, y = universality),
               color = "gray") +
    geom_point(data = plot_df %>% filter(fam1 == "Lachnospiraceae" | fam2 == "Lachnospiraceae"),
               mapping = aes(x = synchrony, y = universality, fill = fam2),
               size = 2,
               shape = 21) +
    scale_fill_manual(values = fpalette) +
    theme_bw() #+
    # theme(legend.position = "none")

  p <- plot_grid(p1, p2, rel_widths = c(1, 1.35))
  ggsave(file.path("output", "figures", "temp.svg"),
         p,
         dpi = 100,
         units = "in",
         height = 4.5,
         width = 12)
}

if(FALSE) {
  # Evaluate which pairs have strong, singificant correlations (synchrony x universality)
  # These are interesting - Lachno-Lachno is at the top - but it's REALLY tough
  # to interpret these. Probably drop this.
  families <- sort(unique(tax$family)) # all families
  # families <- sort(unique(c(plot_df$fam1, plot_df$fam2))) # Johannes' seasonal only

  # Restore family labels
  plot_df$fam1_preserved <- tax[plot_df$idx1,6]
  plot_df$fam2_preserved <- tax[plot_df$idx2,6]

  temp <- combn(1:length(families), m = 2)
  fam_combos <- data.frame(a = 1:length(families), b = 1:length(families), r = NA, p = NA)
  for(i in 1:ncol(temp)) {
    fam_combos <- rbind(fam_combos,
                        data.frame(a = temp[1,i], b = temp[2,i], r = NA, p = NA))
  }
  fam_combos <- fam_combos %>%
    arrange(a, b) %>%
    left_join(data.frame(a = 1:length(families), family_a = families), by = "a") %>%
    left_join(data.frame(b = 1:length(families), family_b = families), by = "b")

  for(i in 1:nrow(fam_combos)) {
    temp <- plot_df %>%
      filter((fam1_preserved == fam_combos$family_a[i] & fam2_preserved == fam_combos$family_b[i]) |
               (fam2_preserved == fam_combos$family_a[i] & fam1_preserved == fam_combos$family_b[i])) %>%
      dplyr::select(c(synchrony, universality))
    if(nrow(temp) >= 5) {
      fit <- lm(scale(temp$universality) ~ scale(temp$synchrony))
      fam_combos$r[i] <- coef(summary(fit))[2,1]
      fam_combos$p[i] <- coef(summary(fit))[2,4]
    }
  }
  fam_combos <- fam_combos %>%
    filter(complete.cases(.))
  fam_combos$q <- p.adjust(fam_combos$p, method = "BH")

  fam_combos <- fam_combos %>%
    filter(q < 0.01) %>%
    mutate(r = round(r, 3),
           q = round(q, 3)) %>%
    dplyr::select(c(family_a, family_b, r, q)) %>%
    arrange(q, desc(r))
  fam_combos

  # Labeled: arbitrary family pairs
  all_family_pair_palette <- readRDS(file.path("output", "family-family_palette.rds"))
  fams <- c("Clostridiaceae 1", "Prevotellaceae")
  temp <- plot_df %>%
    mutate(label = case_when(
      fam1_preserved < fam2_preserved ~ paste0(fam1_preserved, " - ", fam2_preserved),
      TRUE ~ paste0(fam2_preserved, " - ", fam1_preserved)
    ))
  temp$label[temp$label != ifelse(fams[1] < fams[2],
                                  paste0(fams[1], " - ", fams[2]),
                                  paste0(fams[2], " - ", fams[1]))] <- NA
  ggplot() +
    geom_point(data = temp %>% filter(is.na(label)),
               mapping = aes(x = synchrony, y = universality),
               color = "gray") +
    geom_point(data = temp %>% filter(!is.na(label)),
               mapping = aes(x = synchrony, y = universality, fill = label),
               size = 3,
               shape = 21) +
    scale_fill_manual(values = all_family_pair_palette) +
    theme_bw() +
    theme(legend.position = "none")
}

# ------------------------------------------------------------------------------
#   Lachnospiraceae-lachnospiraceae pairs show increased universality with
#   increased synchrony
# ------------------------------------------------------------------------------

if(FALSE) {
  x <- plot_df %>% filter(fam1 == "Lachnospiraceae" & fam2 == "Lachnospiraceae")
  y <- plot_df %>% filter(fam1 != "Lachnospiraceae" & fam2 != "Lachnospiraceae")

  x$universality <- scale(x$universality)
  x$synchrony <- scale(x$synchrony)
  res <- summary(lm(universality ~ synchrony, x))
  cat(paste0("Lachno-lachno beta (synchrony x universality): ", round(coef(res)[2,1], 3), "\n"))
  cat(paste0("\tp-value: ", round(coef(res)[2,4], 3), "\n"))

  y$universality <- scale(y$universality)
  y$synchrony <- scale(y$synchrony)
  res <- summary(lm(universality ~ synchrony, y))
  cat(paste0("Non-lachno 'seasonal' pairs beta (synchrony x universality): ", round(coef(res)[2,1], 3), "\n"))
  cat(paste0("\tp-value: ", round(coef(res)[2,4], 3), "\n"))
}

# ------------------------------------------------------------------------------
#   Barplots
# ------------------------------------------------------------------------------

# D here is really D - 1
D <- n_tax - 1
temp <- data.frame(index = 1:D, fam = tax$family[1:D], synchrony = correlations)
temp$seasonal <- "none"
temp$seasonal[temp$fam %in% seasonal_families$wet] <- "wet"
temp$seasonal[temp$fam %in% seasonal_families$dry] <- "dry"
temp$label <- paste0("taxon ", temp$index, " (family ", temp$fam, ")")

season_palette <- c("gray", "#9B59B6", "#59B690")
names(season_palette) <- c("none", "dry", "wet")

p1 <- ggplot(temp %>% filter(fam == "Lachnospiraceae"),
       aes(x = synchrony, y = reorder(label, synchrony), fill = seasonal)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = season_palette) +
  theme_bw() +
  xlim(c(0,0.5)) +
  labs(y = "",
       fill = "Seasonal enrichment") +
  theme(legend.position = "bottom")
legend <- get_legend(p1)
p1 <- p1 +
  theme(legend.position = "none")
p2 <- ggplot(temp %>% filter(fam != "Lachnospiraceae"),
       aes(x = synchrony, y = reorder(label, synchrony), fill = seasonal)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = season_palette) +
  theme_bw() +
  xlim(c(0,0.5)) +
  labs(y = "",
       fill = "Seasonal enrichment") +
  theme(legend.position = "none")

row1 <- plot_grid(p1, p2,
                  ncol = 2,
                  rel_widths = c(1, 1),
                  labels = c("A", "B"),
                  label_size = 18)

p <- plot_grid(row1, legend,
               ncol = 1,
               rel_heights = c(1, 0.1))

temp$lachno <- !is.na(temp$fam) & temp$fam == "Lachnospiraceae"
p <- ggplot(temp, aes(x = reorder(label, synchrony), y = synchrony, fill = seasonal, color = lachno)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = season_palette) +
  scale_color_manual(values = c("white", "black")) +
  theme_bw() +
  ylim(c(0,0.5)) +
  labs(x = "",
       fill = "Seasonal enrichment") +
  guides(color = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
        legend.position = "bottom")

ggsave(file.path("output", "figures", "season_vs_synchrony.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 14,
       width = 7)

# Test for significant enrichment of labeled things and highly synchronous things
temp$seasonal[temp$fam %in% unlist(seasonal_families)] <- "seasonal"
temp$seasonal <- factor(temp$seasonal)
# base::summary(aov(synchrony ~ seasonal, temp))
# p = 0.358

# ------------------------------------------------------------------------------
#   Boxplots
# ------------------------------------------------------------------------------

plot_df$n_seasonal <- as.numeric(!sapply(plot_df$fam1, is.na)) + as.numeric(!sapply(plot_df$fam2, is.na))
plot_df$n_seasonal_factor <- factor(plot_df$n_seasonal)
levels(plot_df$n_seasonal_factor) <- c("no seasonal partners", "one seasonal partner", "two seasonal partners")
plot_df$both_seasonal <- plot_df$n_seasonal == 2

# Test for differences in average universality scores in seasonal pairs
# t.test(universality ~ n_seasonal_factor, plot_df %>% filter(n_seasonal %in% c(0,1))) # not signif
# t.test(universality ~ n_seasonal_factor, plot_df %>% filter(n_seasonal %in% c(0,2)))
t.test(mcs ~ both_seasonal, plot_df)
cat(paste0("Difference of means: ", round(mean(plot_df$mcs[plot_df$n_seasonal == 2]) -
                                            mean(plot_df$mcs[plot_df$n_seasonal != 0]), 3), "\n"))

p <- ggplot(plot_df, aes(x = n_seasonal_factor, y = mcs)) +
  geom_boxplot() + #width = 0.2) +
  theme_bw() +
  labs(x = "", y = "median correlation strength")

ggsave(file.path("output", "figures", "synchrony_seasonal_boxplots.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 6)

# ------------------------------------------------------------------------------
#   Previous version of "significance" thresholding
# ------------------------------------------------------------------------------

if(FALSE) {
  cat(paste0("Synchrony: median = ",
             round(median(correlations), 3),
             ", min = ",
             round(min(correlations), 3),
             ", max = ",
             round(max(correlations), 3),
             "\n"))

  idx <- which(correlations > max(correlations_permuted))
  cat(paste0("Percent of taxa with 'significant' synchrony: ",
             round((length(idx) / length(correlations))*100, 1),
             " (", length(idx), " of 134)\n"))

  tax_syn <- cbind(idx = 1:134, data$taxonomy[1:134,2:ncol(data$taxonomy)], syn = correlations)

  # Most synchronous
  View(tax_syn %>% filter(syn > 0.4) %>% arrange(family))

  # Lease synchronous
  View(tax_syn %>% arrange(syn, family) %>% slice(1:5))

  # Output Supplemental Table S6
  out_df <- data.frame(Synchrony = correlations,
                       Significant = correlations > max(correlations_permuted),
                       ASV = 1:134,
                       Taxonomy = apply(data$taxonomy[1:134,2:ncol(data$taxonomy)], 1, function(x) paste0(x, collapse = " / ")))
  write.csv(out_df %>% arrange(desc(Synchrony)),
            file.path("output", "Table_S6.csv"))
}

# ------------------------------------------------------------------------------
#   Synchrony not greater than expected by chance
#   Evaluated via p-value w/ FDR <= 0.05
# ------------------------------------------------------------------------------

f <- ecdf(correlations_permuted)

pvalues <- sapply(correlations, function(x) 2*min(f(x), 1-f(x)))
pvalues_adj <- p.adjust(pvalues, method = "BH")

prop <- sum(pvalues_adj < 0.05) / length(pvalues_adj)
cat(paste0("Significant after FDR: ", round(prop, 3), "\n"))

# ------------------------------------------------------------------------------
#   Enrichment of top center and top right-hand parts
# ------------------------------------------------------------------------------

# Score enrichment of family pairs or families themselves?
topcenter_pairs <- which(plot_df$synchrony < 0.3 & plot_df$universality > 0.4)
topright_pairs <- which(plot_df$synchrony > 0.3 & plot_df$universality > 0.4)

# topcenter_pairs <- which(plot_df$synchrony <= 0.25)
# topright_pairs <- which(plot_df$synchrony > 0.25)

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

all_pairs_noNA <- all_pairs %>%
  filter(tax1 != "NA" & tax2 != "NA")

all_pairs_noNA <- all_pairs_noNA %>%
  mutate(taxpair = paste0(tax1, " - ", tax2))

# ------------------------------------------------------------------------------
#   Family enrichment
# ------------------------------------------------------------------------------

all_family_palette <- readRDS(file.path("output", "family_palette.rds"))
all_family_pair_palette <- readRDS(file.path("output", "family-family_palette.rds"))

enrichment <- NULL

# Calculate a baseline frequency for all family-family pairs
frequencies <- table(c(all_pairs_noNA$tax1, all_pairs_noNA$tax2))

# Observed frequency in this subregion
all_pairs_noNA_tc <- all_pairs_noNA %>%
  filter(topcenter == TRUE)
frequencies_subset1 <- table(c(all_pairs_noNA_tc$tax1, all_pairs_noNA_tc$tax2))

e_obj <- plot_enrichment(frequencies,
                         frequencies_subset1,
                         type_label = "family",
                         location_label = "Low synchrony, high universality",
                         pt_sz = 2.5,
                         stroke_sz = 0.8)
p2a <- e_obj$p
enrichment <- e_obj$enrichment

# Observed frequency in this region
all_pairs_noNA_tr <- all_pairs_noNA %>%
  filter(topright == TRUE)
frequencies_subset2 <- table(c(all_pairs_noNA_tr$tax1, all_pairs_noNA_tr$tax2))

e_obj <- plot_enrichment(frequencies,
                         frequencies_subset2,
                         type_label = "family",
                         location_label = "High synchrony, high universality",
                         enrichment = enrichment,
                         pt_sz = 2.5,
                         stroke_sz = 0.8)
p2b <- e_obj$p
enrichment <- e_obj$enrichment

# ------------------------------------------------------------------------------
#   Family-pair enrichment
# ------------------------------------------------------------------------------

frequencies <- table(all_pairs_noNA$taxpair)
frequencies_subset1 <- table(all_pairs_noNA$taxpair[all_pairs_noNA$topcenter == TRUE])

e_obj <- plot_enrichment(frequencies,
                         frequencies_subset1,
                         type_label = "family-pair",
                         location_label = "Low synchrony, high universality",
                         enrichment = enrichment,
                         pt_sz = 2.5,
                         stroke_sz = 0.8)
p3a <- e_obj$p
enrichment <- e_obj$enrichment

frequencies_subset2 <- table(all_pairs_noNA$taxpair[all_pairs_noNA$topright == TRUE])

e_obj <- plot_enrichment(frequencies,
                         frequencies_subset2,
                         type_label = "family-pair",
                         location_label = "High synchrony, high universality",
                         enrichment = enrichment,
                         pt_sz = 2.5,
                         stroke_sz = 0.8)
p3b <- e_obj$p
enrichment <- e_obj$enrichment

enrichment <- enrichment %>%
  arrange(location, type, name)
colnames(enrichment) <- c("ASV family or pair name",
                          "Type",
                          "Enrichment evaluated in",
                          "Odds ratio",
                          "P-value (Fisher's exact test)",
                          "Lower 95% CI",
                          "Upper 95% CI",
                          "Threshold for small counts",
                          "Adj. p-value (Benjamini-Hochberg)")
write.table(enrichment,
            file = file.path("output", "enrichment-synchrony.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

pcol1 <- plot_grid(p2a, p3a, ncol = 1,
                   rel_heights = c(1,2.25),
                   labels = c("B", "C"),
                   label_size = 20,
                   label_x = -0.06,
                   label_y = 1.04,
                   scale = 0.95)
pcol2 <- plot_grid(p2b, p3b, ncol = 1,
                   rel_heights = c(1,2.25),
                   labels = c("D", "E"),
                   label_size = 20,
                   label_x = -0.06,
                   label_y = 1.04,
                   scale = 0.95)

pcol3 <- plot_grid(pcol1, NULL, pcol2, ncol = 3,
                   rel_widths = c(1, 0.05, 1))

p5 <- plot_grid(NULL, p0, ncol = 1, rel_heights = c(0.05, 1),
                labels = c("", "A"),
                label_size = 20,
                label_x = -0.01,
                label_y = 1.07)

p <- plot_grid(p5, pcol3, ncol = 2,
               rel_heights = c(1, 1),
               scale = 0.95)

ggsave(file.path("output", "figures", "synchrony_vs_universality.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 6,
       width = 12)

# ------------------------------------------------------------------------------
#   Synchrony calculation "cartoon" figure
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#   Panel (A): sample alignment "cartoon"
# ------------------------------------------------------------------------------

df <- data.frame(x = c(3, 4, 5, 7, 10,
                       1, 3, 6, 8, 9,
                       5, 6, 9,
                       1, 2, 7, 9, 10,
                       2, 3, 4, 8),
                 y = c(5, 5, 5, 5, 5,
                       4, 4, 4, 4, 4,
                       3, 3, 3,
                       2, 2, 2, 2, 2,
                       1, 1, 1, 1),
                 shape = c(1, 0, 0, 0, 0,
                           0, 2, 1, 0, 0,
                           0, 2, 1,
                           0, 1, 0, 2, 0,
                           2, 0, 0, 0))

df$shape <- factor(df$shape)
levels(df$shape) <- c("Unused sample", "Series 1", "Series 2")

p1 <- ggplot(df, aes(x = x, y = y, fill = shape)) +
  geom_point(size = 6, shape = 21) +
  theme_bw() +
  scale_fill_manual(values = c("#dddddd", "#d99e57", "#9e87c9")) +
  scale_x_discrete(name = "sample day",
                   limits = 1:10) +
  scale_y_discrete(name = "host",
                   limits = c("E", "D", "C", "B", "A")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom") +
  labs(fill = "")

# ------------------------------------------------------------------------------
#   Panel (B): correlation in CLR abundance for the most synchronous pair
# ------------------------------------------------------------------------------

p2 <- ggplot(data.frame(x = paired_samples_x, y = paired_samples_y),
             aes(x = x, y = y)) +
  geom_point(size = 3, shape = 21, fill = "#999999") +
  theme_bw() +
  labs(x = "model-estimated logratio abundance (host 1)",
       y = "model-estimated logratio abundance (host 2)")

# ------------------------------------------------------------------------------
#   Panel (C): aligned trajectories of most-synchronous taxon over 5 hosts
# ------------------------------------------------------------------------------

alpha <- 0.6
tax_idx <- 21
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

p3 <- ggplot()
for(this_host in hosts) {
  p3 <- p3 +
    geom_ribbon(data = plot_df %>% filter(host == this_host & coord == tax_idx),
                mapping = aes(x = day, ymin = p25, ymax = p75),
                fill = "#fdbf6f",
                alpha = alpha) +
    geom_line(data = plot_df %>% filter(host == this_host & coord == tax_idx),
              mapping = aes(x = day, y = mean),
              color = "#ff7f00",
              size = 1,
              alpha = alpha)
}
p3 <- p3 +
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

host_labels <- read.delim(file.path("output", "host_labels.tsv"),
                          header = TRUE)

cat("Order of host series (from top to bottom):\n")
for(this_host in hosts[c(1,5,3,4,2)]) {
  cat(paste0("\t", host_labels %>% filter(sname == this_host) %>% pull(host_label), "\n"))
}

# ------------------------------------------------------------------------------
#   Panel (D): distributions of observed vs. permuted synchrony scores
# ------------------------------------------------------------------------------

p4 <- ggplot(data.frame(x = c(correlations, correlations_permuted),
                        type = c(rep("observed", length(correlations)),
                                 rep("permuted", length(correlations_permuted)))),
             aes(x = x, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#59AAD7", "#aaaaaa")) +
  theme_bw() +
  labs(fill = "Source",
       x = "synchrony score")

# ------------------------------------------------------------------------------
#   How many observed synchronies exceed the 95% interval of the "randomly"
#   estimated synchrony distribution?
# ------------------------------------------------------------------------------

# cutoff <- quantile(correlations_permuted, probs = c(0.025, 0.975))
# cat(paste0(sum(correlations > cutoff[2]), " / ",
#            length(correlations),
#            " taxa exceed the 95% interval (", round(cutoff[1],3), ", ", round(cutoff[2],3), ")\n"))
cutoff <- quantile(correlations_permuted, probs = c(0.95))
cat(paste0(sum(correlations > cutoff[1]), " / ",
           length(correlations),
           " taxa exceed the 95% interval (", round(cutoff[1],3), ")\n"))

# ------------------------------------------------------------------------------
#   Assemble panels
# ------------------------------------------------------------------------------

p2_padded <- plot_grid(p2, NULL, ncol = 1,
                       rel_heights = c(1, 0.14))
p2_padded_padded <- plot_grid(NULL, p2_padded, ncol = 2,
                              rel_widths = c(0.05, 1))
prow1 <- plot_grid(p1, NULL, p2_padded_padded, ncol = 3,
                   labels = c("A", "", "B"),
                   rel_widths = c(1, 0.1, 1),
                   label_size = 18,
                   label_x = -0.04,
                   label_y = 1.01)
p3_padded <- plot_grid(NULL, p3, ncol = 2,
                       rel_widths = c(0.06, 1))
prow2 <- plot_grid(p3_padded, NULL, p4, ncol = 3,
                   rel_widths = c(1, 0.1, 1),
                   labels = c("C", "", "D"),
                   label_size = 18,
                   label_x = -0.02)
p <- plot_grid(prow1, prow2, ncol = 1,
               rel_heights = c(1.1, 1),
               scale = 0.95)

ggsave(file.path("output", "figures", "synchrony_cartoon.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 8,
       width = 10)

