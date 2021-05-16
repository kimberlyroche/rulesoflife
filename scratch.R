source("path_fix.R")

library(rulesoflife)
library(uuid)

if(TRUE) {
  data <- load_data(tax_level = "ASV", alr_median = TRUE)
  #data <- load_data(tax_level = "family")
  #data <- load_data(tax_level = "phylum")
}

if(FALSE) {
  results <- downsample_counts(sample_no = c(10, 20, 30, 40, 50, 75, 100, 200))
  saveRDS(results, paste0("output/simdata_", UUIDgenerate(), ".rds"))
}

if(FALSE) {
  # Phylum-level "rugs"
  sensitivity_sweep(output_dir_list = c("phy_days90_diet25_scale1",
                                        "phy_days90_diet25_scale1_scrambled"))

  # Family-level "rugs"
  sensitivity_sweep(output_dir_list = c("fam_days30_diet0_scale1",
                                        "fam_days90_diet0_scale1",
                                        "fam_days90_diet0_scale2",
                                        "fam_days90_diet25_scale1",
                                        "fam_days90_diet50_scale1",
                                        "fam_days90_diet25_scale1_scrambled",
                                        "fam_days90_diet25_scale1_ALRmedian"))

  # ASV-level "rugs"
  sensitivity_sweep(output_dir_list = c("asv_days90_diet25_scale1",
                                        "asv_days90_diet25_scale1_scrambled"))
}

if(FALSE) {
  #levels <- list(phy = "phylum", fam = "family", asv = "ASV")
  levels <- list(phy = "phylum")
  for(ll in 1:length(levels)) {
    short_name <- names(levels)[ll]
    long_name <- levels[[ll]]
    for(append in c("", "_scrambled")) {
      obj <- summarize_Sigmas(output_dir = paste0(short_name, "_days90_diet25_scale1", append))
      plot_correlation_histogram(rug = obj$rug, save_name = paste0("CLR_correlation_histogram_", long_name, append))
    }
  }
}
