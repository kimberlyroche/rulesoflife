#' Render line plot of a given CLR taxon against diet PC1 for a given host.
#' This plot is meant to demonstrate the largely uncorrelated nature of the diet
#' data with the observed logratio abundances.
#'
#' @param sname host short name indicating which baboon's series to fit
#' @param tax_idx symmetric matrix object
#' @param counts filtered 16S count table (taxa x samples)
#' @param metadata annotations data.frame
#' @return NULL
#' @import driver
#' @import ggplot2
#' @export
plot_clr_vs_diet <- function(sname, tax_idx, counts, metadata) {
  sname_idx <- which(metadata$sname == sname)
  sub_md <- metadata[sname_idx,]
  sub_counts <- counts[,sname_idx]
  clr_counts <- clr_array(sub_counts + 0.5, parts = 1)
  plot_df <- data.frame(time = rep(1:nrow(sub_md), 2),
                        value = c(scale(sub_md$diet_PC1, center = TRUE, scale = FALSE),
                                  scale(clr_counts[tax_idx,], center = TRUE, scale = FALSE)),
                        type = c(rep("diet", nrow(sub_md)), rep("taxon", nrow(sub_md))))
  plot_df$type <- factor(plot_df$type)
  p <- ggplot(plot_df, aes(x = time, y = value, color = type)) +
    geom_line(size = 0.5) +
    geom_point() +
    xlab("sample index")
  filename <- paste0("dietPC1_vs_tax", tax_idx, "_", sname, ".png")
  output_dir <- check_dir(c("output", "figures"))
  ggsave(file.path(output_dir, filename),
         p,
         units = "in",
         height = 3,
         width = 6)
}

#' Use the "rug" plot to visualize the output of multiple model runs. If no
#' column/row orderings are specified, all plots will be clustered according
#' to the column/row orderings derived from the first set of models specified.
#'
#' Note: This ASSUMES method of interest is the Gaussian process.
#'
#' @param output_dir_list specifies the subdirectories of the model fits
#' directory in which to look for fitted model output
#' @param canonical_col_order if not NULL, order of pairs (columns) to use
#' @param canonical_row_order if not NULL, order of rows (hosts) to use
#' @return NULL
#' @import driver
#' @import ggplot2
#' @export
sensitivity_sweep <- function(output_dir_list,
                              canonical_col_order = NULL,
                              canonical_row_order = NULL) {
  if(length(output_dir_list) == 0) {
    stop("No model output directories specified!")
  }
  for(output_dir in output_dir_list) {
    for(prop in c(FALSE, TRUE)) {
      sigma_obj <- summarize_Sigmas(output_dir = output_dir,
                                    use_proportionality = prop)
      save_name <- output_dir
      if(prop) {
        save_name <- paste0(save_name, "_proportionality")
      } else {
        save_name <- paste0(save_name, "_CLR")
      }
      ordering <- plot_rug(sigma_obj$rug,
                           canonical_col_order = canonical_col_order,
                           canonical_row_order = canonical_row_order,
                           save_name = save_name)
      if(is.null(canonical_col_order)) {
        canonical_col_order <- ordering$col_order
      }
      if(is.null(canonical_row_order)) {
        canonical_row_order <- ordering$row_order
      }
    }
  }
}
