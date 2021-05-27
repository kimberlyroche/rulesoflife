source("path_fix.R")

library(rulesoflife)
library(tidyverse)

source("ggplot_fix.R")

# ------------------------------------------------------------------------------
#   Print table of all taxa
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")
tax_df <- bind_cols(id = 1:nrow(data$taxonomy),
                    data$taxonomy[,2:ncol(data$taxonomy)],
                    OTU = data$taxonomy[,1])

write.table(tax_df,
            file = file.path("output", "taxon_IDs.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
