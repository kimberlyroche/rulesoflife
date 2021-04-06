library(rulesoflife)

# Filter data; this primarily agglomerates
filtered_data <- filter_data(tax_level = "family")

# Load filtered data
data <- load_filtered_data(tax_level = "family")
names(data)
