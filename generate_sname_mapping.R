library(rulesoflife)
library(tidyverse)

host_labels <- load_data()$metadata %>%
  select(sname, sex) %>%
  distinct() %>%
  slice(sample(1:n())) %>%
  group_by(sex) %>%
  mutate(host_number = 1:n()) %>%
  ungroup() %>%
  mutate(host_label = paste0(sex, host_number))

write.table(host_labels %>% select(sname, host_label) %>% arrange(sname),
            file.path("output", "host_labels.tsv"),
            row.names = FALSE,
            sep = "\t")
