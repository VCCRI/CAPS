# Merges splicing categories into "essential splice"

library(dplyr)
library(magrittr)

by_csq <- read.table(snakemake@input[["by_csq_unprocessed"]], header = TRUE, sep = "\t")

by_csq %>%
  filter(worst_csq %in% c("splice_acceptor_variant", "splice_donor_variant")) %>%
  group_by_at(setdiff(colnames(by_csq), c("variant_count", "singleton_count"))) %>%
  summarise(
    worst_csq = "essential_splice",
    variant_count = sum(variant_count),
    singleton_count = sum(singleton_count)
  ) %>%
  ungroup() %>%
  rbind(by_csq) %>%
  write.table(snakemake@output[["by_csq_processed"]], sep = "\t", quote = FALSE, row.names = FALSE)
