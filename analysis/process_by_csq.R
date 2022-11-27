library(dplyr)
library(magrittr)

by_csq <- read.table(snakemake@input[["by_csq_unprocessed"]], header = TRUE, sep = "\t")

by_csq %>%
  filter(worst_csq %in% c("splice_acceptor_variant", "splice_donor_variant")) %>%
  group_by(context, ref, alt, methylation_level, mu) %>%
  summarise(worst_csq = "essential_splice", variant_count = sum(variant_count), singleton_count = sum(singleton_count)) %>%
  rbind(by_csq) %>%
  write.table(snakemake@output[["by_csq_processed"]], sep = "\t", quote = FALSE, row.names = FALSE)
