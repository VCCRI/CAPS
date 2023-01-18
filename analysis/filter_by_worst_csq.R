library(dplyr)
library(magrittr)

read.table(snakemake@input[["In"]], sep = "\t", header = TRUE) %>%
  filter(coverage >= 30) %>%
  filter(protein_coding == "true") %>%
  filter(worst_csq == snakemake@params[["worst_csq"]]) %>%
  group_by(context, ref, alt, mu, methylation_level) %>%
  summarise(
    variant_count = sum(variant_count),
    singleton_count = sum(singleton_count)
  ) %>%
  write.table(snakemake@output[["Out"]], sep = "\t", quote = FALSE, row.names = FALSE)
