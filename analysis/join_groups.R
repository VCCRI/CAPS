# This function joins groups of variants split into bins by creating
# nested bins. Each bin is defined by two threshold (i.e. whether the
# score is between two values), whereas nested bins can be interpreted
# as top 5%, top 10%, top 50%, etc.

library(magrittr)
library(dplyr)

# Load groups of variants split into bins
groups <- read.table(snakemake@input[["groups"]], header = TRUE)

# Regroup by variable of interest
groups <- groups %>%
  group_by_at(vars("context", "ref", "alt", "methylation_level", "mu", snakemake@params[["variable"]])) %>%
  summarise(
    singleton_count = sum(singleton_count),
    variant_count = sum(variant_count)
  ) %>%
  ungroup()

# All observed bins in their nesting order
bins <- sort(unique(groups[[snakemake@params[["variable"]]]]), decreasing = ifelse(!is.null(snakemake@params[["desc"]]), snakemake@params[["desc"]], FALSE))

df <- data.frame()
for (i in length(bins):1) {
  # Go from last bin to first (the higher the score, the higher the bin number)
  groups %>%
    # Only include the current bin and all higher scores
    filter(.data[[snakemake@params[["variable"]]]] %in% bins[i:length(bins)]) %>%
    # Drop the bins column
    group_by(context, ref, alt, methylation_level, mu) %>%
    summarise(
      singleton_count = sum(singleton_count),
      variant_count = sum(variant_count)
    ) %>%
    ungroup() %>%
    # The new table is all combinations of "context", "ref", "alt",
    # "methylation_level" and "mu", where "variant_count" and
    # "singleton_count" includes variants in the current bin and all
    # variants with higher scores
    mutate(rep(bins[i:length(bins)][1], dim(.)[1])) %>%
    as.data.frame() %>%
    rbind(df) -> df
}

colnames(df) <- c(
  "context",
  "ref",
  "alt",
  "methylation_level",
  "mu",
  "singleton_count",
  "variant_count",
  snakemake@params[["variable"]]
)
write.table(df, snakemake@output[["joined_groups"]],
            row.names = FALSE, quote = FALSE, sep = "\t")
