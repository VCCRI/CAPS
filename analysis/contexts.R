# Contexts by variant class

library(tidyverse)
library(xtable)

# TODO: finish

vars <- read_tsv(snakemake@input[["variants"]])

vars1 <- vars %>%
  group_by_at(c(
    "context", "ref", "alt", "methylation_level", "mu",
    snakemake@params[["variable"]]
  )) %>%
  summarise() %>%
  ungroup() %>%
  mutate(variant_type = ifelse(
    (ref == "A" & alt == "T") |
      (ref == "A" & alt == "C") |
      (ref == "G" & alt == "T") |
      (ref == "G" & alt == "C") |
      (alt == "A" & ref == "T") |
      (alt == "A" & ref == "C") |
      (alt == "G" & ref == "T") |
      (alt == "G" & ref == "C"),
    "transversion", "transition"
  )) %>%
  group_by_at(c(snakemake@params[["variable"]], "variant_type")) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by_at(c(snakemake@params[["variable"]])) %>%
  summarise(
    # Total number of contexts
    `Total` = sum(n),
    # Proportion of transversion contexts
    `Transversions` = round(n / (lag(n) + n), 3)
  ) %>%
    ungroup() %>%
    #TODO: make this filtering specific: vars[!is.na(...),]
  na.omit()

vars2 <- vars %>%
  group_by(context, ref, alt, methylation_level, mu, worst_csq) %>%
  summarise(variant_count = sum(variant_count)) %>%
  ungroup() %>%
  group_by(worst_csq) %>%
  summarise(mu, w = variant_count / sum(variant_count)) %>%
  ungroup() %>%
  group_by(worst_csq) %>%
  summarise(wmean_mu = sum(mu * w))

vars <- inner_join(vars1, vars2)

if (!is.null(snakemake@params[["labels_set"]])) vars <- filter(vars, .data[[snakemake@params[["variable"]]]] %in% snakemake@params[["labels_set"]])

if (!is.null(snakemake@params[["labels"]]) || !is.null(snakemake@params[["new_labels"]])) {
  if (!is.null(snakemake@params[["labels"]]) && !is.null(snakemake@params[["new_labels"]])) {
    vars[[snakemake@params[["variable"]]]] <- factor(vars[[snakemake@params[["variable"]]]],
      levels = snakemake@params[["labels"]],
      labels = snakemake@params[["new_labels"]]
    )
  } else {
    stop("Either original or new labels were not provided")
  }
}

vars %>%
  xtable(
    caption = snakemake@params[["caption"]], label = snakemake@params[["label"]],
    align = "ll|rrr", display = c("s", "s", "s", "s", "s")
  ) %>%
  print(include.rownames = FALSE) %>%
  write(snakemake@output[["Out"]])
