# Applies QC filters on coverage and transcript consequences and
# merges splicing categories into "essential splice"

library(dplyr)
library(magrittr)

vars <- read.table(snakemake@input[["In"]], sep = "\t", header = TRUE)
if (!is.null(snakemake@params[["min_coverage"]])) vars <- filter(vars, coverage >= snakemake@params[["min_coverage"]])
if (!is.null(snakemake@params[["protein_coding"]]) && (snakemake@params[["protein_coding"]] == TRUE)) vars <- filter(vars, protein_coding == "true")

if (!is.null(snakemake@params[["add_essential_splice_cat"]]) && (snakemake@params[["add_essential_splice_cat"]] == TRUE)) {
  vars <- vars %>%
    filter(worst_csq %in% c("splice_acceptor_variant", "splice_donor_variant")) %>%
    group_by_at(setdiff(colnames(vars), c("variant_count", "singleton_count"))) %>%
    summarise(
      worst_csq = "essential_splice",
      variant_count = sum(variant_count),
      singleton_count = sum(singleton_count)
    ) %>%
    ungroup() %>%
    rbind(vars)
}

if (!is.null(snakemake@params[["clinvar_remove_likely_categories"]]) && (snakemake@params[["clinvar_remove_likely_categories"]] == TRUE)) {
  vars <- vars %>%
    filter(ClinicalSignificance %in% c("Benign", "Benign/Likely benign", "Likely benign")) %>%
    group_by_at(setdiff(colnames(vars), c("variant_count", "singleton_count"))) %>%
    summarise(
      ClinicalSignificance = "All benign",
      variant_count = sum(variant_count),
      singleton_count = sum(singleton_count)
    ) %>%
    ungroup() %>%
    rbind(vars)
  vars <- vars %>%
    filter(ClinicalSignificance %in% c("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic")) %>%
    group_by_at(setdiff(colnames(vars), c("variant_count", "singleton_count"))) %>%
    summarise(
      ClinicalSignificance = "All pathogenic",
      variant_count = sum(variant_count),
      singleton_count = sum(singleton_count)
    ) %>%
    ungroup() %>%
    rbind(vars)
}

if (!is.null(snakemake@params[["variable"]])) vars <- filter(vars, .data[[snakemake@params[["variable"]]]] == snakemake@params[["variable_value"]])

write.table(vars, snakemake@output[["Out"]], sep = "\t", quote = FALSE, row.names = FALSE)
