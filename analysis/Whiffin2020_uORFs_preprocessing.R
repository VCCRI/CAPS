# Preprocessing of Whiffin2020 files. This script filters and renames
# some of the columns and annotates the variants with mutability.

library(dplyr)
library(magrittr)
library(stringi)
library(stringr)

missing_contexts <- c(
  "AGA", "AGC", "AGG", "AGT", "ATA",
  "ATC", "ATG", "ATT", "CTA", "CTC",
  "CTG", "CTT", "GTA", "GTC", "GTG",
  "GTT", "TGA", "TGG", "TTA", "TTC",
  "TTG", "TTT", "CGA", "CGC", "CGG",
  "CGT", "GGA", "GGG", "GGT", "TGT"
)

mut_df <- read.table(snakemake@input[["mu"]], sep = "\t", header = TRUE)
vars <- read.table(snakemake@input[["inFile"]], sep = "\t", header = TRUE)

columns <- setdiff(
  colnames(vars),
  c(
    "chr", "pos", "ref", "alt", "context",
    "gnomAD_AC", "methyl_bin", "GRCh38_pos", "Kozak_sequence"
  )
)

vars[!is.na(vars$gnomAD_AC), ] %>%
  rename(methylation_level = methyl_bin) %>%
  mutate(
    variant_count = 1,
    singleton_count = ifelse(gnomAD_AC == 1, 1, 0)
  ) %>%
  # Convert 7-base context into 3-base context
  mutate(context = str_sub(context, 3, 5)) %>%
  # The mutability table is for one strand only -- need to convert
  # some context into their reverse complementary sequence
  mutate(ref = as.character(ref), alt = as.character(alt)) %>%
  mutate(ref = ifelse(
    context %in% missing_contexts,
    chartr("ATGC", "TACG", ref),
    ref
  )) %>%
  mutate(alt = ifelse(
    context %in% missing_contexts,
    chartr("ATGC", "TACG", alt),
    alt
  )) %>%
  mutate(context = ifelse(
    context %in% missing_contexts,
    stri_reverse(chartr("ATGC", "TACG", context)),
    context
  )) %>%
  # Annotate with mutability
  left_join(select(mut_df, context, ref, alt, methylation_level, mu_snp)) %>%
  rename(mu = mu_snp) %>%
  # Re-group to make the format suitable for CAPS functions
  group_by_at(c("context", "ref", "alt", "methylation_level", "mu", columns)) %>%
  mutate(distanceToCDS = ifelse(
    "distanceToCDS" %in% columns, ifelse(.data[["distanceToCDS"]] >= 50,
      ">=50 bp", "<50 bp"
    ), NA
  )) %>%
  summarise(
    variant_count = sum(variant_count),
    singleton_count = sum(singleton_count)
  ) %>%
  write.table(file = snakemake@output[["outFile"]], quote = FALSE, row.names = FALSE, sep = "\t")
