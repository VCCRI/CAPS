library(dplyr)
library(magrittr)

syn_by_context <- vroom::vroom(snakemake@input[["syn_by_context_and_coverage"]])
syn_by_context %>%
  group_by(context, ref, alt, methylation_level, mu) %>%
  summarise(variant_count = sum(variant_count), singleton_count = sum(singleton_count)) %>%
  mutate(proportion_singletons = singleton_count / variant_count) -> syn_by_context
lm(proportion_singletons ~ mu,
  weights = variant_count,
  data = syn_by_context
) -> model
syn_by_context$residuals <- model$residuals
syn_by_context <- select(ungroup(syn_by_context), mu, residuals)

mutation_ht <- read.table(
  file = snakemake@input[["mutation_ht"]],
  sep = "\t",
  header = TRUE
) %>% select(mu_snp, variant_type)

x <- vroom::vroom(snakemake@input[["variants"]]) %>%
  mutate(mu_snp = mu) %>%
  left_join(mutation_ht, by = "mu_snp") %>%
  left_join(syn_by_context, by = "mu") %>%
  na.omit()

if (snakemake@params[["variant_level"]]) x <- mutate(x, variant_count = 1, singleton_count = ifelse(AC == 1, 1, 0))

if (snakemake@params[["only_deleterious"]]) x <- filter(x, AM_CLASS == "pathogenic" & AM_SCORE >= 0.8)

x %>%
  group_by(context, ref, alt, variant_type, methylation_level, mu, residuals) %>%
  summarise() %>%
  as.data.frame() %>%
  arrange(variant_type, ref, alt, desc(residuals))

x <- x %>%
  filter(
    (mu == snakemake@params[["mu1"]]) |
      (mu == snakemake@params[["mu2"]])
  )
x <- x %>% mutate(full_context = paste0(context, " ", ref, "/", alt, " ", methylation_level, " (", substr(variant_type, 1, 20), ")"))

write.table(x, snakemake@output[["Out"]], sep = "\t", quote = FALSE, row.names = FALSE)
