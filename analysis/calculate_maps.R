library(dplyr)

calibrate_on <- read.table(
  file = snakemake@input[["calibrate_on"]],
  sep = "\t",
  header = TRUE
) %>% mutate(proportion_singletons = singleton_count / variant_count)

calibrate_on <- calibrate_on %>% mutate(transformed_mu = unlist(lapply(mu, snakemake@params[["transformation"]])))

mut_model <- lm(proportion_singletons ~ transformed_mu,
  weights = variant_count,
  data = calibrate_on
)

vars <- read.table(
  file = snakemake@input[["variants"]],
  sep = "\t",
  header = TRUE
)
# NOTE: this removes variants that have NAs in 'extra'
# column(-s). The 'NA' group will contain variants that were not
# included in the dataset from which the annotations were taken.
if (!is.null(snakemake@params[["extra"]])) {
  vars <- vars[!is.na(vars[[snakemake@params[["extra"]]]]), ]
}

vars <- vars %>%
  mutate(transformed_mu = unlist(lapply(mu, snakemake@params[["transformation"]]))) %>%
  mutate(expected_singletons = (transformed_mu * mut_model$coefficients[[2]] +
    mut_model$coefficients[[1]]) * variant_count) %>%
  group_by_at(vars(snakemake@params[["extra"]])) %>%
  summarise(
    singleton_count = sum(singleton_count),
    expected_singletons = sum(expected_singletons),
    variant_count = sum(variant_count)
  ) %>%
  ungroup() %>%
  mutate(proportion_singletons = singleton_count / variant_count) %>%
  mutate(maps = (singleton_count - expected_singletons) / variant_count) %>%
  mutate(
    maps_sem = sqrt(proportion_singletons * (1 - proportion_singletons) / variant_count),
    maps_lconf = maps - 1.96 * maps_sem,
    maps_uconf = maps + 1.96 * maps_sem
  ) -> maps_by_variable

maps_by_variable <- maps_by_variable %>%
  mutate(variable = rep(names(maps_by_variable[, 1]), dim(maps_by_variable)[1])) %>%
  rename("variable_value" = snakemake@params[["extra"]])

if (!is.null(snakemake@params[["old_cols"]]) && !is.null(snakemake@params[["new_cols"]])) {
  names(maps_by_variable)[match(snakemake@params[["old_cols"]], names(maps_by_variable))] <- snakemake@params[["new_cols"]]
}


write.table(maps_by_variable, snakemake@output[["scores"]], sep = "\t", quote = FALSE, row.names = FALSE)
