library(dplyr)
library(magrittr)

get_CAPS <- function(dat, phat, confint_method) {
  # using Gaussian approximation of the binomial
  res <- dat %>%
    inner_join(phat, by = c("context", "ref", "alt", "methylation_level")) %>%
    mutate(
      exp = variant_count * phat, # expected number of singletons
      # variance of the binomial
      exp_var = variant_count * phat * (1 - phat)
    ) %>%
    summarise(
      observed = sum(singleton_count),
      expected = sum(exp),
      variant_count = sum(variant_count),
      # variance of sum is sum of variances
      exp_var = sum(exp_var)
    ) %>%
    mutate(
      caps = (observed - expected) / variant_count,
      obs_prop = observed / variant_count,
      # variance of the binomial
      obs_var = variant_count * obs_prop * (1 - obs_prop),
      caps_se = case_when(
        (confint_method == "MAPS") ~ sqrt(obs_var) / variant_count,
        (confint_method == "CAPS") ~ sqrt(obs_var + exp_var) / variant_count
      ),
      caps_lconf = caps - 1.96 * caps_se,
      caps_uconf = caps + 1.96 * caps_se
    ) %>%
    select(caps, caps_se, caps_uconf, caps_lconf)
}

exp_vars <- read.table(snakemake@input[["exp_variants"]], header = TRUE, sep = "\t")

vars <- read.table(snakemake@input[["variants"]], header = TRUE, sep = "\t")

if (!("phat" %in% colnames(exp_vars))) {
    exp_vars <- mutate(exp_vars, phat = singleton_count / variant_count)
    warning("Using proportions observed as best estimates of per-context probabilities (phat)")
}

df <- data.frame()
for (v in unique(vars[[snakemake@params[["variable"]]]])) {
  {
    get_CAPS(filter(vars, .data[[snakemake@params[["variable"]]]] == v),
      exp_vars,
      confint_method = ifelse(is.null(snakemake@params[["confint_method"]]), "CAPS", snakemake@params[["confint_method"]])
    )
  } %>%
    mutate(variable_value = v) -> x
  df <- rbind(df, x)
}

df %>%
  mutate(variable = rep(snakemake@params[["variable"]], dim(df)[1])) %>%
  write.table(file = snakemake@output[["scores"]], quote = FALSE, row.names = FALSE, sep = "\t")
