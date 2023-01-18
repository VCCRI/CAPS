library(dplyr)
library(magrittr)

set.seed(12345)

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

# calculate the CI on CAPS using simulated
# posterior predictive distribution
get_CAPS_pdd <- function(dat, phat_sim, bPDD = FALSE) {
  # dat: data of counts of singletons per context for the class of variants of interest
  #       should have colums context, alt, methylation_level, singleton_count, variant_count
  # phat_sim: posterior predictive distributions of expected prop. sing. per context

  # number of simulations
  N <- nrow(phat_sim)

  # add context string to match context with PDD matrix
  dat$colname <- paste(dat$context, ".", dat$alt, dat$methylation_level, sep = "")

  # simulate expected counts. Could be done in one matrix operation,
  # but simulate context per context for readabilty
  exp_sim <- sapply(seq_len(nrow(dat)), function(context_index) {
    # get context string
    colname <- dat$colname[context_index]

    # get the PPD corresponding to the context
    PPD_index <- match(colname, colnames(phat_sim))
    if (length(PPD_index) != 1) stop(cat("Context not found:", colname, "\n"))

    # generate binomial simulation for context
    rbinom(N, dat$variant_count[context_index], phat_sim[, PPD_index])
  })

  # get expected prop. of sing.
  exp_ps_sim <- rowSums(exp_sim) / sum(dat$variant_count)

  # get observed proportion of sing.
  obs_count <- dat %>% summarise(
    observed = sum(singleton_count),
    variant_count = sum(variant_count)
  )

  # get PDD of observed prop. of sing.
  pobs_sim <- rbeta(N, obs_count$observed + 1, obs_count$variant_count - obs_count$observed + 1)

  # binomial simulation
  obs_ps_sim <- rbinom(N, obs_count$variant_count, pobs_sim) / obs_count$variant_count

  # CAPS
  caps_sim <- obs_ps_sim - exp_ps_sim

  # return summary or whole PDD
  if (bPDD) {
    caps_sim
  } else {
    res <- data.frame(
      caps_pdd = mean(caps_sim), caps_pdd_se = sd(caps_sim),
      caps_pdd_lconf = unname(quantile(caps_sim, probs = c(0.025))),
      caps_pdd_uconf = unname(quantile(caps_sim, probs = c(0.975)))
    )
  }
}

exp_vars <- read.table(snakemake@input[["exp_variants"]], header = TRUE, sep = "\t")

vars <- read.table(snakemake@input[["variants"]], header = TRUE, sep = "\t")

if (!("phat" %in% colnames(exp_vars)) & snakemake@params[["phat_method"]] != "PDD") {
  exp_vars <- mutate(exp_vars, phat = singleton_count / variant_count)
  warning("Using proportions observed as best estimates of per-context probabilities (phat)")
}

if (snakemake@params[["phat_method"]] == "Var") exp_vars <- select(exp_vars, context, ref, alt, methylation_level, phat)

df <- data.frame()
for (v in na.omit(unique(vars[[snakemake@params[["variable"]]]]))) {
  {    if (snakemake@params[["phat_method"]] == "Var") {
    get_CAPS(filter(vars, .data[[snakemake@params[["variable"]]]] == v),
      exp_vars,
      confint_method = ifelse(is.null(snakemake@params[["confint_method"]]), "CAPS", snakemake@params[["confint_method"]])
    )
  } else if (snakemake@params[["phat_method"]] == "PDD") {
    get_CAPS_pdd(filter(vars, .data[[snakemake@params[["variable"]]]] == v), exp_vars)
  } else {
    stop("Only 'Var' and 'PDD' methods are supported ('phat_method')")
  }  } %>%
    mutate(variable_value = v) -> x
  df <- rbind(df, x)
}

df %>%
  mutate(variable = rep(snakemake@params[["variable"]], dim(df)[1])) %>%
  write.table(file = snakemake@output[["scores"]], quote = FALSE, row.names = FALSE, sep = "\t")
