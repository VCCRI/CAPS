# LT 10/11/2022
# first CAPS prototype

library(tidyverse)

# read proportion of singleton per context
phat <- read_tsv("model/phat.tsv")

# read the data to apply CAPS to
# here we use synonymous variants
dat <- read_tsv("files/syn_by_context.tsv")

# calculate CAPS and its standard error (and MAPS se)
get_CAPS <- function (dat, phat) {
    # using Gaussian approximation of the binomial
    res <- dat %>%
        inner_join(phat, by = c('context', 'ref', 'alt', 'methylation_level')) %>%
        mutate(
            exp = variant_count * phat, # expected number of singletons
            # variance of the binomial
            exp_var = variant_count * phat * (1 - phat)) %>%
        summarise(
            observed = sum(singleton_count),
            expected = sum(exp),
            variant_count = sum(variant_count),
            # variance of sum is sum of variances
            exp_var = sum(exp_var)) %>%
        mutate(
            CAPS = (observed - expected) / variant_count,
            obs_prop = observed / variant_count,
            # variance of the binomial
            obs_var = variant_count * obs_prop * (1 - obs_prop),
            CAPS_se = sqrt(obs_var + exp_var) / variant_count,
            MAPS_se = sqrt(obs_var) / variant_count) %>%
        select(CAPS, CAPS_se, MAPS_se)
}

print(get_CAPS(dat, phat))

# load missense data
mis_dat <- read_tsv("files/missense_by_context.tsv")
print(get_CAPS(mis_dat, phat))

# calculate the CI on CAPS using simulated
# posterior predictive distribution
get_CAPS_pdd <- function (dat, phat_sim, bPDD = FALSE) {
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
        variant_count = sum(variant_count))
            
    # get PDD of observed prop. of sing.
    pobs_sim <- rbeta(N, obs_count$observed + 1, obs_count$variant_count - obs_count$observed + 1)

    # binomial simulation
    obs_ps_sim <- rbinom(N, obs_count$variant_count, pobs_sim) / obs_count$variant_count

    # CAPS
    caps_sim <- obs_ps_sim - exp_ps_sim

    # return summary or whole PDD
    if (bPDD)
        caps_sim
    else
       c(CAPS = mean(caps_sim), CAPS_se = sd(caps_sim), quantile(caps_sim, probs = c(0.025, 0.975)))
}

# load PPD of expected prop. of sing. by context
phat_sim <- read.delim("model/phat_sim.tsv")

# synonymous
print(get_CAPS(dat, phat))
print(get_CAPS_pdd(dat, phat_sim))

# missense
print(get_CAPS(mis_dat, phat))
print(get_CAPS_pdd(mis_dat, phat_sim))

# load new file provided by Mikhail with all consequences
all_csq <- read_tsv("files/by_csq.tsv")
stop_gained_dat <- all_csq %>% filter(worst_csq == 'stop_gained')

# stop gain
print(get_CAPS(stop_gained_dat, phat))
print(get_CAPS_pdd(stop_gained_dat, phat_sim))
