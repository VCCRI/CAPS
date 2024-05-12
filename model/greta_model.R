# LT 24/10/2022
#
# implement model to predict the proportion of singletons in synonymous variants using the proportion
# of singletons in intronic variants for each context

require(tidyverse)
require(greta)

syn <- read_tsv(snakemake@input[["syn_vars"]])
intron <- read_tsv(snakemake@input[["intron_vars"]])

if (snakemake@params[["group_by_context"]]) {
  syn %>%
    group_by(context, ref, alt, methylation_level, mu) %>%
    summarise(variant_count = sum(variant_count), singleton_count = sum(singleton_count)) -> syn
  intron %>%
    group_by(context, ref, alt, methylation_level, mu) %>%
    summarise(variant_count = sum(variant_count), singleton_count = sum(singleton_count)) -> intron
}

if ((dim(syn)[1] != 100) || (dim(intron)[1] != 104)) {
  stop("Input should be grouped by context (100 rows for synonymous, 104 rows for intronic)")
}

obs <- intron %>%
  rename(s_intron = singleton_count, n_intron = variant_count) %>%
  inner_join(syn, by = "mu") %>%
  rename(s_syn = singleton_count, n_syn = variant_count) %>%
  select(s_intron, n_intron, s_syn, n_syn)

n_obs <- nrow(obs)

# data
y <- as_data(c(obs$s_intron, obs$s_syn))
n <- as_data(c(obs$n_intron, obs$n_syn))

# variables and priors
# here I removed all priors
# a <- normal(0, 3)
# b <- normal(0, 3)
a <- variable()
b <- variable()
# sig <- student(3, 0, 1, truncation = c(0, Inf))
sig <- variable(lower = 0)

p <- variable(dim = n_obs)
eps <- normal(0, sig, dim = n_obs)

# operations
q <- a + b * p + eps

# likelihood
distribution(y) <- binomial(n, c(iprobit(p), iprobit(q)))

# defining the model
m <- model(a, b, sig)

gc()
draws <- mcmc(m, n_samples = 1000, warmup = 10000)
gc()
plot(draws)

# save(draws, file='model/greta_model.Rdata')

theta <- summary(draws)$statistics
a <- theta[1, 1]
b <- theta[2, 1]

# get predicted proportion of singletons for the missing contexts
missing_context <- intron %>%
  select(context, ref, alt, methylation_level, singleton_count, variant_count) %>%
  anti_join(syn, by = c("context", "ref", "alt", "methylation_level")) %>%
  mutate(phat = pnorm(a + b * qnorm(singleton_count / variant_count)))

phat <- syn %>%
  select(context, ref, alt, methylation_level, singleton_count, variant_count) %>%
  mutate(phat = singleton_count / variant_count) %>%
  rbind(missing_context)

write_tsv(phat, file = snakemake@output[["phat"]])


#####
# calculating the posterior predictive distribution of CAPS
#####

acf(draws[[1]], 200)
# need to thin the chains to get independent samples
# we will thin by 180, with a target of 1000 samples

# add  draws to the existing chains
large_mcmc <- extra_samples(draws, 179001)
indices <- seq(1, 180001, by = 180)
indep_draws <- lapply(large_mcmc, function(chain) chain[indices, ])

acf(indep_draws[[1]])

# there is still evidence of autocorrelation
# thin again by 4
indices <- seq(1, 1001, by = 4)
indep_draws <- lapply(indep_draws, function(chain) chain[indices, ])

acf(indep_draws[[1]])

# all good, cbind and save
final_draws <- rbind(
  indep_draws[[1]], indep_draws[[2]],
  indep_draws[[3]], indep_draws[[4]]
)
write.table(final_draws[1:1000, ], file = snakemake@output[["theta_sample"]], sep = "\t", row.names = FALSE)

# now simulate from phat
N <- 1000
counts <- cbind(success = phat$singleton_count, failure = phat$variant_count - phat$singleton_count)
p_sim <- apply(counts, 1, function(oneline) rbeta(N, oneline["success"] + 1, oneline["failure"] + 1))

# apply uncertainty of model to missing contexts
missing_p <- apply(
  p_sim[, 101:104], 2,
  function(onecol) pnorm(final_draws[1:N, "a"] + final_draws[1:N, "b"] * qnorm(onecol) + rnorm(N, 0, final_draws[1:N, "sig"]))
)

# check we didn't mess up
plot(colMeans(missing_p), phat$phat[101:104])
abline(a = 0, b = 1, col = 2)

# all good, replace the missing contexts probabilities
p_sim[, 101:104] <- missing_p

hist(apply(p_sim, 2, sd))
apply(missing_p, 2, sd)
# all looking nice

# add context as columns names
colnames(p_sim) <- paste(phat$context, ".", phat$alt, phat$methylation_level, sep = "")

# save simulations
write.table(p_sim, file = snakemake@output[["phat_sim"]], sep = "\t", row.names = FALSE)
