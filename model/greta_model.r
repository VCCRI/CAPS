# LT 24/10/2022
# 
# implement model to predict the proportion of singletons in synonymous variants using the proportion 
# of singletons in intronic variants for each context

require(tidyverse)
require(greta)

syn <- read_tsv("files/syn_by_context.tsv")
intron <- read_tsv("files/intronic_by_context.tsv")

obs <- intron %>%
  rename(s_intron = singleton_count, n_intron = variant_count) %>%
  inner_join(syn, by = 'mu') %>%
  rename(s_syn = singleton_count, n_syn = variant_count) %>%
  select(s_intron, n_intron, s_syn, n_syn)

n_obs = nrow(obs)

# data
y <- as_data(c(obs$s_intron, obs$s_syn))
n <- as_data(c(obs$n_intron, obs$n_syn))

# variables and priors
#a <- normal(0, 3)
#b <- normal(0, 3)
a <- variable()
b <- variable()
sig <- student(3, 0, 1, truncation = c(0, Inf))

p <- variable(dim = n_obs)
eps <- normal(0, sig, dim = n_obs)

# operations
q <- a + b * p + eps

# likelihood
distribution(y) <- binomial(n, c(iprobit(p), iprobit(q)))

# defining the model
m <- model(a, b, sig)

draws <- mcmc(m, n_samples = 10000, warmup = 20000)
plot(draws)