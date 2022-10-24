# LT 13/10/2022

# preliminary CAPS work

# Objective: test whether the proportion of singletons observed per context
# differ between synonymous variants and intronic variants
# using gnomAD v2.1 exome data

require(tidyverse)

syn = read_tsv("files/syn_by_context.tsv")
intron = read_tsv("files/intronic_by_context.tsv")

# join synonymous and intronic data
exomes = intron %>%
  rename(s_intron = singleton_count, n_intron = variant_count) %>%
  inner_join(syn, by='mu') %>%
  rename(s_syn = singleton_count, n_syn = variant_count) %>%
  select(s_intron, n_intron, s_syn, n_syn) %>%
  mutate(ps_intron = s_intron/n_intron, ps_syn = s_syn/n_syn)

plot(ps_syn ~ ps_intron, data=exomes)
abline(a = 0, b = 1, col = 2)

m.0 <- lm(ps_syn ~ 0 + offset(ps_intron), data=exomes)
m.1 <- lm(ps_syn ~ offset(ps_intron) + ps_intron, data=exomes)
anova(m.0, m.1)
abline(lm(ps_syn ~ ps_intron, data=exomes), col=3)

# which observations have high leverage ?
syn[order(hatvalues(m.1), decreasing =TRUE)[1:12],]

# as expected all CpG transitions

# any very small p-values
pvals <- exomes %>% rowwise() %>% mutate(pval = prop.test(c(s_intron, s_syn),c(n_intron, n_syn))$p.value)
sort(pvals$pval) [1:10]
# no way the identify line could be the favored model
# with that many very low p-values

# explore a model on transversions only

transv = intron %>%
  rename(s_intron = singleton_count, n_intron = variant_count) %>%
  filter(variant_type == "transversion") %>%
  inner_join(syn, by='mu') %>%
  rename(s_syn = singleton_count, n_syn = variant_count) %>%
  select(s_intron, n_intron, s_syn, n_syn) %>%
  mutate(ps_intron = s_intron/n_intron, ps_syn = s_syn/n_syn)

plot(ps_syn ~ ps_intron, data=transv)
abline(a = 0, b = 1, col = 2)
m.0 <- lm(ps_syn ~ 0 + offset(ps_intron), data=transv)
m.1 <- lm(ps_syn ~ offset(ps_intron) + ps_intron, data=transv)
anova(m.0, m.1)
abline(lm(ps_syn ~ ps_intron, data=transv), col=3)

# a simplest likelihood ratio test toi confirm intuition
# compare a model with one rate per context with a model with two rates

LL.null = sum(dbinom(exomes$s_intron,
                     exomes$n_intron,
                     (exomes$s_intron+exomes$s_syn)/(exomes$n_intron+exomes$n_syn),
                     log = TRUE),
              dbinom(exomes$s_syn,
                     exomes$n_syn,
                     (exomes$s_intron+exomes$s_syn)/(exomes$n_intron+exomes$n_syn),
                     log = TRUE))
LL.alt = sum(dbinom(exomes$s_intron, exomes$n_intron, exomes$s_intron/exomes$n_intron, log=TRUE),
             dbinom(exomes$s_syn, exomes$n_syn, exomes$s_syn/exomes$n_syn, log=TRUE))
pchisq(2*(LL.alt-LL.null), 100)
# model with 2 rates per context is better

# explore linear model on probit scale
plot(qnorm(ps_syn) ~ qnorm(ps_intron), data=exomes)
abline(a = 0, b = 1, col = 2)
m.0 <- lm(qnorm(ps_syn) ~ 0 + offset(qnorm(ps_intron)), data=exomes)
m.1 <- lm(qnorm(ps_syn) ~ offset(qnorm(ps_intron)) + qnorm(ps_intron), data=exomes)
anova(m.0, m.1)
abline(lm(qnorm(ps_syn) ~ qnorm(ps_intron), data=exomes), col=3)
