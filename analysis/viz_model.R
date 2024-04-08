# Generates a figure showing the fit of the model

library(dplyr)
library(ggplot2)

mutation_ht <- read.table(
  file = snakemake@input[["mutation_ht"]],
  sep = "\t",
  header = TRUE
)

calibrate_on <- read.table(
  file = snakemake@input[["calibrate_on"]],
  sep = "\t",
  header = TRUE
)

# To ensure that the script is usable with different input files
calibrate_on <- group_by(calibrate_on, mu) %>%
  summarise(variant_count = sum(variant_count), singleton_count = sum(singleton_count)) %>%
  mutate(proportion_singletons = singleton_count / variant_count)

calibrate_on <- calibrate_on %>% mutate(transformed_mu = unlist(lapply(mu, snakemake@params[["transformation"]])))

mut_model <- lm(proportion_singletons ~ transformed_mu,
  weights = variant_count,
  data = calibrate_on
)

pdf(snakemake@output[["plot"]])
calibrate_on %>%
  mutate(mu_snp = mu) %>%
  left_join(mutation_ht, by = "mu_snp") %>%
  ggplot() +
  aes(x = transformed_mu, y = proportion_singletons, color = variant_type) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    aspect.ratio = 1,
    text = element_text(size = ifelse(is.null(snakemake@params[["text_size"]]), 24, snakemake@params[["text_size"]]))
  ) +
  scale_color_manual(snakemake@params[["legend_title"]],
    values = snakemake@params[["colors"]]
  ) +
  geom_abline(
    intercept = mut_model$coefficients[[1]],
    slope = mut_model$coefficients[[2]]
  ) +
  xlab("Mutability\n(mutation rate per base pair per generation)") +
  ylab("Proportion of singletons")
dev.off()
