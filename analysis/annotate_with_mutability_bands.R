# Annotates variants with mutability bands based on mutation rates ("mu" column)

library(dplyr)

vars <- read.table(snakemake@input[["variants"]],
  header = TRUE,
  sep = "\t"
)

vars <- mutate(vars, mutability_band = case_when(
  mu >= quantile(mu)[1] & mu < quantile(mu)[2] ~ "Lowest",
  mu >= quantile(mu)[2] & mu < quantile(mu)[3] ~ "Lower",
  mu >= quantile(mu)[3] & mu < quantile(mu)[4] ~ "Higher",
  mu >= quantile(mu)[4] ~ "Highest",
  TRUE ~ ""
))

write.table(vars, snakemake@output[["annotated_variants"]], sep = "\t", quote = FALSE, row.names = FALSE)
