library(dplyr)
library(magrittr)

x <- vroom::vroom(snakemake@input[["variants"]])
x %>%
  group_by(full_context) %>%
    summarise(
        alphamissense = mean(AM_SCORE, na.rm = TRUE),
        alphamissense_lconf = t.test(AM_SCORE)$conf.int[1],
        alphamissense_uconf = t.test(AM_SCORE)$conf.int[2],
        ) %>%
  mutate(variable = "full_context", variable_value = full_context) %>%
  write.table(snakemake@output[["scores"]], sep = "\t", quote = FALSE, row.names = FALSE)
