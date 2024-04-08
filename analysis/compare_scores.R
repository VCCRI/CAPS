# Calculates the increase in the size of confidence intervals by
# variable for all models with the first model specified being the
# reference

library(dplyr)
library(reshape2)
library(stringr)
library(tidyr)
library(xtable)

dfs <- list()
for (i in 1:length(snakemake@input[["scores"]])) {
  dfs[[i]] <- read.table(snakemake@input[["scores"]][i], header = TRUE, sep = "\t")
}

if (length(unique(unlist(lapply(dfs, function(df) {
  df[["variable"]]
})))) != 1) {
  stop("Found inconsistency in the list of scores")
}

dfs %>%
  Reduce((function(df1, df2) inner_join(df1, df2, by = "variable_value")), .) %>%
  select(variable_value, (starts_with("caps") | starts_with("maps")) & !ends_with("_se") & !ends_with("_sem")) %>%
  melt() %>%
  rowwise() %>%
  mutate(model = unlist(str_split(variable, "_[lu]conf"))[1], type = str_extract(variable, "[lu]conf")) %>%
  select(-variable) -> df

df[is.na(df$type), ]$type <- "score"

df <- pivot_wider(df, names_from = type, values_from = value)

# Rename models
if (!is.null(snakemake@params[["model_labels"]]) | !is.null(snakemake@params[["new_model_labels"]])) {
  if (!is.null(snakemake@params[["model_labels"]]) & !is.null(snakemake@params[["new_model_labels"]])) {
    df[["model"]] <- factor(df[["model"]],
      levels = snakemake@params[["model_labels"]],
      labels = snakemake@params[["new_model_labels"]]
    )
  } else {
    stop("Either original or new labels were not provided")
  }
}

# Rename and select a subset of variable values
if (!is.null(snakemake@params[["labels_set"]])) df <- filter(df, variable_value %in% snakemake@params[["labels_set"]])
if (!is.null(snakemake@params[["labels"]]) | !is.null(snakemake@params[["new_labels"]])) {
  if (!is.null(snakemake@params[["labels"]]) & !is.null(snakemake@params[["new_labels"]])) {
    df[["variable_value"]] <- factor(df[["variable_value"]],
      levels = snakemake@params[["labels"]],
      labels = snakemake@params[["new_labels"]]
    )
  } else {
    stop("Either original or new labels were not provided")
  }
}

df <- mutate(df, ci_size = uconf - lconf)

filter(df, model == snakemake@params[["new_model_labels"]][1]) -> df_ref

filter(df, model != snakemake@params[["new_model_labels"]][1]) %>%
  rowwise() %>%
  mutate(perc_increase = 100 * (ci_size - df_ref[df_ref$variable_value == .data[["variable_value"]], ]$ci_size) / df_ref[df_ref$variable_value == .data[["variable_value"]], ]$ci_size) %>%
  mutate(perc_increase = paste0(as.character(round(perc_increase)), "%")) %>%
  select(variable_value, model, perc_increase) %>%
  arrange(variable_value, model) %>%
  rename(`Subset` = variable_value, `Model` = model, `CI increase` = perc_increase) %>%
  na.omit() %>% #TODO: make this filtering specific: df[!is.na(...),]
  xtable(
    caption = snakemake@params[["caption"]], label = snakemake@params[["label"]],
    align = "lll|r", display = c("s", "s", "s", "s")
  ) %>%
  print(include.rownames = FALSE) %>%
  write(snakemake@output[["table"]])
