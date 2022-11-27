# Generates a figure showing the distribution of "adjusted proportion
# of singletons" scores by variable for multiple models

library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(tidyr)

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
  mutate(model = unlist(stringr::str_split(variable, "_[lu]conf"))[1], type = stringr::str_extract(variable, "[lu]conf")) %>%
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
if (!is.null(snakemake@params[["xlab_labels_set"]])) df <- filter(df, variable_value %in% snakemake@params[["xlab_labels_set"]])
if (!is.null(snakemake@params[["xlab_labels"]]) | !is.null(snakemake@params[["new_xlab_labels"]])) {
  if (!is.null(snakemake@params[["xlab_labels"]]) & !is.null(snakemake@params[["new_xlab_labels"]])) {
    df[["variable_value"]] <- factor(df[["variable_value"]],
      levels = snakemake@params[["xlab_labels"]],
      labels = snakemake@params[["new_xlab_labels"]]
    )
  } else {
    stop("Either original or new labels were not provided")
  }
}

pdf(snakemake@output[["plot"]])
ggplot(df) +
  aes(
    x = variable_value,
    y = score,
    color = model
  ) +
  geom_pointrange(
    aes(
      ymin = lconf,
      ymax = uconf
    ),
    size = 1.3, linewidth = 1.8, position = position_dodge(width = ifelse(is.null(snakemake@params[["dodge_width"]]), 0.65, snakemake@params[["dodge_width"]]))
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(
      vjust = snakemake@params[["xlab_vjust"]],
      hjust = snakemake@params[["xlab_hjust"]],
      angle = ifelse(is.null(snakemake@params[["xlab_angle"]]), 0, snakemake@params[["xlab_angle"]])
    ),
    text = element_text(size = 24),
    plot.margin = margin(0, 5.5, 0, 5.5)
  ) +
  scale_color_manual(snakemake@params[["legend_title"]],
    values = snakemake@params[["colors"]]
  ) +
  ylab(snakemake@params[["ylab"]]) +
  xlab(snakemake@params[["xlab"]])
dev.off()
