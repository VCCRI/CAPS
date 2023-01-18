# Counts by variant class

library(ggplot2)
library(dplyr)
library(magrittr)
library(scales)
library(xtable)

vars <- read.table(snakemake@input[["variants"]], sep = "\t", header = TRUE)

if (!is.null(snakemake@params[["coverage_threshold"]])) vars <- filter(vars, coverage >= snakemake@params[["coverage_threshold"]])

if (!is.null(snakemake@params[["labels_set"]])) vars <- filter(vars, .data[[snakemake@params[["variable"]]]] %in% snakemake@params[["labels_set"]])

if (!is.null(snakemake@params[["labels"]]) || !is.null(snakemake@params[["new_labels"]])) {
  if (!is.null(snakemake@params[["labels"]]) && !is.null(snakemake@params[["new_labels"]])) {
    vars[[snakemake@params[["variable"]]]] <- factor(vars[[snakemake@params[["variable"]]]],
      levels = snakemake@params[["labels"]],
      labels = snakemake@params[["new_labels"]]
    )
  } else {
    stop("Either original or new labels were not provided")
  }
}

if (snakemake@params[["mode"]] == "pdf") {
  pdf(snakemake@output[["Out"]])
  vars %>%
    group_by_at(vars(snakemake@params[["variable"]])) %>%
    summarise(variant_count = sum(variant_count)) %>%
    ggplot() +
    aes(x = !!sym(snakemake@params[["variable"]]), y = variant_count) +
    geom_bar(stat = "identity") +
    xlab(snakemake@params[["xlab"]]) +
    ylab(snakemake@params[["ylab"]]) +
    scale_y_continuous(labels = label_number()) +
    theme_classic() +
    theme(
      aspect.ratio = ifelse(!is.null(snakemake@params[["aspect_ratio"]]), snakemake@params[["aspect_ratio"]], 1),
      legend.direction = "horizontal",
      text = element_text(size = 24),
      axis.text.x = element_text(
        vjust = snakemake@params[["xlab_vjust"]],
        hjust = snakemake@params[["xlab_hjust"]],
        angle = ifelse(is.null(snakemake@params[["xlab_angle"]]), 0, snakemake@params[["xlab_angle"]])
      ),
      plot.margin = margin(0, 5.5, 0, 5.5)
    )
  dev.off()
}

if (snakemake@params[["mode"]] == "tex") {
  vars %>%
    group_by_at(vars(snakemake@params[["variable"]])) %>%
    summarise(variant_count = sum(variant_count)) %>%
    xtable(
      caption = snakemake@params[["caption"]], label = snakemake@params[["label"]],
      align = "ll|r", display = c("s", "s", "s")
    ) %>%
  print(include.rownames = FALSE) %>%
  write(snakemake@output[["Out"]])
}
