# Generates a figure showing the distribution of "adjusted proportion
# of singletons" scores by variable

library(dplyr)
library(ggplot2)

scores <- read.table(snakemake@input[["scores"]],
  header = TRUE,
  sep = "\t"
)

if (!is.null(snakemake@params[["xlab_labels_set"]])) scores <- filter(scores, variable_value %in% snakemake@params[["xlab_labels_set"]])

if (!is.null(snakemake@params[["xlab_labels"]]) || !is.null(snakemake@params[["new_xlab_labels"]])) {
  if (!is.null(snakemake@params[["xlab_labels"]]) && !is.null(snakemake@params[["new_xlab_labels"]])) {
    scores[["variable_value"]] <- factor(scores[["variable_value"]],
      levels = snakemake@params[["xlab_labels"]],
      labels = snakemake@params[["new_xlab_labels"]]
    )
  } else {
    stop("Either original or new labels were not provided")
  }
}

pdf(snakemake@output[["plot"]])
ggplot(scores) +
  {
    if (snakemake@params[["reorder_xlab_by_score"]]) {
      aes(x = reorder(factor(variable_value), !!sym(snakemake@params[["score_name"]])), y = !!sym(snakemake@params[["score_name"]]))
    } else {
      aes(x = factor(variable_value), y = !!sym(snakemake@params[["score_name"]]))
    }
  } +
  ylab(ifelse(snakemake@params[["score_name"]] %in% c("maps", "caps", "caps_pdd"), case_when(
    (snakemake@params[["score_name"]] == "maps") ~ "MAPS",
    (snakemake@params[["score_name"]] == "caps") ~ "CAPS",
    (snakemake@params[["score_name"]] == "caps_pdd") ~ "CAPS-PDD"
  ), stop("Score name error"))) +
  xlab(snakemake@params[["xlab"]]) +
  geom_pointrange(
    aes(
      ymin = !!sym(snakemake@params[["lconf"]]),
      ymax = !!sym(snakemake@params[["uconf"]])
    ),
    size = 1.3, linewidth = 1.8
  ) +
  {
    if (!is.null(snakemake@params[["ylim_min"]]) &&
      !is.null(snakemake@params[["ylim_max"]])) {
      ylim(
        snakemake@params[["ylim_min"]],
        snakemake@params[["ylim_max"]]
      )
    }
  } +
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
