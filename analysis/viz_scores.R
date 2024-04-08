# Generates a figure showing the distribution of "adjusted proportion
# of singletons" scores by variable

library(dplyr)
library(ggplot2)
library(ggprism)
library(BSDA)

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

# TODO: make this filtering specific: scores[!is.na(...),]
if (!is.null(snakemake@params[["NA_omit"]]) && snakemake@params[["NA_omit"]]) {
  scores <- na.omit(scores)
}

if (!is.null(snakemake@params[["add_pvalues"]]) && snakemake@params[["add_pvalues"]] == TRUE) {
  pvals <- combn(m = 2, x = scores$variable_value, simplify = FALSE) %>%
    as.data.frame() %>%
    t() %>%
    as_tibble() %>%
    rename(group1 = V1, group2 = V2) %>%
    filter(group1 != group2)
  if (snakemake@params[["pvalue_test"]] == "Two-sample t-test") {
    pvals <- pvals %>%
      rowwise() %>%
      mutate(
        SE = sqrt(
          (scores[scores$variable_value == group1, ][[paste0(snakemake@params[["score_name"]], "_sem")]])^2 +
            (scores[scores$variable_value == group2, ][[paste0(snakemake@params[["score_name"]], "_sem")]])^2
        ),
        z = (
          scores[scores$variable_value == group1, ][[snakemake@params[["score_name"]]]] -
            scores[scores$variable_value == group2, ][[snakemake@params[["score_name"]]]]) / SE,
        p = 2 * (1 - pnorm(abs(z)))
      )
  } else if (snakemake@params[["pvalue_test"]] == "Welch") {
    pvals <- pvals %>%
      rowwise() %>%
      mutate(
        p = tsum.test(
          mean.x = scores[scores$variable_value == group1, ][[snakemake@params[["score_name"]]]],
          s.x = scores[scores$variable_value == group1, ][[paste0(snakemake@params[["score_name"]], "_sem")]] * sqrt(scores[scores$variable_value == group1, ]$variant_count),
          n.x = scores[scores$variable_value == group1, ]$variant_count,
          mean.y = scores[scores$variable_value == group2, ][[snakemake@params[["score_name"]]]],
          s.y = scores[scores$variable_value == group2, ][[paste0(snakemake@params[["score_name"]], "_sem")]] * sqrt(scores[scores$variable_value == group2, ]$variant_count),
          n.y = scores[scores$variable_value == group2, ]$variant_count,
        )$p.value
      )
  }
  if (snakemake@params[["bonferroni"]] == TRUE) {
    pvals <- as.data.frame(pvals) %>% mutate(p = round(p.adjust(p, "bonferroni"), 5))
  } else {
    pvals <- as.data.frame(pvals) %>% mutate(p = round(p, 5))
  }
  if (snakemake@params[["pvalue_stars"]] == TRUE) {
    pvals <- mutate(pvals, p = case_when(
      (p < 0.0001) ~ "****",
      (p < 0.001) ~ "***",
      (p < 0.01) ~ "**",
      (p < 0.05) ~ "*",
      TRUE ~ "ns"
    ))
  }
  pvals <- pvals %>% mutate(y.position = max(scores[[paste0(snakemake@params[["score_name"]], "_uconf")]]) + 0.01 * as.integer(rownames(pvals)))
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
  ylab(ifelse(snakemake@params[["score_name"]] %in% c("maps", "caps", "caps_ppd"), case_when(
    (snakemake@params[["score_name"]] == "maps") ~ "MAPS",
    (snakemake@params[["score_name"]] == "caps") ~ "CAPS",
    (snakemake@params[["score_name"]] == "caps_ppd") ~ "CAPS-PPD"
  ), stop("Score name error"))) +
  {
    if (!is.null(snakemake@params[["ylab"]])) {
      ylab(
        snakemake@params[["ylab"]]
      )
    }
  } +
  xlab(snakemake@params[["xlab"]]) +
  geom_pointrange(
    aes(
      ymin = !!sym(paste0(snakemake@params[["score_name"]], "_lconf")),
      ymax = !!sym(paste0(snakemake@params[["score_name"]], "_uconf"))
    ),
    size = 1.3, linewidth = 1.8,
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
      size = ifelse(is.null(snakemake@params[["xlab_size"]]), 24, snakemake@params[["xlab_size"]]),
      vjust = snakemake@params[["xlab_vjust"]],
      hjust = snakemake@params[["xlab_hjust"]],
      angle = ifelse(is.null(snakemake@params[["xlab_angle"]]), 0, snakemake@params[["xlab_angle"]])
    ),
    plot.margin = margin(0, 5.5, 0, 5.5)
  ) +
  {
    if (!is.null(snakemake@params[["add_pvalues"]]) && snakemake@params[["add_pvalues"]] == TRUE) {
      add_pvalue(pvals, label.size = 4.5)
    }
  }
dev.off()
