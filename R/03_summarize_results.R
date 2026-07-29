#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]))
} else {
  normalizePath("work/real_data_analysis/R/03_summarize_results.R")
}
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

required <- c("data.table", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
library(data.table)
library(ggplot2)

result_dir <- file.path(project_dir, "results")
figure_dir <- file.path(project_dir, "figures")
table_dir <- file.path(project_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

metrics <- fread(file.path(result_dir, "split_metrics.csv"))
groups <- fread(file.path(result_dir, "scale_group_metrics.csv"))
scores <- fread(file.path(result_dir, "score_pivotality.csv"))
descriptor <- fread(file.path(project_dir, "data", "processed", "dataset_descriptives.csv"))

method_order <- c("Parametric-TI", "SR-TI", "ASR-TI", "CQR-TI")
metrics[, method := factor(method, levels = method_order)]
groups[, method := factor(method, levels = method_order)]
scores[, method := factor(method, levels = method_order)]

summary_table <- metrics[, .(
  repetitions = .N,
  coverage_mean = mean(coverage),
  coverage_sd = sd(coverage),
  coverage_min = min(coverage),
  pac_success_rate = mean(coverage >= 0.90),
  mean_width = mean(mean_width),
  median_width = mean(median_width),
  q90_width = mean(q90_width),
  lower_miss = mean(lower_miss),
  upper_miss = mean(upper_miss)
), by = .(dataset, dataset_label, method)]

conditional_table <- groups[, .(
  coverage_mean = mean(coverage),
  coverage_sd = sd(coverage),
  mean_width = mean(mean_width)
), by = .(dataset, method, group)]

conditional_worst <- groups[, {
  split_summary <- .(
    min_group_coverage = min(coverage),
    max_group_coverage = max(coverage),
    max_abs_gap = max(abs(coverage - 0.90))
  )
  split_summary
}, by = .(dataset, method, repetition)][, .(
  mean_min_group_coverage = mean(min_group_coverage),
  worst_group_coverage = min(min_group_coverage),
  mean_max_abs_gap = mean(max_abs_gap)
), by = .(dataset, method)]

pivotality_table <- scores[, .(
  mean_ks = mean(ks_to_marginal),
  q90_ks = unname(quantile(ks_to_marginal, 0.90)),
  max_ks = max(ks_to_marginal)
), by = .(dataset, method)]

fwrite(summary_table, file.path(table_dir, "method_summary.csv"))
fwrite(conditional_table, file.path(table_dir, "scale_quintile_summary.csv"))
fwrite(conditional_worst, file.path(table_dir, "conditional_worst_case.csv"))
fwrite(pivotality_table, file.path(table_dir, "pivotality_summary.csv"))

label_map <- unique(metrics[, .(dataset, dataset_label)])
short_labels <- c(
  cps_wage = "CPS wages",
  superconductivity = "Superconductivity",
  meps_expenditure = "MEPS expenditure",
  meps_utilization = "MEPS utilization"
)
primary_order <- c("cps_wage", "superconductivity", "meps_expenditure")
summary_table <- merge(summary_table, conditional_worst, by = c("dataset", "method"))
summary_table <- merge(summary_table, pivotality_table, by = c("dataset", "method"), all.x = TRUE)
summary_table[, dataset_short := unname(short_labels[dataset])]
setcolorder(
  summary_table,
  c(
    "dataset", "dataset_label", "method", "repetitions",
    "coverage_mean", "coverage_sd", "coverage_min", "pac_success_rate",
    "mean_width", "median_width", "q90_width",
    "lower_miss", "upper_miss",
    "mean_min_group_coverage", "worst_group_coverage", "mean_max_abs_gap",
    "mean_ks", "q90_ks", "max_ks"
  )
)
fwrite(summary_table, file.path(table_dir, "manuscript_summary.csv"))

theme_set(theme_minimal(base_size = 12))
method_colors <- c(
  "Parametric-TI" = "#6B7280",
  "SR-TI" = "#0072B2",
  "ASR-TI" = "#D55E00",
  "CQR-TI" = "#009E73"
)

primary_summary <- summary_table[dataset %in% primary_order]
primary_summary[, dataset_short := factor(dataset_short, levels = short_labels[primary_order])]

p_coverage <- ggplot(
  primary_summary,
  aes(x = method, y = coverage_mean, color = method)
) +
  geom_hline(yintercept = 0.90, linetype = 2, color = "black") +
  geom_errorbar(
    aes(
      ymin = pmax(0, coverage_mean - coverage_sd),
      ymax = pmin(1, coverage_mean + coverage_sd)
    ),
    width = 0.15
  ) +
  geom_point(size = 2.6) +
  facet_wrap(~ dataset_short, nrow = 1) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  coord_cartesian(ylim = c(0.86, 0.96)) +
  labs(
    x = NULL,
    y = "Held-out empirical coverage (mean ± SD)",
    color = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1),
    legend.position = "none"
  )
ggsave(
  file.path(figure_dir, "coverage_by_method.png"),
  p_coverage, width = 12.5, height = 4.3, dpi = 220
)
ggsave(
  file.path(figure_dir, "coverage_by_method.pdf"),
  p_coverage, width = 12.5, height = 4.3, device = "pdf"
)

p_width <- ggplot(
  primary_summary,
  aes(x = method, y = mean_width, fill = method)
) +
  geom_col(width = 0.72) +
  facet_wrap(~ dataset_short, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = NULL, y = "Mean interval width (original response scale)", fill = NULL) +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1),
    legend.position = "none"
  )
ggsave(
  file.path(figure_dir, "width_by_method.png"),
  p_width, width = 12.5, height = 4.3, dpi = 220
)
ggsave(
  file.path(figure_dir, "width_by_method.pdf"),
  p_width, width = 12.5, height = 4.3, device = "pdf"
)

conditional_plot_data <- merge(conditional_table, label_map, by = "dataset")
conditional_plot_data <- conditional_plot_data[dataset %in% primary_order]
conditional_plot_data[, dataset_short := factor(
  unname(short_labels[dataset]), levels = short_labels[primary_order]
)]
p_conditional <- ggplot(
  conditional_plot_data,
  aes(x = group, y = coverage_mean, color = method, group = method)
) +
  geom_hline(yintercept = 0.90, linetype = 2, color = "black") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ dataset_short, nrow = 1) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_x_continuous(breaks = 1:5) +
  coord_cartesian(ylim = c(0.68, 1.00)) +
  labs(
    x = "Predicted-scale quintile (low to high)",
    y = "Held-out coverage",
    color = NULL
  ) +
  theme(legend.position = "bottom")
ggsave(
  file.path(figure_dir, "coverage_by_predicted_scale.png"),
  p_conditional, width = 12.5, height = 4.5, dpi = 220
)
ggsave(
  file.path(figure_dir, "coverage_by_predicted_scale.pdf"),
  p_conditional, width = 12.5, height = 4.5, device = "pdf"
)

pivotality_plot_data <- merge(pivotality_table, label_map, by = "dataset")
pivotality_plot_data <- pivotality_plot_data[dataset %in% primary_order]
pivotality_plot_data[, dataset_short := factor(
  unname(short_labels[dataset]), levels = short_labels[primary_order]
)]
p_pivotality <- ggplot(
  pivotality_plot_data,
  aes(x = method, y = q90_ks, fill = method)
) +
  geom_col(width = 0.72) +
  facet_wrap(~ dataset_short, nrow = 1) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(
    x = NULL,
    y = "90th percentile KS distance to marginal score CDF",
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none"
  )
ggsave(
  file.path(figure_dir, "score_pivotality.png"),
  p_pivotality, width = 11.5, height = 4.1, dpi = 220
)
ggsave(
  file.path(figure_dir, "score_pivotality.pdf"),
  p_pivotality, width = 11.5, height = 4.1, device = "pdf"
)

sensitivity <- summary_table[dataset == "meps_utilization"]
fwrite(sensitivity, file.path(table_dir, "meps_utilization_sensitivity.csv"))

format_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
latex_rows <- summary_table[
  dataset != "meps_utilization",
  paste(
    gsub("&", "\\\\&", dataset_label),
    as.character(method),
    format_num(coverage_mean),
    format_num(pac_success_rate),
    format_num(mean_width, 2),
    format_num(q90_width, 2),
    format_num(lower_miss),
    format_num(upper_miss),
    sep = " & "
  )
]
latex <- c(
  "\\begin{table}[t]",
  "\\centering",
  "\\caption{Real-data empirical coverage and interval width over repeated random splits.}",
  "\\label{tab:real-data-new}",
  "\\begin{tabular}{llrrrrrr}",
  "\\toprule",
  "Dataset & Method & Coverage & PAC success & Mean width & Q90 width & Lower miss & Upper miss \\\\",
  "\\midrule",
  paste0(latex_rows, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(latex, file.path(table_dir, "real_data_table.tex"))

cat("Summary:\n")
print(summary_table[, .(
  dataset,
  method,
  coverage = round(coverage_mean, 3),
  pac_success = round(pac_success_rate, 3),
  mean_width = round(mean_width, 2),
  min_group = round(mean_min_group_coverage, 3),
  mean_ks = round(mean_ks, 3)
)])
