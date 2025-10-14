# Overlap of significant DEG sets between cell lines

library(tidyverse)
library(data.table)
library(patchwork)
library(UpSetR)
library(grid)
library(ggvenn)

pairs_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs"
QC_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/data/perturb_QC"
plots_dir <- "/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_exploratory/plots/cross_cell_line"
cell_lines <- c("K562_essential", "Jurkat", "RPE1", "HepG2")

# Construct paths
pairs_paths <- file.path(pairs_dir, cell_lines, "perturbation_pairs_eff.csv")
base_expression_paths <- file.path(QC_dir, cell_lines, "base_expression.csv")

# Read all CSVs into a list
pairs_list <- lapply(pairs_paths, fread)
base_expression_list <- lapply(base_expression_paths, fread)

# Optionally, add a column indicating the cell line
pairs_list <- Map(function(df, cl) {
  df$cell_line <- cl
  df
}, pairs_list, cell_lines)
base_expression_list <- Map(function(df, cl) {
  df$cell_line <- cl
  df
}, base_expression_list, cell_lines)

# Combine all into a single data frame
pairs <- rbindlist(pairs_list, fill=TRUE)
base_expression <- rbindlist(base_expression_list)

deg_sets <- pairs[, .(deg = list(unique(effect))), by = .(perturbation, cell_line)]

perturb_cell_counts <- deg_sets[, .N, by = perturbation]

perturb_in_all_CLs <- deg_sets[, uniqueN(cell_line), by = perturbation][V1 == 4, perturbation]

pairs <- pairs[perturbation %in% perturb_in_all_CLs]
deg_sets <- deg_sets[perturbation %in% perturb_in_all_CLs]

# Function to compute pairwise Jaccard overlaps
jaccard_index <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) return(NA_real_)
  inter / union
}

# Function for pairwise Jaccard on one perturb
compute_jaccard_mat <- function(df) {
  cell_lines <- df$cell_line
  degs <- df$deg
  m <- matrix(NA, nrow = length(cell_lines), ncol = length(cell_lines),
              dimnames = list(cell_lines, cell_lines))
  for (i in seq_along(cell_lines)) {
    for (j in seq_along(cell_lines)) {
      m[i, j] <- jaccard_index(degs[[i]], degs[[j]])
    }
  }
  m
}

# Apply to all perturbations
jaccard_list <- deg_sets[, .(jaccard_mat = list(compute_jaccard_mat(.SD))), by = perturbation]

# function to extract mean pairwise Jaccard (excluding diagonal)
mean_jaccard <- function(x) {
  mat <- matrix(as.numeric(x), nrow = 4, ncol = 4)
  upper_vals <- mat[upper.tri(mat)]
  mean(upper_vals)
}

# compute mean Jaccard for each perturbation
jaccard_summary <- jaccard_list[, .(
  mean_jaccard = mean_jaccard(jaccard_mat[[1]])
), by = perturbation]

mean_jaccard.p <- ggplot(jaccard_summary, aes(x = mean_jaccard)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of cross-cell-line Jaccard similarity",
    x = "Mean Jaccard index (across 4 cell lines)",
    y = "Number of perturbations"
  )

ggsave(
    filename = file.path(plots_dir, "mean_jaccard_all.png"),
    plot = mean_jaccard.p,
    width = 8, height = 8, dpi = 300
  )

plots <- lapply(cell_lines, function(cl) {
  
  dt <- base_expression[cell_line == cl]
  
  # Compute number and % of lowly expressed genes
  n_total <- nrow(dt)
  n_low <- sum(dt$base_mean_per_cell < 0.1)
  pct_low <- round(100 * n_low / n_total, 1)
  
  ggplot(dt, aes(x = base_mean_per_cell)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    scale_x_log10() +
    geom_vline(xintercept = 0.1, color = "red", linetype = "dashed") +
    labs(
      title = glue("{cl} ({n_low} / {n_total} ({pct_low}%) genes < 0.1)"),
      x = "Base mean per cell (log scale)",
      y = "Number of genes"
    ) +
    theme_bw()
})

# Combine plots with patchwork
base_expression.p <- wrap_plots(plots, ncol = 2)

ggsave(
  filename = file.path(plots_dir, "base_expression_dist.png"),
  plot = base_expression.p,
  width = 12, height = 12, dpi = 300
)

expressed_genes_list <- list(
  K562   = base_expression[cell_line == "K562_essential" & base_mean_per_cell > 0.1, gene_id],
  HepG2  = base_expression[cell_line == "HepG2" & base_mean_per_cell > 0.1, gene_id],
  Jurkat = base_expression[cell_line == "Jurkat" & base_mean_per_cell > 0.1, gene_id],
  RPE1   = base_expression[cell_line == "RPE1" & base_mean_per_cell > 0.1, gene_id]
)

png(file.path(plots_dir, "upset_expressed_genes.png"), width = 1200, height = 800, res = 150)
upset(
  fromList(expressed_genes_list),
  mainbar.y.label = "# genes in intersection (mean expression per cell > 0.1)",
  sets.x.label = "# expressed genes per cell line"
)
dev.off()

expressed_genes.venn <- ggvenn(
  expressed_genes_list,
  fill_color = c("red", "green", "blue", "yellow"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 3,
  show_percentage = FALSE
) + 
  ggtitle("Overlap of expressed genes across 4 cell lines (mean expression per cell > 0.1)")

# Save to file
ggsave(
  filename = file.path(plots_dir, "venn_expressed_genes.png"),
  plot = expressed_genes.venn,
  width = 8, height = 6, dpi = 300
)


common_expressed_genes <- Reduce(intersect, expressed_genes_list)

pairs <- pairs[effect_ensg %in% common_expressed_genes]
deg_sets <- pairs[, .(deg = list(unique(effect))), by = .(perturbation, cell_line)]

jaccard_list <- deg_sets[, .(jaccard_mat = list(compute_jaccard_mat(.SD))), by = perturbation]

# function to extract mean pairwise Jaccard (excluding diagonal)
mean_jaccard <- function(x) {
  mat <- matrix(as.numeric(x), nrow = 4, ncol = 4)
  upper_vals <- mat[upper.tri(mat)]
  mean(upper_vals)
}

# compute mean Jaccard for each perturbation
jaccard_summary <- jaccard_list[, .(
  mean_jaccard = mean_jaccard(jaccard_mat[[1]])
), by = perturbation]

mean_jaccard.p <- ggplot(jaccard_summary, aes(x = mean_jaccard)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of cross-cell-line Jaccard similarity (expressed genes)",
    x = "Mean Jaccard index (across 4 cell lines)",
    y = "Number of perturbations"
  )

ggsave(
    filename = file.path(plots_dir, "mean_jaccard_expressed.png"),
    plot = mean_jaccard.p,
    width = 8, height = 8, dpi = 300
  )

