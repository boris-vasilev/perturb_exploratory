library(Seurat)
library(sceasy)
library(reticulate)
library(tidyverse)
library(patchwork)
library(here)
source(here("src/functions_ashr.R"))
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflicts_prefer(base::setdiff)

use_condaenv('env_nf')

print("Read Jurkat-Essential")
seurat_obj.J <- sceasy::convertFormat(here("data/perturb/seurat/GSE264667_jurkat_raw_singlecell_01.h5ad"), from="anndata", to="seurat")
print("Read K562-Essential")
seurat_obj.K <- sceasy::convertFormat(here("data/perturb/seurat/K562_essential_raw_singlecell_01.h5ad"), from="anndata", to="seurat")


get_baseline_table <- function(seurat_obj, assay = "RNA", use_log = FALSE, mc.cores = 2) {
  all_perturbations <- unique(seurat_obj$gene)
  all_perturbations <- setdiff(all_perturbations, "non-targeting")
  n_perturbations <- length(all_perturbations)
  
  message("Starting baseline expression calculation for ", n_perturbations, " perturbations.")
  expr_mat <- if (use_log) {
    message("Using log-normalized expression values from assay: ", assay)
    seurat_obj[[assay]]@data
  } else {
    message("Using raw counts from assay: ", assay)
    seurat_obj[[assay]]@counts
  }
  
  # Define function to run in parallel
  calc_baseline <- function(pert) {
    message("Processing: ", pert)
    cells_subset <- WhichCells(seurat_obj, expression = gene %in% c(pert, "non-targeting"))
    sub_expr <- expr_mat[, cells_subset, drop = FALSE]
    baseline_expr <- rowMeans(sub_expr)
    data.frame(
      effect = rownames(sub_expr),
      perturbation = pert,
      baseMean = baseline_expr
    )
  }
  
  # Run in parallel
  result_list <- mclapply(all_perturbations, calc_baseline, mc.cores = mc.cores)
  
  baseline_table <- do.call(rbind, result_list)
  rownames(baseline_table) <- NULL
  baseline_table <- map_ensg_to_symbol(baseline_table)
  message("Baseline expression table complete.")
  return(baseline_table)
}

base_exp.J <- get_baseline_table(seurat_obj.J, mc.cores = 16)
base_exp.K <- get_baseline_table(seurat_obj.K, mc.cores = 16)

write_csv(base_exp.J, here("data/perturb/QC/Jurkat/baseline_expression.csv"))
write_csv(base_exp.K, here("data/perturb/QC/K562_essential/baseline_expression.csv"))


plot_cell_QC_metrics <- function(seurat_object) {
  df <- seurat_object@meta.data
  
  # Gene detection rate
  raw_counts <- GetAssayData(seurat_object, slot = "counts")
  gene_detection_rate <- Matrix::rowMeans(raw_counts > 0)
  gene_min <- min(gene_detection_rate, na.rm = TRUE)
  
  p.gene <- ggplot(data.frame(rate = gene_detection_rate), aes(x = rate)) +
    geom_histogram(binwidth = 0.01, fill = "#009E73", color = "black", alpha = 0.7) +
    geom_vline(xintercept = gene_min, color = "red", linetype = "dashed") +
    annotate("text", x = gene_min, y = Inf, 
             label = paste0("min = ", signif(gene_min, 3)), 
             vjust = 2, hjust = -0.1, color = "red", size = 6) +
    theme_minimal() +
    ggtitle("Gene Detection Rate") +
    xlab("Fraction of Cells per Gene") +
    ylab("Gene Count")
  
  # UMI per cell
  umi_min <- min(df$UMI_count, na.rm = TRUE)
  
  p.UMI <- ggplot(df, aes(x = UMI_count)) +
    geom_histogram(bins = 100, fill = "#56B4E9", color = "black", alpha = 0.7) +
    scale_x_log10() +
    geom_vline(xintercept = umi_min, color = "red", linetype = "dashed") +
    annotate("text", x = umi_min, y = Inf, 
             label = paste0("min = ", signif(umi_min, 3)), 
             vjust = 2, hjust = -0.1, color = "red", size = 6) +
    theme_minimal() +
    ggtitle("Unique UMIs per Cell (library size)") +
    xlab("Unique UMI Count") +
    ylab("Cell Count")
  
  # Compute number of detected genes per cell (non-zero genes)
  raw_counts <- GetAssayData(seurat_object, slot = "counts")
  df$nFeature_RNA <- Matrix::colSums(raw_counts > 0)
  
  nFeature_min <- min(df$nFeature_RNA, na.rm = TRUE)
  
  p.nFeature <- ggplot(df, aes(x = nFeature_RNA)) +
    geom_histogram(fill = "#CCCCCC", color="black", bins = 100) +
    geom_vline(xintercept = nFeature_min, color = "red", linetype = "dashed") +
    annotate("text", x = nFeature_min, y = 0, 
             label = paste0("min = ", signif(nFeature_min, 3)), 
             vjust = -0.5, hjust = -0.1, color = "red", size = 6) +
    theme_minimal() +
    ggtitle("Detected Genes per Cell") +
    xlab("Number of Expressed Genes") +
    ylab("Cell Count")
  
  # Mitochondrial % per cell
  mito_max <- max(df$mitopercent, na.rm = TRUE)
  
  p.mito <- ggplot(df, aes(x = mitopercent)) +
    geom_histogram(binwidth = 0.005, fill = "#E69F00", color = "black", alpha = 0.7) +
    geom_vline(xintercept = mito_max, color = "red", linetype = "dashed") +
    annotate("text", x = mito_max, y = Inf, 
             label = paste0("max = ", signif(mito_max, 3)), 
             vjust = 2, hjust = 1.1, color = "red", size = 6) +
    theme_minimal() +
    ggtitle("Mitochondrial % per Cell") +
    xlab("Mitochondrial RNA %") +
    ylab("Cell Count")
  
  return(p.UMI + p.nFeature + p.mito + plot_layout(ncol = 3))
}

# Generate QC panels
p.J <- plot_cell_QC_metrics(seurat_obj.J)
p.K <- plot_cell_QC_metrics(seurat_obj.K)

p.J_wrapped <- wrap_elements(p.J + plot_annotation(title = "Jurkat-Essential"))
p.K_wrapped <- wrap_elements(p.K + plot_annotation(title = "K562-Essential"))

# Stack vertically and tag only the two panels A and B
final_plot <- (p.J_wrapped / p.K_wrapped) +
  plot_annotation(
    title = "Quality Control Metrics for Jurkat and K562 Perturb-seq",
    tag_levels = "A",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )

print("Saving QC plots")
ggsave(here("plots/exploratory/cell_QC.png"), plot = final_plot, width = 14, height = 8)

DEGs.J.bulk <- list.files(here("data/perturb/DEGs/bulk_Jurkat"), full.names = T)
DEGs.K.bulk <- list.files(here("data/perturb/DEGs/bulk_K562_essential"), full.names = T)
pairs.J.bulk <- extract_DEGs(DEGs.J.bulk, cores = 32, bulk=T)$perturbation_effect_df %>%
  select(perturbation, effect_ensg)
pairs.K.bulk <- extract_DEGs(DEGs.K.bulk, cores = 32, bulk = T)$perturbation_effect_df %>%
  select(perturbation, effect_ensg)

DEGs.J.sc <- list.files(here("data/perturb/DEGs/Jurkat"), full.names = T)
DEGs.K.sc <- list.files(here("data/perturb/DEGs/K562_essential"), full.names = T)
pairs.J.sc <- extract_DEGs(DEGs.J.sc, cores = 32, bulk=F)$perturbation_effect_df %>%
  select(perturbation, effect_ensg)
pairs.K.sc <- extract_DEGs(DEGs.K.sc, cores = 32, bulk = F)$perturbation_effect_df %>%
  select(perturbation, effect_ensg)

gene_QC_metrics <- function(seurat_object, pairs, expression_threshold = 0, min_cells = 7) {
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  
  expr_mat <- GetAssayData(seurat_object, slot = "data")
  meta <- seurat_object@meta.data
  
  # Filter pairs to valid ones
  valid_genes <- rownames(expr_mat)
  valid_perts <- unique(meta$gene)
  pairs <- pairs %>%
    filter(effect_ensg %in% valid_genes, perturbation %in% valid_perts)
  
  cat("Total valid pairs:", nrow(pairs), "\n")
  
  # Get control expression
  ctrl_cells <- rownames(meta)[meta$gene == "non-targeting"]
  ctrl_expr_mat <- expr_mat[, ctrl_cells, drop = FALSE]
  ctrl_counts <- Matrix::rowSums(ctrl_expr_mat > expression_threshold)
  ctrl_df <- tibble(gene = names(ctrl_counts), n_expr_ctrl = ctrl_counts)
  
  # Group treatment cells by perturbation
  pert_cell_index <- split(rownames(meta), meta$gene)
  pert_cell_index <- pert_cell_index[names(pert_cell_index) != "non-targeting"]
  
  trt_names <- names(pert_cell_index)
  total_perts <- length(trt_names)
  
  # Process each perturbation
  trt_counts_list <- mclapply(seq_along(trt_names), function(i) {
    pert <- trt_names[i]
    cat(sprintf("Processing %d / %d: %s\n", i, total_perts, pert))
    
    pert_cells <- pert_cell_index[[pert]]
    trt_expr_mat <- expr_mat[, pert_cells, drop = FALSE]
    counts <- Matrix::rowSums(trt_expr_mat > expression_threshold)
    
    tibble(gene = names(counts), perturbation = pert, n_expr_trt = counts)
  }, mc.cores = 32, mc.preschedule = FALSE)
  
  trt_counts_df <- bind_rows(trt_counts_list)
  
  # Join with pair list and control counts
  df <- pairs %>%
    rename(gene = effect_ensg) %>%
    left_join(trt_counts_df, by = c("perturbation", "gene")) %>%
    left_join(ctrl_df, by = "gene") %>%
    filter(!is.na(n_expr_trt), !is.na(n_expr_ctrl)) %>%
    mutate(QC_pass = ifelse(n_expr_trt >= min_cells & n_expr_ctrl >= min_cells, "Pass", "Fail"))
  
  return(df)
}

plot_gene_QC_metrics <- function(gene_QC, min_cells = 7, title = "Pairwise QC") {
  
  gene_QC <- gene_QC %>%
    mutate(
      QC_pass = ifelse(n_expr_trt >= min_cells & n_expr_ctrl >= min_cells, "Pass", "Fail"),
      is_high_logFC = abs(avg_log2FC) > 15
    )
  
  # Calculate % of failed QC pairs
  qc_fail_pct <- mean(gene_QC$QC_pass == "Fail") * 100
  qc_fail_label <- sprintf("Excluded by Pairwise QC: %.2f%%", qc_fail_pct)
  
  p <- ggplot(gene_QC, aes(x = n_expr_ctrl, y = n_expr_trt)) +
    geom_point(aes(color = abs(avg_log2FC)), alpha = 0.5, size = 1) +
    geom_point(
      data = filter(gene_QC, is_high_logFC),
      aes(color = abs(avg_log2FC)),
      shape = 21, stroke = 0.8, size = 2.2, fill = NA, color = "red"
    ) +
    scale_x_log10() +
    scale_y_log10() +
    geom_vline(xintercept = min_cells, linetype = "solid") +
    geom_hline(yintercept = min_cells, linetype = "solid") +
    labs(
      title = title,
      x = "Number of non-zero NT control cells",
      y = "Number of non-zero gene perturbation cells",
      color = "abs(logFC)"
    ) +
    annotate("text", 
             x = 10, 
             y = max(gene_QC$n_expr_trt, na.rm = TRUE),
             label = qc_fail_label, 
             hjust = 0, vjust = 1, size = 6, fontface = "italic") +
    theme_minimal()
  
  return(p)
}


gene_QC.J <- gene_QC_metrics(seurat_obj.J, pairs.J.bulk)
write_csv(gene_QC.J, here("plots/exploratory/J.B.gene_QC.csv"))
gene_QC.K <- gene_QC_metrics(seurat_obj.K, pairs.K.bulk)
write_csv(gene_QC.K, here("plots/exploratory/K.B.gene_QC.csv"))

gene_QC.J.S <- gene_QC_metrics(seurat_obj.J, pairs.J.sc)
write_csv(gene_QC.J.S, here("plots/exploratory/J.S.gene_QC.csv"))
gene_QC.K.S <- gene_QC_metrics(seurat_obj.K, pairs.K.sc)
write_csv(gene_QC.K.S, here("plots/exploratory/K.S.gene_QC.csv"))

gene_QC.J.B <- read_csv(here("plots/exploratory/J.B.gene_QC.csv"))
gene_QC.K.B <- read_csv(here("plots/exploratory/K.B.gene_QC.csv"))
gene_QC.J.S <- read_csv(here("plots/exploratory/J.S.gene_QC.csv"))
gene_QC.K.S <- read_csv(here("plots/exploratory/K.S.gene_QC.csv"))



# Generate QC panels
p.J.B <- plot_gene_QC_metrics(gene_QC.J.B, title = "Jurkat-Essential (DESeq2 LRT)")
p.K.B <- plot_gene_QC_metrics(gene_QC.K.B, title = "K562-Essential (DESeq2 LRT)")
p.J.S <- plot_gene_QC_metrics(gene_QC.J.S, title = "Jurkat-Essential (Seurat-Wilcox)")
p.K.S <- plot_gene_QC_metrics(gene_QC.K.S, title = "K562-Essential (Seurat-Wilcox)")

# Stack vertically and tag only the two panels A and B
final_plot <- wrap_elements(p.J.B + p.K.B) / wrap_elements(p.J.S + p.K.S) +
  plot_annotation(
    title = "Perturbation-expressed gene pair QC (differential expression FDR < 0.05)",
    tag_levels = "A",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )

print("Saving QC plots")
ggsave(here("plots/exploratory/gene_QC.png"), plot = final_plot, width = 14, height = 14)


DRGs.J.B <- read_csv(here("data/perturb/pairs/bulk_Jurkat/bulk_perturbation_pairs_eff.csv"))
DRGs.K.B <- read_csv(here("data/perturb/pairs/bulk_K562_essential/bulk_perturbation_pairs_eff.csv"))
DRGs.J.S <- read_csv(here("data/perturb/pairs/Jurkat/perturbation_pairs_eff.csv"))
DRGs.K.S <- read_csv(here("data/perturb/pairs/K562_essential/perturbation_pairs_eff.csv"))


map_ensg_to_symbol <- function(gene_QC) {
  gene.ensg <- gene_QC$gene %>% unique
  
  gene.symbols <- mapIds(org.Hs.eg.db, keys = gene.ensg,
                            column = "SYMBOL", keytype = "ENSEMBL",
                            multiVals = "first")
  
  gene.symbols <- gene.symbols %>% na.omit
  
  gene_QC <- gene_QC %>%
    mutate(effect = gene.symbols[gene])
  
  return(gene_QC)
}

gene_QC.J.B <- gene_QC.J.B %>% map_ensg_to_symbol %>% filter(perturbation != effect)
gene_QC.K.B <- gene_QC.K.B %>% map_ensg_to_symbol %>% filter(perturbation != effect)
gene_QC.J.S <- gene_QC.J.S %>% map_ensg_to_symbol %>% filter(perturbation != effect)
gene_QC.K.S <- gene_QC.K.S %>% map_ensg_to_symbol %>% filter(perturbation != effect)


gene_QC.J.B <- left_join(
  gene_QC.J.B,
  DRGs.J.B,
  by = c("perturbation", "effect")
) %>%
  select(perturbation, effect, n_expr_trt, n_expr_ctrl, avg_log2FC, p_val_adj)
gene_QC.K.B <- left_join(
  gene_QC.K.B,
  DRGs.K.B,
  by = c("perturbation", "effect")
) %>%
  select(perturbation, effect, n_expr_trt, n_expr_ctrl, avg_log2FC, p_val_adj)
gene_QC.J.S <- left_join(
  gene_QC.J.S,
  DRGs.J.S,
  by = c("perturbation", "effect")
) %>%
  select(perturbation, effect, n_expr_trt, n_expr_ctrl, avg_log2FC, p_val_adj)
gene_QC.K.S <- left_join(
  gene_QC.K.S,
  DRGs.K.S,
  by = c("perturbation", "effect")
) %>%
  select(perturbation, effect, n_expr_trt, n_expr_ctrl, avg_log2FC, p_val_adj)

plot_gene_QC_metrics <- function(gene_QC, title) {
  # Plot
  p <- ggplot(gene_QC, aes(x = n_expr_ctrl, y = n_expr_trt, color = abs(avg_log2FC))) +
    geom_point(alpha = 0.5, size = 1) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_gradient(
      low = "blue", high = "red",
      oob = squish,  # caps values outside limits
      limits = c(0, 10)  # cap color scale at 10
    ) +
    labs(
      title = title,
      x = "N nonzero ctrl. cells",
      y = "N nonzero trt. cells",
      color = "Average LogFC"
    ) +
    theme_minimal()
  return(p)
}

# Generate QC panels
p.J.B <- plot_gene_QC_metrics(gene_QC.J.B, title = "Jurkat-Essential (DESeq2 LRT)")
p.K.B <- plot_gene_QC_metrics(gene_QC.K.B, title = "K562-Essential (DESeq2 LRT)")
p.J.S <- plot_gene_QC_metrics(gene_QC.J.S, title = "Jurkat-Essential (Seurat Wilcox)")
p.K.S <- plot_gene_QC_metrics(gene_QC.K.S, title = "K562-Essential (Seurat Wilcox)")

# Stack vertically and tag only the two panels A and B
final_plot <- wrap_elements(p.J.B + p.K.B) / wrap_elements(p.J.S + p.K.S) +
  plot_annotation(
    title = "Perturbation-expressed gene pair QC (differential expression FDR < 0.05)",
    tag_levels = "A",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )

print("Saving QC plots")
ggsave(here("plots/exploratory/gene_QC.png"), plot = final_plot, width = 14, height = 14)


write_csv(gene_QC.J.B, here("data/perturb/QC/bulk_Jurkat/gene_QC.csv"))
write_csv(gene_QC.K.B, here("data/perturb/QC/bulk_K562_essential/gene_QC.csv"))
write_csv(gene_QC.J.S, here("data/perturb/QC/Jurkat/gene_QC.csv"))
write_csv(gene_QC.K.S, here("data/perturb/QC/K562_essential/gene_QC.csv"))


gene_QC_metrics.calibration <- function(seurat_object, pairs, expression_threshold = 0, min_cells = 7) {
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  
  expr_mat <- GetAssayData(seurat_object, slot = "data")
  meta <- seurat_object@meta.data %>%
    filter(gene == "non-targeting") %>%
    droplevels
  
  # Filter pairs to valid ones
  valid_genes <- rownames(expr_mat)
  valid_perts <- unique(meta$sgID_AB) %>% as.character
  pairs <- pairs %>%
    filter(effect_ensg %in% valid_genes, perturbation %in% valid_perts)
  
  cat("Total valid pairs:", nrow(pairs), "\n")
  
  # Get control expression
  ctrl_cells <- rownames(meta)[meta$gene == "non-targeting"]
  ctrl_expr_mat <- expr_mat[, ctrl_cells, drop = FALSE]
  ctrl_counts <- Matrix::rowSums(ctrl_expr_mat > expression_threshold)
  ctrl_df <- tibble(gene = names(ctrl_counts), n_expr_ctrl = ctrl_counts)
  
  # Group treatment cells by perturbation
  pert_cell_index <- split(rownames(meta), meta$sgID_AB)
  
  trt_names <- names(pert_cell_index)
  total_perts <- length(trt_names)
  
  # Process each perturbation
  trt_counts_list <- mclapply(seq_along(trt_names), function(i) {
    pert <- trt_names[i]
    cat(sprintf("Processing %d / %d: %s\n", i, total_perts, pert))
    
    pert_cells <- pert_cell_index[[pert]]
    trt_expr_mat <- expr_mat[, pert_cells, drop = FALSE]
    counts <- Matrix::rowSums(trt_expr_mat > expression_threshold)
    
    tibble(gene = names(counts), perturbation = pert, n_expr_trt = counts)
  }, mc.cores = 32, mc.preschedule = FALSE)
  
  trt_counts_df <- bind_rows(trt_counts_list)
  
  # Join with pair list and control counts
  df <- pairs %>%
    rename(gene = effect_ensg) %>%
    left_join(trt_counts_df, by = c("perturbation", "gene")) %>%
    left_join(ctrl_df, by = "gene") %>%
    filter(!is.na(n_expr_trt), !is.na(n_expr_ctrl)) %>%
    mutate(QC_pass = ifelse(n_expr_trt >= min_cells & n_expr_ctrl >= min_cells, "Pass", "Fail"))
  
  return(df)
}

DEGs.J.bulk <- list.files(here("data/perturb/calibration/bulk_Jurkat"), full.names = T)
DEGs.K.bulk <- list.files(here("data/perturb/calibration/bulk_K562_essential"), full.names = T)
pairs.J.bulk <- extract_DEGs(DEGs.J.bulk, cores = 32, bulk=T, significant=F, efficient = F)$perturbation_effect_df %>%
  select(perturbation, effect_ensg) %>%
  mutate(perturbation = sub("(^[^_]+_[^_]+)_(.+$)", "\\1|\\2", perturbation))
pairs.K.bulk <- extract_DEGs(DEGs.K.bulk, cores = 32, bulk = T, significant=F, efficient = F)$perturbation_effect_df %>%
  select(perturbation, effect_ensg) %>%
  mutate(perturbation = sub("(^[^_]+_[^_]+)_(.+$)", "\\1|\\2", perturbation))

DEGs.J.sc <- list.files(here("data/perturb/calibration/Jurkat"), full.names = T)
DEGs.K.sc <- list.files(here("data/perturb/calibration/K562_essential"), full.names = T)
pairs.J.sc <- extract_DEGs(DEGs.J.sc, cores = 32, bulk=F, significant=F, efficient = F)$perturbation_effect_df %>%
  select(perturbation, effect_ensg)
pairs.K.sc <- extract_DEGs(DEGs.K.sc, cores = 32, bulk = F, significant=F, efficient = F)$perturbation_effect_df %>%
  select(perturbation, effect_ensg)


gene_QC.J <- gene_QC_metrics.calibration(seurat_obj.J, pairs.J.bulk)
write_csv(gene_QC.J, here("data/perturb/QC/bulk_Jurkat/calibration_gene_QC.csv"))
gene_QC.K <- gene_QC_metrics.calibration(seurat_obj.K, pairs.K.bulk)
write_csv(gene_QC.K, here("data/perturb/QC/bulk_K562_essential/calibration_gene_QC.csv"))

gene_QC.J.S <- gene_QC_metrics.calibration(seurat_obj.J, pairs.J.sc)
write_csv(gene_QC.J.S, here("data/perturb/QC/Jurkat/calibration_gene_QC.csv"))
gene_QC.K.S <- gene_QC_metrics.calibration(seurat_obj.K, pairs.K.sc)
write_csv(gene_QC.K.S, here("data/perturb/QC/K562_essential/calibration_gene_QC.csv"))