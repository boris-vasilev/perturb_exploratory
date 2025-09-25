# Sign concordance between two Perturb-seq datasets (e.g. K562-Essential and K562-GenomeWide)

library(tidyverse)
library(patchwork)
library(glue)
library(igraph)
library(ggraph)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--cell1", type = "character", help = "Cell type/line 1", required = TRUE)
parser$add_argument("--cell2", type = "character", help = "Cell type/line 2", required = TRUE)

args <- parser$parse_args()
cell1 <- args$cell1
cell2 <- args$cell2

plot_logFC_correlation.cells <- function(cell1, cell2) {
  perturb_dat1 <- read_csv(glue("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/{cell1}/perturbation_pairs_eff.csv"))
  perturb_dat2 <- read_csv(glue("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/{cell2}/perturbation_pairs_eff.csv"))

  perturb_dat <- perturb_dat1 %>% inner_join(perturb_dat2, by=c("perturbation", "effect"))
  correlation <- cor(perturb_dat$avg_log2FC.x, perturb_dat$avg_log2FC.y, method = "pearson")

  p <- ggplot(perturb_dat, aes(x = avg_log2FC.x, y = avg_log2FC.y)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "red") +
  labs(
    x = glue("avg_log2FC ({cell1} cell line)"),
    y = glue("avg_log2FC ({cell2} cell line)"),
    title = paste0("Correlation = ", round(correlation, 3))
  ) +
  ggtitle(glue("logFC correlation {cell1} vs {cell2} (R={round(correlation, 3)})")) +
  theme_minimal()

  return(p)
}

correlation <- cor(no_degs$no_DEGs.x, no_degs$no_DEGs.y, method = "pearson")
p <- ggplot(no_degs, aes(x = no_DEGs.x, y = no_DEGs.y)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "red") +
  labs(
    x = glue("Number of DEGs (K562-E)"),
    y = glue("Number of DEGs (K562-GW)"),
    title = paste0("Correlation = ", round(correlation, 3))
  ) +
  ggtitle(glue("Number of DEGs K562-E vs K562-GW (R={round(correlation, 3)})")) +
  theme_minimal()

cell_lines <- c("Jurkat", "K562_essential", "RPE1", "HepG2")

# Generate all pairwise combinations
pairwise_plots <- combn(cell_lines, 2, function(cells) {
  plot_logFC_correlation.cells(cells[1], cells[2])
}, simplify = FALSE)

# Combine with patchwork
combined_plot <- wrap_plots(pairwise_plots, ncol = 2)

plot_sign_concordance_graph.cells <- function(cell1, cell2) {

  perturb_dat1 <- read_csv(glue("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/{cell1}/perturbation_pairs_eff.csv"))
  perturb_dat2 <- read_csv(glue("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/{cell2}/perturbation_pairs_eff.csv"))

  perturb_dat <- perturb_dat1 %>% inner_join(perturb_dat2, by=c("perturbation", "effect"))
  # Add column for sign comparison
  perturb_dat <- perturb_dat %>% 
    mutate(Sign = ifelse(sign(avg_log2FC.x) == sign(avg_log2FC.y), "Same effect", "Opposite effect"))
  
  
  # Create edge list for igraph
  edges <- perturb_dat %>%
    select(from = perturbation, to = effect, Sign)
  
  # Create graph
  g <- graph_from_data_frame(edges, directed = TRUE)
  
  # Plot using ggraph
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(color = Sign), arrow = arrow(length = unit(3, "mm")), end_cap = circle(3, 'mm')) +
    geom_node_point(size = 4) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_color_manual(values = c("Same effect" = "red", "Opposite effect" = "black")) +
    theme_void() +
    ggtitle(glue("{cell1} vs {cell2}")) +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  
  return(p)
  
}


plot_sign_concordance_graph.cells.all <- function() {

  cell_lines <- c(
    "Jurkat",
    "K562_essential",
    "RPE1",
    "HepG2"
  )

  perturb_dats <- lapply(cell_lines, function(cell_line){
    read_csv(glue("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/{cell_line}/perturbation_pairs_eff.csv")) %>%
      mutate(Sign = sign(avg_log2FC))
  })
  
  edges <- lapply(perturb_dats, function(dat) {
    dat %>% select(from = perturbation, to = effect, Sign)
  })
  
  # Find intersection of edges across all cell lines
  common_edges <- reduce(edges, function(x, y) inner_join(x, y, by = c("from", "to"))) %>%
    select(from, to) %>% distinct()
  
  # Restrict each edge set to only common edges
  edges_common <- lapply(edges, function(ed) {
    ed %>% inner_join(common_edges, by = c("from", "to"))
  })
  
  # Create graphs
  graphs <- lapply(edges_common, function(graph_edges) {
    graph_from_data_frame(graph_edges, directed = TRUE)
  })
  
  # Generate plots
  plots <- lapply(seq_along(graphs), function(i) {
    g <- graphs[[i]]
    ggraph(g, layout = "fr") +
      geom_edge_link(aes(color = as.factor(Sign)), 
                     arrow = arrow(length = unit(3, "mm")), 
                     end_cap = circle(3, 'mm')) +
      geom_node_point(size = 3) +
      geom_node_text(aes(label = name), repel = TRUE, size = 3) +
      ggtitle(cell_lines[i]) +
      theme_void() +
      scale_edge_color_manual(values = c("-1" = "red", "0" = "grey70", "1" = "green"))
  })
  
  return(wrap_plots(plots, ncol = 2))
  
}