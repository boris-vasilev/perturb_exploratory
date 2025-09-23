# SCAP has a significant slope in both K562 and Jurkat (logistic VIVS)
library(tidyverse)
library(here)
library(patchwork)
library(glue)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(igraph)
library(ggraph)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("intersect", "base")
conflicted::conflicts_prefer(base::`%in%`)



source(here("src/bayesian/functions_plot_patchwork.R"))

jurkat <- read_csv(here("data/perturb/pairs/bulk_Jurkat/bulk_logit_dat_eff.csv"))
k562 <- read_csv(here("data/perturb/pairs/bulk_K562_essential/bulk_logit_dat_eff.csv"))

get_shared_perturb_edges <- function(gene_name) {
  gene.jurkat <- jurkat %>%
    filter(perturb == gene_name, x == 1) %>%
    transmute(pair = paste(perturb, effect, sep = "___"))
  
  gene.k562 <- k562 %>%
    filter(perturb == gene_name, x == 1) %>%
    transmute(pair = paste(perturb, effect, sep = "___"))
  
  shared <- intersect(gene.jurkat$pair, gene.k562$pair)
  return(shared)
}


plot_gene_eQTL_perturb_network.single <- function(gene_name, dataset, shared_edges) {
  df <- if(dataset == "bulk_Jurkat") jurkat else k562
  title <- if(dataset == "bulk_Jurkat") "Jurkat-essential" else "K562-essential"
  
  gene.df <- df %>%
    filter(perturb == gene_name,
           (x == 1 | y == 1)) %>%
    mutate(Effect = case_when(
      x == 1 & y == 1 ~ "Perturbation and trans-eQTL effect",
      x == 1 & y == 0 ~ "Only perturbation effect",
      x == 0 & y == 1 ~ "Only trans-eQTL effect"
    ),
    pair = paste(perturb, effect, sep = "___"),
    is_shared = pair %in% shared_edges,
    edge_linetype = if_else(is_shared, "dotted", "solid"))
  
  edges <- gene.df %>%
    select(from = perturb, to = effect, Effect, edge_linetype)
  
  g <- graph_from_data_frame(edges, directed = TRUE)
  
  ggraph(g, layout = "kk") +
    geom_edge_link(
      aes(color = Effect, linetype = edge_linetype),
      arrow = arrow(length = unit(4, 'mm')),
      end_cap = circle(3, 'mm'),
      lineend = "round"
    ) +
    scale_edge_linetype_identity() +
    scale_edge_width(range = c(0.5, 1.5), guide = "none") +
    geom_node_point(size = 5, color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(color = "Effect", title=title)
}

plot_gene_eQTL_perturb_network <- function(gene_name) {
  shared_edges <- get_shared_perturb_edges(gene_name)
  
  p.J <- plot_gene_eQTL_perturb_network.single(gene_name, "bulk_Jurkat", shared_edges)
  p.K <- plot_gene_eQTL_perturb_network.single(gene_name, "bulk_K562_essential", shared_edges)
  
  combined_plot <- p.J + p.K +
    plot_layout(guides = "collect") +
    plot_annotation(title = glue("{gene_name} regulatory network - eQTL and Perturb-seq evidence")) &
    theme(legend.position = "bottom")
  
  # ggsave(
  #   filename = here(glue("plots/genes/reg_network_{gene_name}.png")),
  #   plot = combined_plot,
  #   width = 12, height = 6, dpi = 300
  # )
  
  return(combined_plot)
}

reactome_enrich_shared_effects <- function(gene_name) {
# Get genes from both Jurkat and K562 with shared effects
  genes.jurkat <- jurkat %>%
    filter(perturb == gene_name, x == 1) %>%
    pull(effect)

  genes.k562 <- k562 %>%
    filter(perturb == gene_name, x == 1) %>%
    pull(effect)

  # Union of both effects plus the original gene
  genelist <- unique(c(gene_name, intersect(genes.jurkat, genes.k562)))

  # Convert gene symbols to Entrez IDs
  entrez_df <- bitr(genelist, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
  
  # Reactome enrichment
  reactome_res <- enrichPathway(gene = entrez_df$ENTREZID,
                                organism = "human",
                                pAdjustMethod = "fdr",
                                readable = TRUE)

  # Add a `highlight` flag and a `label` column to use in the plot
  reactome_res@result <- reactome_res@result %>%
    mutate(
      highlight = grepl(gene_name, geneID),
      label = ifelse(highlight, "*", "")
    )
  
  # Generate the dotplot
  p <- dotplot(reactome_res, showCategory = 20)
  
  # Add asterisks (*) next to terms containing the gene
  # Merge the result into the plotting data
  highlight_df <- reactome_res@result %>%
    filter(ID %in% p$data$ID) %>%
    select(ID, Description, label)
  
  # Join to original plot data to match coordinates
  plot_data <- p$data %>%
    left_join(highlight_df, by = c("ID", "Description"))
  
  # Add asterisk labels on the geneset
  p <- p + geom_text(
    data = plot_data %>% filter(label.x != ""),
    aes(x = GeneRatio, y = Description, label = label.x),
    inherit.aes = FALSE,
    color = "black", size = 10,
    vjust = 0.1  # places label above the dot
  )+ labs(
    title=glue("Reactome pathway enrichment of {gene_name} perturbation DEGs shared between K562 and Jurkat")
  ) + theme(
    plot.title = element_text(hjust = 1, size = 12, face = "bold", margin = margin(b = 10))
  )
  
  return(p)
  
  ggsave(
    filename = here(glue("plots/genes/reactome_{gene_name}.png")),
    plot = p,
    width = 6, height = 6, dpi = 300
  )
}


plot_gene_network_and_reactome <- function(gene_name) {
  networks <- plot_gene_eQTL_perturb_network(gene_name)
  reactome <- reactome_enrich_shared_effects(gene_name)
  
  combined_plot <- (networks / reactome) + plot_annotation(tag_levels = 'A')
  
  
  ggsave(
    filename = here(glue("plots/genes/{gene_name}.png")),
    plot = combined_plot,
    width = 12, height = 9, dpi = 300
  )
}







