library(tidyverse)
library(patchwork)
library(glue)
library(igraph)
library(ggraph)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE)

args <- parser$parse_args()
cells <- args$cells


# sign_table <- table(sign(perturb_dat$x), sign(perturb_dat$y))

plot_sign_concordance_graph <- function(cells) {

  perturb_dat <- read_csv(glue("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/{cells}/gaussian_ratio_dat_eff.csv"))

  # Add column for sign comparison
  perturb_dat <- perturb_dat %>% 
    mutate(Sign = ifelse(sign(x) == sign(y), "Same (discordant effect)", "Opposite (concordant effect)"))
  
  
  # Create edge list for igraph
  edges <- perturb_dat %>%
    select(from = perturb, to = effect, Sign)
  
  # Create graph
  g <- graph_from_data_frame(edges, directed = TRUE)
  
  # Plot using ggraph
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(color = Sign), arrow = arrow(length = unit(3, "mm")), end_cap = circle(3, 'mm')) +
    geom_node_point(size = 4) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_color_manual(values = c("Same (discordant effect)" = "black", "Opposite (concordant effect)" = "red")) +
    theme_void() +
    ggtitle(glue("{cells}")) +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  
  return(p)
  
}

plot_sign_concordance_graph.multi <- function(cells1, cells2) {
  p1 <- plot_sign_concordance_graph(cells1)
  p2 <- plot_sign_concordance_graph(cells2)
  
  combined_plot <- p1 + p2 +
    plot_annotation(tag_levels = 'A',
                    title="Sign Concordance Networks",
                    theme = theme(
                      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
                    )) +
    plot_layout(guides = "collect") 
  ggsave(
    filename = glue("/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_exploratory/plots/sign_networks_{cells1}_{cells2}.png"),
    plot = combined_plot,
    width = 16, height = 8, dpi = 300
  )
}
