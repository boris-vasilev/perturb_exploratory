library(tidyverse)
library(here)
library(glue)
library(igraph)
library(ggraph)


jurkat <- read_csv(here("data/perturb/pairs/bulk_Jurkat/bulk_gaussian_ratio_dat_eff0.05.csv"))
k562 <- read_csv(here("data/perturb/pairs/bulk_K562_essential/bulk_gaussian_ratio_dat_eff0.05.csv"))

jurkat.table <- table(sign(jurkat$x), sign(jurkat$y))
k562.table <- table(sign(k562$x), sign(k562$y))


plot_sign_concordance_graph <- function(dataset) {
  df <- if(dataset == "Jurkat") jurkat else k562
  
  # Add column for sign comparison
  df <- df %>% 
    mutate(Sign = ifelse(sign(x) == sign(y), "Same (discordant effect)", "Opposite (concordant effect)"))
  
  
  # Create edge list for igraph
  edges <- df %>%
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
    ggtitle(glue("{dataset}")) +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  
  return(p)
  
}

plot_sign_concordance_graph.multi <- function() {
  p.J <- plot_sign_concordance_graph("Jurkat")
  p.K <- plot_sign_concordance_graph("K562")
  
  combined_plot <- p.J + p.K +
    plot_annotation(tag_levels = 'A',
                    title="Sign Concordance Networks",
                    theme = theme(
                      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
                    )) +
    plot_layout(guides = "collect") 
  ggsave(
    filename = here(glue("plots/genes/sign_networks.png")),
    plot = combined_plot,
    width = 16, height = 8, dpi = 300
  )
}
