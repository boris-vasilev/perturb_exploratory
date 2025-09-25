library(tidyverse)
library(glue)
library(argparse)
library(data.table)
library(parallel)

parser <- ArgumentParser()
parser$add_argument("--cells", type = "character", help = "Cell type/line", required = TRUE)

args <- parser$parse_args()
cells <- args$cells

DEG_files <- list.files(glue("/rds/project/rds-csoP2nj6Y6Y/biv22/data/perturb/{cells}"), full.names = TRUE)

n_cores <- 32

DEGs <- mclapply(DEG_files, function(f) {
  fread(f) %>% select(gene, log2FoldChange)
}, mc.cores = n_cores)

names(DEGs) <- sub("\\.csv$", "", basename(DEG_files))
logFC_dt <- rbindlist(DEGs, idcol = "perturbation")
fwrite(logFC_dt, glue("/rds/project/rds-csoP2nj6Y6Y/biv22/data/pairs/{cells}/logFC_dat.csv"))