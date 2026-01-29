library(tidyverse)
library(data.table)
library(glue)
library(argparse)
library(susieR)
library(ggrepel)
library(here)

parser <- ArgumentParser()
parser$add_argument("--gene", type = "character", help = "Gene name", required = TRUE)

args <- parser$parse_args()
gene <- args$gene


LD_loci_dir <- here("data/LD_loci")

ld <- fread(Sys.glob(glue("{LD_loci_dir}/{gene}_*.bgz")), header = FALSE) %>% as.matrix

snps <- fread(glue("{LD_loci_dir}/{gene}_SNPs.tsv"))
# Extract chromosome from lead SNP locus
chrom <- strsplit(snps$locus[1], ":") %>% unlist %>% first

# Read phenotype summary from INTERVAL to map ENSG to gene name
phenotypes <- fread("/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/eqtl/INTERVAL_eQTL_summary_statistics/INTERVAL_eQTL_phenotype_summary.tsv")

cis_ensg <- phenotypes[gene_name == gene, phenotype_id]

# Obtain INTERVAL summary stats for the same chromosome
cis_eqtl <- fread(glue("/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/eqtl/INTERVAL_eQTL_summary_statistics/INTERVAL_eQTL_nominal_chr{chrom}.tsv")) %>%
  filter(phenotype_id == cis_ensg)


### WARNING:
# INTERVAL eQTL effect allele assignment seems to be flipped compared to UKB and OpenTargets
# effect_allele is NOT the minor allele as one might expect
# alleleB in the snps file we derive in create_LD_matrices.R is the minor (effect)allele
# matching OpenTargets. INTERVAL has those flipped and explicitly called effect/reference
# which doesn't match the minor/major allele definition.
# Therefore majority of slopes need to be flipped and af adjusted accordingly.
align_eqtl <- function(cis_eqtl, snps) {
  merged <- cis_eqtl %>% inner_join(snps, by = c("variant_id" = "rsid")) %>%
    mutate(
      slope = case_when(
          effect_allele == alleleB ~ slope,
          effect_allele == alleleA ~ -slope,
          TRUE ~ NA_real_
    ),
    af = case_when(
          effect_allele == alleleB ~ af,
          effect_allele == alleleA ~ 1 - af,
          TRUE ~ NA_real_
    )
  ) %>% arrange(idx) %>%
  filter(!is.na(slope))
}

cis_eqtl_aligned <- align_eqtl(cis_eqtl, snps)

ld_idx <- which(snps$idx %in% cis_eqtl_aligned$idx)

LD_sub <- ld[ld_idx, ld_idx]

susie_fit <- susie_rss(
  bhat = cis_eqtl_aligned$slope,
  shat = cis_eqtl_aligned$slope_se,
  maf = cis_eqtl_aligned$af,
  R = LD_sub,
  n = 4732,
  # L = 10,
  unmappable_effects = "inf",
  convergence_method = "pip"
)

# # HBS1L
# lead_snp <- "rs12526055"
# selected_snp <- "rs9399135"
# # PXK
# lead_snp <- "rs7633553"
# selected_snp <- "rs9311676"

# # Find index of selected SNP in cis_eqtl_aligned
# lead_snp_idx <- match(lead_snp, cis_eqtl_aligned$SNP)
# selected_snp_idx <- match(selected_snp, cis_eqtl_aligned$SNP)

# print(paste("Lead SNP PIP:", susie_fit$pip[lead_snp_idx]))
# print(paste("Selected SNP PIP:", susie_fit$pip[selected_snp_idx]))

# susie_plot(susie_fit, y = "PIP", add_legend = TRUE)


eqtl_susie <- cbind(cis_eqtl_aligned, summary(susie_fit)$vars %>% arrange(variable))

manhattan_plot <- ggplot(eqtl_susie, aes(x=pos_b37, y=-log10(pval_nominal))) +
  geom_point(alpha=0.8, size=1.3) +
  geom_point(data=subset(eqtl_susie, eqtl_susie$cs!=-1), aes(color=factor(cs)), size=4) +
  geom_label_repel( data=subset(eqtl_susie, eqtl_susie$cs!=-1), aes(label=variant_id), size=4, max.overlaps = 20) +
  labs(title=gene, color="Credible Set")

saveRDS(susie_fit, glue(here("finemap_fit/∫{gene}_INTERVAL_finemap_susieR.rds")))
ggsave(here(glue("plots/{gene}_INTERVAL_finemap_manhattan.png")), manhattan_plot, width=10, height=6)