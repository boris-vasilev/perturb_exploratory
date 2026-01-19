library(here)
library(data.table)
library(tidyverse)
library(ComplexUpset)
library(LDlinkR)


jurkat.eQTL <- fread(here("data/perturb/pairs/bulk_Jurkat/bulk_eQTL_pairs_eff.csv"))
jurkat.perturb <- fread(here("data/perturb/pairs/bulk_Jurkat/bulk_perturbation_pairs_eff.csv"))

jurkat.eQTL %>% filter(FDR.effect < 0.05) %>%
  group_by(SNP, perturbation) %>%
  summarise(n_trans_genes = n_distinct(effect), .groups = "drop") %>%
  arrange(desc(n_trans_genes))

eqtl <- jurkat.eQTL %>% filter(FDR.effect < 0.05) %>%
  group_by(perturbation) %>%
  summarise(n_trans_genes = n_distinct(effect), trans_genes=list(unique(effect)), .groups = "drop") %>%
  arrange(desc(n_trans_genes))

perturb <- jurkat.perturb %>% filter(p_val_adj < 0.05) %>%
  group_by(perturbation) %>%
  summarise(n_DEGs = n_distinct(effect), DEGs=list(unique(effect)), .groups = "drop") %>%
  arrange(desc(n_DEGs))

merged <- eqtl %>% inner_join(perturb, by="perturbation") %>%
  rowwise() %>% 
  mutate(overlap = list(intersect(trans_genes, DEGs)),
         n_overlap = length(overlap),
         perc_overlap = n_overlap/n_trans_genes) %>%
  select(perturbation, n_trans_genes, n_DEGs, n_overlap, overlap, perc_overlap)

# Step 1: Create a summary of trans gene counts per SNP
trans_counts <- jurkat.eQTL %>%
  filter(perturbation == "NAA25", FDR.effect < 0.05) %>%
  group_by(SNP) %>%
  summarise(n_trans_genes = n_distinct(effect), .groups = "drop")

# Step 2: Create the cis p-value table
cis_pvals <- jurkat.eQTL %>%
  filter(perturbation == "NAA25") %>%
  select(SNP, Pvalue.perturb, FDR.perturb) %>%
  distinct() %>%
  arrange(Pvalue.perturb) %>%
  rename(cis_pval = Pvalue.perturb,
         cis_fdr = FDR.perturb)

# Step 3: Join both tables
final <- cis_pvals %>%
  left_join(trans_counts, by = "SNP") %>%
  arrange(cis_pval) %>%
  replace_na(list(n_trans_genes = 0)) %>%
  mutate(
    cis_pval = signif(cis_pval, 3),
    cis_fdr = signif(cis_fdr, 3)
  )

trans_sig <- jurkat.eQTL %>%
  filter(perturbation == "NAA25", FDR.effect < 0.05)

# Get list of trans-eGenes per SNP
snp_to_trans <- trans_sig %>%
  group_by(SNP) %>%
  summarise(trans_genes = list(unique(effect)), .groups = "drop")

# Step 1: Convert tibble to named list
snp_sets <- setNames(snp_to_trans$trans_genes, snp_to_trans$SNP)
snps <- names(snp_sets)

# Step 2: Compute % of each SNP's genes found in other SNPs
percent_shared <- sapply(snps, function(snp) {
  my_genes <- snp_sets[[snp]]
  other_genes <- unique(unlist(snp_sets[setdiff(snps, snp)]))
  shared_genes <- length(intersect(my_genes, other_genes))
  percent <- shared_genes / length(my_genes) * 100
  return(round(percent, 1))
})

# Step 3: Make a summary table
shared_summary <- data.frame(
  SNP = snps,
  Percent_shared = percent_shared,
  Total_trans_genes = sapply(snp_sets, length)
)

# Step 4: View
print(shared_summary)

snp_list <- c("rs10774610", "rs10774624", "rs10849915", "rs1265564",  "rs3184504",  "rs3809272",  "rs7310615", "rs3026445", "rs4766428")

# Get LD matrix for EUR population using r²
ld_matrix <- LDmatrix(snp_list, pop = "EUR", r2d = "r2", token = LDlink_token)

rownames(ld_matrix) <- ld_matrix$RS_number
ld_mat <- as.matrix(ld_matrix[, -1])

pheatmap(ld_mat,
         cluster_rows = T,
         cluster_cols = T,
         display_numbers = TRUE,
         number_format = "%.2f",
         main = "Pairwise LD (r²) between NAA25 cis-eSNPs")


# Overlap of trans-effects of all cis-eSNPs with Perturb-seq
trans_genes <- snp_to_trans$trans_genes %>% unlist %>% unique

perturb_DEGs <- jurkat.perturb %>% filter(perturbation == "NAA25") %>% pull(effect)

which(perturb_DEGs %in% trans_genes) %>% length()
length(perturb_DEGs)

cis_eQTL <- fread(here("data/summary_stats/cis_eQTLs_eQTLgen.tsv"),
                  sep = "\t")