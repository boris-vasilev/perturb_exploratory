library(here)
library(data.table)
library(tidyverse)
library(ComplexUpset)
library(LDlinkR)


# Selecting only lead cis-eSNP
jurkat.eQTL <- fread(here("data/perturb/pairs/bulk_Jurkat/bulk_eQTL_pairs_eff.csv")) %>%
  filter(FDR.perturb < 0.05) %>%
  group_by(perturbation) %>%
  filter(Pvalue.perturb == min(Pvalue.perturb)) %>%
  filter(abs(Beta.perturb) == max(abs(Beta.perturb))) %>%
  ungroup()

jurkat.perturb <- fread(here("data/perturb/pairs/bulk_Jurkat/bulk_perturbation_pairs_eff.csv"))

jurkat.eQTL %>%
  filter(FDR.effect < 0.05) %>%
  group_by(SNP, perturbation) %>%
  summarise(n_trans_genes = n_distinct(effect), .groups = "drop") %>%
  arrange(desc(n_trans_genes)) %>% rename(`cis-eGene`="perturbation")

trans_sig <- jurkat.eQTL %>%
  filter(FDR.effect < 0.05)

snp_to_trans <- trans_sig %>%
  group_by(perturbation) %>%
  summarise(trans_genes = list(unique(effect)), .groups = "drop",
            n_trans_genes = uniqueN(effect)) %>%
  arrange(desc(n_trans_genes))

perturb_to_effect <- jurkat.perturb %>% filter(p_val_adj < 0.05) %>%
  group_by(perturbation) %>%
  summarise(DEGs = list(unique(effect)), .groups = "drop",
            n_DEGs = uniqueN(effect)) %>%
  arrange(desc(n_DEGs))

merged <- snp_to_trans %>% inner_join(perturb_to_effect, by="perturbation")

merged.overlap <- merged %>%
  rowwise() %>% 
  mutate(overlap = list(intersect(trans_genes, DEGs)),
         n_overlap = length(overlap)) %>%
  select(perturbation, n_trans_genes, n_DEGs, n_overlap, overlap)


merged %>% filter(perturbation == "SMG5")