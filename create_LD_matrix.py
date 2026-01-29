import os

# Increase Java Heap Size for PySpark
os.environ["PYSPARK_SUBMIT_ARGS"] = (
    "--driver-memory 6g "
    "--executor-memory 32g "
    "--conf spark.executor.memoryOverhead=10g "
    "pyspark-shell"
)

import polars as pl
import argparse
import hail as hl


def create_ld_matrix(gene_name, lead_snp, chrom, start_pos, end_pos):
    # Configure Hail for S3 access
    hl.init(
        spark_conf={
            "spark.jars.packages": "org.apache.hadoop:hadoop-aws:3.3.4",
            "spark.hadoop.fs.s3a.impl": "org.apache.hadoop.fs.s3a.S3AFileSystem",
            # Use ~/.aws/credentials
            "spark.hadoop.fs.s3a.aws.credentials.provider": "com.amazonaws.auth.DefaultAWSCredentialsProviderChain",
            # REQUIRED for pan-ukb
            "spark.hadoop.fs.s3a.requester.pays.enabled": "true",
        },
        idempotent=True,  # Run init only once
    )

    # Read variant indices table
    ht_idx = hl.read_table(
        "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht"
    )

    # Load LD matrix
    bm = hl.linalg.BlockMatrix.read(
        "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm"
    )

    # Create interval for region
    interval = hl.parse_locus_interval(f"{chrom}:{start_pos}-{end_pos}")

    # Filter variants in the region
    ht_idx_region = ht_idx.filter(interval.contains(ht_idx.locus))

    ht_idx_region.select(
        alleleA=ht_idx_region.alleles[0],
        alleleB=ht_idx_region.alleles[1],
        AF=ht_idx_region.AF,
        idx=ht_idx_region.idx,
        rsid=ht_idx_region.rsid,
    ).export(f"{gene_name}_SNPs.tsv")

    idx = ht_idx_region.idx.collect()

    # Filter LD matrix to region
    bm = bm.filter(idx, idx)

    # Export filtered LD matrix to flat file
    bm.write(
        f"/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/perturb_exploratory/data/LD_loci/{gene_name}_{lead_snp}",
        force_row_major=True,
    )
    hl.linalg.BlockMatrix.export(
        f"/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/perturb_exploratory/data/LD_loci/{gene_name}_{lead_snp}",
        f"/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/perturb_exploratory/data/LD_loci/{gene_name}_{lead_snp}.bgz",
        delimiter=" ",
    )


parser = argparse.ArgumentParser(
    description="Generate LD matrix for locus around lead cis-eSNP"
)

parser.add_argument("--gene", type=str, required=True, help="Gene name")
parser.add_argument(
    "--window",
    type=int,
    required=False,
    help="Window size around the lead cis-eSNP (in Kbs)",
    default=1000,
)

args = parser.parse_args()

# Example args for testing
# args = argparse.Namespace(gene="PXK", window=1000)

args = argparse.Namespace(gene="HBS1L", window=1000)

# Read cis-eQTL data for the specified gene

selected_cols = [
    "SNP",
    "GeneSymbol",
    "Zscore",
    "Pvalue",
    "FDR",
    "SNPPos",
    "SNPChr",
    "GeneChr",
    "GenePos",
    "AssessedAllele",
    "OtherAllele",
    "NrSamples",
]

cis_eQTL = (
    pl.scan_csv(
        f"/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/cis_eQTLs_eQTLgen.tsv",
        separator="\t",
    )
    .select(selected_cols)
    .filter(pl.col("GeneSymbol") == args.gene)
    .collect()
)

# Save cis-eQTL summary stats for the gene
cis_eQTL.write_csv(
    f"/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/perturb_exploratory/data/LD_loci/{args.gene}_cis_eQTLs.tsv",
)

lead_snp, lead_snp_chr, lead_snp_pos = (
    cis_eQTL.sort("Pvalue").select("SNP", "SNPChr", "SNPPos").row(0)
)

# Read LD matrix for the locus around the lead SNP
window_size = args.window * 1000
start_pos = lead_snp_pos - window_size
end_pos = lead_snp_pos + window_size

create_ld_matrix(args.gene, lead_snp, lead_snp_chr, start_pos, end_pos)
print(
    f"LD matrix for {args.gene} around lead SNP {lead_snp} ({lead_snp_chr}:{lead_snp_pos}) created."
)
