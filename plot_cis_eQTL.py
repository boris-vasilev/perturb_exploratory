import polars as pl
import argparse
import numpy as np
import matplotlib.pyplot as plt
import hail as hl

# from hail.linalg import BlockMatrix


# Configure Hail for S3 access
hl.init(
    spark_conf={
        "spark.jars.packages": "org.apache.hadoop:hadoop-aws:3.3.4",
        "spark.hadoop.fs.s3a.impl": "org.apache.hadoop.fs.s3a.S3AFileSystem",
        # Use ~/.aws/credentials
        "spark.hadoop.fs.s3a.aws.credentials.provider": "com.amazonaws.auth.DefaultAWSCredentialsProviderChain",
        # REQUIRED for pan-ukb
        "spark.hadoop.fs.s3a.requester.pays.enabled": "true",
    }
)

# def read_ld_matrix(chrom, start_pos, end_pos, cis_eQTL_snps):
#     # Read variant indices table
#     ht_idx = hl.read_table(
#         "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht"
#     )

#     # Load LD matrix
#     bm = BlockMatrix.read("s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm")

#     # Create interval for region
#     interval = hl.parse_locus_interval(f"{chrom}:{start_pos}-{end_pos}")

#     # Filter variants in the region
#     ht_idx_region = ht_idx.filter(interval.contains(ht_idx.locus))

#     # Filter to only variants tested for cis-eQTLs for the gene of interest
#     rsids_eqtl = cis_eQTL_snps.to_pandas()
#     ht_rsids_eqtl = hl.Table.from_pandas(rsids_eqtl)
#     ht_rsids_eqtl = ht_rsids_eqtl.rename({'SNP': 'rsid'}).key_by('rsid')

#     ht_idx_by_rsid = ht_idx_region.key_by('rsid')
#     ht_idx_region_tested = ht_idx_by_rsid.join(ht_rsids_eqtl, how='inner')
#     ht_idx_region_tested = ht_idx_region_tested.key_by('locus', 'alleles')

#     idx = ht_idx_region_tested.idx.collect()

#     # Filter LD matrix to region
#     bm = bm.filter(idx, idx)

#     # Export filtered LD matrix to flat file
#     bm.write('/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/perturb_exploratory/data/LD_loci/CSNK2B', force_row_major=True)
#     BlockMatrix.export(
#         '/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/perturb_exploratory/data/LD_loci/CSNK2B',
#         '/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/perturb_exploratory/data/LD_loci/CSNK2B.bgz',
#         delimiter=' '
#     )


def check_snps_in_ld_matrix(chrom, start_pos, end_pos, snps):
    # Read variant indices table
    ht_idx = hl.read_table(
        "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht"
    )

    # Create interval for region
    interval = hl.parse_locus_interval(f"{chrom}:{start_pos}-{end_pos}")

    # Filter variants in the region
    ht_idx_region = ht_idx.filter(interval.contains(ht_idx.locus))

    # eQTL SNP set -> broadcast
    rsid_set = hl.literal(set(snps.to_list()))

    # Filter LD index only by rsIDs of interest
    ht_idx_subset = ht_idx_region.filter(rsid_set.contains(ht_idx_region.rsid))

    # Return list of rsIDs present in LD matrix
    return ht_idx_subset.rsid.collect()


def get_r2_between_snps(snp1, snp2):
    # Read variant indices table
    ht_idx = hl.read_table(
        "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht"
    )

    # Load LD matrix
    bm = hl.linalg.BlockMatrix.read(
        "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm"
    )

    def flip_variant(v):
        chrom, pos, a1, a2 = v.split(":")
        return f"{chrom}:{pos}:{a2}:{a1}"

    flipped_snp1 = flip_variant(snp1)
    flipped_snp2 = flip_variant(snp2)

    # Variants provided as Chr:Pos:Ref:Alt strings which get parsed in the
    # appropriate locus alleles struct to reference the Hail table
    variants = [
        hl.parse_variant(snp1),
        hl.parse_variant(snp2),
        # To resolve issues with allele ordering (flipped ref/alt in LD matrix)
        # Because only one will exist in the LD matrix this will end up being
        # only 2 variants in the output
        hl.parse_variant(flipped_snp1),
        hl.parse_variant(flipped_snp2),
        # TODO: resolve ambiguities such as strand flips
        # TODO: drop those that cannot be resolved
    ]

    # Filter the variant index table to only the two variants of interest
    # and select their indices (to match in the LD BlockMatrix)
    ht_idx = ht_idx.filter(
        hl.literal(variants).contains(
            hl.struct(locus=ht_idx.locus, alleles=ht_idx.alleles)
        )
    )
    idx = ht_idx.idx.collect()

    # Get correlation from LD matrix.
    # Sort indices first because the LD matrix is symmetric but only the lower triangle is stored.
    i, j = sorted(idx)
    r = bm[i, j]

    return r**2


parser = argparse.ArgumentParser(description="Load cis-eQTL data")

parser.add_argument("--gene", type=str, required=True, help="Gene name")
parser.add_argument(
    "--mark_lead", action="store_true", help="Mark lead SNP", default=True
)
parser.add_argument(
    "--selected_snp", type=str, required=True, help="Selected SNP to highlight"
)
parser.add_argument(
    "--data",
    type=str,
    required=True,
    help="Source of cis-eQTL data",
    choices=["INTERVAL", "eQTLgen"],
)
parser.add_argument(
    "--out",
    type=str,
    required=True,
    help="Output path for locuszoom plot",
)

args = parser.parse_args()
## Example usage:
# args = argparse.Namespace(
#     gene="PXK",
#     mark_lead=True,
#     selected_snp="rs9311676",
#     data="eQTLgen",
#     out="locuszoom_PXK_eQTLgen.png",
# )

if args.data == "eQTLgen":
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
    ]

    cis_eQTL = (
        pl.scan_csv(
            f"/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/cis_eQTLs_{args.data}.tsv",
            separator="\t",
        )
        .select(selected_cols)
        .filter(pl.col("GeneSymbol") == args.gene)
        .collect()
    )
else:
    # TODO: INTERVAL currently not working because of missing SNP and gene positions
    selected_cols = [
        "variant_rsid",
        "pvalue_nominal",
        "qval",
        "gene_name",
        "slope",
        "slope_se",
    ]

    cis_eQTL = (
        pl.scan_csv(
            f"/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/cis_eQTLs_{args.data}.tsv",
            separator="\t",
        )
        .select(selected_cols)
        .filter(pl.col("GeneSymbol") == args.gene)
        .with_columns((pl.col("slope") / pl.col("slope_se")).alias("Zscore"))
        .drop("slope", "slope_se")
        .rename({"variant_rsid": "SNP", "qval": "FDR", "gene_name": "GeneSymbol"})
        .collect()
    )


def locuszoom_plot(
    df,
    gene,
    selected_snp,
    out,
    mark_lead=True,
    fdr_thresh=0.05,
    window_kb=500,
    dpi=300,
):
    # Convert to pandas for matplotlib
    pdf = df.to_pandas()

    # Compute -log10(P)
    pdf["neglog10p"] = -np.log10(pdf["Pvalue"])

    # Identify lead SNP (smallest p-value)
    lead_snp = df.sort("Pvalue").select("SNP").row(0)[0]

    # Select chromosome of selected SNP (same chromosome as gene and lead SNP)
    chrom = df.filter(pl.col("SNP") == selected_snp).select("SNPChr").item()

    lead_snp_id = (
        df.filter(pl.col("SNP") == lead_snp)
        .select(
            pl.format(
                "{}:{}:{}:{}",
                pl.col("SNPChr"),
                pl.col("SNPPos"),
                pl.col("OtherAllele"),
                pl.col("AssessedAllele"),
            ).alias("SNP_id")
        )
        .item()
    )

    selected_snp_id = (
        df.filter(pl.col("SNP") == selected_snp)
        .select(
            pl.format(
                "{}:{}:{}:{}",
                pl.col("SNPChr"),
                pl.col("SNPPos"),
                pl.col("OtherAllele"),
                pl.col("AssessedAllele"),
            ).alias("SNP_id")
        )
        .item()
    )

    lead_selected_r2 = get_r2_between_snps(lead_snp_id, selected_snp_id)

    # Gene TSS
    gene_pos = pdf["GenePos"].iloc[0]

    # Optional zoom window
    if window_kb is not None:
        lo = gene_pos - window_kb * 1_000
        hi = gene_pos + window_kb * 1_000
        pdf = pdf.query("SNPPos >= @lo and SNPPos <= @hi")

    # Significance mask
    sig = pdf["FDR"] < fdr_thresh

    SNPs_in_LD_matrix = check_snps_in_ld_matrix(chrom, lo, hi, pdf["SNP"])

    pdf["in_ld"] = pdf["SNP"].isin(SNPs_in_LD_matrix)

    fig, ax = plt.subplots(figsize=(10, 5))

    # Non-significant SNPs
    ax.scatter(
        pdf.loc[~sig, "SNPPos"],
        pdf.loc[~sig, "neglog10p"],
        c="lightgray",
        s=10,
        label="FDR ≥ 0.05",
        zorder=1,
    )

    # Significant SNPs
    ax.scatter(
        pdf.loc[sig, "SNPPos"],
        pdf.loc[sig, "neglog10p"],
        c="steelblue",
        s=15,
        label="FDR < 0.05",
        zorder=2,
    )

    # SNPs in LD matrix
    ax.scatter(
        pdf.loc[pdf["in_ld"], "SNPPos"],
        pdf.loc[pdf["in_ld"], "neglog10p"],
        facecolors="none",
        edgecolors="purple",
        s=40,
        linewidths=0.8,
        label="In LD matrix?",
        zorder=2.5,
    )

    ax.text(
        0.02,
        0.98,
        rf"$r^2$(lead, selected) = {lead_selected_r2:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
    )

    # Lead SNP
    if mark_lead:
        lead = pdf[pdf["SNP"] == lead_snp]
        ax.scatter(
            lead["SNPPos"],
            lead["neglog10p"],
            c="red",
            s=70,
            marker="*",
            label="Lead SNP",
            zorder=3,
        )

        ax.annotate(
            lead_snp,
            (lead["SNPPos"].values[0], lead["neglog10p"].values[0]),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            color="red",
            weight="bold",
        )

    # Selected SNP
    if selected_snp in set(pdf["SNP"]):
        sel = pdf[pdf["SNP"] == selected_snp]
        ax.scatter(
            sel["SNPPos"],
            sel["neglog10p"],
            c="orange",
            s=60,
            marker="D",
            label="Selected SNP",
            zorder=3,
        )

        ax.annotate(
            selected_snp,
            (sel["SNPPos"].values[0], sel["neglog10p"].values[0]),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            color="darkorange",
            weight="bold",
        )

    # Gene TSS
    ax.axvline(
        gene_pos,
        color="black",
        linestyle="--",
        linewidth=1,
        label=f"{gene} TSS",
        zorder=0,
    )

    # Labels and title
    ax.set_xlabel(f"Position on chr{chrom}")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_title(f"cis-eQTL locus: {gene}")

    ax.legend(frameon=False)
    plt.tight_layout()

    # Save and close
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


locuszoom_plot(
    cis_eQTL,
    gene=args.gene,
    selected_snp=args.selected_snp,
    mark_lead=args.mark_lead,
    out=args.out,
)
