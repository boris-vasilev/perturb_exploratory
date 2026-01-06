import polars as pl
import argparse
import numpy as np
import matplotlib.pyplot as plt


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
#     gene="CSNK2B",
#     mark_lead=True,
#     selected_snp="rs4151657",
#     data="eQTLgen",
#     out="locuszoom_CSNK2B_eQTLgen.png",
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
    lead_idx = pdf["Pvalue"].idxmin()
    lead_snp = pdf.loc[lead_idx, "SNP"]

    # Gene TSS
    gene_pos = pdf["GenePos"].iloc[0]
    chrom = pdf["SNPChr"].iloc[0]

    # Optional zoom window
    if window_kb is not None:
        lo = gene_pos - window_kb * 1_000
        hi = gene_pos + window_kb * 1_000
        pdf = pdf.query("SNPPos >= @lo and SNPPos <= @hi")

    # Significance mask
    sig = pdf["FDR"] < fdr_thresh

    fig, ax = plt.subplots(figsize=(10, 5))

    # Non-significant
    ax.scatter(
        pdf.loc[~sig, "SNPPos"],
        pdf.loc[~sig, "neglog10p"],
        c="lightgray",
        s=10,
        label="FDR ≥ 0.05",
        zorder=1,
    )

    # Significant
    ax.scatter(
        pdf.loc[sig, "SNPPos"],
        pdf.loc[sig, "neglog10p"],
        c="steelblue",
        s=15,
        label="FDR < 0.05",
        zorder=2,
    )

    # Lead SNP
    if mark_lead:
        lead = pdf[pdf["SNP"] == lead_snp]
        ax.scatter(
            lead["SNPPos"],
            lead["neglog10p"],
            c="red",
            s=60,
            marker="*",
            label="Lead SNP",
            zorder=3,
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

    # Gene TSS
    ax.axvline(
        gene_pos,
        color="black",
        linestyle="--",
        linewidth=1,
        label=f"{gene} TSS",
        zorder=0,
    )

    # Labels
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
