import polars as pl
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import gzip


def load_interval_eqtl_for_gene(gene_id, chromosome, snp_rsids, eqtl_dir):
    """
    Load INTERVAL eQTL data for a specific gene from chromosome-specific files.

    Parameters
    ----------
    gene_id : str
        Ensembl gene ID (e.g., ENSG00000120645)
    chromosome : str or int
        Chromosome number
    snp_rsids : list
        List of rsIDs to filter for (from SNP file)
    eqtl_dir : str
        Directory containing INTERVAL eQTL summary statistics

    Returns
    -------
    polars.DataFrame
        Filtered eQTL data for the gene
    """
    eqtl_file = Path(eqtl_dir) / f"INTERVAL_eQTL_nominal_chr{chromosome}.tsv"

    if not eqtl_file.exists():
        raise FileNotFoundError(f"INTERVAL eQTL file not found: {eqtl_file}")

    # Read only the gene of interest
    eqtl_data = pl.read_csv(eqtl_file, separator="\t").filter(
        pl.col("phenotype_id") == gene_id
    )

    # Filter to SNPs that are in our LD matrix
    rsid_set = set(snp_rsids)
    eqtl_data = eqtl_data.filter(pl.col("variant_id").is_in(rsid_set))

    # Calculate Z-score from slope and slope_se
    eqtl_data = eqtl_data.with_columns(
        (pl.col("slope") / pl.col("slope_se")).alias("Zscore")
    )

    # Calculate p-value based FDR (simple Bonferroni for now, or use actual FDR if available)
    # Note: You might want to use a proper FDR calculation here
    eqtl_data = eqtl_data.with_columns(
        (pl.col("pval_nominal") * pl.len()).alias(
            "FDR"
        )  # Simple Bonferroni approximation
    )

    return eqtl_data


def align_eqtl_to_snps(eqtl_data, snps):
    """
    Merge eQTL data with SNP info and align effect directions.

    Parameters
    ----------
    eqtl_data : polars.DataFrame
        INTERVAL eQTL data
    snps : polars.DataFrame
        SNP information with alleles and indices

    Returns
    -------
    polars.DataFrame
        Merged and aligned data
    """
    # Merge on rsID
    merged = eqtl_data.join(
        snps.rename({"rsid": "variant_id"}), on="variant_id", how="inner"
    )

    # Align Z-scores: if effect_allele matches alleleB, keep sign; if matches alleleA, flip
    merged = merged.with_columns(
        pl.when(pl.col("effect_allele") == pl.col("alleleB"))
        .then(pl.col("Zscore"))
        .when(pl.col("effect_allele") == pl.col("alleleA"))
        .then(-pl.col("Zscore"))
        .otherwise(None)
        .alias("Zscore_aligned")
    )

    # Filter out SNPs where alignment couldn't be determined
    merged = merged.filter(pl.col("Zscore_aligned").is_not_null())

    # Sort by index
    merged = merged.sort("idx")

    return merged


def locuszoom_plot_interval(
    gene,
    gene_id,
    selected_snp,
    out,
    snp_dir="/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_exploratory/data/LD_loci",
    eqtl_dir="/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/INTERVAL_eQTL_summary_statistics",
    ld_dir="/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_exploratory/data/LD_loci",
    mark_lead=True,
    pval_thresh=None,
    window_kb=None,
    dpi=300,
):
    """
    Create a locuszoom plot for INTERVAL cis-eQTL data.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g., 'COPZ1')
    gene_id : str
        Ensembl gene ID (e.g., 'ENSG00000120645')
    selected_snp : str
        SNP to highlight (rsID)
    out : str
        Output file path for plot
    snp_dir : str
        Directory containing SNP info files
    eqtl_dir : str
        Directory containing INTERVAL eQTL summary statistics
    ld_dir : str
        Directory containing LD matrices
    mark_lead : bool
        Whether to mark the lead SNP
    pval_thresh : float or None
        P-value threshold for significance (None = use gene-specific threshold from phenotype summary)
    window_kb : int or None
        Window size around gene TSS in kb (None = use all data)
    dpi : int
        DPI for output figure
    """

    # Load SNP info
    snp_file = Path(snp_dir) / f"{gene}_SNPs.tsv"
    if not snp_file.exists():
        raise FileNotFoundError(f"SNP file not found: {snp_file}")

    snps = pl.read_csv(snp_file, separator="\t")

    # Load gene-specific p-value threshold from phenotype summary if not provided
    if pval_thresh is None:
        phenotype_summary_file = Path(eqtl_dir) / "INTERVAL_eQTL_phenotype_summary.tsv"
        if phenotype_summary_file.exists():
            phenotype_summary = pl.read_csv(phenotype_summary_file, separator="\t")
            gene_info = phenotype_summary.filter(pl.col("phenotype_id") == gene_id)

            if len(gene_info) > 0:
                pval_thresh = gene_info.select("pval_nominal_threshold").item(
                    row=0, column=0
                )
                print(f"Using gene-specific p-value threshold: {pval_thresh:.2e}")
            else:
                print(
                    f"Warning: Gene {gene_id} not found in phenotype summary, using default threshold 5e-8"
                )
                pval_thresh = 5e-8
        else:
            print(
                f"Warning: Phenotype summary file not found, using default threshold 5e-8"
            )
            pval_thresh = 5e-8
    else:
        print(f"Using provided p-value threshold: {pval_thresh:.2e}")

    # Extract chromosome from first locus
    first_locus = snps.select("locus").item(row=0, column=0)
    chrom = first_locus.split(":")[0]

    # Get gene position (use median of SNP positions as proxy for TSS)
    snps = snps.with_columns(
        pl.col("locus").str.split(":").list.get(1).cast(pl.Int64).alias("pos_b37")
    )
    gene_pos = snps.select(pl.col("pos_b37").median()).item(row=0, column=0)

    # Load INTERVAL eQTL data for this gene
    eqtl_data = load_interval_eqtl_for_gene(
        gene_id=gene_id,
        chromosome=chrom,
        snp_rsids=snps.select("rsid").to_series().to_list(),
        eqtl_dir=eqtl_dir,
    )

    if len(eqtl_data) == 0:
        raise ValueError(f"No eQTL data found for gene {gene_id} on chromosome {chrom}")

    # Align eQTL data with SNP info
    df = align_eqtl_to_snps(eqtl_data, snps)

    # Convert to pandas for plotting
    pdf = df.to_pandas()

    # Compute -log10(P)
    pdf["neglog10p"] = -np.log10(pdf["pval_nominal"])

    # Identify lead SNP (smallest p-value)
    lead_idx = pdf["pval_nominal"].idxmin()
    lead_snp = pdf.loc[lead_idx, "variant_id"]

    # Optional zoom window
    if window_kb is not None:
        lo = gene_pos - window_kb * 1_000
        hi = gene_pos + window_kb * 1_000
        pdf = pdf.query("pos_b37 >= @lo and pos_b37 <= @hi").copy()

    # Significance mask
    sig = pdf["pval_nominal"] < pval_thresh

    # Check which SNPs are in the LD matrix
    # The LD matrix contains SNPs that were used in the finemapping
    # These correspond to the indices in the SNP file
    ld_files = glob.glob(str(Path(ld_dir) / f"{gene}_*.bgz"))
    snps_in_ld = set()
    ld_matrix = None

    if ld_files:
        try:
            # Load LD matrix (handle gzip compression)
            with gzip.open(ld_files[0], "rb") as f:
                ld_matrix = np.loadtxt(f)

            # Get the indices that are in the LD matrix
            # These should be the sorted unique indices from the filtered data
            all_indices_in_data = sorted(pdf["idx"].unique())

            # The LD matrix dimensions tell us how many SNPs are included
            n_snps_in_ld = ld_matrix.shape[0]

            # The first n_snps_in_ld indices are in the LD matrix
            if n_snps_in_ld <= len(all_indices_in_data):
                snps_in_ld = set(all_indices_in_data[:n_snps_in_ld])
        except Exception as e:
            print(f"Could not load LD matrix: {e}")

    # Mark SNPs that are in LD matrix
    pdf["in_ld"] = pdf["idx"].isin(snps_in_ld)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 5))

    # Non-significant SNPs
    ax.scatter(
        pdf.loc[~sig, "pos_b37"],
        pdf.loc[~sig, "neglog10p"],
        c="lightgray",
        s=10,
        label=f"p ≥ {pval_thresh:.0e}",
        zorder=1,
    )

    # Significant SNPs
    ax.scatter(
        pdf.loc[sig, "pos_b37"],
        pdf.loc[sig, "neglog10p"],
        c="steelblue",
        s=15,
        label=f"p < {pval_thresh:.0e}",
        zorder=2,
    )

    # SNPs in LD matrix (hollow purple circles)
    if snps_in_ld:
        ax.scatter(
            pdf.loc[pdf["in_ld"], "pos_b37"],
            pdf.loc[pdf["in_ld"], "neglog10p"],
            facecolors="none",
            edgecolors="purple",
            s=40,
            linewidths=0.8,
            label="In LD matrix",
            zorder=2.5,
        )

    # Lead SNP
    if mark_lead:
        lead = pdf[pdf["variant_id"] == lead_snp]
        if len(lead) > 0:
            ax.scatter(
                lead["pos_b37"],
                lead["neglog10p"],
                c="red",
                s=70,
                marker="*",
                label="Lead SNP",
                zorder=3,
            )

            ax.annotate(
                lead_snp,
                (lead["pos_b37"].values[0], lead["neglog10p"].values[0]),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                fontsize=9,
                color="red",
                weight="bold",
            )

    # Selected SNP
    if selected_snp in set(pdf["variant_id"]):
        sel = pdf[pdf["variant_id"] == selected_snp]
        ax.scatter(
            sel["pos_b37"],
            sel["neglog10p"],
            c="orange",
            s=60,
            marker="D",
            label="Selected SNP",
            zorder=3,
        )

        ax.annotate(
            selected_snp,
            (sel["pos_b37"].values[0], sel["neglog10p"].values[0]),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            color="darkorange",
            weight="bold",
        )

        # Calculate r² with lead SNP if LD matrix is available
        if ld_matrix is not None:
            try:
                # Get indices
                lead_idx_ld = pdf[pdf["variant_id"] == lead_snp]["idx"].values[0]
                sel_idx_ld = sel["idx"].values[0]

                # Check if both SNPs are in the LD matrix
                if lead_idx_ld in snps_in_ld and sel_idx_ld in snps_in_ld:
                    # Find positions in the LD matrix
                    all_indices = sorted(snps_in_ld)
                    lead_pos = all_indices.index(lead_idx_ld)
                    sel_pos = all_indices.index(sel_idx_ld)

                    # Get r² from LD matrix
                    if lead_pos < ld_matrix.shape[0] and sel_pos < ld_matrix.shape[1]:
                        r2 = ld_matrix[lead_pos, sel_pos] ** 2

                        # Add r² annotation
                        ax.text(
                            0.02,
                            0.98,
                            rf"$r^2$(lead, selected) = {r2:.3f}",
                            transform=ax.transAxes,
                            ha="left",
                            va="top",
                            fontsize=10,
                            bbox=dict(
                                boxstyle="round,pad=0.3",
                                fc="white",
                                ec="gray",
                                alpha=0.8,
                            ),
                        )
                else:
                    print(f"Warning: Lead or selected SNP not in LD matrix")
            except Exception as e:
                print(f"Could not calculate LD: {e}")

    # Gene TSS (approximate)
    ax.axvline(
        gene_pos,
        color="black",
        linestyle="--",
        linewidth=1,
        label=f"{gene} (approx. TSS)",
        zorder=0,
    )

    # Labels and title
    ax.set_xlabel(f"Position on chr{chrom} (GRCh37)")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_title(f"cis-eQTL locus: {gene} (INTERVAL)")

    ax.legend(frameon=False, fontsize=8, loc="best")
    plt.tight_layout()

    # Save and close
    fig.savefig(out, dpi=dpi)
    plt.close(fig)

    print(f"Plot saved to: {out}")
    print(f"Lead SNP: {lead_snp} (p = {pdf.loc[lead_idx, 'pval_nominal']:.2e})")
    print(f"Total SNPs plotted: {len(pdf)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create locuszoom plot for INTERVAL cis-eQTL data"
    )

    parser.add_argument(
        "--gene", type=str, required=True, help="Gene symbol (e.g., COPZ1)"
    )
    parser.add_argument(
        "--gene_id",
        type=str,
        required=True,
        help="Ensembl gene ID (e.g., ENSG00000120645)",
    )
    parser.add_argument(
        "--selected_snp",
        type=str,
        required=True,
        help="Selected SNP to highlight (rsID)",
    )
    parser.add_argument(
        "--out", type=str, required=True, help="Output path for locuszoom plot"
    )
    parser.add_argument(
        "--snp_dir",
        type=str,
        default="/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_exploratory/data/LD_loci",
        help="Directory containing SNP files",
    )
    parser.add_argument(
        "--eqtl_dir",
        type=str,
        default="/rds/project/rds-csoP2nj6Y6Y/biv22/data/eqtl/INTERVAL_eQTL_summary_statistics",
        help="Directory containing INTERVAL eQTL summary statistics",
    )
    parser.add_argument(
        "--ld_dir",
        type=str,
        default="/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_exploratory/data/LD_loci",
        help="Directory containing LD matrices",
    )
    parser.add_argument(
        "--mark_lead",
        action="store_true",
        default=True,
        help="Mark lead SNP (default: True)",
    )
    parser.add_argument(
        "--no_mark_lead",
        action="store_false",
        dest="mark_lead",
        help="Do not mark lead SNP",
    )
    parser.add_argument(
        "--pval_thresh",
        type=float,
        default=None,
        help="P-value threshold for significance (default: None, use gene-specific threshold from phenotype summary)",
    )
    parser.add_argument(
        "--window_kb",
        type=int,
        default=None,
        help="Window size around gene TSS in kb (default: None, use all data)",
    )
    parser.add_argument(
        "--dpi", type=int, default=300, help="DPI for output figure (default: 300)"
    )

    args = parser.parse_args()

    # Example usage:
    # python plot_cis_eQTL_INTERVAL.py --gene COPZ1 --gene_id ENSG00000120645 --selected_snp rs11170507 --out locuszoom_COPZ1_INTERVAL.png

    locuszoom_plot_interval(
        gene=args.gene,
        gene_id=args.gene_id,
        selected_snp=args.selected_snp,
        out=args.out,
        snp_dir=args.snp_dir,
        eqtl_dir=args.eqtl_dir,
        ld_dir=args.ld_dir,
        mark_lead=args.mark_lead,
        pval_thresh=args.pval_thresh,
        window_kb=args.window_kb,
        dpi=args.dpi,
    )
