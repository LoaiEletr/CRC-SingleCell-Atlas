"""
Module: Visualization Core
Project: Colorectal Cancer (CRC) & Liver Metastasis Single-Cell Atlas
Purpose: Standardized plotting engine for scRNA-seq quality control,
         integration diagnostics, and differential expression analysis.

Author: Loai Eletr
"""

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import scvi
import numpy as np
from pathlib import Path
import decoupler as dc
import pandas as pd
import seaborn as sns

# =============================================================================
# 1. QUALITY CONTROL & PREPROCESSING PLOTS
# =============================================================================

def plot_qc_metrics(adata: ad.AnnData, file_prefix: str, plot_label: str, save_dir: str = "figures/qc"):
    """
    Generates and saves diagnostic violin and scatter plots for QC.

    This function produces two standard diagnostic outputs:
    1. Violin Plots: Distribution of library size, gene richness, and MT content.
    2. Scatter Plots: Correlation between library size and gene richness,
       colored by mitochondrial percentage to identify apoptotic cells.

    Args:
        adata: AnnData object containing QC metrics in .obs.
        file_prefix: Unique identifier for the sample (e.g., 'CRC_Patient1_Tumor').
        plot_label: Descriptive stage for the title (e.g., 'Post-Filtering').
        save_dir: Path to save the resulting PDF files.
    """

    # --- 1. Directory Setup ---
    # Ensure the output path exists and configure Scanpy to save figures there.
    out_path = Path(save_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = out_path

    # Normalize labels for consistent filesystem naming
    file_label = plot_label.lower().replace(" ", "_")

    # --- 2. Violin Plot Generation ---
    # Used to see the 'spread' of your data. High MT counts are a sign of stress.
    print(f"{'-' * 60}")
    print(f"📊 Saving QC Violins: violin_{file_prefix}_{file_label}.pdf")
    sc.pl.violin(
        adata,
        ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
        jitter=0.4,  # Adds visibility to individual cell density
        multi_panel=True,  # Keeps metrics grouped for easier comparison
        show=False,  # Prevents GUI popup in headless mode
        save=f"_{file_prefix}_{file_label}.pdf"
    )

    # --- 3. Scatter Plot Generation ---
    # Vital for identifying the 'Empty Droplet' vs 'Healthy Cell' threshold.
    # Cells in the top left (high MT, low genes) are typically discarded.
    print(f"{'-' * 60}")
    print(f"📈 Saving QC Scatters: scatter_{file_prefix}_{file_label}.pdf")
    sc.pl.scatter(
        adata,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        title=f"{file_prefix} - {plot_label}",
        show=False,
        save=f"_{file_prefix}_{file_label}.pdf"
    )
    print(f"{'-' * 60}\n")

def plot_pca_variance(adata: ad.AnnData, file_prefix: str, n_pcs: int = 50, save_dir: str = "figures/qc"):
    """
    Plots the PCA variance ratio to identify the 'elbow' point.

    The 'elbow' represents the point where the addition of more PCs
    results in diminishing returns of explained variance. Choosing
    the correct number of PCs ensures that biological signal is
    preserved while technical noise is discarded.

    Args:
        adata: AnnData object after PCA has been computed.
        file_prefix: Unique sample identifier (e.g., 'Tumor_Sample_01').
        n_pcs: Number of PCs to display (default: 50).
        save_dir: Path to save the resulting PDF.
    """
    print(f"{'-' * 60}")
    print(f"📉 Plotting PCA Variance Ratio: {file_prefix}")
    print(f"💾 Saving to: pca_{file_prefix}_pca_variance.pdf")

    # --- 1. Path Management ---
    out_path = Path(save_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = out_path

    # --- 2. Variance Ratio Visualization ---
    # Log scale is used here (log=True) to make the 'elbow' more visible,
    # especially when the first few PCs explain a massive amount of variance.
    sc.pl.pca_variance_ratio(
        adata,
        n_pcs=n_pcs,
        log=True,
        show=False,
        save=f"_{file_prefix}_pca_variance.pdf"
    )
    print(f"{'-' * 60}\n")

# =============================================================================
# 2. CLUSTERING & ANNOTATION PLOTS
# =============================================================================

def plot_cluster_comparison(
        adata: ad.AnnData,
        color_keys: list[str],
        file_prefix: str,
        suptitle: str = "Cluster Comparison",
        save_dir: str = "figures/clustering",
        **kwargs
):
    """
    Generates a multi-panel UMAP plot for any specified observation keys.

    Args:
        adata: The AnnData object with 'X_umap' coordinates.
        color_keys: List of adata.obs columns to visualize (e.g., ['batch', 'leiden_0.5']).
        file_prefix: Sample/Project name for the filename.
        suptitle: Main title for the entire figure.
        save_dir: Path to store the PDF output.
        **kwargs: Passed to scanpy.pl.umap (e.g., size, palette, legend_loc).
    """
    print(f"{'-' * 60}")
    print(f"🎨 Visualizing UMAP layers: {', '.join(color_keys)}")

    # --- 1. Infrastructure Setup ---
    out_path = Path(save_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = out_path

    # --- 2. UMAP Generation ---
    # We use dynamic column adjustment (ncols) to ensure the panels are readable.
    # 'frameon=False' removes the bounding box for a cleaner, high-impact aesthetic.
    sc.pl.umap(
        adata,
        color=color_keys,
        ncols=min(3, len(color_keys)),
        title=[f"{k}" for k in color_keys],
        frameon=False,
        show=False,  # Important for manual saving/headlining
        **kwargs
    )

    # --- 3. Figure Customization & Persistence ---
    # If a suptitle is provided, we add it and use a specific naming convention.
    # 'bbox_inches=tight' is crucial to prevent legends from being cut off.
    if suptitle:
        filename = f"umap_{file_prefix}_{suptitle.lower().replace(' ', '_')}.pdf"
        print(f"💾 Saving annotated plot: {filename}")
        plt.suptitle(f"{suptitle}", fontsize=16)
        plt.savefig(out_path / filename, bbox_inches='tight')
    else:
        filename = f"umap_{file_prefix}.pdf"
        print(f"💾 Saving plot: {filename}")
        plt.savefig(out_path / filename, bbox_inches='tight')

    plt.close()  # Prevents memory leaks during large batch runs
    print(f"{'-' * 60}\n")


def plot_marker_dotplot(
        adata: ad.AnnData,
        cluster_key: str,
        marker_key: str,
        n_genes: int = 10,
        file_prefix: str = "",
        save_dir: str = "figures/clustering"
):
    """
    Generates a formatted dotplot for identified marker genes.

    This visualization is essential for Stage 03 (Characterization) to
    validate that identified markers are specific to their clusters.

    Args:
        adata: AnnData object with rank_genes_groups results.
        cluster_key: The obs column used for grouping (e.g., 'leiden').
        marker_key: The key in .uns where marker results are stored.
        n_genes: Number of top markers to show per cluster (default: 10).
        file_prefix: Unique sample identifier for naming.
        save_dir: Output directory for the PDF file.
    """
    print(f"{'-' * 60}")
    print(f"🧬 Generating Marker Dotplot for: {cluster_key}")
    print(f"📊 Using marker set: {marker_key}")

    # --- 1. Path Management ---
    out_path = Path(save_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = out_path

    # --- 2. Dotplot Execution ---
    # 'standard_scale="var"' scales each gene from 0 to 1 across clusters,
    # making it easier to compare genes with vastly different baseline
    # expression levels (e.g., GAPDH vs a specific transcription factor).
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=cluster_key,
        standard_scale="var",
        n_genes=n_genes,
        key=marker_key,
        save=f"_{file_prefix}_{cluster_key}.pdf",
        show=False  # Prevents window popups in automated runs
    )

    print(f"💾 Plot saved to: {save_dir}/dotplot_{file_prefix}_{cluster_key}.pdf")
    print(f"{'-' * 60}\n")


# =============================================================================
# 3. INTEGRATION (scVI) DIAGNOSTIC PLOTS
# =============================================================================

def plot_scvi_training_history(model: scvi.model.SCVI, prefix: str = "", save_dir: str = "figures/integration"):
    """
    Plots and saves ELBO training/validation curves.

    The Negative ELBO should show a sharp decline and eventually plateau.
    A significant gap between Train and Validation curves can indicate
    overfitting, while a curve that is still dropping at the final epoch
    suggests that 'n_epochs' should be increased.

    Args:
        model: The trained scVI model object.
        prefix: Unique identifier for the run (e.g., 'CRC_Atlas_v1').
        save_dir: Path to save the resulting PDF.
    """
    print(f"{'-' * 60}")
    print(f"📈 Plotting scVI Training History: {prefix}")

    # --- 1. Plot Initialization ---
    fig, ax = plt.subplots(figsize=(10, 6))

    # --- 2. History Extraction ---
    # We plot the 'Negative ELBO'. Lower is better.
    if 'elbo_train' in model.history:
        model.history['elbo_train'].plot(ax=ax, label='Train')
    if 'elbo_validation' in model.history:
        model.history['elbo_validation'].plot(ax=ax, label='Validation')

    # --- 3. Formatting ---
    ax.set_title(f'scVI Training Progress (ELBO) - {prefix}')
    ax.set_xlabel('Epochs')
    ax.set_ylabel('Loss (Negative ELBO)')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.6)  # Added for better readability of the plateau

    # --- 4. Persistence ---
    out_path = Path(save_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    save_file = out_path / f"{prefix}_scvi_training_history.pdf"
    plt.savefig(save_file, bbox_inches='tight')
    print(f"💾 Diagnostic plot saved to: {save_file}")

    plt.show()
    print(f"{'-' * 60}\n")


def plot_latent_variance(adata: ad.AnnData, prefix: str = "", save_dir: str = "figures/integration"):
    """
    Creates and saves a bar plot of scVI latent dimension variances.

    A high variance in a latent dimension suggests it is capturing a
    significant biological or technical signal. If only 2-3 dimensions
    have high variance and the rest are zero, the model may be under-utilizing
    its capacity. Conversely, if all 30 have high variance, the latent
    space is dense with information.

    Args:
        adata: AnnData object containing the 'X_scVI' representation in .obsm.
        prefix: Unique identifier for the integration run.
        save_dir: Path to save the resulting PDF.
    """
    print(f"{'-' * 60}")
    print(f"📊 Analyzing scVI Latent Variance: {prefix}")

    # --- 1. Variance Calculation ---
    # Retrieve the scVI embedding and calculate the empirical variance
    # for each of the n_latent (30) dimensions.
    if "X_scVI" not in adata.obsm:
        raise KeyError("X_scVI not found in obsm. Ensure scVI integration has run.")

    latent_representation = adata.obsm["X_scVI"]
    variances = np.var(latent_representation, axis=0)

    # --- 2. Sorting ---
    # Sort dimensions by variance in descending order for a cleaner "scree-like" plot.
    sorted_idx = np.argsort(variances)[::-1]
    sorted_vars = variances[sorted_idx]

    # --- 3. Visualization ---
    plt.figure(figsize=(8, 4))
    plt.bar(range(len(sorted_vars)), sorted_vars, color='skyblue', edgecolor='navy', alpha=0.7)
    plt.title(f"Empirical Variance of scVI Latent Dimensions - {prefix}")
    plt.ylabel("Variance")
    plt.xlabel("Latent Dimension (Sorted)")
    plt.grid(axis='y', linestyle='--', alpha=0.3)

    # --- 4. Persistence ---
    out_path = Path(save_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    save_file = out_path / f"{prefix}_latent_variance.pdf"
    plt.savefig(save_file, bbox_inches='tight')
    print(f"💾 Latent variance plot saved to: {save_file}")

    plt.show()
    print(f"{'-' * 60}\n")

# =============================================================================
# 4. DIFFERENTIAL EXPRESSION & PATHWAY PLOTS
# =============================================================================

def plot_pseudobulk_filtering(
        pdata: ad.AnnData,
        prefix: str = "pseudobulk",
        min_count: int = 10,
        min_total_count: int = 15,
        large_n: int = 10,
        min_prop_expr: float = 0.7,
        min_prop_cell: float = 0.1,
        min_smpls: int = 2,
        save_dir: str = "figures/pseudobulk"
):
    """
    Diagnostic plots specifically for assessing pseudobulk filtering thresholds.

    Args:
        pdata: The pseudobulk AnnData object (cells aggregated by sample).
        prefix: Filename prefix for the saved PDFs.
        min_count: Minimum count for a gene in at least 'large_n' samples.
        min_total_count: Minimum total count across all samples.
        large_n: Minimum number of samples to consider for expression filtering.
        min_prop_expr: Minimum proportion of samples expressing the gene.
        min_prop_cell: Minimum proportion of cells expressing the gene (at sc level).
        min_smpls: Minimum number of samples requiring 'min_prop_cell'.
        save_dir: Output directory for the diagnostic plots.
    """
    print(f"{'-' * 60}")
    print(f"🧪 Assessing Pseudobulk Filtering for: {prefix}")

    # Ensure output directory exists
    Path(save_dir).mkdir(parents=True, exist_ok=True)

    # --- 1. Filter by Expression (Counts) ---
    # This plot identifies genes with sufficient library size.
    # It helps you see if your threshold is too aggressive (losing real signal)
    # or too lax (keeping background noise).
    print("[*] Generating Filter-by-Expression diagnostic...")
    dc.pl.filter_by_expr(
        pdata,
        group="condition",
        min_count=min_count,
        min_total_count=min_total_count,
        large_n=large_n,
        min_prop=min_prop_expr
    )
    plt.savefig(f"{save_dir}/{prefix}_filter_expr.pdf", bbox_inches='tight')
    plt.show()

    # --- 2. Filter by Proportion (Prevalence) ---
    # This ensures genes aren't just highly expressed in a single 'outlier' sample.
    # Essential for CRC where patient-to-patient variability is high.
    print("[*] Generating Filter-by-Proportion diagnostic...")
    dc.pl.filter_by_prop(pdata, min_prop=min_prop_cell, min_smpls=min_smpls)
    plt.savefig(f"{save_dir}/{prefix}_filter_prop.pdf", bbox_inches='tight')
    plt.show()

    print(f"✅ Pseudobulk diagnostics saved to: {save_dir}")
    print(f"{'-' * 60}\n")


def plot_pseudobulk_pca(
        pdata: ad.AnnData,
        prefix: str,
        color_keys: list = ["condition", "batch"],
        save_dir: str = "figures/pseudobulk"
):
    """
    Flexible PCA plotting for aggregated pseudobulk samples.

    Args:
        pdata: The pseudobulk AnnData object (obs = samples, var = genes).
        prefix: Unique identifier for the cell type or group being analyzed.
        color_keys: Metadata columns to color the PCA points (e.g., condition).
        save_dir: Path to store the resulting PDF.
    """
    print(f"{'-' * 60}")
    print(f"🧬 Generating Pseudobulk PCA for: {prefix}")

    # --- 1. Infrastructure Setup ---
    out_path = Path(save_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = out_path

    # --- 2. PCA Visualization ---
    # We use a larger point size (300) because pseudobulk datasets have
    # fewer observations (samples) than single-cell datasets (cells).
    # 'frameon=True' is used here to provide a formal axis for sample distribution.
    sc.pl.pca(
        pdata,
        color=color_keys,
        ncols=1,
        size=300,
        frameon=True,
        show=False,
        save=f"_{prefix}_pca.pdf"
    )

    print(f"💾 Sample-level PCA saved to: {save_dir}/pca_{prefix}_pca.pdf")
    print(f"{'-' * 60}\n")


def plot_volcano(
        results_df,
        prefix="",
        lfc_thr=1.0,
        padj_thr=0.05,
        save_dir: str = "figures/differential_expression"
):
    """
    Generates and saves a volcano plot for differential expression results.

    Genes in the top-right are significantly upregulated in the Tumor,
    while genes in the top-left are significantly downregulated (or
    upregulated in the Normal condition).

    Args:
        results_df: DataFrame containing 'log2FoldChange' and 'padj' columns.
        prefix: Unique identifier (e.g., the cell type name like 'C1QA_Macro').
        lfc_thr: Log2 Fold Change threshold (default 1.0 = 2-fold change).
        padj_thr: Adjusted p-value significance threshold (default 0.05).
        save_dir: Path to store the PDF output.
    """
    print(f"{'-' * 60}")
    print(f"🌋 Generating Volcano Plot for: {prefix}")

    # Ensure output directory exists
    Path(save_dir).mkdir(parents=True, exist_ok=True)

    # --- 1. Volcano Visualization ---
    # Decoupler's volcano plot automatically handles the -log10 transformation
    # for the p-values and colors genes based on the provided thresholds.
    dc.pl.volcano(
        results_df,
        x="log2FoldChange",
        y="padj",
        thr_stat=lfc_thr,
        thr_sign=padj_thr
    )

    # --- 2. Persistence ---
    save_file = f"{save_dir}/{prefix}_volcano.pdf"
    plt.savefig(save_file, bbox_inches='tight')
    print(f"💾 Volcano plot saved to: {save_file}")

    plt.show()
    print(f"{'-' * 60}\n")


def plot_gsea_results(
        gsea_results: pd.DataFrame,
        title: str = "GSEA",
        top_n: int = 20,
        prefix="",
        save_dir: str = "figures/gsea"
):
    """
    Plots a bar chart of the top enriched pathways based on significance.

    Args:
        gsea_results: DataFrame containing 'source' (pathway name) and 'pval'.
        title: Description of the comparison (e.g., 'Tumor vs Normal').
        top_n: Number of top-scoring pathways to display.
        prefix: Filename prefix (e.g., 'C1QA_Macro').
        save_dir: Path to store the resulting PDF.
    """
    print(f"{'-' * 60}")
    print(f"🧬 Visualizing Functional Enrichment: {prefix}")

    # Ensure output directory exists
    Path(save_dir).mkdir(parents=True, exist_ok=True)

    # 1. Figure Initialization
    # 10x8 is the 'sweet spot' for horizontal bar charts in publications.
    plt.figure(figsize=(10, 8))

    # 2. Data Transformation
    # We use -log10(p-value) because it turns small fractions into
    # large, intuitive "significance scores."
    data = gsea_results.head(top_n).copy()
    data["-log10_pval"] = -np.log10(data["pval"])

    # 3. Visualization
    # 'viridis' palette ensures the most significant hits are visually distinct.
    sns.barplot(data=data, x="-log10_pval", y="source", palette="viridis")

    # 4. Formatting
    plt.title(f"Top {top_n} Enriched Pathways: {title}")
    plt.xlabel("-log10(p-value)")
    plt.ylabel("Biological Pathway")
    plt.grid(axis='x', linestyle='--', alpha=0.4)

    # 5. Persistence
    save_file = f"{save_dir}/{prefix}_gsea.pdf"
    plt.savefig(save_file, bbox_inches='tight')
    print(f"💾 GSEA bar plot saved to: {save_file}")

    plt.show()
    print(f"{'-' * 60}\n")