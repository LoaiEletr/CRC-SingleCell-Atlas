"""
Pseudobulk Orchestration Module: Differential Expression (DE) Pipeline.

This script manages the transition from single-cell resolution to patient-level
statistical inference. It leverages pseudobulk aggregation to resolve
pseudoreplication issues and applies rigorous filtering to ensure robust
DE results across Normal, Tumor, and Metastatic conditions.

Key Workflow:
-------------
1. Aggregation: Summing counts per cell-type per sample using Decoupler.
2. Statistical Filtering: edgeR-style expression and occupancy checks.
3. Modeling: Negative Binomial regression via DESeq2.
4. Export: Automated saving of plots, AnnData checkpoints, and CSV results.

Author: Loai Eletr
"""

# --- Standard Library & Data Science ---
import pandas as pd

# --- Bioinformatics Frameworks ---
import anndata as ad
import decoupler as dc
import scanpy as sc

# --- Custom Project Utilities ---
from preprocessing import normalize_adata, run_pca
from visualization import plot_pseudobulk_filtering, plot_pseudobulk_pca, plot_volcano
from deg import run_deseq2_analysis, get_significant_genes
from utils import save_adata

# ==============================================================================
# 1. PSEUDOBULK GENERATION
# ==============================================================================

def create_pseudobulk(
        adata: ad.AnnData,
        sample_col: str = "sample_id",
        groups_col: str = "manual_celltype_annotation",
        layer: str = "soupX_counts"
) -> ad.AnnData:
    """
    Aggregates single-cell data into pseudobulk profiles.

    Using the 'sum' mode, this function collapses thousands of individual cells
    into a single expression profile per cell type per sample. This is essential
    for multi-factor experimental designs (e.g., Patient + Tissue Type).

    Parameters
    ----------
    adata : ad.AnnData
        The integrated and annotated AnnData object.
    sample_col : str, default "sample_id"
        The column in .obs identifying biological replicates (e.g., Patient_N1).
    groups_col : str, default "manual_celltype_annotation"
        The column in .obs identifying the cell populations (e.g., T-cells, B-cells).
    layer : str, default "soupX_counts"
        The count layer to aggregate. 'soupX_counts' is recommended as it
        removes technical ambient noise before aggregation.

    Returns
    -------
    ad.AnnData
        A new AnnData object where each 'observation' is a unique
        Sample-Group combination.
    """
    print("-" * 64)
    print(f"Generating pseudobulk using {groups_col} grouped by {sample_col}...")
    print(f"Input Layer: {layer}")

    # Utilize Decoupler's high-performance pseudobulk utility
    # 'mode="sum"' is standard for count-based differential expression tools.
    pdata = dc.pp.pseudobulk(
        adata=adata,
        sample_col=sample_col,
        groups_col=groups_col,
        layer=layer,
        mode="sum"
    )

    # Note: pdata.obs will now contain the cross-product of samples and groups.
    print(f"✅ Pseudobulk complete. Profiles generated: {pdata.n_obs}")
    print("-" * 64)

    return pdata


# ==============================================================================
# 2. STATISTICAL PREPROCESSING
# ==============================================================================

def preprocess_pseudobulk(
        pdata: ad.AnnData,
        min_count: int = 10,
        min_total_count: int = 15,
        large_n: int = 10,
        min_prop_expr: float = 0.7,  # Argument for filter_by_expr
        min_prop_cell: float = 0.1,  # Argument for filter_by_prop (e.g., 10% of cells)
        min_smpls: int = 2  # Min samples that must meet min_prop_cell
):
    """
    Handles the math: filtering, normalization, and PCA for pseudobulk.

    Workflow:
    1. edgeR-style filtering based on library size and count depth.
    2. Occupancy filtering to ensure genes are present in enough samples.
    3. Normalization and PCA for exploratory visualization (PCoA/PCA).
    4. Reset to raw counts to satisfy DESeq2/edgeR input requirements.
    """
    print("-" * 64)
    print("🚀 Starting Pseudobulk Preprocessing...")

    # 1. Filter by expression (edgeR-style)
    # This keeps genes that have a minimum count in a minimum number of samples,
    # adjusted for the total library size of the samples.
    print(f"Applying count-based filtering (min_count={min_count})...")
    dc.pp.filter_by_expr(
        pdata, group="condition", min_count=min_count,
        min_total_count=min_total_count, large_n=large_n, min_prop=min_prop_expr
    )

    # 2. Filter by proportion (decoupler-style)
    # Ensures the gene isn't just a high-expression outlier in one single cell/sample.
    # It must be expressed in at least 'min_prop_cell' of cells in 'min_smpls' samples.
    print(f"Applying occupancy filtering (min_prop_cell={min_prop_cell})...")
    dc.pp.filter_by_prop(
        pdata, min_prop=min_prop_cell, min_smpls=min_smpls
    )

    # 3. Layer & Metadata Management
    # We must keep raw counts because DE tools (DESeq2) model the raw data.
    pdata.layers["counts"] = pdata.X.copy()
    pdata.obs["lib_size"] = pdata.X.sum(1)

    # 4. Exploratory Dimensionality Reduction
    # We normalize and run PCA purely to check if samples group by 'condition' (N vs T).
    print("Performing technical normalization and PCA for QC visualization...")
    normalize_adata(pdata)
    run_pca(pdata)

    # 5. Reset to Raw Counts
    # Crucial step: swap back to the 'counts' layer before handing off to DE functions.
    dc.pp.swap_layer(adata=pdata, key="counts", inplace=True)

    print(f"✅ Preprocessing complete. Final gene count: {pdata.n_vars}")
    print("Layer 'counts' is currently active for downstream DE analysis.")
    print("-" * 64)

# ==============================================================================
# 3. FULL PIPELINE ORCHESTRATION
# ==============================================================================

def run_pseudobulk_pipeline(
        adata: ad.AnnData,
        prefix: str,
        design_factors: str,
        contrast: list,
        lfc_thr: float = 1.0,
        padj_thr: float = 0.05,
        # --- Metadata & Layer Inputs ---
        sample_col: str = "sample_id",
        groups_col: str = "manual_celltype_annotation",
        layer: str = "soupX_counts",
        # --- Threshold Inputs ---
        min_count: int = 10,
        min_total_count: int = 15,
        large_n: int = 10,
        min_prop_expr: float = 0.7,
        min_prop_cell: float = 0.1,
        min_smpls: int = 2,
        color_keys: list = ["condition", "batch"],
        save_dir_plots: str = "figures/pseudobulk_analysis",
        save_adata_folder: str = "results/06_pseudobulk_analysis"
):
    """
    Orchestrates the entire workflow from raw adata to DE results.

    Workflow Stages:
    1. Aggregation: Collapses cells into pseudobulk samples (Decoupler).
    2. Diagnostics: Visualizes the impact of filtering thresholds.
    3. Math: Performs count-based filtering and library-size normalization.
    4. Dimensionality Reduction: PCA to verify biological separation (e.g., N vs T).
    5. DE Testing: Executes DESeq2 and identifies significant genes.
    6. Persistence: Exports H5AD objects, full results, and significant gene lists.
    """
    print(f"----------------------------------------------------------------")
    print(f"🚀 Initializing Pseudobulk Pipeline for: {prefix}")

    # --- 1. PSEUDOBULK GENERATION ---
    # Aggregates single-cell clusters into "Bulk-like" profiles
    pdata = create_pseudobulk(
        adata,
        sample_col=sample_col,
        groups_col=groups_col,
        layer=layer
    )

    # --- 2. QC & FILTERING DIAGNOSTICS ---
    # Crucial for verifying that thresholds aren't too aggressive for rare cell types
    plot_pseudobulk_filtering(
        pdata, prefix=prefix, min_count=min_count, min_total_count=min_total_count, large_n=large_n,
        min_prop_expr=min_prop_expr, min_prop_cell=min_prop_cell, min_smpls=min_smpls,
        save_dir=f"{save_dir_plots}/{prefix}"
    )

    # --- 3. STATISTICAL PREPROCESSING ---
    # Filters low-signal genes and normalizes for library size differences
    preprocess_pseudobulk(
        pdata, min_count=min_count, min_total_count=min_total_count, large_n=large_n,
        min_prop_expr=min_prop_expr, min_prop_cell=min_prop_cell, min_smpls=min_smpls
    )

    # --- 4. EXPLORATORY PCA ---
    # If samples don't separate by 'condition' here, DE results will likely be weak
    plot_pseudobulk_pca(pdata, prefix=prefix, color_keys=color_keys, save_dir=f"{save_dir_plots}/{prefix}")

    # --- 5. DIFFERENTIAL EXPRESSION (DESeq2) ---
    # The heavy lifting: Modeling counts with a Negative Binomial distribution
    results_df = run_deseq2_analysis(pdata, design_factors, contrast)

    up, down = get_significant_genes(results_df)
    print(f"✅ Analysis complete. Found {len(up)} UP and {len(down)} DOWN genes.")

    # Visualization of global expression shifts
    plot_volcano(results_df, prefix, lfc_thr, padj_thr, f"{save_dir_plots}/{prefix}")

    # --- 6. DATA PERSISTENCE ---
    # Ensure H5AD serialization doesn't fail due to complex design objects
    if 'design_matrix' in pdata.obsm:
        pdata.obsm['design_matrix'] = pd.DataFrame(pdata.obsm['design_matrix'])

    # Save outputs for the specific cell type/subset defined by 'prefix'
    save_adata(pdata, f"{save_adata_folder}/{prefix}/{prefix}_pseudobulk.h5ad")
    results_df.to_csv(f"{save_adata_folder}/{prefix}/{prefix}_all_genes_deseq2.csv")

    # Filter significant hits for quick downstream review (Pathway analysis/GSEA)
    sig_genes = results_df[
        (results_df['padj'] < padj_thr) &
        (results_df['log2FoldChange'].abs() > lfc_thr)
        ].sort_values('padj')

    sig_genes.to_csv(f"{save_adata_folder}/{prefix}/{prefix}_significant_DE_genes.csv")

    print(f"🏁 Files successfully exported to {save_adata_folder}/{prefix}/")
    print("-" * 64)

    return pdata, results_df