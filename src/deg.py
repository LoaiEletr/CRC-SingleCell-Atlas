"""
Differential Expression & Marker Identification Module.

This script manages two primary workflows:
1. Marker Identification: Finding cluster-specific genes using Wilcoxon tests.
2. Pseudobulk DE: Comparing biological conditions (e.g., Tumor vs Normal)
   using the PyDESeq2 framework for statistically robust results.

Author: Loai Eletr
"""

import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats
from utils import swap_adata_layers

# ==============================================================================
# 1. CLUSTER MARKER DISCOVERY (Single-Cell Level)
# ==============================================================================

def identify_cluster_markers(
        adata: ad.AnnData,
        cluster_key: str,
        min_in_group_fraction: float = 0.25,
        max_out_group_fraction: float = 0.50
) -> tuple[str, str]:
    """
    Runs Wilcoxon rank-sum test and applies fraction-based filtering.

    This function finds genes differentially expressed in one cluster vs. the rest.
    The filtering step is crucial for "clean" annotations, removing genes that
    might be statistically significant but are expressed in too many other
    cell types to be useful markers.

    Parameters
    ----------
    adata : ad.AnnData
        The processed AnnData object.
    cluster_key : str
        The observation (.obs) column containing the clusters (e.g., 'leiden_res_0.5').
    min_in_group_fraction : float, default 0.25
        Minimum proportion of cells in the cluster that must express the gene.
    max_out_group_fraction : float, default 0.50
        Maximum proportion of cells outside the cluster allowed to express the gene.

    Returns
    -------
    tuple of (str, str)
        The keys added to .uns containing the raw and filtered DEA results.
    """
    dea_key = f"dea_{cluster_key}"
    dea_filtered_key = f"{dea_key}_filtered"

    print("-" * 64)
    print(f"Calculating marker genes for: {cluster_key}")

    # 1. Wilcoxon Rank-Sum Test
    # A non-parametric test ideal for single-cell data, as it doesn't
    # assume a specific distribution for gene expression.
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
        key_added=dea_key
    )

    # 2. Occupancy Filtering
    # This removes genes like ribosomal/mitochondrial genes that might be
    # 'significant' but are not cell-type specific.
    print(f"Applying filters (min_in: {min_in_group_fraction}, max_out: {max_out_group_fraction})...")
    sc.tl.filter_rank_genes_groups(
        adata,
        groupby=cluster_key,
        min_in_group_fraction=min_in_group_fraction,
        max_out_group_fraction=max_out_group_fraction,
        key=dea_key,
        key_added=dea_filtered_key
    )

    print(f"✅ Marker identification complete. Keys: {dea_key}, {dea_filtered_key}")
    print("-" * 64)

    return dea_key, dea_filtered_key


def export_marker_results(
        adata: ad.AnnData,
        dea_key: str,
        output_path: Path,
        prefix: str,
        n_top: int = 20
):
    """
    Exports both the full marker list and a top-N summary per cluster.

    This function bridges the gap between Scanpy's internal storage and
    external reporting, ensuring results are saved in a tidy format
    compatible with spreadsheet software and downstream enrichment tools.

    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object containing the results in .uns.
    dea_key : str
        The key where rank_genes_groups results are stored (e.g., 'dea_leiden').
    output_path : Path
        A pathlib.Path object pointing to the target directory.
    prefix : str
        The identifier for the run (e.g., 'CRC_Atlas_v1').
    n_top : int, default 20
        The number of top genes to filter for the summary table.

    Returns
    -------
    tuple of (pd.DataFrame, pd.DataFrame)
        The full marker DataFrame and the subsetted Top N DataFrame.
    """
    print("-" * 64)
    print(f"📦 Exporting marker results for key: {dea_key}")

    # 1. Data Retrieval
    # sc.get.rank_genes_groups_df converts complex dict structures in .uns
    # into a tidy-format DataFrame with 'group', 'names', 'logfoldchanges', etc.
    full_df = sc.get.rank_genes_groups_df(adata, group=None, key=dea_key)

    # 2. Curated Summary Extraction
    # We group by cluster and take the leading entries (lowest p-values/highest scores).
    top_n_df = full_df.groupby('group').head(n_top).copy()

    # Ensure numerical stability by removing empty entries (if any)
    top_n_df.dropna(subset=['names'], inplace=True)

    # 3. Persistence
    # Constructing paths using pathlib ensures OS compatibility (Linux/Windows)
    full_csv = output_path / f"{prefix}_all_markers_full.csv"
    top_csv = output_path / f"{prefix}_top{n_top}_markers.csv"

    full_df.to_csv(full_csv, index=False)
    top_n_df.to_csv(top_csv, index=False)

    print(f"✅ Successfully exported markers to: {output_path}")
    print(f"   ↳ Full Results: {full_csv.name}")
    print(f"   ↳ Summary (Top {n_top}): {top_csv.name}")
    print("-" * 64)

    return full_df, top_n_df

# ==============================================================================
# 2. COMPARATIVE TRANSCRIPTOMICS (Pseudobulk Level)
# ==============================================================================

def run_deseq2_analysis(pdata, design_factors, contrast_list, n_cpus=2):
    """
    Executes the DESeq2 pipeline and returns the statistics results.

    This function performs the three core DESeq2 steps:
    1. Estimation of size factors (normalization for sequencing depth).
    2. Estimation of dispersion (measuring gene-wise noise).
    3. GLM fitting and Wald testing (identifying significant changes).

    Parameters
    ----------
    pdata : ad.AnnData
        The aggregated pseudobulk object.
    design_factors : str
        The experimental design formula (e.g., 'condition' or 'batch + condition').
    contrast_list : list
        The specific comparison to make (e.g., ['condition', 'Tumor', 'Normal']).
    n_cpus : int, default 2
        Number of cores for parallel processing.

    Returns
    -------
    pd.DataFrame
        A table of statistics including log2FoldChange, p-values, and padj.
    """
    print("-" * 64)
    print(f"Running DESeq2 analysis with design: {design_factors}")

    # 1. LAYER MANAGEMENT
    # DESeq2 requires raw, non-normalized integers.
    # We swap to our 'counts' layer to ensure normalization hasn't skewed the model.
    pdata = swap_adata_layers(pdata, "counts")
    pdata.X = pdata.layers["counts"].copy()

    # 2. MODEL INITIALIZATION
    # We use DefaultInference for standard Wald tests and Cooks distance filtering.
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        adata=pdata,
        design=design_factors,
        refit_cooks=True,  # Handles outliers that might skew results
        inference=inference
    )

    # 3. EXECUTION
    # This runs the normalization, dispersion estimation, and model fitting.
    dds.deseq2()

    # 4. CONTRAST EXTRACTION
    # We extract the specific comparison defined in contrast_list.
    print(f"Extracting contrast: {contrast_list}")
    stat_res = DeseqStats(dds, contrast=contrast_list, inference=inference)
    stat_res.summary()

    print("✅ DESeq2 analysis complete.")
    print("-" * 64)

    return stat_res.results_df

def get_significant_genes(results_df, padj_thr=0.05, lfc_thr=1.0):
    """
    Filters the DESeq2 results table for biologically and statistically
    significant genes.

    Parameters
    ----------
    results_df : pd.DataFrame
        The output from run_deseq2_analysis.
    padj_thr : float, default 0.05
        The False Discovery Rate (FDR) threshold.
    lfc_thr : float, default 1.0
        The minimum Log2 Fold Change (2-fold difference).

    Returns
    -------
    tuple of (list, list)
        Lists of gene symbols for UP and DOWN regulated genes.
    """
    # UP: Signficant p-value AND positive change (Tumor > Normal)
    up = results_df[
        (results_df["padj"] < padj_thr) &
        (results_df["log2FoldChange"] > lfc_thr)
        ].index.tolist()

    # DOWN: Significant p-value AND negative change (Tumor < Normal)
    down = results_df[
        (results_df["padj"] < padj_thr) &
        (results_df["log2FoldChange"] < -lfc_thr)
        ].index.tolist()

    return up, down