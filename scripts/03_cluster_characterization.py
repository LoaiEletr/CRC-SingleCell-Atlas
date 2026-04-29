"""
Orchestration Script: 03_Characterization
Project: Colorectal Cancer (CRC) Single-Cell Atlas
Stage: Cluster Marker Discovery & Cell-Type Validation

This script executes the cluster characterization workflow. It identifies
statistically significant marker genes for each integrated cluster,
applies stringent occupancy filters to ensure marker specificity, and
generates the visual evidence (DotPlots) required for manual cell-type
annotation (e.g., identifying T-cell subsets, TAMs, or Malignant Epithelium).

Author: Loai Eletr
"""

from src.deg import identify_cluster_markers, export_marker_results
from src.visualization import plot_marker_dotplot
import scanpy as sc
import config
from pathlib import Path

# ==============================================================================
# MAIN EXECUTION BLOCK
# ==============================================================================

if __name__ == "__main__":
    """
    Entry point for the cluster characterization workflow.

    The workflow follows these steps:
    1. Loads markers parameters (cluster_key, thresholds) from config.
    2. Ingests the integrated atlas (post-batch correction).
    3. Performs Wilcoxon-based Marker Discovery.
    4. Exports 'Raw' vs 'Clean' marker lists for supplementary data.
    5. Visualizes marker specificity using a high-resolution DotPlot.
    """

    # --- 1. Configuration & Environment Setup ---
    # Accesses parameters like cluster_key (e.g., 'leiden_res_0.8') and filtering thresholds.
    cluster_cfg = config.ClusterCharacterizationConfig()

    # --- 2. Data Ingestion ---
    # Automatically locate the most recent integrated dataset from Stage 02.
    # Uses rglob to find files containing 'post' in the name (e.g., integrated_post_scVI.h5ad).
    try:
        SAMPLE = next(config.INTEGR_OUT.rglob("*post*.h5ad"))
        print(f"📂 Loading integrated atlas: {SAMPLE.name}")
        adata = sc.read_h5ad(SAMPLE)
    except StopIteration:
        raise FileNotFoundError("❌ No integrated h5ad files found in the Stage 02 output directory.")

    # --- 3. Differential Expression Analysis (Discovery) ---
    # Wilcoxon rank-sum test determines cluster uniqueness.
    # Returns keys for both the 'raw' statistical hits and the 'filtered' biological hits.
    print(f"🧬 Running marker discovery on: {cluster_cfg.MARKER_PARAMS['cluster_key']}")
    dea_key, dea_filtered_key = identify_cluster_markers(
        adata=adata,
        **cluster_cfg.MARKER_PARAMS
    )

    # --- 4. Data Export (Documentation) ---
    # We export two versions:
    # 1. Full results (raw) - necessary for volcano plots or exhaustive searches.
    # 2. Filtered results - high-confidence markers for cell type naming.

    print(f"📦 Exporting marker CSVs to: {cluster_cfg.IO_PARAMS['marker_file']}")

    # Export raw marker results
    export_marker_results(
        adata=adata,
        dea_key=dea_key,
        output_path=cluster_cfg.IO_PARAMS["marker_file"],
        prefix=SAMPLE.stem,
        n_top=cluster_cfg.N_GENES_PARAMS
    )

    # Export filtered marker results (The 'clean' set)
    export_marker_results(
        adata=adata,
        dea_key=dea_filtered_key,
        output_path=cluster_cfg.IO_PARAMS["marker_file"],
        prefix=f"{SAMPLE.stem}_filtered",
        n_top=cluster_cfg.N_GENES_PARAMS
    )

    # --- 5. Visualization (Evidence) ---
    # Generates a DotPlot: The X-axis is the markers, the Y-axis is the clusters.
    # Dot size = fraction of cells; Dot color = average expression.
    print(f"🖼️  Generating Marker DotPlot (Top {cluster_cfg.N_GENES_PARAMS} genes)...")
    plot_marker_dotplot(
        adata=adata,
        cluster_key=cluster_cfg.MARKER_PARAMS["cluster_key"],
        marker_key=dea_filtered_key,
        n_genes=cluster_cfg.N_GENES_PARAMS,
        file_prefix=f"{SAMPLE.stem}_filtered",
        save_dir=cluster_cfg.IO_PARAMS["save_dir_plots"]
    )

    print(f"✅ Characterization complete for {SAMPLE.name}.")
    print(f"📊 Results available in: {cluster_cfg.IO_PARAMS['marker_file']}")