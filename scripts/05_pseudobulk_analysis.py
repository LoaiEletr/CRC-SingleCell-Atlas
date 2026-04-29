"""
Orchestration Script: 05_Pseudobulk_DE
Project: Colorectal Cancer (CRC) Single-Cell Atlas
Stage: Cell-Type Specific Comparative Transcriptomics

This script executes the pseudobulk differential expression pipeline.
By aggregating single cells into biological replicates (pseudobulk),
it identifies robust gene expression shifts between Tumor and Normal
states. We focus on specific 'stories': epithelial transformation
(enterocytes) and immune modulation (macrophages).

Author: Loai Eletr
"""

from src.pseudobulk import run_pseudobulk_pipeline
import config
import scanpy as sc
from config import PseudobulkConfig

# ==============================================================================
# MAIN EXECUTION BLOCK
# ==============================================================================

if __name__ == "__main__":
    # --- 1. Configuration & Data Ingestion ---
    pseudobulk_cfg = PseudobulkConfig()

    # Locate the final annotated pre-integration file (Stage 04 output)
    # We use 'pre' integration counts to ensure we use the most 'raw' biological data
    SAMPLE = next(config.ANNOTATION_OUT.glob("*pre*.h5ad"))
    print(f"📦 Loading annotated atlas for Pseudobulk: {SAMPLE.name}")
    adata = sc.read_h5ad(SAMPLE)

    # --- 2. Cohort Filtering ---
    # Focusing the analysis on a clear 'Tumor vs Normal' contrast.
    # Metastasis is excluded here to maintain a clean binary comparison.
    adata = adata[adata.obs["condition"] != "Metastasis"].copy()

    # --- 3. Biological Replicate Definition ---
    # Creating a unique identifier for each pseudobulk 'bin'.
    # This ensures that cells from the same patient/replicate are grouped together.
    adata.obs['sample_id'] = (
            adata.obs['condition'].astype(str) + "_" +
            adata.obs['replicate'].astype(str)
    )

    # ==========================================================================
    # STORY SELECTION: Cell-Type Specific Subsetting
    # ==========================================================================

    # STORY B: Epithelial Transformation
    # Analyzing how Malignant Enterocytes differ from their Normal counterparts.
    print("🔬 Isolating Malignant Enterocytes...")
    adata_epi = adata[adata.obs['manual_celltype_annotation'] == 'Malig_Entero'].copy()

    # STORY C: Monocyte/Macrophage Activation
    # Analyzing the polarization of C1QA+ Macrophages in the tumor niche.
    print("🔬 Isolating C1QA+ Macrophages...")
    adata_mono = adata[adata.obs['manual_celltype_annotation'] == 'C1QA_Macro'].copy()

    # ==========================================================================
    # EXECUTION: Pseudobulk Pipeline
    # ==========================================================================

    # Each cell type is processed through the PyDESeq2-backed pipeline.
    # The thresholds (min_count, min_prop_expr) ensure only high-quality,
    # consistently expressed genes are modeled.

    print(f"🚀 Running Pseudobulk for Epithelial Cells (Tumor vs Normal)")
    run_pseudobulk_pipeline(
        adata=adata_epi,
        prefix="Tumor_vs_Normal_in_enterocytes",
        min_count=10,  # Minimum counts for a gene to be kept
        min_total_count=15,  # Minimum total counts across samples
        min_prop_expr=0.7,  # Gene must be in this % of samples
        min_prop_cell=0.1,  # Gene must be in this % of cells
        min_smpls=3,  # Minimum number of samples per condition
        large_n=50,  # Threshold for large sample size approximations
        **pseudobulk_cfg.PSEUDOBULK_MODEL_PARAMS,
        **pseudobulk_cfg.PSEUDOBULK_IO_PARAMS
    )

    print(f"🚀 Running Pseudobulk for Macrophages (Tumor vs Normal)")
    run_pseudobulk_pipeline(
        adata=adata_mono,
        prefix="Tumor_vs_Normal_in_macrophages",
        min_count=10,
        min_total_count=15,
        min_prop_expr=0.7,
        min_prop_cell=0.1,
        min_smpls=3,
        large_n=50,
        **pseudobulk_cfg.PSEUDOBULK_MODEL_PARAMS,
        **pseudobulk_cfg.PSEUDOBULK_IO_PARAMS
    )

    print("✅ Pseudobulk DE Analysis Complete.")