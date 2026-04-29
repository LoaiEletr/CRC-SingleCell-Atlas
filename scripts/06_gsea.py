"""
Orchestration Script: 06_GSEA
Project: Colorectal Cancer (CRC) Single-Cell Atlas
Stage: Gene Set Enrichment (GSEA)

This script performs functional enrichment analysis on the results of the
pseudobulk DE stage. It maps gene-level changes to biological pathways
(MSigDB, Hallmark, KEGG), allowing for the identification of systemic
biological shifts in the tumor microenvironment.

Author: Loai Eletr
"""

from src.gsea import run_functional_enrichment_pipeline
import config
import scanpy as sc
from config import GSEAConfig
import pandas as pd

# ==============================================================================
# MAIN EXECUTION BLOCK
# ==============================================================================

if __name__ == "__main__":
    # --- 1. Configuration & Data Discovery ---
    gsea_cfg = GSEAConfig()

    # Dynamically locate all DESeq2 results from the Stage 05 output directory.
    # This automatically picks up 'enterocyte' and 'macrophage' results.
    SAMPLES = list(config.PSEUDOBULK_OUT.glob("**/*all_genes_deseq2.csv"))

    if not SAMPLES:
        print("⚠️ No DESeq2 results found. Ensure Stage 05 has completed successfully.")

    # --- 2. ENRICHMENT LOOP ---
    for sample in SAMPLES:
        # Load the DEG table (index_col=0 ensures gene symbols are the identifiers)
        # BUG-FIX NOTE: Using sample instead of SAMPLES[0] inside the loop
        # to ensure each specific cell-type story is processed.
        df = pd.read_csv(sample, index_col=0)

        # Clean the prefix for clean file naming (e.g., 'Tumor_vs_Normal_in_enterocytes')
        analysis_prefix = f"{sample.stem.split('_all_genes')[0]}_gsea"

        print(f"--- 🧬 Running GSEA for: {analysis_prefix} ---")

        # --- 3. PIPELINE EXECUTION ---
        run_functional_enrichment_pipeline(
            deg_df=df,
            prefix=analysis_prefix,
            **GSEAConfig.GSEA_ANALYSIS_PARAMS,  # e.g., collection=['Hallmark']
            **GSEAConfig.GSEA_IO_PARAMS  # e.g., save_dir_plots,
        )

    print("--- ✅ Functional Enrichment Complete ---")