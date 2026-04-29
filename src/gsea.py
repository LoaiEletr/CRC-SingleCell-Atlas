"""
Functional Enrichment Module: GSEA Orchestration.

This module automates the Gene Set Enrichment Analysis (GSEA) workflow. It
connects Differential Gene Expression (DGE) results with curated biological
knowledge bases (MSigDB) to identify coordinated pathway shifts in colorectal
cancer progression and liver metastasis.

Author: Loai Eletr
"""

import pandas as pd
import decoupler as dc
import numpy as np
from visualization import plot_gsea_results
from pathlib import Path

# ==============================================================================
# 1. KNOWLEDGE BASE RETRIEVAL (MSigDB)
# ==============================================================================

def get_msigdb_resource(collection: str = "hallmark"):
    """
    Downloads and formats MSigDB gene sets for decoupler compatibility.

    This function fetches specific gene set collections (like 'hallmark' or 'cp_kegg')
    and converts them into a long-format DataFrame. This format is the standard
    input for enrichment methods like ORA (Over-Representation Analysis) or
    GSEA (Gene Set Enrichment Analysis).

    Parameters
    ----------
    collection : str, default "hallmark"
        The MSigDB collection to retrieve.
        Common options: 'hallmark', 'go_biological_process', 'cp_kegg'.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns ['source', 'target'] mapping gene sets
        to their constituent gene symbols.
    """
    print("-" * 64)
    print(f"🌐 Fetching MSigDB collection: {collection}...")

    # 1. Connect to Decoupler's OmniPath/MSigDB mirror
    msigdb = dc.op.resource("MSigDB")

    # 2. Subset to the requested collection
    resource = msigdb[msigdb["collection"] == collection]

    # 3. Data Cleaning
    # Remove duplicates to ensure statistical tests aren't biased by redundant entries
    resource = resource[~resource.duplicated(("geneset", "genesymbol"))]

    # 4. Format for Decoupler API
    # Decoupler expects 'source' (the pathway) and 'target' (the gene)
    formatted_resource = resource[["geneset", "genesymbol"]].rename(
        columns={"geneset": "source", "genesymbol": "target"}
    )

    print(f"✅ Resource loaded. {len(formatted_resource['source'].unique())} gene sets found.")
    print("-" * 64)

    return formatted_resource


# ==============================================================================
# 2. CORE ENRICHMENT LOGIC (GSEA)
# ==============================================================================

def run_gsea(deg_df: pd.DataFrame, gene_sets: pd.DataFrame):
    """
    Calculates GSEA scores and p-values using decoupler.

    Unlike ORA (which only looks at 'significant' genes), GSEA considers the
    entire transcriptomic profile, making it much more sensitive to subtle
    but coordinated changes in pathway activity.

    Parameters
    ----------
    deg_df : pd.DataFrame
        The output from DESeq2/run_pseudobulk_pipeline.
        Must contain a 'log2FoldChange' column with gene symbols as the index.
    gene_sets : pd.DataFrame
        A formatted MSigDB resource (from get_msigdb_resource)
        containing 'source' and 'target' columns.

    Returns
    -------
    pd.DataFrame
        A table of pathways (source) with their corresponding
        enrichment scores and p-values, sorted by significance.
    """
    print("-" * 64)
    print("🧬 Running GSEA (Gene Set Enrichment Analysis)...")

    # 1. Rank Generation
    # GSEA requires genes to be ranked by a metric of 'association'.
    # We use LFC to identify which pathways are moving toward the 'Tumor'
    # (positive LFC) vs. the 'Normal' (negative LFC).
    ranking = deg_df[["log2FoldChange"]].sort_values("log2FoldChange", ascending=False).T

    # 2. Score Calculation
    # We use decoupler's GSEA implementation.
    # The 'seed' ensures the permutation-based p-values are reproducible.
    scores, pvals = dc.mt.gsea(ranking, gene_sets, seed=123)

    # 3. Result Formatting
    # Combine scores and p-values into a single tidy DataFrame.
    res = pd.concat({"score": scores.T, "pval": pvals.T}, axis=1).droplevel(1, axis=1)
    res = res.sort_values("pval").reset_index().rename(columns={'index': 'source'})

    # 4. Numerical Stability
    # Clipping p-values prevents math errors (log(0)) in downstream volcano
    # or bubble plot visualizations.
    res["pval"] = res["pval"].clip(lower=np.finfo(float).eps)

    print(f"✅ GSEA complete. Found {len(res[res['pval'] < 0.05])} significant pathways (p < 0.05).")
    print("-" * 64)

    return res


# ==============================================================================
# 3. PIPELINE ORCHESTRATOR
# ==============================================================================

def run_functional_enrichment_pipeline(
        deg_df: pd.DataFrame,
        collection: str = "hallmark",
        prefix: str = "GSEA",
        top_n: int = 20,
        save_dir_plots: str = "figures/gsea",
        save_dir_csv: str = "results/gsea"
):
    """
    Master Wrapper: Downloads gene sets, runs GSEA on DGE results,
    and generates visualization.

    Workflow:
    ---------
    1. Retrieval: Fetches a specific MSigDB collection (e.g., Hallmark, GO).
    2. Calculation: Executes GSEA to rank pathways by enrichment score.
    3. Visualization: Produces a bar/dot plot of the most significant pathways.
    4. Persistence: Saves the full enrichment table for downstream reporting.

    Parameters
    ----------
    deg_df : pd.DataFrame
        The DGE results table (from DESeq2). Must contain 'log2FoldChange'.
    collection : str, default "hallmark"
        The MSigDB category to test against.
    prefix : str, default "GSEA"
        A naming identifier for saving files (e.g., 'CD8_Tcells_Tumor_vs_Normal').
    top_n : int, default 20
        The number of top pathways to show in the summary plot.
    """
    print(f"----------------------------------------------------------------")
    print(f"🚀 Initializing Functional Enrichment Pipeline: {prefix}")

    # Ensure output directories exist
    Path(save_dir_plots).mkdir(parents=True, exist_ok=True)
    Path(save_dir_csv).mkdir(parents=True, exist_ok=True)

    # --- 1. RESOURCE ACQUISITION ---
    # Maps gene symbols to curated biological categories
    print(f"Fetching MSigDB resource: {collection}...")
    gene_sets = get_msigdb_resource(collection=collection)

    # --- 2. ENRICHMENT ANALYSIS ---
    # Ranks pathways based on the coordinated shift of their member genes
    print(f"Running GSEA for {prefix}...")
    gsea_results = run_gsea(deg_df, gene_sets)

    # --- 3. VISUALIZATION ---
    # Generates a visual summary of pathway activation/inhibition
    print(f"Generating summary plots (Top {top_n} pathways)...")
    plot_gsea_results(
        gsea_results,
        title=prefix,
        top_n=top_n,
        save_dir=save_dir_plots,
        prefix=prefix
    )

    # --- 4. DATA EXPORT ---
    # Save the full results for supplemental tables or cross-study comparison
    save_path = f"{save_dir_csv}/{prefix}_{collection}_results.csv"
    gsea_results.to_csv(save_path, index=False)

    print(f"✅ GSEA Complete.")
    print(f"📊 Results: {save_path}")
    print(f"🖼️  Plots: {save_dir_plots}")
    print("-" * 64)

    return gsea_results