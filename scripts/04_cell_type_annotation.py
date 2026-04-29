"""
Orchestration Script: 04_Annotation
Project: Colorectal Cancer (CRC) Single-Cell Atlas
Stage: Cluster-to-CellType Mapping & Label Synchronization

This script performs the final mapping of Leiden clusters to biological
identities. It utilizes a curated annotation dictionary derived from
marker gene analysis (Stage 03). Crucially, it processes the integrated
'post' dataset first to establish a master label set, which is then
synchronized across related AnnData objects to ensure consistency.

Author: Loai Eletr
"""

import scanpy as sc
import config
from src.visualization import plot_cluster_comparison
from config import ClusterCharacterizationConfig
from utils import save_adata

import matplotlib

matplotlib.use('Agg')  # Non-interactive backend for headless/server execution


# ==============================================================================
# 1. ANNOTATION DICTIONARY
# ==============================================================================

def apply_manual_annotation(adata, cluster_key="leiden_res_0.5"):
    """
    Assigns high-resolution biological cell type labels to numerical clusters.

    Annotation Rationale (CRC-Specific):
    - C1QA_Macro: Immunosuppressive myeloid niche.
    - Malig_Entero: Malignant transition of the intestinal epithelium.
    - iCAF: Inflammatory CAFs often associated with immune exclusion.
    - IgG/kappa/lambda PC: Resolving the plasma cell repertoire.
    """

    cl_annotation = {
        "0": "C1QA_Macro",  # Resident-like macrophages
        "1": "Malig_Entero",  # Enterocyte-like colorectal tumor cells
        "2": "Treg",  # Regulatory T cells
        "3": "CTL",  # Cytotoxic T lymphocytes
        "4": "Mural_Act",  # Activated tumor-associated mural cells
        "5": "Neutro",  # Tumor-associated neutrophils
        "6": "Endo",  # Blood endothelial cells
        "7": "iCAF",  # Inflammatory CAFs
        "8": "IgG_PC",  # IgG plasma cells
        "9": "kappa_PC",  # κ-restricted plasma cells
        "10": "Naive_B",  # Naive / early B cells
        "11": "lambda_PC",  # λ-restricted plasma cells
        "12": "LEC",  # Lymphatic endothelial cells
        "13": "pDC",  # Plasmacytoid dendritic cells
        "14": "Mast",  # Mast cells
        "15": "Schwann",  # Schwann cells
        "16": "Fibro_quies",  # Adventitial / Quiescent fibroblasts
        "17": "Tuft"  # Tuft cells
    }

    # Map clusters to annotations safely, handling any unassigned clusters as 'Unknown'
    adata.obs["manual_celltype_annotation"] = (
        adata.obs[cluster_key]
        .astype(str)
        .map(cl_annotation)
        .fillna("Unknown")
    )

    return adata.obs["manual_celltype_annotation"]


# ==============================================================================
# 2. ORCHESTRATION LOGIC
# ==============================================================================

if __name__ == "__main__":
    cluster_cfg = ClusterCharacterizationConfig()

    # 1. FILE DISCOVERY & PRIORITIZATION
    # We sort to ensure the integrated 'post-integration' file is processed first.
    # This establishes the 'annotation_cache' needed for the 'pre-integration' files.
    SAMPLES = sorted(
        config.INTEGR_OUT.glob("*.h5ad"),
        key=lambda x: "post" in x.stem,
        reverse=True
    )

    if not SAMPLES:
        raise FileNotFoundError("❌ No .h5ad files found in integration output directory.")

    annotation_cache = None

    # 2. PROCESSING LOOP
    for sample in SAMPLES:
        print(f"--- 🔄 Processing: {sample.name} ---")
        adata = sc.read_h5ad(sample)

        if "post" in sample.stem:
            # MASTER ANNOTATION: Generate labels from the integrated coordinates
            print(f"--- 🏷️  Annotating INTEGRATED dataset: {sample.name} ---")
            adata.obs["manual_celltype_annotation"] = apply_manual_annotation(adata)

            # Cache the labels to maintain identity across file versions
            annotation_cache = adata.obs["manual_celltype_annotation"]

            # VISUAL VERIFICATION: UMAP with 'on data' legends for final review
            plot_cluster_comparison(
                adata=adata,
                color_keys=["manual_celltype_annotation"],
                file_prefix=f"{sample.stem}_final_cell_type_annotations",
                suptitle=None,
                save_dir=cluster_cfg.IO_PARAMS["save_dir_annotation_plots"],
                legend_loc="on data"
            )

        elif annotation_cache is not None:
            # LABEL SYNC: Apply the cached labels to the pre-integration version
            print(f"--- 🔗 Syncing labels to PRE-INTEGRATION dataset: {sample.name} ---")
            adata.obs["manual_celltype_annotation"] = annotation_cache
        else:
            # SAFETY CHECK: Prevent processing if the master list hasn't been created
            raise ValueError(f"❌ Aborted: Integrated 'post' file must be processed before '{sample.name}'.")

        # 3. DATA PERSISTENCE
        # Save to the final '04_annotation' directory
        out_path = config.ANNOTATION_OUT / f"{sample.stem}_final_annotated.h5ad"
        save_adata(adata, out_path)

    print("--- ✅ Annotation Pipeline Complete ---")