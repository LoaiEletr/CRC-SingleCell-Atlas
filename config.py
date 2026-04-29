import os
import sys
import json
import warnings
import matplotlib
import scanpy as sc
from pathlib import Path

# =============================================================================
# 1. SYSTEM & PATHWAY INFRASTRUCTURE
# =============================================================================
# Silence TensorFlow warnings to keep the console clean during scVI integration
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

# Define the root of the project relative to this config file
BASE_DIR = Path(__file__).resolve().parent

# Define a standard directory tree for raw data, results, and figures
DATA_RAW      = BASE_DIR / "data" / "raw"
RESULTS_ROOT  = BASE_DIR / "results"
FIGURES_ROOT  = BASE_DIR / "figures"

# Stage-specific output directories for modular pipeline execution
PREPROC_OUT   = RESULTS_ROOT / "01_preprocessing"
INTEGR_OUT    = RESULTS_ROOT / "02_integration"
CLUSTER_CHARACTERIZATION_OUT   = RESULTS_ROOT / "03_cluster_characterization"
ANNOTATION_OUT = RESULTS_ROOT / "04_cell_type_annotation"
PSEUDOBULK_OUT = RESULTS_ROOT / "05_pseudobulk_analysis"
GSEA_OUT = RESULTS_ROOT / "06_gsea"
INDIVIDUAL_CELLCHAT_DIR = RESULTS_ROOT / "07_cell_cell_communication_analysis"
DIFF_CELLCHAT_DIR = RESULTS_ROOT / "08_cell_cell_communication_differential"

# Auto-generate the folder tree to prevent "FileNotFound" errors on new systems
for folder in [DATA_RAW, PREPROC_OUT, INTEGR_OUT, CLUSTER_CHARACTERIZATION_OUT,
               ANNOTATION_OUT, PSEUDOBULK_OUT, GSEA_OUT, FIGURES_ROOT]:
    folder.mkdir(parents=True, exist_ok=True)

# =============================================================================
# 2. SOFTWARE ENGINE SETTINGS
# =============================================================================
# Use 'Agg' for headless environments (servers/VMs) to prevent GUI crashes
matplotlib.use('Agg')
sc.settings.autoshow = False
sc.settings.verbosity = 1
warnings.simplefilter(action='ignore', category=FutureWarning)

# Point to the specific R installation for rpy2/reticulate interop
R_HOME = r'C:\Program Files\R\R-4.5.2'
SAVE_PLOTS = True

# =============================================================================
# 3. STAGE 01: PREPROCESSING (Single-Sample)
# =============================================================================
class PreprocessingConfig:
    # Orchestration metadata
    GENERAL_PARAMS = {
        "prefix": None,
        "batch_key": None,
        "save_plots": SAVE_PLOTS,
    }

    # Quality Control: Using MAD (Median Absolute Deviation) for adaptive filtering
    QC_PARAMS = {
        "mad_usage": True,
        "total_counts_bounds": None,
        "n_genes_by_counts_bounds": None,
        "pct_counts_mt_upper": 12, # Mitochondrial threshold for CRC samples
        "mad_cutoff": 3,
    }

    # Normalization & Dimensionality Reduction (PCA)
    NORM_PARAMS = {
        "target_sum": 1e4,
        "n_top_genes": 2000,
        "n_pcs": 30
    }

    # Initial clustering for SoupX ambient RNA correction
    CLUSTERING_PARAMS = {
        "n_neighbors": 30,
        "use_rep": "X_pca",
        "resolutions": [0.25, 0.5, 0.75, 1],
        "resolution_for_soupx": 0.5
    }

# =============================================================================
# 4. STAGE 02: INTEGRATION (scVI Batch Correction)
# =============================================================================
class IntegrationConfig:
    # Generative model parameters for batch-effect removal
    MODEL_PARAMS = {
        "n_layers": 2,
        "n_latent": 30,
        "n_epochs": 150,
        "n_top_genes": 2000,
        "batch_key": "batch",
    }

    # Evaluation metrics and UMAP visualization keys
    EVAL_PARAMS = {
        "reference_batch": "S1",
        "resolutions": [0.25, 0.5, 0.75, 1],
        "n_neighbors": 50,
        "use_rep": "X_scVI",
        "color_keys": ["condition", "batch"],
    }

    # Input/Output paths for scVI models and integrated AnnData
    IO_PARAMS = {
        "prefix": "merged_dataset",
        "model_path": INTEGR_OUT / "scvi_model",
        "save_dir_adata": INTEGR_OUT,
        "save_dir_plots": FIGURES_ROOT / "integration",
        "save_plots": SAVE_PLOTS
    }

# =============================================================================
# 5. STAGE 03 & 04: CHARACTERIZATION & ANNOTATION
# =============================================================================
class ClusterCharacterizationConfig:
    # Thresholds for Wilcoxon rank-sum marker discovery
    MARKER_PARAMS = {
        "cluster_key": "leiden_res_0.5",
        "min_in_group_fraction": 0.25,     # Sensitivity: gene must be in 25% of cells
        "max_out_group_fraction": 0.50,    # Specificity: gene must be restricted
    }

    N_GENES_PARAMS = 20

    IO_PARAMS = {
        "marker_file": CLUSTER_CHARACTERIZATION_OUT,
        "save_dir_plots": FIGURES_ROOT / "cluster_characterization",
        "save_dir_annotation_plots": FIGURES_ROOT / "cell_type_annotation",
    }

# =============================================================================
# 6. STAGE 05: PSEUDOBULK DIFFERENTIAL EXPRESSION
# =============================================================================
class PseudobulkConfig:
    PSEUDOBULK_FILTER_PARAMS = {
        "layer": "soupX_counts", # Use raw-ish counts for DESeq2 statistical modeling
    }

    # Contrast definitions (Tumor vs Normal) for PyDESeq2
    PSEUDOBULK_MODEL_PARAMS = {
        "layer": "soupX_counts",
        "groups_col": "manual_celltype_annotation",
        "sample_col": "sample_id",
        "contrast": ["condition", "Tumor", "Normal"],
        "design_factors": "~ batch + condition", # Modeling batch effects in DE
        "lfc_thr": 1.0,                # 2-fold change threshold
        "padj_thr": 0.05,              # 5% FDR significance
        "color_keys": ["condition", "batch"]
    }

    PSEUDOBULK_IO_PARAMS = {
        "save_dir_plots": FIGURES_ROOT / "pseudobulk",
        "save_adata_folder": PSEUDOBULK_OUT,
    }

# =============================================================================
# 7. STAGE 06: GSEA (Pathway Analysis)
# =============================================================================
class GSEAConfig:
    GSEA_ANALYSIS_PARAMS = {
        "collection": "hallmark", # Molecular Signatures Database (MSigDB) targets
        "top_n": 20,
    }

    GSEA_IO_PARAMS = {
        "save_dir_plots": FIGURES_ROOT / "gsea",
        "save_dir_csv": GSEA_OUT,
    }

# =============================================================================
# 8. STAGE 07 & 08: CELLCHAT (Communication Analysis)
# =============================================================================
class CellChatConfig:
    # Optimized communication pattern counts (from diagnostic K-selection)
    K_SELECTION = {
        "Normal_samples_annotated": {"k_out": 3, "k_in": 4},
        "Tumor_samples_annotated": {"k_out": 4, "k_in": 4}
    }

    CELLCHAT_PARAMS = {
        "species": "human",
        "min_cells": 10 # Minimum cells per group to be included in signaling
    }

    CELLCHAT_IO = {
        "input_dir": INDIVIDUAL_CELLCHAT_DIR,
        "individual_cellchat_dir": INDIVIDUAL_CELLCHAT_DIR,
        "diff_cellchat_dir": DIFF_CELLCHAT_DIR,
        "fig_dir": FIGURES_ROOT / "cellchat_communication",
        "k_json": INDIVIDUAL_CELLCHAT_DIR / "k_config.json"
    }

    @classmethod
    def prepare_environment(cls):
        """
        Initializes the cross-language (Python/R) bridge.
        Writes K-selection to JSON and cleans the environment for reticulate.
        """
        cls.CELLCHAT_IO["individual_cellchat_dir"].mkdir(parents=True, exist_ok=True)
        cls.CELLCHAT_IO["diff_cellchat_dir"].mkdir(parents=True, exist_ok=True)
        cls.CELLCHAT_IO["fig_dir"].mkdir(parents=True, exist_ok=True)

        # Serialize K-values for the R worker script
        with open(cls.CELLCHAT_IO["k_json"], "w") as f:
            json.dump(cls.K_SELECTION, f, indent=4)

        # Environment cleanup to prevent reticulate/Python path collisions
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env.pop("PYTHON_SESSION_INITIALIZED", None)
        env["RETICULATE_PYTHON"] = sys.executable
        return env

    @classmethod
    def split_data_by_condition(cls, save_fn):
        """
        Partitions the atlas into condition-specific objects (Normal vs Tumor).
        Required for comparative CellChat analysis.
        """
        input_path = next(ANNOTATION_OUT.rglob("*pre*.h5ad"))
        adata = sc.read_h5ad(input_path)

        for condition in ["Normal", "Tumor"]:
            filename = f"{condition}_samples_annotated.h5ad"
            sub_adata = adata[adata.obs["condition"] == condition].copy()
            save_fn(sub_adata, str(cls.CELLCHAT_IO["input_dir"] / filename))