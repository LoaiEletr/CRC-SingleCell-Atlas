"""
Atlas Integration Module: scVI Workflow for CRC Single-Cell Data.

This script orchestrates the merging of Normal (N), Primary Tumor (T),
and Liver Metastasis (LM) samples. It leverages scVI's deep learning
framework to regress technical batch effects while preserving
the complex biological variance inherent in the tumor microenvironment.

Key Features:
-------------
1. Automated Metadata Mapping: Extracts condition/replicate info from file headers.
2. Latent Space Integration: Learns a non-linear batch-corrected representation (scVI).
3. Counterfactual Normalization: Projects expression into a unified reference batch.
4. Integrated Clustering: Optimized Leiden and UMAP on corrected latent coordinates.

Author: Loai Eletr
"""

# --- Standard Library & Data Science ---
import numpy as np
from pathlib import Path

# --- Bioinformatics Frameworks ---
import anndata as ad
import scanpy as sc
import scvi

# --- Custom Project Utilities ---
from utils import swap_adata_layers, save_adata, save_scvi_model
from preprocessing import normalize_adata, identify_highly_variable_genes, run_pca
from clustering import run_umap_clustering
from visualization import plot_cluster_comparison, plot_scvi_training_history, plot_latent_variance


# ==============================================================================
# 1. DATA CONCATENATION & ANNOTATION
# ==============================================================================

def merge_and_annotate_samples(processed_paths: list[str | Path]) -> ad.AnnData:
    """
    Aggregate individual h5ad files and automatically assign biological metadata.

    This utility reads preprocessed samples, extracts experimental metadata
    (Condition, Replicate, Batch) from the filename structure, and
    merges them using an outer join to preserve the union of all genes.

    Parameters
    ----------
    processed_paths : list of (str or Path)
        A list of paths to the preprocessed .h5ad files.

    Returns
    -------
    ad.AnnData
        A unified AnnData object with unique observation names and
        populated .obs columns for experimental design.

    Notes
    -----
    The function expects filenames to start with specific prefixes:
    - 'N': Normal Tissue
    - 'T': Primary Tumor
    - 'LM': Liver Metastasis
    """
    adatas = []

    for path in processed_paths:
        path_obj = Path(path)
        adata = sc.read_h5ad(path_obj)

        # 1. Filename Parsing
        # Expected format: "Sample_N1_preprocessed.h5ad" -> Stem: "Sample_N1_preprocessed"
        # We split by '_' and find the segment containing the ID
        filename_parts = path_obj.stem.split('_')

        # Logic to find the tag (N1, T1, etc.) regardless of filename length
        tag = next((p for p in filename_parts if p.startswith(('N', 'T', 'LM'))), "Unknown")

        # 2. Assign Condition and Replicate
        if tag.startswith("N"):
            condition, replicate = "Normal", tag[1:]
        elif tag.startswith("T"):
            condition, replicate = "Tumor", tag[1:]
        elif tag.startswith("LM"):
            condition, replicate = "Metastasis", tag[2:]
        else:
            condition, replicate = "Unknown", "0"

        # 3. Metadata Assignment
        # Using .obs.assign or direct assignment is standard
        adata.obs["condition"] = condition
        adata.obs["replicate"] = replicate
        adata.obs["batch"] = f"Patient_{replicate}"

        # Store for merging
        adatas.append(adata)

    print("-" * 64)
    print(f"📊 Merging {len(adatas)} samples into a unified atlas...")

    # 4. Memory-Efficient Concatenation
    # 'join="outer"' ensures no genes are lost if samples have slightly different features
    merged_adata = sc.concat(adatas, join="outer", label="sample_id", index_unique="-")

    # Ensure cell names are unique across the entire study
    merged_adata.obs_names_make_unique()

    print(f"✅ Integration complete. Final shape: {merged_adata.shape}")
    print("-" * 64)

    return merged_adata

# ==============================================================================
# 2. scVI DEEP LEARNING INTEGRATION
# ==============================================================================

def run_scvi_integration(
        adata: ad.AnnData,
        n_layers: int = 2,
        batch_key: str = "batch",
        n_latent: int = 30,
        n_epochs: int = 150,
        model_path: str = "results/integration/scvi_model"
):
    """
    Setup, train, and extract latent representations using the scVI model.

    This function implements a variational autoencoder (VAE) to learn a
    non-linear, low-dimensional representation of the data. It is specifically
    configured to use SoupX-corrected counts to ensure the latent space
    is not biased by ambient RNA contamination.

    Parameters
    ----------
    adata : ad.AnnData
        The integrated AnnData object. Must contain the 'soupX_counts' layer.
    n_layers : int, default 2
        Number of hidden layers in the encoder and decoder networks.
    batch_key : str, default "batch"
        The .obs column identifying different sequencing batches or patients.
    n_latent : int, default 30
        The dimensionality of the latent space (similar to number of PCs).
    n_epochs : int, default 150
        Maximum number of training iterations.
    model_path : str, default "results/integration/scvi_model"
        The directory path for saving or loading the trained model.

    Returns
    -------
    tuple (ad.AnnData, scvi.model.SCVI)
        The AnnData object with 'X_scVI' in .obsm and the trained model instance.
    """
    # Check for existing model to avoid redundant heavy computation
    if Path(f"{model_path}/model.pt").is_file():
        print("-" * 64)
        print(f"📂 Loading existing scVI model from: {model_path}")
        model = scvi.model.SCVI.load(model_path, adata=adata)
        print("Model loaded successfully.")
    else:
        print("-" * 64)
        print("🚀 Initializing scVI Training...")

        # 1. Setup: Define the data architecture for the VAE
        # We use 'soupX_counts' as the raw input for the model's likelihood
        scvi.model.SCVI.setup_anndata(adata, layer="soupX_counts", batch_key=batch_key)

        # 2. Build and Train
        model = scvi.model.SCVI(adata, n_latent=n_latent, n_layers=n_layers)

        print(f"Training for up to {n_epochs} epochs (with early stopping)...")
        # Early stopping prevents overfitting and saves time if convergence is reached
        model.train(max_epochs=n_epochs, early_stopping=True)

        print("✅ Training complete.")

    # 3. Extract Latent Space
    # This matrix is the batch-corrected 'coordinates' for UMAP and Clustering
    adata.obsm["X_scVI"] = model.get_latent_representation()

    print("📊 Latent representations stored in adata.obsm['X_scVI']")
    print("-" * 64)

    return adata, model


def extract_scvi_normalized_counts(
        adata: ad.AnnData,
        model: scvi.model.SCVI,
        reference_batch: str = "S1",
        target_sum: float = 1e4
) -> ad.AnnData:
    """
    Extract batch-corrected, normalized expression from the trained scVI model.

    This function uses the VAE's decoder to project cells back into gene
    expression space while 'standardizing' them to a reference batch. This
    effectively removes technical variation while preserving biological signal.

    Parameters
    ----------
    adata : ad.AnnData
        The integrated AnnData object used to train the scVI model.
    model : scvi.model.SCVI
        The trained scVI model instance.
    reference_batch : str, default "S1"
        The batch ID used as the 'standard' for correction. All cells
        will be projected as if they belonged to this batch.
    target_sum : float, default 1e4
        The library size to normalize to (e.g., 1e4 for CP10k).

    Returns
    -------
    ad.AnnData
        The AnnData object with the 'scvi_normalized' layer added.

    Notes
    -----
    ⚠️ Memory Warning: This operation creates a dense-like matrix in a new
    layer. For very large datasets (e.g., >100k cells), ensure your system
    has sufficient RAM to hold a second full-sized data matrix.
    """
    print("-" * 64)
    print(f"🧪 Extracting scVI-normalized expression...")
    print(f"Target reference batch: {reference_batch}")
    print(f"Library size scaling: {target_sum}")

    # 1. Denoising and Batch Correction
    # 'transform_batch' is the critical argument that performs the
    # counterfactual projection to remove batch effects in gene space.
    try:
        adata.layers["scvi_normalized"] = model.get_normalized_expression(
            transform_batch=reference_batch,
            library_size=target_sum
        )
        print("✅ Correction complete. Layer 'scvi_normalized' is ready.")

    except Exception as e:
        print(f"❌ Error during normalization: {e}")
        raise

    print("-" * 64)

    return adata


# ==============================================================================
# 3. FULL PIPELINE ORCHESTRATION
# ==============================================================================

def run_full_integration_pipeline(
        processed_h5ad_paths: list[str],
        prefix: str,
        n_layers: int = 2,
        reference_batch: str = "S1",
        batch_key: str = "batch",
        n_latent: int = 30,
        n_epochs: int = 150,
        target_sum: float = 1e4,
        n_top_genes: int = 2000,
        resolutions: list[float] = [0.25, 0.5, 0.75, 1],
        color_keys: list[str] = ["condition", "batch"],
        n_neighbors: int = 15,
        model_path: str = "results/integration/scvi_model",
        save_dir_adata: str = "results/integration/",
        save_dir_plots: str = "figures/integration/",
        use_rep: str = "X_scVI",
        save_plots: bool = True
):
    """
    Execute the end-to-end integration workflow for the CRC Atlas.

    Workflow Stages:
    1. Metadata Annotation: Merges sample-specific h5ads and assigns patient/condition tags.
    2. Technical QC: Performs PCA and identifies highly variable genes (HVGs) on raw data.
    3. Deep Learning Integration: Trains scVI to learn a batch-corrected latent space.
    4. Normalization: Extracts denoised, batch-corrected expression layers.
    5. Final Clustering: Computes UMAP and multiple Leiden resolutions on the latent space.
    6. Persistence: Saves the integrated AnnData and the trained scVI model.

    Returns
    -------
    tuple (ad.AnnData, scvi.model.SCVI)
        The finalized, integrated AnnData object and the trained scVI model.
    """

    # --------------------------------------------------------------------------
    # STAGE 1: DATA CONCATENATION & METADATA ANNOTATION
    # --------------------------------------------------------------------------
    # Aggregates N, T, and LM samples into a single object
    merged_adata = merge_and_annotate_samples(processed_h5ad_paths)

    # Use log-normalized data as the standard input for initial PCA
    print("🔄 Swapping to log1p layer for initial dimensionality reduction...")
    merged_adata = swap_adata_layers(merged_adata, "log1p_counts")

    # --------------------------------------------------------------------------
    # STAGE 2: PRE-INTEGRATION ANALYSIS (Technical QC)
    # --------------------------------------------------------------------------
    print(f"🔍 Performing pre-integration check for: {prefix}...")

    # These functions are external to this script but essential for the workflow
    identify_highly_variable_genes(merged_adata, n_top_genes=n_top_genes, batch_key=batch_key)
    run_pca(merged_adata)

    # Proof of Batch Effect: Visual check before correction
    if save_plots:
        plot_cluster_comparison(
            merged_adata,
            color_keys,
            prefix,
            suptitle="Pre-Integration Composition (Batch Effect Presence)",
            save_dir=save_dir_plots
        )

    # Type Safety: Ensure categorical columns from R (scDblFinder) are string-cast
    # to prevent H5AD serialization errors in Python.
    if 'scDblFinder_class' in merged_adata.obs.columns:
        merged_adata.obs['scDblFinder_class'] = merged_adata.obs['scDblFinder_class'].astype(str)

    # Save checkpoint
    save_adata(merged_adata, save_dir_adata / f"{prefix}_pre_integration.h5ad")

    # Subset to Highly Variable Genes for scVI efficiency (improves speed and signal)
    adata_scvi = merged_adata[:, merged_adata.var["highly_variable"]].copy()

    # --------------------------------------------------------------------------
    # STAGE 3: scVI TRAINING (Batch Correction)
    # --------------------------------------------------------------------------
    # Map cells from different patients into a unified latent space (X_scVI)
    integrated_adata, scvi_model = run_scvi_integration(
        adata_scvi,
        n_layers=n_layers,
        batch_key=batch_key,
        n_latent=n_latent,
        n_epochs=n_epochs,
        model_path=model_path
    )

    if save_plots:
        plot_scvi_training_history(scvi_model, prefix, save_dir_plots)
        plot_latent_variance(integrated_adata, prefix, save_dir_plots)

    # --------------------------------------------------------------------------
    # STAGE 4: BATCH-CORRECTED NORMALIZATION
    # --------------------------------------------------------------------------
    # Corrects gene expression to a reference (e.g., 'S1') for clean visualization
    integrated_adata = extract_scvi_normalized_counts(
        integrated_adata,
        scvi_model,
        reference_batch=reference_batch,
        target_sum=target_sum
    )

    # --------------------------------------------------------------------------
    # STAGE 5: FINAL INTEGRATED CLUSTERING
    # --------------------------------------------------------------------------
    print("🗺️ Computing final UMAP based on scVI latent representation...")
    # NOTE: use_rep MUST be 'X_scVI' to utilize the batch-corrected coordinates
    res_list = run_umap_clustering(
        integrated_adata,
        resolutions,
        n_neighbors,
        n_latent,
        use_rep=use_rep
    )

    if save_plots:
        plot_cluster_comparison(
            integrated_adata,
            color_keys,
            prefix,
            suptitle="Post-Integration Composition (scVI Corrected)",
            save_dir=save_dir_plots
        )
        plot_cluster_comparison(
            integrated_adata,
            res_list,
            f"{prefix}_leiden_resolution_optimization",
            "Leiden Resolution Optimization for Cell-Type Identification",
            save_dir_plots
        )

    # --------------------------------------------------------------------------
    # STAGE 6: PERSISTENCE
    # --------------------------------------------------------------------------
    # Save the final Atlas
    save_adata(integrated_adata, save_dir_adata / f"{prefix}_post_integration.h5ad")

    # Persist the model if it was just trained
    if not Path(f"{model_path}/model.pt").is_file():
        save_scvi_model(scvi_model, model_path)

    print("-" * 64)
    print(f"🏁 Full Integration Pipeline Complete for {prefix}")
    print("-" * 64)

    return integrated_adata, scvi_model