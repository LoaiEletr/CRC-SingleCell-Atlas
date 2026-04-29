"""
Utility Module for scRNA-seq Data Management and Environment Setup.

This module provides helper functions for reproducibility, outlier detection,
R-Python interoperability, and AnnData state management.

Author: Loai Eletr
"""

import os
import random
import numpy as np
import pandas as pd
import torch
import scanpy as sc
import anndata as ad
from pathlib import Path
from scipy.stats import median_abs_deviation

# ==============================================================================
# CONFIGURATION & SETUP
# ==============================================================================

# Global plotting settings
sc.settings.set_figure_params(figsize=(8, 8), dpi=80)

# ==============================================================================
# 1. REPRODUCIBILITY & STATISTICS
# ==============================================================================

def fix_seeds(seed: int = 42):
    """
    Ensure global reproducibility across stochastic libraries.

    This function synchronizes the random state for Python's built-in random
    module, NumPy, PyTorch, and scVI-tools. It is critical for ensuring that
    dimensionality reduction (PCA/UMAP) and deep learning models (scVI)
    produce consistent results across different runs.

    Parameters
    ----------
    seed : int, default 42
        The seed value used to initialize the random number generators.
    """
    # Standard Python and NumPy seeding
    random.seed(seed)
    np.random.seed(seed)

    # PyTorch seeding (CPU and GPU)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    # scVI-tools specific seeding
    scvi.settings.seed = seed

    # Force PyTorch to use deterministic algorithms where possible.
    # 'warn_only=True' prevents crashes if a deterministic version
    # of a specific operation isn't available on the current hardware.
    torch.use_deterministic_algorithms(True, warn_only=True)

    print(f"🔒 Reproducibility lock: Global seed set to {seed}")


def is_outlier(adata: ad.AnnData, metric: str, nmads: int) -> pd.Series:
    """
    Detect outliers using Median Absolute Deviation (MAD).

    This function identifies data points that fall outside a specified
    number of MADs from the median. MAD is used instead of standard
    deviation because it is more robust to extreme outliers, ensuring
    that the threshold is not biased by the very cells we aim to remove.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix.
    metric : str
        The column name in `adata.obs` to be tested for outliers
        (e.g., 'total_counts' or 'pct_counts_mt').
    nmads : int
        The number of Median Absolute Deviations used to set the
        threshold. Typically, a value of 3-5 is used.

    Returns
    -------
    pd.Series
        A boolean mask (True for outliers, False for keep) mapped
        to the original index of `adata.obs`.

    Notes
    -----
    The outlier bounds are calculated as:
    $$[ \text{median} - (\text{nmads} \cdot \text{MAD}), \text{median} + (\text{nmads} \cdot \text{MAD}) ]$$
    """
    # Extract the metric values
    M = adata.obs[metric]

    # Calculate robust statistics
    median_val = np.median(M)
    mad_val = median_abs_deviation(M)

    # Define thresholds
    lower_bound = median_val - (nmads * mad_val)
    upper_bound = median_val + (nmads * mad_val)

    # Flag cells outside the bounds
    outlier_mask = (M < lower_bound) | (M > upper_bound)

    return outlier_mask

# ==============================================================================
# 2. R INTEROPERABILITY (rpy2)
# ==============================================================================

def load_r_environment(r_home_path: str):
    """
    Configure the system environment for R-Python interoperability.

    This function sets the 'R_HOME' environment variable and updates the
    system 'PATH' to include the R binary directory. This is required
    for rpy2 to locate the R executable and shared libraries,
    particularly on Windows systems.

    Parameters
    ----------
    r_home_path : str
        The absolute path to the R installation directory
        (e.g., 'C:/Program Files/R/R-4.x.x').

    Notes
    -----
    This function should be called at the very beginning of the pipeline,
    before any R-dependent packages (like SoupX or scDblFinder) are invoked.
    """
    # 1. Set R_HOME so rpy2 knows where the base installation lives
    os.environ['R_HOME'] = r_home_path

    # 2. Derive the binary path. Using Path handles OS-specific slashes.
    r_bin_path = str(Path(r_home_path) / "bin" / "x64")

    # 3. Prepend to PATH to ensure this R version takes priority
    if r_bin_path not in os.environ['PATH']:
        os.environ['PATH'] = f"{r_bin_path}{os.pathsep}{os.environ['PATH']}"

    print("-" * 64)
    print(f"R environment successfully configured at: {r_home_path}")
    print("-" * 64)

def r_write_console(message: str):
    """
    Redirect R console output to the Python standard output.

    By default, rpy2 does not print R-side console messages (warnings,
    progress bars, or standard print statements) to the Python terminal.
    This function acts as a callback to bridge the two streams.

    Parameters
    ----------
    message : str
        The message string sent from the R interpreter.

    Notes
    -----
    To activate this redirect, you must register it with rpy2:
    >>> from rpy2.rinterface_lib import callbacks
    >>> callbacks.consolewrite_print = r_write_console
    """
    # Using end="" because R's output often includes its own
    # newline characters. Adding another would cause double-spacing.
    print(message, end="")

# ==============================================================================
# 3. ANNDATA LAYER & STATE MANAGEMENT
# ==============================================================================

def swap_adata_layers(adata: ad.AnnData, layer_name: str) -> ad.AnnData:
    """
    Swap the active data matrix (.X) with a specified data layer.

    This utility allows the user to switch the primary data matrix used for
    downstream analysis (e.g., swapping raw counts for SoupX-corrected counts).
    The current content of .X is overwritten by a copy of the target layer.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix containing the target layer.
    layer_name : str
        The key of the layer to be moved into the active .X slot.

    Returns
    -------
    ad.AnnData
        The AnnData object with the updated .X matrix.

    Raises
    ------
    KeyError
        If the specified `layer_name` is not present in `adata.layers`.
    """
    if layer_name not in adata.layers:
        raise KeyError(
            f"❌ Layer '{layer_name}' not found. "
            f"Available layers: {list(adata.layers.keys())}"
        )

    print("-" * 64)
    print(f"🔄 Swapping active .X with layer: '{layer_name}'")

    # We use .copy() to ensure .X is an independent matrix.
    # This prevents accidental 'inplace' modifications to the
    # backup stored in .layers.
    adata.X = adata.layers[layer_name].copy()

    print(f"Swap complete. .X now contains {layer_name} data.")
    print("-" * 64)

    return adata

def add_counts_layer(
    adata: ad.AnnData,
    counts_matrix,
    layer_key: str = "unknown_layer"
) -> ad.AnnData:
    """
    Incorporate a new matrix into the adata.layers slot.

    This utility safely adds a data matrix (e.g., SoupX-corrected counts or
    raw counts) to the AnnData layers. It ensures that the input matrix
    dimensions match the existing AnnData object before assignment.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix receiving the new layer.
    counts_matrix : array_like or sparse matrix
        The matrix to be stored. Must have dimensions (n_obs, n_vars).
    layer_key : str, default "unknown_layer"
        The key under which the matrix will be stored in `adata.layers`.

    Returns
    -------
    ad.AnnData
        The AnnData object with the newly incorporated layer.

    Raises
    ------
    ValueError
        If the shape of the `counts_matrix` does not match `adata.shape`.
    """
    # Defensive check: Matrix dimensions must match AnnData (Cells x Genes)
    if counts_matrix.shape != adata.shape:
        raise ValueError(
            f"❌ Dimension mismatch: Matrix shape {counts_matrix.shape} "
            f"does not match AnnData shape {adata.shape}."
        )

    print("-" * 64)
    print(f"📥 Adding new data layer: '{layer_key}'")

    # Store a copy to prevent downstream 'inplace' modifications
    # to the source matrix from affecting the AnnData layer.
    adata.layers[layer_key] = counts_matrix.copy()

    print(f"Layer '{layer_key}' successfully added.")
    print("-" * 64)

    return adata


# ==============================================================================
# 4. DATA SERIALIZATION (IO)
# ==============================================================================

def save_adata(adata: ad.AnnData, file_path: str | Path):
    """
    Serialize the AnnData object to an HDF5-based .h5ad file.

    This utility ensures the target directory exists before writing and
    applies GZIP compression to minimize disk space usage. It is designed
    to handle the high-dimensional sparse matrices common in single-cell
    projects like the CRC atlas.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data object to be saved.
    file_path : str or Path
        The destination path. Should end with the '.h5ad' extension.

    Notes
    -----
    - Uses 'gzip' compression by default for a balance between speed and size.
    - Automatically creates parent directories if they do not exist.
    """
    path = Path(file_path)

    # 1. Ensure the directory structure exists to prevent FileNotFoundError
    path.parent.mkdir(parents=True, exist_ok=True)

    print("-" * 64)
    print(f"📝 Writing AnnData object to: {path.name}")

    # 2. Perform the write operation
    # Gzip is highly recommended for scRNA-seq matrices to save ~60-80% space
    adata.write(path, compression="gzip")

    # 3. Log success and file size for verification
    if path.exists():
        file_size_mb = path.stat().st_size / (1024 * 1024)
        print(f"✅ Save complete. Final file size: {file_size_mb:.2f} MB")

    print("-" * 64)


def save_scvi_model(
        model: scvi.model.SCVI,
        model_dir: str | Path,
        overwrite: bool = True
):
    """
    Serialize a trained scVI model to a directory.

    This function saves the model weights, latent representation registry,
    and training hyperparameters. This is essential for reproducibility,
    allowing the model to be reloaded later for integration or
    differential expression without re-training.

    Parameters
    ----------
    model : scvi.model.SCVI
        The trained scVI model instance.
    model_dir : str or Path
        The directory path where the model will be stored.
    overwrite : bool, default True
        If True, replaces any existing model at the target directory.

    Notes
    -----
    The output is a directory, not a single file. To reload the model, use:
    >>> model = scvi.model.SCVI.load(model_dir, adata=adata)
    """
    path = Path(model_dir)

    # Ensure the parent directory exists, but scvi manages the
    # creation of the final 'path' folder itself.
    path.parent.mkdir(parents=True, exist_ok=True)

    print("-" * 64)
    print(f"💾 Saving scVI model to directory: {path}")

    try:
        # scVI's internal save method
        model.save(path, overwrite=overwrite)

        if path.exists():
            print(f"✅ Model successfully saved at: {path}")
        else:
            print(f"⚠️ Warning: Model save completed but directory not found.")

    except Exception as e:
        print(f"❌ Error saving scVI model: {e}")
        raise

    print("-" * 64)