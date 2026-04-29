"""
Preprocessing and Quality Control Module for Single-Cell RNA-Seq.

This module provides a unified pipeline for the initial processing of 10x Genomics
data. It integrates standard Python-based analysis (Scanpy/AnnData) with
specialized R-based biological corrections (SoupX/scDblFinder) to ensure
high-quality downstream analysis.

Key Features:
-------------
1. Automated QC: Calculation of mitochondrial, ribosomal, and hemoglobin metrics.
2. Adaptive Filtering: Data-driven cell selection using Median Absolute Deviation (MAD).
3. Ambient RNA Correction: Leveraging SoupX to remove background mRNA contamination.
4. Doublet Detection: Identifying heterotypic doublets via the scDblFinder algorithm.
5. Systematic Dimensionality Reduction: Optimized PCA and HVG selection.

Workflow:
---------
The primary entry point is `run_single_sample_preprocessing`, which orchestrates
the entire pipeline from raw MTX files to a cleaned, normalized AnnData object.

Dependencies:
-------------
- Python: scanpy, anndata, pandas, rpy2
- R: SoupX, scDblFinder, BiocParallel, Seurat, scater

Author: Loai Eletr
"""

# --- Standard Library & Data Science ---
import pandas as pd
from pathlib import Path

# --- Bioinformatics Frameworks ---
import scanpy as sc
import anndata as ad

# --- Custom Project Utilities ---
from src.utils import is_outlier, swap_adata_layers, add_counts_layer
from src.visualization import plot_qc_metrics, plot_pca_variance, plot_cluster_comparison
from src.clustering import run_umap_clustering

# --- R-Python Interoperability (rpy2) ---
import anndata2ri
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.conversion import localconverter


# ==============================================================================
# 1. QUALITY CONTROL & FILTERING
# ==============================================================================
def load_and_annotate_adata(adata_path: str, sample_name: str, prefix: str) -> ad.AnnData:
    """
    Load 10x Genomics MTX data and annotate biological gene groups.

    This function reads raw count data, ensures unique identifiers, and creates
    boolean masks for mitochondrial, ribosomal, and hemoglobin genes.

    Parameters
    ----------
    adata_path : str
        Path to the directory containing the matrix, genes, and barcodes files.
    sample_name : str
        The name/ID associated with this specific sample.
    prefix : str
        Prefix string for the 10x files (e.g., 'filtered_feature_bc_matrix_').

    Returns
    -------
    ad.AnnData
        Annotated AnnData object with gene categories in `.var` and raw
        counts preserved in the 'counts' layer.
    """
    # Load the 10x sparse matrix into an AnnData object
    adata = sc.read_10x_mtx(adata_path, var_names='gene_symbols', cache=True, prefix=prefix)

    # Ensure gene names are unique to avoid downstream indexing issues
    adata.var_names_make_unique()

    # Preserve a copy of raw counts before normalization/transformation
    adata.layers["counts"] = adata.X.copy()

    # --- Gene Group Annotation ---
    # Mitochondrial genes: Typically start with 'MT-' in human gene symbols
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    # Ribosomal genes: RPS (Small) and RPL (Large) subunit proteins
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    # Hemoglobin genes: Matches 'HB' but excludes pseudogenes (ending in P)
    # to avoid noise in red blood cell contamination checks
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    return adata


def calculate_sample_qc(adata: ad.AnnData) -> ad.AnnData:
    """
    Calculate comprehensive Quality Control (QC) metrics for the dataset.

    Computes statistics including total counts per cell, genes per cell,
    and the percentage of reads mapped to MT, Ribosomal, and HB genes.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix. Must contain 'mt', 'ribo', and 'hb'
        columns in `.var`.

    Returns
    -------
    ad.AnnData
        The input object with calculated metrics populated in `.obs` and `.var`.

    See Also
    --------
    scanpy.pp.calculate_qc_metrics : The underlying function used for computation.
    """
    print("-" * 64)
    print("Calculating standard QC metrics for MT, Ribo, and HB genes...")

    # Calculate metrics.
    # percent_top=[20] identifies the fraction of total counts
    # occupied by the top 20 most highly expressed genes (complexity check).
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        inplace=True,
        percent_top=[20],
        log1p=True
    )

    return adata

def filter_adata(
        adata: ad.AnnData,
        mad_usage: bool = True,
        mad_cutoff: int = 3,
        total_counts_bounds: tuple = None,
        n_genes_by_counts_bounds: tuple = None,
        pct_counts_mt_upper: float = None
) -> ad.AnnData:
    """
    Filter cells based on Median Absolute Deviation (MAD) or manual thresholds.

    This function identifies and removes low-quality cells by evaluating library
    size (total counts), gene richness (n_genes), and mitochondrial content.
    It supports both data-adaptive filtering (MAD) and fixed hard cutoffs.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix containing QC metrics in `.obs`.
    mad_usage : bool, default True
        If True, use Median Absolute Deviation to dynamically find outliers based
        on the sample's distribution. If False, use fixed thresholds.
    mad_cutoff : int, default 3
        The number of MADs away from the median to consider a cell an outlier.
        Standard practice is 3-5.
    total_counts_bounds : tuple, optional
        Manual (lower, upper) bounds for total UMI counts. Required if mad_usage=False.
    n_genes_by_counts_bounds : tuple, optional
        Manual (lower, upper) bounds for number of detected genes. Required if mad_usage=False.
    pct_counts_mt_upper : float, optional
        Maximum allowable percentage of mitochondrial counts. Required if mad_usage=False.
        If mad_usage=True, this acts as an additional hard ceiling.

    Returns
    -------
    ad.AnnData
        A filtered copy of the original AnnData object.

    Raises
    ------
    ValueError
        If manual filtering is selected but required bounds are missing.
    """

    if mad_usage:
        print("-" * 64)
        print(f"Filtering via MAD (cutoff={mad_cutoff})...")

        # Flag outliers for technical metrics.
        # Log-transformation is used here because sequencing depth and gene counts
        # typically follow a log-normal distribution.
        adata.obs["outlier"] = (
                is_outlier(adata, "log1p_total_counts", mad_cutoff) |
                is_outlier(adata, "log1p_n_genes_by_counts", mad_cutoff) |
                is_outlier(adata, "pct_counts_in_top_20_genes", mad_cutoff)
        )

        # MT filtering: Combine adaptive MAD filtering with a hard ceiling if provided.
        # This prevents 'high quality' cells from being kept if they have biologically
        # impossible MT percentages.
        mt_mad = is_outlier(adata, "pct_counts_mt", mad_cutoff)
        if pct_counts_mt_upper is not None:
            adata.obs["mt_outlier"] = mt_mad | (adata.obs["pct_counts_mt"] > pct_counts_mt_upper)
        else:
            adata.obs["mt_outlier"] = mt_mad

    else:
        # --- Manual Threshold Logic ---
        if not (total_counts_bounds and n_genes_by_counts_bounds and pct_counts_mt_upper is not None):
            raise ValueError(
                "Manual filtering mode requires: total_counts_bounds, "
                "n_genes_by_counts_bounds, and pct_counts_mt_upper."
            )

        print("-" * 64)
        print(f"Filtering via Manual Thresholds: Counts={total_counts_bounds}, "
              f"Genes={n_genes_by_counts_bounds}, MT% < {pct_counts_mt_upper}")

        tc_low, tc_high = total_counts_bounds
        ng_low, ng_high = n_genes_by_counts_bounds

        adata.obs["outlier"] = (
                (adata.obs["total_counts"] < tc_low) | (adata.obs["total_counts"] > tc_high) |
                (adata.obs["n_genes_by_counts"] < ng_low) | (adata.obs["n_genes_by_counts"] > ng_high)
        )
        adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > pct_counts_mt_upper

    # --- Filtering Execution ---
    n_cells_before = adata.n_obs

    # Subset to keep only non-outliers
    filtered_adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    n_cells_after = filtered_adata.n_obs

    print(f"Cells before filtering: {n_cells_before}")
    print(f"Cells after filtering:  {n_cells_after}")
    print(f"Removed {n_cells_before - n_cells_after} cells ({(1 - n_cells_after / n_cells_before):.1%})")
    print("-" * 64)

    return filtered_adata

# ==============================================================================
# 2. NORMALIZATION & DIMENSIONALITY REDUCTION
# ==============================================================================

def normalize_adata(adata: ad.AnnData, target_sum: float = 1e4) -> ad.AnnData:
    """
    Perform total count normalization and log1p transformation.

    This function scales the counts in each cell so that the total counts
    per cell sum to a specific target value (defaulting to 10,000, or 'CP10k').
    It then applies a natural logarithm transformation to stabilize variance.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix. The active layer (.X) should contain
        non-normalized counts (e.g., ambient-removed counts).
    target_sum : float, default 1e4
        The number each cell's total counts will sum to after normalization.
        1e4 is the industry standard for single-cell data.

    Returns
    -------
    ad.AnnData
        The AnnData object with normalized and log-transformed data in `.X`.
        The transformed data is also backed up in the 'log1p_counts' layer.

    Notes
    -----
    The calculation performed is:
    $$X_{norm} = \ln\left(\frac{X}{\sum X} \cdot \text{target\_sum} + 1\right)$$
    """
    # Clean up old log1p metadata if it exists to prevent Scanpy
    # from assuming the data is already transformed.
    if 'log1p' in adata.uns:
        del adata.uns['log1p']

    print("-" * 64)
    print(f"Normalizing total counts to {target_sum}...")

    # Scale counts to target_sum
    sc.pp.normalize_total(adata, target_sum=target_sum)

    print("Applying log1p transformation...")
    # Log transformation: x = log(x + 1)
    sc.pp.log1p(adata)

    # Store the log1p data in a dedicated layer.
    # This ensures that even if .X is later scaled (mean=0, std=1),
    # the unscaled log-normalized counts remain accessible.
    adata.layers["log1p_counts"] = adata.X.copy()

    print("Normalization complete: .X and layers['log1p_counts'] now contain log1p data.")
    print("-" * 64)

    return adata



def identify_highly_variable_genes(
    adata: ad.AnnData,
    n_top_genes: int = 2000,
    batch_key: str = None
) -> ad.AnnData:
    """
    Identify highly variable genes (HVGs) for dimensionality reduction.

    This function flags genes that contribute the most to the biological
    variance of the dataset. By filtering for HVGs, we reduce noise from
    constitutively expressed genes and improve the performance of PCA.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix. .X should contain log-normalized counts.
    n_top_genes : int, default 2000
        The number of genes to be flagged as highly variable.
    batch_key : str, optional
        If provided, HVGs are selected within each batch and merged. This
        prevents batch-specific effects from dominating the gene selection.

    Returns
    -------
    ad.AnnData
        The input object with HVG metadata added to `.var` and `.uns`.

    Notes
    -----
    The function uses the 'seurat' flavor by default (expecting log-normalized
    data). If `batch_key` is used, the function ensures that the top genes
    are represented across all batches to avoid picking up technical artifacts.
    """
    print("-" * 64)
    print(f"Identifying top {n_top_genes} HVGs...")
    if batch_key:
        print(f"Applying batch-aware selection using key: '{batch_key}'")

    # Run Scanpy's HVG selection
    # 'seurat' flavor is standard for log-normalized data.
    # If flavor='seurat_v3' is desired, it would require raw counts in .X.
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        batch_key=batch_key,
        subset=False # We keep all genes, only flagging HVGs in .var['highly_variable']
    )

    # Log summary statistics
    hvg_count = adata.var['highly_variable'].sum()
    print(f"Successfully flagged {hvg_count} highly variable genes.")
    print("-" * 64)

    return adata


def run_pca(adata: ad.AnnData) -> ad.AnnData:
    """
    Compute Principal Component Analysis (PCA) on highly variable genes.

    PCA is used to reduce the dimensionality of the gene expression space
    into a set of linear combinations (Principal Components) that capture
    the maximum variance in the data.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix. Expected to have log-normalized counts
        and highly variable genes (HVGs) already identified.

    Returns
    -------
    ad.AnnData
        The input object with PCA results stored in `.obsm['X_pca']`
        and loadings in `.varm['PCs']`.
    """
    print("-" * 64)
    print("Computing Principal Component Analysis (PCA)...")

    # Run PCA using the ARPACK solver for sparse matrix efficiency.
    # This automatically uses genes flagged in .var['highly_variable'].
    sc.pp.pca(adata, svd_solver="arpack")

    print(f"PCA computation complete. Results stored in .obsm['X_pca'].")
    print("-" * 64)

    return adata


# ==============================================================================
# 3. BIOLOGICAL CORRECTION (R-INTEGRATION)
# ==============================================================================

def ambient_removal(adata: ad.AnnData, soupx_groups: pd.Series) -> ad.AnnData:
    """
    Remove ambient mRNA contamination using the SoupX R package via rpy2.

    This function serves as a Python wrapper for the SoupX algorithm. It
    estimates the background "soup" of mRNA and subtracts it from the
    raw counts to account for ambient contamination in droplet-based
    single-cell RNA sequencing.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix. Should contain raw UMI counts.
    soupx_groups : pd.Series
        Cluster or group assignments (e.g., from Leiden clustering).
        SoupX uses these groups to estimate the contamination fraction.

    Returns
    -------
    ad.AnnData
        An AnnData object where the active `.X` layer is now the
        background-corrected "soupX_counts".

    Notes
    -----
    - This function requires an active R installation with the 'SoupX'
      library installed.
    - Data is automatically transposed during the transfer to R, as R
      expects a [genes x cells] matrix format.
    """
    # SoupX must be run on raw counts. We swap to the raw 'counts' layer
    # to ensure normalized data isn't accidentally processed.
    print("Swapping to 'counts' layer for ambient estimation...")
    adata = swap_adata_layers(adata, "counts")

    print("-" * 64)
    print("Running SoupX ambient RNA removal...")

    # Utilize localconverter for seamless AnnData/Pandas to R conversion
    with localconverter(anndata2ri.converter + pandas2ri.converter):
        # Assign Python variables to the R global environment
        r.assign("genes", adata.var_names)
        r.assign("cells", adata.obs_names)
        r.assign("data_tod", adata.X.T)  # Transpose for R (Genes x Cells)
        r.assign("soupx_groups", soupx_groups)

        r('''
            suppressPackageStartupMessages(library("SoupX"))

            # Reconstruct the sparse matrix in R with proper identifiers
            rownames(data_tod) = genes
            colnames(data_tod) = cells
            data_tod <- as(data_tod, "sparseMatrix")

            # Initialize SoupChannel
            # We use data_tod for both parameters in this context, 
            # effectively assuming the 'table of droplets' is the input.
            sc = SoupChannel(data_tod, data_tod, calcSoupProfile = FALSE)

            # Manually calculate and set the soup profile
            soupProf = data.frame(
                row.names = rownames(data_tod), 
                est = rowSums(data_tod)/sum(data_tod), 
                counts = rowSums(data_tod)
            )
            sc = setSoupProfile(sc, soupProf)
            sc = setClusters(sc, soupx_groups)

            # Automatically estimate contamination fraction
            sc = autoEstCont(sc, doPlot=FALSE)

            # Adjust the counts and round to integers for downstream tools
            corrected_matrix = adjustCounts(sc, roundToInt = TRUE)
        ''')
        # Retrieve the corrected matrix back into the Python environment
        corrected_matrix = r['corrected_matrix']

    # Add the result as a new layer. Note: matrix must be transposed back to (Cells x Genes)
    print("Adding soupX_counts layer to AnnData...")
    adata = add_counts_layer(adata, corrected_matrix.T, "soupX_counts")

    # Set the corrected counts as the active X for subsequent analysis
    print("Activating soupX_counts layer...")
    adata = swap_adata_layers(adata, "soupX_counts")

    print("Ambient RNA removal complete.")
    print("-" * 64)

    return adata


def doublet_identification(adata: ad.AnnData) -> ad.AnnData:
    """
    Identify droplet doublets using scDblFinder via rpy2.

    This function wraps the Bioconductor package `scDblFinder` to detect
    heterotypic doublets (two different cell types in one droplet). It
    converts the AnnData object to an R SingleCellExperiment object,
    runs the detection algorithm, and maps the results back to `.obs`.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix. Ensure the active layer (.X) contains
        raw or soup-corrected counts (non-log transformed).

    Returns
    -------
    ad.AnnData
        The input AnnData object with two new columns in `.obs`:
        - 'scDblFinder_score': Continuous probability score of being a doublet.
        - 'scDblFinder_class': Categorical classification ('singlet' or 'doublet').

    Notes
    -----
    - This function requires R packages: `scDblFinder`, `scater`, and `Seurat`.
    - A fixed seed (123) is set within the R environment for reproducibility,
      as scDblFinder involves random artificial doublet generation.
    """
    print("-" * 64)
    print("Detection of droplet doublets using scDblFinder...")

    # Use localconverter for AnnData -> SingleCellExperiment transition
    with localconverter(anndata2ri.converter):
        # Transpose to (Genes x Cells) for R/Bioconductor compatibility
        r.assign("data_mat", adata.X.T)

        r('''
            suppressPackageStartupMessages(library(Seurat))
            suppressPackageStartupMessages(library(scater))
            suppressPackageStartupMessages(library(scDblFinder))
            suppressPackageStartupMessages(library(BiocParallel))

            # Set seed for reproducible doublet simulation
            set.seed(123)

            # Initialize SCE object and run scDblFinder
            # By default, this uses the 'cluster-based' approach to simulate doublets
            sce = scDblFinder(SingleCellExperiment(list(counts=data_mat)))

            # Extract metrics
            doublet_score = sce$scDblFinder.score
            doublet_class = sce$scDblFinder.class
        ''')
        # Pull R vectors back into Python
        doublet_score = r['doublet_score']
        doublet_class = r['doublet_class']

    # Assign results to metadata
    adata.obs['scDblFinder_score'] = doublet_score
    adata.obs['scDblFinder_class'] = doublet_class

    # Quick summary of results for the user
    num_doublets = (adata.obs['scDblFinder_class'] == 2).sum()
    print(f"Identification complete: Found {num_doublets} doublets.")
    print("-" * 64)

    return adata


def doublet_removal(adata: ad.AnnData) -> ad.AnnData:
    """
    Remove identified droplet doublets from the AnnData object.

    This function filters the dataset to retain only cells classified as
    'singlets' (integer 1) by the scDblFinder R factor conversion.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix. Must contain the 'scDblFinder_class'
        column in `.obs` populated by the identification step.

    Returns
    -------
    ad.AnnData
        A filtered, memory-independent copy of the AnnData object containing
        only singlets.

    Notes
    -----
    In the R-to-Python conversion of scDblFinder factors:
    - 1 = 'singlet'
    - 2 = 'doublet'
    """
    # Defensive check to ensure column exists
    if "scDblFinder_class" not in adata.obs.columns:
        raise KeyError(
            "Column 'scDblFinder_class' not found in adata.obs. "
            "Run 'doublet_identification' first."
        )

    print("-" * 64)
    print("Removing droplet doublets (class == 1)...")

    # Record metrics for logging
    n_before = adata.n_obs

    # Subset to singlets (Factor level 1)
    # Using .copy() is essential to break the view and prevent
    # ImplicitModificationWarnings later.
    adata = adata[adata.obs["scDblFinder_class"] == 1].copy()

    n_after = adata.n_obs
    print(f"Removed {n_before - n_after} doublets. {n_after} singlets remaining.")

    # Ensure we are working with corrected counts for downstream steps
    print("Swapping to 'soupX_counts' layer...")
    adata = swap_adata_layers(adata, "soupX_counts")

    print("-" * 64)

    return adata

# ==============================================================================
# 4. PIPELINE ORCHESTRATION (WRAPPER)
# ==============================================================================

def run_single_sample_preprocessing(
        adata_path: str,
        sample_name: str,
        prefix: str,
        mad_usage: bool = True,
        total_counts_bounds: tuple = None,
        n_genes_by_counts_bounds: tuple = None,
        pct_counts_mt_upper: float = None,
        mad_cutoff: int = 3,
        target_sum: float = 1e4,
        n_top_genes: int = 2000,
        n_pcs: int = 50,
        resolutions: list[float] = [0.25, 0.5, 0.75, 1],
        resolution_for_soupx: float = 0.5,
        batch_key: str = None,
        save_plots: bool = True,
        n_neighbors: int = 15,
        use_rep: str = "X_pca",
        save_dir_qc: str = "figures/qc",
        save_dir_clustering: str = "figures/clustering",
) -> ad.AnnData:
    """
    Orchestrate the end-to-end preprocessing pipeline for a single 10x sample.

    This wrapper automates the standard single-cell preprocessing trajectory:
    1. Data loading and QC metric calculation.
    2. Adaptive (MAD) or manual cell filtering.
    3. Technical normalization and dimensionality reduction.
    4. Ambient RNA estimation and removal (SoupX).
    5. Doublet detection and exclusion (scDblFinder).

    Parameters
    ----------
    adata_path : str
        Directory containing the 10x MTX files.
    sample_name : str
        Identifier for the sample, used in plot titles and metadata.
    prefix : str
        File prefix for the 10x matrix files.
    mad_usage : bool, default True
        Whether to use Median Absolute Deviation for outlier detection.
    total_counts_bounds : tuple, optional
        (min, max) UMI counts for manual filtering.
    n_genes_by_counts_bounds : tuple, optional
        (min, max) genes per cell for manual filtering.
    pct_counts_mt_upper : float, optional
        Upper threshold for mitochondrial percentage.
    mad_cutoff : int, default 3
        Number of MADs for adaptive filtering.
    target_sum : float, default 1e4
        Scaling factor for library size normalization.
    n_top_genes : int, default 2000
        Number of highly variable genes to select.
    n_pcs : int, default 50
        Number of principal components to calculate.
    resolutions : list of float, default [0.25, 0.5, 0.75, 1]
        Leiden clustering resolutions to evaluate.
    resolution_for_soupx : float, default 0.5
        The specific resolution used to define clusters for SoupX estimation.
    batch_key : str, optional
        Column in `.var` used for batch-aware HVG selection.
    save_plots : bool, default True
        If True, exports QC and clustering visualizations to disk.
    n_neighbors : int, default 15
        Size of the local neighborhood for graph construction.
    use_rep : str, default "X_pca"
        Representation used for UMAP and clustering.
    save_dir_qc : str, default "figures/qc"
        Output directory for quality control plots.
    save_dir_clustering : str, default "figures/clustering"
        Output directory for clustering and doublet plots.

    Returns
    -------
    ad.AnnData
        A fully processed and cleaned AnnData object.
    """

    # --- 1. LOAD & ANNOTATE ---
    # Initial ingestion and basic biological group flagging (MT, Ribo, HB)
    adata = load_and_annotate_adata(adata_path, sample_name, prefix)
    adata = calculate_sample_qc(adata)

    # --- 2. VISUALIZE: Raw Data State ---
    if save_plots:
        plot_qc_metrics(adata, sample_name, "Initial QC", save_dir_qc)

    # --- 3. FILTERING (MAD-based or Manual) ---
    # Removes low-quality droplets and potential debris
    adata = filter_adata(
        adata,
        mad_usage=mad_usage,
        mad_cutoff=mad_cutoff,
        total_counts_bounds=total_counts_bounds,
        n_genes_by_counts_bounds=n_genes_by_counts_bounds,
        pct_counts_mt_upper=pct_counts_mt_upper
    )

    # --- 4. VISUALIZE: Post-Cell Filtering ---
    if save_plots:
        plot_qc_metrics(adata, sample_name, "Post-Filtering", save_dir_qc)

    # --- 5. TECHNICAL PCA & CLUSTERING ---
    # Essential for SoupX/scDblFinder which require coarse cluster assignments
    print(f"Running technical dimension reduction for {sample_name}...")
    normalize_adata(adata, target_sum)
    identify_highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
    run_pca(adata)

    if save_plots:
        plot_pca_variance(adata, sample_name, n_pcs, save_dir_qc)

    # Generate UMAP and multiple Leiden resolutions
    res_list = run_umap_clustering(adata, resolutions, n_neighbors, n_pcs, use_rep)

    if save_plots:
        plot_cluster_comparison(
            adata,
            res_list,
            f"{sample_name}_leiden_resolution_optimization",
            "Optimization of Leiden Clustering Resolutions",
            save_dir_clustering
        )

    # --- 6 & 7. BIOLOGICAL CORRECTIONS ---
    # Statistical corrections (SoupX and scDblFinder) require a minimum cell count for stability
    if adata.n_obs >= 100:
        # AMBIENT RNA REMOVAL (SoupX)
        soupx_key = f"leiden_res_{resolution_for_soupx}"
        adata = ambient_removal(adata, adata.obs[soupx_key])

        # DOUBLET IDENTIFICATION & REMOVAL (scDblFinder)
        adata = doublet_identification(adata)

        if save_plots:
            plot_cluster_comparison(
                adata,
                ["scDblFinder_class", soupx_key],
                f"{sample_name}_doublet_qc",
                "Doublet Distribution and Cluster Enrichment Analysis",
                save_dir_clustering
            )

        adata = doublet_removal(adata)
    else:
        print(f"⚠️ SKIPPING SoupX/Doublet removal for {sample_name}: Only {adata.n_obs} cells.")

        # Ensure metadata consistency for samples with low cell counts
        adata.layers["soupX_counts"] = adata.X.copy()
        adata.obs["scDblFinder_class"] = "skipped"
        adata.obs["scDblFinder_score"] = 0.0
        adata.uns["preprocessing_note"] = "Skipped SoupX/Doublet detection due to low cell count."

    # --- 8. VISUALIZE: Final Processed State ---
    if save_plots:
        plot_qc_metrics(adata, sample_name, "Final Cleaned", save_dir_qc)

    print("-" * 64)
    print(f"Successfully processed {sample_name}. Final cell count: {adata.n_obs}")
    print("-" * 64)

    return adata
