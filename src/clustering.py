"""
Clustering and Visualization Module: UMAP & Leiden Workflow.

This module handles the manifold learning and community detection phases of the
pipeline. It constructs a k-nearest neighbor (kNN) graph and optimizes the
Leiden algorithm across multiple resolutions to allow for multi-scale
biological discovery (e.g., broad cell lineages vs. specific sub-states).

Author: Loai Eletr
"""

# --- Bioinformatics Frameworks ---
import scanpy as sc
import anndata as ad


# ==============================================================================
# 1. GRAPH CONSTRUCTION & COMMUNITY DETECTION
# ==============================================================================

def run_umap_clustering(
        adata: ad.AnnData,
        resolutions: list[float],
        n_neighbors: int = 15,
        n_pcs: int = 30,
        use_rep: str = "X_pca"
):
    """
    Computes neighbors, UMAP, and Leiden clustering for multiple resolutions.

    This function transforms high-dimensional data (PCA or scVI latent space)
    into a neighborhood graph, projects it into 2D space via UMAP, and
    segments the graph into clusters using the Leiden algorithm.

    Parameters
    ----------
    adata : ad.AnnData
        The annotated data matrix.
    resolutions : list of float
        A list of scales for the Leiden algorithm (e.g., [0.1, 0.5, 1.0]).
        Lower = fewer/larger clusters; Higher = more/smaller clusters.
    n_neighbors : int, default 15
        Size of the local neighborhood. Lower values preserve more local
        structure, higher values provide a more global view.
    n_pcs : int, default 30
        Number of dimensions (PCs or latent dims) to use for graph construction.
    use_rep : str, default "X_pca"
        The representation to use (e.g., 'X_pca' for single samples,
        'X_scVI' for integrated atlases).

    Returns
    -------
    list of str
        Names of the observation columns (.obs) where clustering results
        were stored.
    """
    print(f"----------------------------------------------------------------")
    print(f"Computing neighbors and UMAP using representation: {use_rep}...")

    # 1. Neighbor Graph Construction
    # This identifies the 'k' most similar cells for every cell in the dataset.
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)

    # 2. UMAP Projection
    # Learning a 2D manifold for visualization. Fixed random_state ensures
    # the 'map' looks the same if re-run.
    sc.tl.umap(adata, random_state=123)

    res_list = []

    # 3. Multi-resolution Clustering
    # We iterate through resolutions to find the 'Goldilocks' zone for
    # cell-type annotation (e.g., separating CD4+ from CD8+ T cells).
    for res in resolutions:
        print(f"  > Clustering at resolution: {res}")

        # 'flavor="igraph"' is the high-performance implementation
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=f"leiden_res_{res}",
            flavor="igraph",
            n_iterations=-1,
            directed=False,
            random_state=123
        )
        res_list.append(f"leiden_res_{res}")

    print(f"✅ Clustering complete. Resolutions stored: {', '.join(res_list)}")
    print(f"----------------------------------------------------------------")

    return res_list