# Complete Computational Methods: CRC-SingleCell-Atlas

---

## Table of Contents

1. [Dataset Description](#1-dataset-description)
2. [Per-Sample Preprocessing](#2-per-sample-preprocessing) 
3. [Ambient RNA Removal (SoupX)](#3-ambient-rna-removal-soupx)
4. [Doublet Detection and Removal](#4-doublet-detection-and-removal)
5. [Batch Integration with scVI](#5-batch-integration-with-scvi)
6. [Clustering and Marker Identification](#6-clustering-and-marker-identification)
7. [Cell Type Annotation](#7-cell-type-annotation)
8. [Pseudobulk Differential Expression](#8-pseudobulk-differential-expression)
9. [Gene Set Enrichment Analysis](#9-gene-set-enrichment-analysis)
10. [Cell-Cell Communication Analysis](#10-cell-cell-communication-analysis)
11. [Differential Communication Analysis](#11-differential-communication-analysis)
12. [Reproducibility](#12-reproducibility)

---

## 1. Dataset Description

### Study Design
Single-cell RNA sequencing was performed on 15 samples from 6 colorectal cancer patients (Sun et al., 2026, GEO: [GSE315534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315534)):



| Tissue Type | Count | Sample IDs |
|-------------|-------|-------------|
| Adjacent normal colon | 6 | N1-N6 |
| Primary colorectal tumor | 6 | T1-T6 |
| Paired liver metastasis | 3 | LM1-LM3 (patients 1-3) |

**Patient matching:** LM1 corresponds to (T1,N1), LM2 to (T2,N2), LM3 to (T3,N3).

### Platform
- **Instrument:** Illumina NovaSeq 6000
- **Platform ID:** GPL24676
- **Chemistry:** 10x Genomics

---

## 2. Per-Sample Preprocessing

**Code:** `scripts/01_preprocessing.py` → `src/preprocessing.py`  
**Output:** `results/01_preprocessing/{sample}/{sample}_cleaned.h5ad`  
**Logs:** `results/01_preprocessing/{sample}/{sample}_processing_log.txt`  
**Summary:** `results/01_preprocessing/Summary.csv`

### 2.1 Data Loading

```python
adata = sc.read_10x_mtx(adata_path, var_names='gene_symbols', cache=True)
adata.layers["counts"] = adata.X.copy()  # Store raw counts
```

Gene annotations are added:

- **Mitochondrial:** `var_names.startswith("MT-")`
- **Ribosomal:** `var_names.startswith(("RPS", "RPL"))`
- **Hemoglobin:** `var_names.contains("^HB[^(P)]")`

### 2.2 Quality Control Metrics

For each cell, the following metrics are calculated (via `sc.pp.calculate_qc_metrics`):

- `total_counts` (nUMI)
- `n_genes_by_counts` (nFeature)
- `pct_counts_mt` (mitochondrial percentage)
- `pct_counts_ribo`, `pct_counts_hb`

### 2.3 Dual Filtering Strategy

Users select strategy in `config.py` via `PreprocessingConfig.QC_PARAMS['mad_usage']`:

**Strategy A: MAD-based (Median Absolute Deviation)**

```python
def is_outlier(adata, metric, nmads=3):
    M = adata.obs[metric]
    median_val = np.median(M)
    mad_val = median_abs_deviation(M)
    return (M < median_val - nmads * mad_val) | (M > median_val + nmads * mad_val)

# Outlier criteria
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 3) |
    is_outlier(adata, "log1p_n_genes_by_counts", 3) |
    is_outlier(adata, "pct_counts_in_top_20_genes", 3)
)
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 12)
```

**Strategy B: Threshold-based**

```python
adata.obs["outlier"] = (
    (adata.obs["total_counts"] < 500) | (adata.obs["total_counts"] > 50000) |
    (adata.obs["n_genes_by_counts"] < 200) | (adata.obs["n_genes_by_counts"] > 6000)
)
adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > 20
```

**Cells retained:** `(~adata.obs.outlier) & (~adata.obs.mt_outlier)`

### 2.4 Normalization and Technical PCA

After filtering:

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p_counts"] = adata.X.copy()  # Store for later
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, svd_solver="arpack")
```

### 2.5 Clustering for Downstream QC

To identify the optimal clustering granularity, we evaluated the Leiden algorithm at incremental resolutions (0.25 to 1.0). Selecting the most representative partition allowed us to:

1. Estimate cluster structure for SoupX ambient RNA estimation
2. Identify potential doublet-enriched clusters

```python
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca")
sc.tl.leiden(adata, resolution=res, key_added=f"leiden_res_{res}")
```

**Output figures:** `figures/clustering/{sample}/umap_{sample}_leiden_resolution_optimization_*.pdf`

---

## 3. Ambient RNA Removal (SoupX)

**Code:** `src/preprocessing.py` → `ambient_removal()`  
**R dependency:** SoupX (via rpy2)

### 3.1 Rationale

Ambient RNA from lysed cells contaminates droplet-based scRNA-seq, creating false positive signals. SoupX estimates contamination fraction using empty droplets or, in our implementation, clusters of low-quality cells [^1].

### 3.2 Implementation (R via rpy2)

```r
# Create SoupChannel
sc = SoupChannel(data_tod, data_tod, calcSoupProfile = FALSE)

# Set soup profile from total counts
soupProf = data.frame(row.names = rownames(data_tod), 
                      est = rowSums(data_tod)/sum(data_tod), 
                      counts = rowSums(data_tod))
sc = setSoupProfile(sc, soupProf)

# Assign clusters (from user-selected Leiden resolution optimal for the dataset)
sc = setClusters(sc, soupx_groups)

# Automatically estimate contamination fraction
sc = autoEstCont(sc, doPlot=FALSE)

# Apply correction
corrected_matrix = adjustCounts(sc, roundToInt = TRUE)
```

### 3.3 Output

- Corrected counts stored in `adata.layers["soupX_counts"]`
- Estimator `rho` (contamination fraction) logged in processing log
- For samples with <100 cells, SoupX is skipped (logged as "SKIPPED")

**QC visualization:** `figures/qc/{sample}/violin_{sample}_final_cleaned.pdf` (post-correction)

---

## 4. Doublet Detection and Removal

**Code:** `src/preprocessing.py` → `doublet_identification()`, `doublet_removal()`  
**R dependency:** scDblFinder (via rpy2)

### 4.1 Rationale

Doublets (two cells captured in one droplet) create artificial hybrid transcriptomes. We used scDblFinder to identify and remove doublets [^2].

### 4.2 Implementation

```r
library(scDblFinder)
set.seed(123)
sce = scDblFinder(SingleCellExperiment(list(counts=data_mat)))
doublet_score = sce$scDblFinder.score
doublet_class = sce$scDblFinder.class  # 1 -> "singlet" or  2 -> "doublet"
```

### 4.3 Quality Control

Doublet enrichment per cluster is visualized:

```python
plot_cluster_comparison(adata, ["scDblFinder_class", soupx_key], 
                        f"{sample_name}_doublet_qc",
                        "Doublet Distribution and Cluster Enrichment Analysis")
```

**Output:** `figures/clustering/{sample}/umap_{sample}_doublet_qc_*.pdf`

### 4.4 Doublet Removal

Cells classified as `"doublet"` are removed:

```python
adata = adata[adata.obs["scDblFinder_class"] == 1].copy()
```

### 4.5 Low-Cell-Count Samples

For samples with <100 cells, doublet detection is skipped:

- Placeholder values added: `scDblFinder_class = "skipped"`, `scDblFinder_score = 0.0`
- Note stored in `adata.uns["preprocessing_note"]`

---

## 5. Batch Integration with scVI

**Code:** `scripts/02_integration.py` → `src/integration.py`  
 **Output:**

- `results/02_integration/merged_dataset_pre_integration.h5ad`
- `results/02_integration/merged_dataset_post_integration.h5ad`
- `results/02_integration/scvi_model/` (saved model)

### 5.1 Data Merging

```python
def merge_and_annotate_samples(processed_paths):
    # Extract metadata from filename (e.g., "N1_cleaned.h5ad" → tag="N1")
    if tag.startswith("N"): condition, replicate = "Normal", tag[1:]
    elif tag.startswith("T"): condition, replicate = "Tumor", tag[1:]
    elif tag.startswith("LM"): condition, replicate = "Metastasis", tag[2:]
    
    adata.obs["condition"] = condition
    adata.obs["replicate"] = replicate
    adata.obs["batch"] = f"S{replicate}"  # Patient ID as batch
```

**Key design:** `batch_key = "batch"` (patient ID), NOT condition, ensuring biological differences are preserved.

### 5.2 Preprocessing for scVI

```python
sc.pp.highly_variable_genes(
    adata, 
    n_top_genes=2000, 
    batch_key="batch"
)

# Subset to highly variable genes for scVI efficiency
adata = adata[:, adata.var["highly_variable"]].copy()
```

### 5.3 scVI Model Configuration

Batch integration was performed using scVI [^3], a deep generative model that learns a low-dimensional latent representation while removing batch effects.

```python
scvi.model.SCVI.setup_anndata(adata, layer="soupX_counts", batch_key="batch")

model = scvi.model.SCVI(
    adata,
    n_latent=30,      # Latent dimensionality
    n_layers=2,       # Neural network depth
)
```

### 5.4 Model Training

```python
model.train(
    max_epochs=150,           # Configurable
    early_stopping=True
)

# Extract corrected latent representation
adata.obsm["X_scVI"] = model.get_latent_representation()

# Extract batch-corrected normalized counts (reference batch = "S1")
adata.layers["scvi_normalized"] = model.get_normalized_expression(
    transform_batch="S1", library_size=1e4
)
```

### 5.5 Quality Control Visualizations

| Figure | Location | Purpose |
| --- | --- | --- |
| Pre-integration UMAP | `figures/integration/umap_*_Pre-Integration Composition.pdf` | Visualize batch effects |
| Post-integration UMAP | `figures/integration/umap_*_Post-Integration Composition.pdf` | Assess correction |
| Training history | `figures/integration/*_scvi_training_history.pdf` | Check convergence |
| Latent variance | `figures/integration/*_latent_variance.pdf` | Inspect latent dimensions |
| Leiden optimization | `figures/integration/umap_*_leiden_resolution_optimization_*.pdf` | Choose resolution |

### 5.6 Model Persistence

The trained scVI model is saved to `results/02_integration/scvi_model/` for:

- Reproducibility (reuse without retraining)
- Downstream analyses requiring batch-corrected counts

---

## 6. Clustering and Marker Identification

**Code:** `scripts/03_cluster_characterization.py` → `src/deg.py`, `src/clustering.py`  
 **Output:** `results/03_cluster_characterization/`

### 6.1 Final Clustering (Post-Integration)

Using the scVI latent representation:

```python
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=50)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_res_0.5")
sc.tl.umap(adata)
```

### 6.2 Marker Gene Identification

**Method:** Wilcoxon rank-sum test (Scanpy default)

```python
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden_res_0.5",
    method="wilcoxon",
    key_added="dea_leiden_res_0.5"
)
```

**Filtering:** Genes meeting criteria (min_in_group_fraction ≥ 0.25, max_out_group_fraction ≤ 0.50):

```python
sc.tl.filter_rank_genes_groups(
    adata,
    groupby="leiden_res_0.5",
    min_in_group_fraction=0.25,
    max_out_group_fraction=0.50,
    key="dea_leiden_res_0.5",
    key_added="dea_leiden_res_0.5_filtered"
)
```

### 6.3 Output Files

| File | Content |
| --- | --- |
| `*_all_markers_full.csv` | Complete Wilcoxon results (all genes) |
| `*_top20_markers.csv` | Top 20 genes per cluster (unfiltered) |
| `*_filtered_all_markers_full.csv` | Filtered Wilcoxon results |
| `*_filtered_top20_markers.csv` | Top 20 genes per cluster (filtered) |

**Figure:** `figures/cluster_characterization/dotplot_*_filtered_leiden_res_0.5.pdf` (dot plot of top markers)

---

## 7. Cell Type Annotation

**Code:** `scripts/04_cell_type_annotation.py`  
 **Output:** `results/04_cell_type_annotation/*_final_annotated.h5ad`

### 7.1 Manual Annotation Strategy

Clusters were assigned cell identities based on:

1. Interpretation of multiple top marker genes from filtered Wilcoxon results
2. Canonical marker literature review

### 7.2 Cell Type Mapping

| Cluster | Cell Type | Key Markers 
| --- | --- | --- |
| 0 | C1QA_Macro | $C1QA$, $CD14$, $AIF1$, $HLA-DRA$, $MS4A7$ |
| 1 | Malig_Entero | $CEACAM5$, $KRT20$, $EPCAM$, $CLDN4$, $LGALS4$ |
| 2 | Treg | $FOXP3$, $CTLA4$, $IL2RA$, $TIGIT$, $BATF$, $TNFRSF18$ | 
| 3 | CTL | $GZMA$, $NKG7$, $IFNG$, $CCL5$, $CD3E$, $KLRD1$ | 
| 4 | Mural_Act | $ACTA2$, $TAGLN$, $MYL9$, $TPM2$, $MYLK$, $RGS5$, $CSRP2$, $ADIRF$, $NDUFA4L2$ | 
| 5 | Neutro | $S100A8$, $S100A9$, $CSF3R$, $FCGR3B$, $CXCL8$ | 
| 6 | Endo | $PECAM1$, $VWF$, $PLVAP$, $RAMP2$, $ACKR1$ | 
| 7 | iCAF | $COL1A1$, $CCDC80$, $COL6A3$, $FN1$, $C1R$, $C1S$ | 
| 8 | IgG_PC | $IGHG1-4$, $MZB1$, $JCHAIN$, $CD79A$, $DERL3$ | 
| 9 | kappa_PC | $IGKV4-1$, $IGKV3-20$, $IGKV3-11$, $MZB1$, $JCHAIN$, $CD79A$, $DERL3$ | 
| 10 | Naive_B | $MS4A1$, $VPREB3$, $CD79A/B$, $HLA-DRA$ | 
| 11 | lambda_PC | $IGLV2-14$, $IGLV6-57$, $IGLV2-8$, $MZB1$, $JCHAIN$, $CD79A$, $DERL3$ | 
| 12 | LEC | $PROX1$, $LYVE1$, $CCL21$, $MMRN1$, $PECAM1$ |
| 13 | pDC | $SPIB$, $PLAC8$, $PTGDS$, $GZMB$, $HLA-DRA$ | 
| 14 | Mast | $TPSB2$, $TPSAB1$, $CPA3$, $KIT$, $GATA2$, $FCER1A$ |
| 15 | Schwann | $PLP1$, $S100B$, $CDH19$, $MAL$, $CRYAB$ |  
| 16 | Fibro_quies | $TCF21$, $SCARA5$, $CFD$ | 
| 17 | Tuft | $POU2F3$, $TRPM5$, $BMX$, $AVIL$, $IL17RB$ |

### 7.3 Annotation Application

```python
cl_annotation = {
    "0": "C1QA_Macro", "1": "Malig_Entero", "2": "Treg", "3": "CTL",
    "4": "Mural_Act", "5": "Neutro", "6": "Endo", "7": "iCAF",
    "8": "IgG_PC", "9": "kappa_PC", "10": "Naive_B", "11": "lambda_PC",
    "12": "LEC", "13": "pDC", "14": "Mast", "15": "Schwann",
    "16": "Fibro_quies", "17": "Tuft"
}
adata.obs["manual_celltype_annotation"] = adata.obs["leiden_res_0.5"].astype(str).map(cl_annotation)
```

### 7.4 Output

**Figure:** `figures/cell_type_annotation/umap_*_final_cell_type_annotations.pdf` (UMAP colored by cell type, legend on data points)

**Objects:** Both pre-integration and post-integration objects are annotated and saved separately for compatibility with downstream tools.

---

## 8. Pseudobulk Differential Expression

**Code:** `scripts/05_pseudobulk_analysis.py` → `src/pseudobulk.py`  
 **R dependency:** DESeq2 (via pydeseq2)  
**Output:** `results/05_pseudobulk_analysis/Tumor_vs_Normal_in_{cell_type}/`

### 8.1 Rationale

**Metastasis samples excluded** (n=3) due to insufficient statistical power. Analysis compares  **Tumor (n=6)**  vs  **Normal (n=6)**  only.

### 8.2 Cell Type Subsetting

Two cell types analyzed:

- **Malig_Entero** (malignant epithelial cells)
- **C1QA_Macro** (macrophages)

```python
adata = adata[adata.obs["condition"] != "Metastasis"].copy()
adata_epi = adata[adata.obs["manual_celltype_annotation"] == "Malig_Entero"].copy()
adata_macro = adata[adata.obs["manual_celltype_annotation"] == "C1QA_Macro"].copy()
```

### 8.3 Pseudobulk Aggregation (Decoupler)

```python
pdata = dc.pp.pseudobulk(
    adata,
    sample_col="sample_id",  # e.g., "Normal_1", "Tumor_3"
    groups_col="manual_celltype_annotation",
    layer="soupX_counts",    # Ambient-corrected counts
    mode="sum"               # Sum counts per gene per sample
)
```

### 8.4 Gene Filtering (Pre-DESeq2)

Two-step filtering to remove low-quality genes:

**Filter by expression (edgeR-style):**

```python
dc.pp.filter_by_expr(
    pdata, group="condition",
    min_count=10,          # Minimum counts in a sample
    min_total_count=15,    # Minimum total counts across all samples
    large_n=10,            # Number of top samples for filtering
    min_prop=0.7           # Must pass in ≥70% of samples per condition
)
```

**Filter by proportion (cell occupancy):**

```python
dc.pp.filter_by_prop(
    pdata,
    min_prop=0.1,          # Gene expressed in ≥10% of cells per sample
    min_smpls=2            # Must pass in ≥2 samples per condition
)
```

### 8.5 DESeq2 Analysis (via pydeseq2)

```python
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

inference = DefaultInference(n_cpus=2)
dds = DeseqDataSet(
    adata=pdata,
    design="~ batch + condition",  # condition = Normal or Tumor
    refit_cooks=True,
    inference=inference
)
dds.deseq2()

stat_res = DeseqStats(dds, contrast=["condition", "Tumor", "Normal"])
stat_res.summary()
results_df = stat_res.results_df
```

### 8.6 Significance Thresholds

- **Adjusted p-value (padj):** < 0.05 (Benjamini-Hochberg)
- **Log2 fold change:** |log2FC| > 1.0 (2-fold change)

### 8.7 Output Files

| File | Description |
| --- | --- |
| `*_pseudobulk.h5ad` | Pseudobulk AnnData object |
| `*_all_genes_deseq2.csv` | Complete DESeq2 results |
| `*_significant_DE_genes.csv` | Filtered (padj < 0.05, \|log2FC\| > 1) |

### 8.8 Figures

| Figure | Location |
| --- | --- |
| Filter expression diagnostics | `figures/pseudobulk/{prefix}/{prefix}_filter_expr.pdf` |
| Filter proportion diagnostics | `figures/pseudobulk/{prefix}/{prefix}_filter_prop.pdf` |
| PCA (colored by condition) | `figures/pseudobulk/{prefix}/pca_{prefix}_pca.pdf` |
| Volcano plot | `figures/pseudobulk/{prefix}/{prefix}_volcano.pdf` |

---

## 9. Gene Set Enrichment Analysis

**Code:** `scripts/06_gsea.py` → `src/gsea.py`  
 **Method:** Pre-ranked GSEA (decoupler)  
**Output:** `results/06_gsea/`, `figures/gsea/`

### 9.1 Gene Ranking

Genes are ranked by log2FoldChange:

```python
ranking = deg_df[["log2FoldChange"]].sort_values("log2FoldChange", ascending=False).T
```

### 9.2 Gene Set Databases

GSEA was performed using MSigDB hallmark gene sets [^4] via the decoupler framework.

Two collections from MSigDB (via decoupler):

| Collection | Description |
| --- | --- | 
| `hallmark` | hallmark pathways |
| `immunesigdb` | Immune-related signatures |

```python
msigdb = dc.op.resource("MSigDB")
gene_sets = msigdb[msigdb["collection"] == collection][["geneset", "genesymbol"]]
```

### 9.3 GSEA Implementation

```python
scores, pvals = dc.mt.gsea(ranking, gene_sets, seed=123)
results = pd.concat({"score": scores.T, "pval": pvals.T}, axis=1)
```

### 9.4 Output Files

| File | Content |
| --- | --- |
| `{prefix}_hallmark_results.csv` | Hallmark GSEA results (score, pval) |
| `{prefix}_immunesigdb_results.csv` | ImmuneSigDB GSEA results |

**Figure:** `figures/gsea/{prefix}_gsea.pdf` (Bar plot of top 20 enriched pathways, -log10 p-value)

### 9.5 Analysis Performed

Two comparisons × two gene set collections = 4 total GSEA runs:

1. Tumor_vs_Normal_in_enterocytes × hallmark
2. Tumor_vs_Normal_in_enterocytes × immunesigdb
3. Tumor_vs_Normal_in_macrophages × hallmark
4. Tumor_vs_Normal_in_macrophages × immunesigdb

---

## 10. Cell-Cell Communication Analysis

**Code:** `scripts/07_cell_cell_communication.py` → `scripts/07_cell_cell_communication.R` → `src/cellchats.R`  
 **R dependency:** CellChat (v1.6+)  
**Output:** `results/07_cell_cell_communication_analysis/`

### 10.1 Important Note

Cell-cell communication was inferred using CellChat [^5], which quantifies ligand-receptor interactions from single-cell data.

**Metastasis samples (n=3) excluded**  from CellChat analysis due to:

- Insufficient statistical power (3 vs 6 for normal/tumor)
- Lower total cell counts per sample
- Risk of spurious ligand-receptor inferences

Analysis performed  **separately**  on:

- **Normal samples only** (N1-N6) → `Normal_samples_annotated`
- **Tumor samples only** (T1-T6) → `Tumor_samples_annotated`

### 10.2 Data Preparation

The Python wrapper (`08_cellchat_pipeline.py`) performs:

1. Splits the annotated AnnData by condition (Normal, Tumor)
2. Creates per-condition .h5ad files
3. Generates `k_config.json` for K-selection parameters
4. Launches R script with appropriate arguments

### 10.3 CellChat Object Creation

The annotated AnnData objects were converted to Seurat format using anndataR [^6] for compatibility with CellChat.

```r
create_cc_obj <- function(path, assay = "RNA", slot = "log1p_counts") {
    seurat_obj <- anndataR::read_h5ad(path, as = "Seurat")
    data.input <- GetAssayData(seurat_obj, assay = assay, layer = slot)
    meta <- seurat_obj@meta.data
    cc_obj <- createCellChat(object = data.input, meta = meta, 
                             group.by = "manual_celltype_annotation")
    return(cc_obj)
}
```

### 10.4 Database Configuration

```r
db_use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
ppi <- PPI.human
cc_obj@DB <- db_use
```

### 10.5 Core Inference Pipeline

```r
# Step 1: Subset and identify overexpressed genes
cc_obj <- subsetData(cc_obj)
cc_obj <- identifyOverExpressedGenes(cc_obj)
cc_obj <- identifyOverExpressedInteractions(cc_obj)

# Step 2: Project to PPI network
cc_obj <- projectData(cc_obj, PPI.human)

# Step 3: Compute communication probabilities
cc_obj <- computeCommunProb(cc_obj, raw.use = FALSE)
cc_obj <- filterCommunication(cc_obj, min.cells = 10)

# Step 4: Aggregate pathways
cc_obj <- computeCommunProbPathway(cc_obj)
cc_obj <- aggregateNet(cc_obj)
```

### 10.6 K-Selection (Pattern Number)

Optimal number of outgoing/incoming patterns determined via:

```r
selectK(cc_obj, pattern = "outgoing", k.range = seq(2, 10))
selectK(cc_obj, pattern = "incoming", k.range = seq(2, 10))
```

**Output:** `figures/cellchat_communication/{condition}_samples_annotated/k_diagnostics/*_selectK_diagnostics.pdf`

### 10.7 Systems Analysis

```r
# Compute centrality
cc_obj <- netAnalysis_computeCentrality(cc_obj, slot.name = "netP")

# Identify patterns (heatmaps + river plots)
cc_obj <- identifyCommunicationPatterns(cc_obj, pattern = "outgoing", k = k_out)
cc_obj <- identifyCommunicationPatterns(cc_obj, pattern = "incoming", k = k_in)

# Similarity and clustering
for (type in c("functional", "structural")) {
    cc_obj <- computeNetSimilarity(cc_obj, type = type)
    cc_obj <- netEmbedding(cc_obj, type = type)
    cc_obj <- netClustering(cc_obj, type = type, do.parallel = FALSE)
}
```

**Output figures:**

- `figures/cellchat_communication/{condition}_samples_annotated/individual_patterns/estimationNumCluster__functional_dataset_single.pdf`
- `figures/cellchat_communication/{condition}_samples_annotated/individual_patterns/estimationNumCluster__structural_dataset_single.pdf`
- `figures/cellchat_communication/{condition}_samples_annotated/individual_patterns/{condition}_comprehensive_patterns.pdf`

### 10.8 Saved Objects

- RDS: `results/08_cell_cell_communication_analysis/{condition}_samples_annotated_cellchat.rds`
- AnnData: `results/08_cell_cell_communication_analysis/{condition}_samples_annotated.h5ad`

---

## 11. Differential Communication Analysis

**Code:** `scripts/08_differential_communication.py` → `scripts/08_differential_communication.R` → `src/cellchats.R`  
 **Output:** `results/08_cell_cell_communication_differential/merged_cellchat_final.rds`

### 11.1 Merging CellChat Objects

```r
cc_list <- list(Normal = cellchat_normal, Tumor = cellchat_tumor)

# Align cell identities
group_names <- unique(unlist(lapply(cc_list, function(x) levels(x@idents))))
cc_list <- lapply(cc_list, function(x) liftCellChat(x, group.new = group_names))

merged_cc <- mergeCellChat(cc_list, add.names = names(cc_list))
```

### 11.2 Differential Analysis Functions

**Global statistics:**

```r
compareInteractions(merged_cc, measure = "count")
compareInteractions(merged_cc, measure = "weight")
netVisual_heatmap(merged_cc, measure = "count")
netVisual_heatmap(merged_cc, measure = "weight")
```

**Signaling roles:**

```r
rankNet(merged_cc, mode = "comparison", stacked = TRUE, do.stat = TRUE)
netAnalysis_signalingRole_scatter(cc_obj, title = condition)
netAnalysis_signalingRole_heatmap(cc_obj, pattern = "all")
```

### 11.3 Targeted Pathway Analysis

Three pathways of interest were analyzed with specific sender/receiver configurations:

| Pathway | Senders | Targets |
| --- | --- | --- |
| **SPP1** | C1QA_Macro | Endo, LEC, Mural_Act, Fibro_quies, iCAF  |
| **MIF** | iCAF | C1QA_Macro, CTL, Naive_B, Treg, pDC |
| **MK** | iCAF | Endo, Mural_Act, Malig_Entero |

```r
# Bubble plot for specific senders/receivers
netVisual_bubble(merged_cc, sources.use = sources, targets.use = targets, 
                 signaling = pathway, comparison = c(1,2))

# Role change scatter
netAnalysis_signalingChanges_scatter(merged_cc, signaling.use = pathway)
```

### 11.4 Cell-Type Specific Differential Analysis

```r
# For each cell type of interest
netAnalysis_signalingChanges_scatter(merged_cc, idents.use = cell_type, comparison = c(1,2))
```

### 11.5 Output Figures

| Figure | Location | Content |
| --- | --- | --- |
| `01_global_interactions.pdf` | `figures/cellchat_communication/differential_communication/` | Interaction count & strength comparison |
| `02_signaling_roles.pdf` | Same | RankNet, role scatter plots, pathway heatmaps |
| `03_network_topology.pdf` | Same | Global and pathway-specific circle plots |
| `04_pathway_deep_dive_targeted.pdf` | Same | SPP1, MIF, MK bubble plots + heatmaps |
| `05_cell_specific_diff_scatters.pdf` | Same | Cell-type specific signaling changes |

### 11.6 Saved Object

`results/08_cell_cell_communication_differential/merged_cellchat_final.rds` — Contains merged CellChat object with both Normal and Tumor conditions, enabling direct comparison of signaling networks.

---

## 12. Reproducibility

### 12.1 Random Seeds

All stochastic processes fixed in `src/utils.py`:

```python
def fix_seeds(seed: int = 0):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    scvi.settings.seed = seed
    torch.use_deterministic_algorithms(True)
```

**Applied in:** scVI training.

### 12.2 Configuration Management

All parameters centralized in `config.py` with dataclass-based configuration. For example:

```python
@dataclass
class PreprocessingConfig:
    # Grouping the individual settings into one dictionary
    GENERAL_PARAMS = {
        "prefix": None,
        "batch_key": None,
        "save_plots": True,
    }

    QC_PARAMS = {
        "mad_usage": True,
        "total_counts_bounds": None,
        "n_genes_by_counts_bounds": None,
        "pct_counts_mt_upper": 12,
        "mad_cutoff": 3,
    }

    NORM_PARAMS = {
        "target_sum": 1e4,
        "n_top_genes": 2000,
        "n_pcs": 30
    }

    CLUSTERING_PARAMS = {
        "n_neighbors": 30,
        "use_rep": "X_pca",
        "resolutions": [0.25, 0.5, 0.75, 1],
        "resolution_for_soupx": 0.5
    }
```

### 12.3 Logging

Each preprocessing run generates a detailed log:

- Timestamps (start/end)
- Cell counts before/after filtering
- SoupX rho estimate
- Doublet rate
- Duration

**Summary report:** `results/01_preprocessing/Summary.csv` (generated by `src/qc_utils.py`)

### 12.4 Version Tracking

Key package versions (as used in development):

#### Python: 3.11
- scanpy: 1.10+
- scvi-tools: 1.1+
- anndata: 0.10+
- decoupler: 1.5+
- pydeseq2: 0.4+
- pandas: 2.0+
- numpy: 1.24+
- matplotlib: 3.7+
- seaborn: 0.12+
- scipy: 1.10+
- rpy2: 3.5+
- anndata2ri: 1.2+
- torch: 2.0+

#### R: 4.5.2
- Seurat: 4.0+
- CellChat: 1.6+
- SoupX: 1.6+
- scDblFinder: 1.18+
- ComplexHeatmap: 2.20+
- patchwork: 1.2+
- ggalluvial: 0.12+
- NMF: 0.28+
- reticulate: 1.40+
- anndataR: 0.1+
- SingleCellExperiment: 1.26+
- scater: 1.32+

**Reproducible environment:** Full environment specification is available in `environment.yml` (included in the repository).

### 12.5 Data Availability

Due to patient privacy constraints:

- **Raw FASTQ:** Not publicly available
- **Count matrices:** Available in [GSE315534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315534).

---

## References

[^1]: Young, M. D., & Behjati, S. (2020). SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. *GigaScience*, 9(12), giaa151.

[^2]: Germain, P. L., Lun, A., Garcia Meixide, C., Macnair, W., & Robinson, M. D. (2021). Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research*, 10, 979.

[^3]: Lopez, R., Regier, J., Cole, M. B., Jordan, M. I., & Yosef, N. (2018). Deep generative modeling for single-cell transcriptomics. *Nature Methods*, 15(12), 1053-1058.

[^4]: Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., & Tamayo, P. (2015). The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Systems*, 1(6), 417-425.

[^5]: Jin, S., Guerrero-Juarez, C. F., Zhang, L., et al. (2021). Inference and analysis of cell-cell communication using CellChat. *Nature Communications*, 12(1), 1088.

[^6]: Hao, Y., Hao, S., Andersen-Nissen, E., et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13), 3573-3587.