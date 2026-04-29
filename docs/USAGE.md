# Usage Guide: CRC-SingleCell-Atlas Pipeline

This guide provides **step-by-step instructions** for executing each pipeline component, from raw data to final figures.

---

## Prerequisites

### System Requirements
- **RAM:** 8GB minimum (16GB recommended for integration)
- **Storage:** 10GB free (for intermediate objects)
- **OS:** Windows 10/11 (tested) or Linux/MacOS (should work with path adjustments)
- **CPU:** 2+ cores recommended

### Software
- Python 3.11+
- R 4.5.2+ (with libraries: CellChat, SoupX, scDblFinder, Seurat, anndataR)

---

## Step 0: Initial Setup

### 0.1 Clone and Navigate
```bash
git clone https://github.com/yourlab/CRC_Atlas.git
cd CRC_Atlas
```

### 0.2 Set up Conda Environment

```bash
# Create environment from provided file
conda env create -f environment.yml

# Activate environment
conda activate scrna_pipeline
```

### 0.3 Verify Installation

```bash
# Activate environment first
conda activate scrna_pipeline

# Check Python version (should be 3.11)
python --version

# Verify Python packages
python -c "import scanpy, scvi, anndata, decoupler, pydeseq2, rpy2, anndata2ri, torch; print('✅ Python packages OK')"

# Check Python package versions
python -c "import scanpy; print(f'scanpy: {scanpy.__version__}')"
python -c "import scvi; print(f'scvi-tools: {scvi.__version__}')"
python -c "import anndata; print(f'anndata: {anndata.__version__}')"

# Verify R packages (run in R session)
R -e "cat('R version:', as.character(getRversion()), '\n'); library(Seurat); library(CellChat); library(SoupX); library(scDblFinder); library(ComplexHeatmap); library(patchwork); library(ggalluvial); library(NMF); library(reticulate); library(anndataR); library(SingleCellExperiment); library(scater); cat('✅ R packages OK\n')"
```

### 0.4 Configure Paths in `config.py`

Edit the following lines:

```python
DATA_RAW = Path("C:/your/path/to/data/raw")
R_HOME = r'C:\Program Files\R\R-4.5.2'  # Windows example
```

### 0.5 Verify Data Structure

Your `data/raw/` folder should contain:

```text
data/raw/
├── GSM9430352_CRC_N1/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── GSM9430353_CRC_T1/
│   ├── ...
└── ... (15 total folders)
```

---

## Step 1: Per-Sample Preprocessing

**Script:** `scripts/01_preprocessing.py`  
 **Duration:** 10-30 minutes per sample (depending on cell count)

### Execute

```bash
python scripts/01_preprocessing.py
```

### What Happens

1. Loads each 10x dataset
2. Calculates QC metrics (total_counts, n_genes, mt%, ribo%, hb%)
3. Filters cells using MAD or thresholds
4. Normalizes (1e4 target sum) and log1p transforms
5. Identifies HVGs (2000 genes)
6. Runs PCA and Leiden clustering (resolutions 0.25-1.0)
7. **If ≥100 cells:** Runs SoupX ambient RNA correction
8. **If ≥100 cells:** Runs scDblFinder doublet detection and removal
9. Saves cleaned object to `results/01_preprocessing/{sample}/{sample}_cleaned.h5ad`

### Output Files

```text
results/01_preprocessing/
├── N1/
│   ├── N1_cleaned.h5ad
│   └── N1_processing_log.txt
├── T1/
│   ├── T1_cleaned.h5ad
│   └── T1_processing_log.txt
├── LM1/
│   ├── LM1_cleaned.h5ad
│   └── LM1_processing_log.txt
├── ... (15 samples)
└── Summary.csv
```

### Check Success

- Open `results/01_preprocessing/Summary.csv` - verify cell counts are reasonable
- Check logs for errors: `grep "ERROR" results/01_preprocessing/*/*.txt`
- Visualize QC: `figures/qc/{sample}/*.pdf`

---

## Step 2: Batch Integration (scVI)

**Script:** `scripts/02_integration.py`  
 **Duration:** 1-3 hours (depending on cell count, GPU available)

### Execute

```bash
python scripts/02_integration.py
```

### What Happens

1. Merges all 15 cleaned samples into one AnnData
2. Adds metadata: `condition` (Normal/Tumor/Metastasis), `batch` (patient ID)
3. Runs scVI model (30 latent dimensions, 2 layers, 150 epochs)
4. Extracts batch-corrected latent representation (`X_scVI`)
5. Extracts batch-corrected normalized counts (`scvi_normalized` layer)
6. Saves pre-integration and post-integration objects

### Output Files

```text
results/02_integration/
├── merged_dataset_pre_integration.h5ad
├── merged_dataset_post_integration.h5ad
└── scvi_model/
    ├── model.pt
    └── ...
```

### Figures

```text
figures/integration/
├── merged_dataset_Pre-Integration Composition.pdf
├── merged_dataset_Post-Integration Composition.pdf
├── merged_dataset_scvi_training_history.pdf
├── merged_dataset_latent_variance.pdf
└── umap_merged_dataset_leiden_resolution_optimization_*.pdf
```

### Validation

- **Check mixing:** Compare pre/post integration UMAPs - patients should intermix
- **Check convergence:** Training history should show decreasing ELBO
- **Check latent variance:** ~80% variance in first 10 dimensions is typical

---

## Step 3: Cluster Characterization (Marker Identification)

**Script:** `scripts/03_cluster_characterization.py`  
 **Duration:** 10-20 minutes

### Execute

```bash
python scripts/03_cluster_characterization.py
```

### What Happens

1. Loads post-integration object
2. Computes UMAP on `X_scVI`
3. Runs Leiden clustering (resolution 0.5)
4. Identifies marker genes (Wilcoxon rank-sum test)
5. Applies filters (min_in_group_fraction=0.25, max_out_group_fraction=0.50)
6. Exports marker CSVs
7. Generates dot plot of top markers

### Output Files

```text
results/03_cluster_characterization/
├── merged_dataset_post_integration_all_markers_full.csv
├── merged_dataset_post_integration_top20_markers.csv
├── merged_dataset_post_integration_filtered_all_markers_full.csv
└── merged_dataset_post_integration_filtered_top20_markers.csv
```

### Figure

`figures/cluster_characterization/dotplot__merged_dataset_post_integration_filtered_leiden_res_0.5.pdf`

### Next Step

Review the dot plot and marker CSVs to prepare for manual annotation. Identify which clusters correspond to which cell types.

---

## Step 4: Cell Type Annotation

**Script:** `scripts/04_cell_type_annotation.py`  
 **Duration:** 5 minutes

### Execute

```bash
python scripts/04_cell_type_annotation.py
```

### What Happens

1. Loads post-integration object
2. Applies manual annotation mapping (clusters → cell types)
3. Saves annotated object (both pre and post integration)
4. Generates UMAP with cell type labels

### Output Files

```text
results/04_cell_type_annotation/
├── merged_dataset_pre_integration_final_annotated.h5ad
└── merged_dataset_post_integration_final_annotated.h5ad
```

### Figure

`figures/cell_type_annotation/umap_merged_dataset_post_integration_final_cell_type_annotations.pdf`

### Customizing Annotation

Edit `apply_manual_annotation()` in `scripts/04_cell_type_annotation.py` to modify cell type mapping.

---

## Step 5: Pseudobulk Differential Expression

**Script:** `scripts/05_pseudobulk_analysis.py`  
 **Duration:** 15 minutes per cell type

### Execute

```bash
python scripts/05_pseudobulk_analysis.py
```

### What Happens

1. Excludes metastasis samples (keeps Normal + Tumor only)
2. Subsets Malig_Entero and C1QA_Macro cells
3. Aggregates counts per patient (pseudobulk)
4. Filters low-expression genes
5. Runs DESeq2 (Tumor vs Normal)
6. Saves results and generates figures

### Output Files

```text
results/05_pseudobulk_analysis/
├── Tumor_vs_Normal_in_enterocytes/
│   ├── Tumor_vs_Normal_in_enterocytes_pseudobulk.h5ad
│   ├── Tumor_vs_Normal_in_enterocytes_all_genes_deseq2.csv
│   └── Tumor_vs_Normal_in_enterocytes_significant_DE_genes.csv
└── Tumor_vs_Normal_in_macrophages/
    ├── ... (same structure)
```

### Figures

```text
figures/pseudobulk/
├── Tumor_vs_Normal_in_enterocytes/
│   ├── Tumor_vs_Normal_in_enterocytes_filter_expr.pdf
│   ├── Tumor_vs_Normal_in_enterocytes_filter_prop.pdf
│   ├── pca_Tumor_vs_Normal_in_enterocytes_pca.pdf
│   └── Tumor_vs_Normal_in_enterocytes_volcano.pdf
└── Tumor_vs_Normal_in_macrophages/
    └── ... (same structure)
```

### Check Results

- **Volcano plot:** Should show clear up/down regulated genes
- **PCA:** Tumor and Normal samples should separate
- **Significant genes count:** Expect 100-1000 DE genes per comparison

---

## Step 6: GSEA

**Script:** `scripts/06_gsea.py`  
 **Duration:** 5-10 minutes per comparison

### Execute

```bash
python scripts/06_gsea.py
```

### What Happens

1. Loads DESeq2 results for each comparison
2. Ranks genes by log2FoldChange
3. Runs pre-ranked GSEA (Hallmark + ImmuneSigDB)
4. Saves results and generates bar plots

### Output Files

```text
results/06_gsea/
├── Tumor_vs_Normal_in_enterocytes_gsea_hallmark_results.csv
├── Tumor_vs_Normal_in_enterocytes_gsea_immunesigdb_results.csv
├── Tumor_vs_Normal_in_macrophages_gsea_hallmark_results.csv
└── Tumor_vs_Normal_in_macrophages_gsea_immunesigdb_results.csv
```

### Figures

```text
figures/gsea/
├── Tumor_vs_Normal_in_enterocytes_gsea_gsea.pdf
├── Tumor_vs_Normal_in_macrophages_gsea_gsea.pdf
└── (PDFs for ImmuneSigDB if generated)
```

### Interpretation

- **Positive enrichment score:** Pathway upregulated in Tumor
- **Negative enrichment score:** Pathway downregulated in Tumor
- **p < 0.05:** Statistically significant enrichment

---

## Step 7: Cell-Cell Communication (CellChat)

**Script:** `scripts/07_cell_cell_communication.py` (calls R script)  
**Duration:** 30-60 minutes per condition

### Execute

```bash
python scripts/07_cell_cell_communication.py
```

### What Happens

1. Splits annotated object into Normal and Tumor subsets
2. Launches R script for each condition
3. Creates CellChat objects
4. Computes communication probabilities
5. Identifies outgoing/incoming patterns (K-selection)
6. Saves RDS objects and generates figures

### Output Files

```text
results/08_cell_cell_communication_analysis/
├── Normal_samples_annotated.h5ad
├── Normal_samples_annotated_cellchat.rds
├── Tumor_samples_annotated.h5ad
├── Tumor_samples_annotated_cellchat.rds
└── k_config.json
```

### Figures

```text
figures/cellchat_communication/
├── Normal_samples_annotated/
│   ├── individual_patterns/
│   │   ├── estimationNumCluster__functional_dataset_single.pdf
│   │   ├── estimationNumCluster__structural_dataset_single.pdf
│   │   └── Normal_samples_annotated_comprehensive_patterns.pdf
│   └── k_diagnostics/
│       └── Normal_samples_annotated_selectK_diagnostics.pdf
└── Tumor_samples_annotated/
    └── ... (same structure)
```

### Validation

- Check K-selection plots: Identify "elbow" for optimal k
- Comprehensive patterns PDF should show heatmaps and river plots

---

## Step 8: Differential Communication

**Script:** `scripts/08_differential_communication.py` (calls R script)  
**Duration:** 15-20 minutes

### Execute

```bash
python scripts/08_differential_communication.py
```

### What Happens

1. Loads Normal and Tumor CellChat RDS objects
2. Merges them for comparison
3. Compares global interaction counts and weights
4. Compares signaling roles (outgoing/incoming)
5. Analyzes targeted pathways (SPP1, MIF, MK)
6. Generates differential figures

### Output Files

```text
results/08_cell_cell_communication_differential/
└── merged_cellchat_final.rds
```

### Figures

```text
figures/cellchat_communication/differential_communication/
├── 01_global_interactions.pdf
├── 02_signaling_roles.pdf
├── 03_network_topology.pdf
├── 04_pathway_deep_dive_targeted.pdf
└── 05_cell_specific_diff_scatters.pdf
```

### Key Interpretations

- **01_global_interactions.pdf:** Interaction count/strength higher in Tumor?
- **02_signaling_roles.pdf:** Which pathways gain/lose outgoing/incoming signaling?
- **03_network_topology.pdf:** Circle plots showing SPP1, MIF, MK networks
- **04_pathway_deep_dive_targeted.pdf:** Bubble plots for specific senders/receivers
- **05_cell_specific_diff_scatters.pdf:** How individual cell types' signaling changes

---

## Step 9: Generate Final Report

### Summary Statistics

Run the QC report generator:

```python
from src.qc_utils import run_report_generation
run_report_generation("results/01_preprocessing", "Summary.csv")
```

This updates `results/01_preprocessing/Summary.csv` with all sample metrics.

---

## Common Issues and Solutions

| Issue | Likely Cause | Solution |
| --- | --- | --- |
| R script fails | Missing R packages | Run `install.packages(c(...))` in R |
| scVI won't train | GPU memory | Reduce `batch_size` in `integration.py` |
| SoupX fails | <100 cells | Normal - script auto-skips |
| Memory error | Too many cells | Use subsampling or increase RAM |
| Path errors | Windows vs Linux | Use `Path()` objects, not raw strings |
| doublet detection fails | scDblFinder version | Update: `remotes::install_github('plger/scDblFinder')` |
