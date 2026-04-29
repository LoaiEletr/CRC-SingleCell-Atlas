"""
Worker Script: 08_differential_communication.R
Project: Colorectal Cancer (CRC) Single-Cell Atlas
Stage: Comparative Signaling Analysis & Pathway Deep-Dive

This script merges condition-specific CellChat objects to identify 
differential signaling. It focuses on specific pathways critical to CRC 
pathology (SPP1, MIF, MK) and evaluates how cell-type roles shift between 
Normal and Tumor states.

Input: Per-condition .rds files
Output: Merged .rds object and comparative visualization suite.

Author: Loai Eletr
"""

library(CellChat)

# --- 1. ARGUMENT PARSING & SETUP ---
args <- commandArgs(trailingOnly = TRUE)
project_dir   <- args[1] # Root directory
rds_indiv_dir <- args[2] # Where Stage 08 saved individual RDS files
rds_diff_dir  <- args[3] # Where the merged object will be saved
fig_dir       <- args[4] # Output directory for differential plots

source(file.path(project_dir, "src/cellchats.R"))

# --- 2. DATA MERGING & ALIGNMENT ---
# Loading all condition RDS files (Normal, Tumor, etc.)
rds_files <- list.files(rds_indiv_dir, pattern = "_cellchat.rds$", full.names = TRUE)
cc_list   <- lapply(rds_files, readRDS)
names(cc_list) <- gsub("_cellchat.rds", "", basename(rds_files))

# IMPORTANT: 'liftCellChat' ensures all objects share the same cell type labels.
# This is required for the matrices to align during subtraction/comparison.
group_names <- unique(unlist(lapply(cc_list, function(x) levels(x@idents))))
cc_list     <- lapply(cc_list, function(x) liftCellChat(x, group.new = group_names))
merged_cc   <- mergeCellChat(cc_list, add.names = names(cc_list))

# --- 3. BIOLOGICAL GROUPING (Storytelling) ---
# Mapping fine-grained clusters to broad lineages for cleaner visualization.
broad_groups <- c(
  "Endo" = "Endothelial", "LEC" = "Endothelial",
  "iCAF" = "Fibroblasts", "Fibro_quies" = "Fibroblasts", "Mural_Act" = "Fibroblasts",
  "Treg" = "T_Cells", "CTL" = "T_Cells",
  "Malig_Entero" = "Malignant",
  "C1QA_Macro" = "Myeloid", "Mono" = "Myeloid", "Neutro" = "Myeloid"
)

# --- 4. TARGETED PATHWAY CONFIGURATIONS ---
# These pathways represent the "Malignant-Myeloid-Stroma" axis in CRC.
pathway_configs <- list(
  "SPP1" = list(
    sources = c("C1QA_Macro"), 
    targets = c("Endo", "LEC", "Mural_Act", "Fibro_quies", "iCAF")
  ),
  "MIF" = list(
    sources = "iCAF", 
    targets = c("C1QA_Macro", "CTL", "Naive_B", "Treg", "pDC")
  ),
  "MK" = list(
    sources = "iCAF", 
    targets = c("Endo","Mural_Act", "Malig_Entero")
  )
)

# --- 5. VISUALIZATION SUITE ---

# 📊 GLOBAL STATS: Number of interactions and interaction strength comparisons.
plot_global_statistics(merged_cc, fig_dir)

# 🗺️ SIGNALING ROLES: Who are the dominant senders and receivers in each condition?
plot_signaling_roles(cc_list, merged_cc, fig_dir)

# ⭕ TOPOLOGY: Differential circle plots for the key pathways (SPP1, MIF, MK).
plot_topology_comparison(cc_list, fig_dir, pathways = names(pathway_configs), group_map = broad_groups)

# 🔍 DEEP DIVE: Bubble plots and chord diagrams for specific source-target pairs.
plot_pathway_details(merged_cc, cc_list, pathway_configs, fig_dir)

# 📈 SCATTER PLOTS: Identify which incoming/outgoing signals changed most for specific cells.
all_focus_cells <- c("Malig_Entero","C1QA_Macro", "iCAF")
plot_cell_specific_changes(merged_cc, cell_types = all_focus_cells, fig_dir = fig_dir)

# --- 6. DATA PERSISTENCE ---
saveRDS(merged_cc, file = file.path(rds_diff_dir, "merged_cellchat_final.rds"))
cat("\n✅ Differential Analysis Complete. Merged object saved.\n")