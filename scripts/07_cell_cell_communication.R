"""
Worker Script: 07_cell_cell_communication.R
Project: Colorectal Cancer (CRC) Single-Cell Atlas
Stage: Per-Condition Signaling Inference & Pattern Discovery

This script iterates through condition-specific h5ad files (e.g., Normal, Tumor),
performs ligand-receptor interaction inference using CellChat, and executes
manifold learning to identify global communication patterns.

Input: Condition-specific .h5ad files
Output: .rds objects for each condition and diagnostic figures.

Author: Loai Eletr
"""

library(jsonlite)
library(reticulate)

# --- 1. ENVIRONMENT VALIDATION ---
# CellChat uses Python's UMAP for manifold learning (Pattern Analysis).
if (!py_module_available("umap")) {
  stop("❌ UMAP not found. Ensure your Python environment is correctly linked to R via reticulate.")
}

# --- 2. ARGUMENT PARSING ---
# These arguments are passed from the 07_cell_cell_communication.py orchestrator.
args <- commandArgs(trailingOnly = TRUE)
project_dir  <- args[1] # Root project directory
input_dir    <- args[2] # Directory containing split h5ad files
rds_dir      <- args[3] # Destination for processed R objects
fig_dir      <- args[4] # Destination for individual signaling plots
species      <- args[5] # 'human' or 'mouse'
min_cells    <- as.numeric(args[6]) # Threshold for cell group inclusion
k_json_path  <- args[7] # Path to JSON containing optimized K values

# --- 3. RESOURCE INITIALIZATION ---
source(file.path(project_dir, "src/cellchats.R")) # Custom wrapper functions
k_map <- fromJSON(k_json_path)                    # Load k_out and k_in per sample

# Create output directories if they don't exist
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Discover the data files to be processed
files <- list.files(input_dir, pattern = "*.h5ad", full.names = TRUE)

# --- 4. MAIN PROCESSING LOOP ---
for (f in files) {
  name <- tools::file_path_sans_ext(basename(f))

  # Retrieve optimized K values from JSON, fallback to 3 if not found
  k_out <- if (!is.null(k_map[[name]])) k_map[[name]]$k_out else 3
  k_in  <- if (!is.null(k_map[[name]])) k_map[[name]]$k_in else 3

  cat("\n", rep("=", 50), "\n", sep="")
  cat("🚀 PROCESSING SAMPLE: ", name, "\n")
  cat(rep("=", 50), "\n", sep="")

  # --- STEP 1: INFERENCE ---
  # Initialize CellChat object and set the L-R database (e.g., CellChatDB.human)
  obj <- create_cc_obj(f)
  setup <- set_cc_db(obj, species)

  # Calculate communication probabilities and filter out rare cell groups
  obj <- run_cc_inference(setup$obj, setup$ppi, min_cells = min_cells)

  # --- STEP 2: DIAGNOSTIC K-SELECTION ---
  # Generates plots to verify if the number of patterns (k_out/k_in) is appropriate
  cat("[*] Generating K-selection diagnostics...\n")
  run_cc_k_selection(obj, name = name, plots_dir = file.path(fig_dir, name, "k_diagnostics"))

  # --- STEP 3: SYSTEMS ANALYSIS (PATTERNS) ---
  # Identifies which cell types are 'outgoing' (secretors) and 'incoming' (receivers)
  # across coordinated signaling pathways.
  cat("[*] Identifying communication patterns (k_out=", k_out, ", k_in=", k_in, ")...\n", sep="")
  obj <- run_cc_systems_analysis(
    obj, k_out = k_out, k_in = k_in, name = name,
    plots_dir = file.path(fig_dir, name, "individual_patterns")
  )

  # --- STEP 4: DATA PERSISTENCE ---
  # Save as RDS for use in the Stage 08 Differential Analysis
  saveRDS(obj, file = file.path(rds_dir, paste0(name, "_cellchat.rds")))
  cat("✅ Saved: ", paste0(name, "_cellchat.rds"), "\n")
}

cat("\n🎉 Individual processing complete. Ready for differential analysis.\n")