# ==============================================================================
# Cell-Cell Communication (CCC) Analysis Pipeline: CellChat Integration
#
# This script facilitates the inference, visualization, and comparison of
# ligand-receptor interactions across biological conditions (e.g., Normal vs Tumor).
#
# Key Features:
# 1. Automated CellChat object creation from Scanpy/H5AD.
# 2. Signaling pattern identification (Outgoing vs Incoming).
# 3. Comparative topology and pathway-specific "deep-dives".
#
# Author: Loai Eletr
# ==============================================================================

library(anndataR)
library(CellChat)
library(Seurat)
library(reticulate)
library(ComplexHeatmap)
library(patchwork)
library(NMF)
library(ggalluvial)

# --- 1. Object Creation ---
create_cc_obj <- function(path, assay = "RNA", slot = "log1p_counts") {
  cat("\n[STEP 1/5] Loading Data...\n")
  message("Reading h5ad: ", basename(path))
  
  # Load H5AD and convert to Seurat for CellChat compatibility
  seurat_obj <- anndataR::read_h5ad(path, as = "Seurat")
  data.input <- GetAssayData(seurat_obj, assay = assay, layer = slot)
  meta <- seurat_obj@meta.data

  # Initialize CellChat object using your manual cell type annotations
  cc_obj <- createCellChat(object = data.input, meta = meta, group.by = "manual_celltype_annotation")
  
  message("✅ CellChat object created with ", ncol(data.input), " cells.")
  return(cc_obj)
}

# --- 2. Database Selection ---
# Selects the signaling molecules to test based on the species and pathway type.
set_cc_db <- function(cc_obj, species = "human", db_type = "Secreted Signaling") {
  cat("[STEP 2/5] Setting up Database...\n")
  
  if (species == "human") {
    db_use <- subsetDB(CellChatDB.human, search = db_type)
    ppi <- PPI.human
  } else {
    db_use <- subsetDB(CellChatDB.mouse, search = db_type)
    ppi <- PPI.mouse
  }
  
  cc_obj@DB <- db_use
  message("✅ Database set to: ", species, " (", db_type, ")")
  return(list(obj = cc_obj, ppi = ppi))
}

# --- 3. Core Inference Pipeline ---
# Executes the statistical modeling of interaction probabilities.
run_cc_inference <- function(cc_obj, ppi, min_cells = 10) {
  cat("[STEP 3/5] Inferring Communication Probabilities...\n")

  # Identify differentially expressed ligands/receptors
  cc_obj <- subsetData(cc_obj)
  cc_obj <- identifyOverExpressedGenes(cc_obj)
  cc_obj <- identifyOverExpressedInteractions(cc_obj)

  # Project to PPI to smooth and validate potential interactions
  cc_obj <- projectData(cc_obj, ppi)

  # Compute probabilities using the Law of Mass Action
  cc_obj <- computeCommunProb(cc_obj, raw.use = FALSE)
  cc_obj <- filterCommunication(cc_obj, min.cells = min_cells)

  # Consolidate individual L-R pairs into cohesive pathways (e.g., SPP1, WNT)
  cc_obj <- computeCommunProbPathway(cc_obj)
  cc_obj <- aggregateNet(cc_obj)
  
  message("✅ Core inference complete.")
  return(cc_obj)
}

# --- 4. Diagnostic: K-Selection ---
# Determines the optimal number of patterns (K) for matrix factorization (NMF).
run_cc_k_selection <- function(cc_obj, name = "sample", plots_dir = "plots", max_k = 10) {
  cat("[DIAGNOSTIC] Running K-selection...\n")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  plot_path <- file.path(plots_dir, paste0(name, "_selectK_diagnostics.pdf"))
  pdf(plot_path, width = 12, height = 6)

  # Metrics look for the "elbow" where cophenetic correlation drops
  print(selectK(cc_obj, pattern = "outgoing", k.range = seq(2, max_k)))
  print(selectK(cc_obj, pattern = "incoming", k.range = seq(2, max_k)))

  dev.off()
  message("✅ K-selection plots saved to: ", plot_path)
}

# --- 5. Systems Analysis ---
# Identifies communication "programs" or themes across cell types.
run_cc_systems_analysis <- function(cc_obj, k_out = 3, k_in = 5, name = "sample", plots_dir = "plots") {
  cat("[STEP 4/5] Running Systems Analysis...\n")

  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  # 1. Temporary WD to catch internal clustering PDFs
  old_wd <- getwd()
  setwd(plots_dir)

  message("--> Computing network centrality...")
  cc_obj <- netAnalysis_computeCentrality(cc_obj, slot.name = "netP")

  # 2. Open the Comprehensive PDF Device
  # This file will now contain BOTH Heatmaps and River plots
  final_pdf_path <- file.path(getwd(), paste0(name, "_comprehensive_patterns.pdf"))
  pdf(final_pdf_path, width = 12, height = 10)

  message("--> Identifying patterns (Saving Heatmaps to PDF)...")
  # Page 1 & 2: Heatmaps
  # We set heatmap.show = TRUE here because the PDF device is open and waiting
  cc_obj <- identifyCommunicationPatterns(cc_obj, pattern = "outgoing", k = k_out, heatmap.show = TRUE)
  cc_obj <- identifyCommunicationPatterns(cc_obj, pattern = "incoming", k = k_in, heatmap.show = TRUE)

  message("--> Generating River plots (Adding to PDF)...")
  # Page 3 & 4: River Plots
  try({
    print(netAnalysis_river(cc_obj, pattern = "outgoing", font.size = 2))
    print(netAnalysis_river(cc_obj, pattern = "incoming", font.size = 2))
  })

  # Close the PDF device
  dev.off()

  # 3. Handle Similarity & Clustering
  # Note: These still generate their own "estimationNumCluster" PDFs in the current WD (plots_dir)
  message("--> Calculating similarity and embeddings...")
  for (type in c("functional", "structural")) {
    cc_obj <- computeNetSimilarity(cc_obj, type = type)
    cc_obj <- netEmbedding(cc_obj, type = type)
    cc_obj <- netClustering(cc_obj, type = type, do.parallel = FALSE)
  }

  # Return to original directory
  setwd(old_wd)

  message("✅ Analysis complete. Comprehensive PDF saved: ", basename(final_pdf_path))
  return(cc_obj)
}

# --- 6. Global Interaction Stats ---
# Compares overall "talkativeness" between datasets (e.g., Normal vs Tumor).
plot_global_statistics <- function(merged_cc, fig_dir) {
  message("--> Plotting Global Interaction Stats...")
  pdf(file.path(fig_dir, "01_global_interactions.pdf"), width = 12, height = 7)
  
  p1 <- compareInteractions(merged_cc, show.legend = FALSE, measure = "count")
  p2 <- compareInteractions(merged_cc, show.legend = FALSE, measure = "weight")
  print(p1 + p2)
  
  print(netVisual_heatmap(merged_cc, measure = "count", title.name = "Differential Interaction Count"))
  print(netVisual_heatmap(merged_cc, measure = "weight", title.name = "Differential Interaction Strength"))
  
  dev.off()
}

# --- 7. Signaling Role & Pathway Alignment ---
# Identifies which cell types act as hubs for specific pathways.
plot_signaling_roles <- function(cc_list, merged_cc, fig_dir) {
  message("--> Plotting Signaling Roles...")
  pdf(file.path(fig_dir, "02_signaling_roles.pdf"), width = 14, height = 12)
  
  print(rankNet(merged_cc, mode = "comparison", stacked = TRUE, do.stat = TRUE))
  
  # Role Identification Scatters (All cell types)
  all_weights <- sapply(cc_list, function(x) { 
    rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count) 
  })
  weight.MinMax <- c(min(all_weights), max(all_weights))
  
  gg <- list()
  for (i in 1:length(cc_list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(cc_list[[i]], 
                                                 title = names(cc_list)[i], 
                                                 weight.MinMax = weight.MinMax)
  }
  print(patchwork::wrap_plots(plots = gg))
  
  # Aligned Pathway Heatmaps
  all_pathways <- unique(unlist(lapply(cc_list, function(x) x@netP$pathways)))
  ht1 <- netAnalysis_signalingRole_heatmap(cc_list[[1]], pattern = "all", signaling = all_pathways, 
                                           title = names(cc_list)[1], width = 5, height = 12)
  ht2 <- netAnalysis_signalingRole_heatmap(cc_list[[2]], pattern = "all", signaling = all_pathways, 
                                           title = names(cc_list)[2], width = 5, height = 12)
  draw(ht1 + ht2, ht_gap = unit(5, "cm"))
  
  dev.off()
}


# --- 8. Topology Functions (Global & Pathway-Specific Circles) ---
# Visualizes "who talks to whom" for global networks and specific pathways.
plot_topology_comparison <- function(cc_list, fig_dir, pathways = c("SPP1", "MIF"), group_map = NULL) {
  message("--> Plotting Network Topology and Pathway Aggregates...")
  pdf(file.path(fig_dir, "03_network_topology.pdf"), width = 14, height = 10)
  
  # 1. GLOBAL INTERACTION COUNT (Number of interactions)
  # weight.max[2] corresponds to the "count" attribute
  weight.max_global <- getMaxWeight(cc_list, attribute = c("idents", "count"))
  
  par(mfrow = c(1, length(cc_list)), xpd = TRUE)
  for (i in 1:length(cc_list)) {
    netVisual_circle(cc_list[[i]]@net$count, 
                     weight.scale = TRUE, 
                     label.edge = FALSE,			 
                     edge.weight.max = weight.max_global[2], 
                     edge.width.max = 10, 
                     title.name = paste0("Global Count: ", names(cc_list)[i]))
  }

  # 2. PATHWAY-SPECIFIC AGGREGATE CIRCLES (e.g., SPP1, MIF)
  # This matches the specific "netVisual_aggregate" logic you provided
  for (pw in pathways) {
    message(paste("   - Plotting Aggregate Circle for:", pw))
    
    # Calculate max weight for THIS specific pathway to keep scales consistent
    weight.max_pw <- getMaxWeight(cc_list, slot.name = c("netP"), attribute = pw)
    
    par(mfrow = c(1, length(cc_list)), xpd = TRUE)
    for (i in 1:length(cc_list)) {
      if (pw %in% cc_list[[i]]@netP$pathways) {
        netVisual_aggregate(cc_list[[i]], 
                            signaling = pw, 
                            layout = "circle",
                            edge.weight.max = weight.max_pw[1], 
                            edge.width.max = 10, 
                            arrow.size = 0.05, 
                            signaling.name = paste(pw, "-", names(cc_list)[i]))
      }
    }
  }
  
  dev.off()
}

# --- 9. Targeted Pathway Deep-Dives (Pathway-Specific Senders/Receivers) ---
# Analyzes specific cell-cell pairs for a given pathway.
plot_pathway_details <- function(merged_cc, cc_list, pathway_configs, fig_dir) {
  message("--> Plotting Targeted Pathway Deep-Dives with custom Senders/Receivers...")
  pdf(file.path(fig_dir, "04_pathway_deep_dive_targeted.pdf"), width = 14, height = 11)
  
  # pathway_configs is a named list: e.g., list("SPP1" = list(sources = ..., targets = ...))
  pathways <- names(pathway_configs)

  for (pw in pathways) {
    try({
      # Extract pathway-specific filters
      current_sources <- pathway_configs[[pw]]$sources
      current_targets <- pathway_configs[[pw]]$targets
      
      message(paste("   - Processing:", pw))
      
      # 1. Bubble Plot with Specific Senders/Receivers for THIS pathway
      print(netVisual_bubble(merged_cc, 
                             sources.use = current_sources, 
                             targets.use = current_targets,
                             signaling = pw, 
                             comparison = c(1, 2), 
                             max.dataset = 2,
                             title.name = paste0(pw, " Signaling: Targeted Senders/Receivers"), 
                             angle.x = 45, remove.isolate = TRUE))
      
      # 2. Individual Heatmaps (Structural view)
      ht <- list()
      for (i in 1:length(cc_list)) {
        if (pw %in% cc_list[[i]]@netP$pathways) {
           ht[[i]] <- netVisual_heatmap(cc_list[[i]], signaling = pw, 
                                        title.name = paste(pw, "Signaling -", names(cc_list)[i]),
                                        color.heatmap = "Reds")
        }
      }
      if(length(ht) == 2) draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
      
      # 3. Role Change Scatter (Global shift for this pathway)
      print(netAnalysis_signalingChanges_scatter(merged_cc, signaling.use = pw))
      
    }, silent = TRUE)
  }
  dev.off()
}

# --- 10. Cell-Type Specific Differential Scatters ---
# Focused on how specific cells (e.g., Treg, Mono) shift across all signaling
plot_cell_specific_changes <- function(merged_cc, cell_types = c("Treg", "Mono"), fig_dir) {
  message("--> Plotting Cell-Type Specific Role Shifts...")
  pdf(file.path(fig_dir, "05_cell_specific_diff_scatters.pdf"), width = 10, height = 8)
  
  for (ct in cell_types) {
    try({
      # Global shift for this cell type
      print(netAnalysis_signalingChanges_scatter(merged_cc, idents.use = ct, comparison = c(1,2)))
    }, silent = TRUE)
  }
  dev.off()
}