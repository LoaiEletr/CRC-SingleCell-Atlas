"""
Orchestration Script: 08_Differential_Communication
Project: Colorectal Cancer (CRC) Single-Cell Atlas
Stage: Comparative Cell-Cell Communication (CCC) Analysis

This script manages the high-level comparison between inferred signaling
networks (e.g., Normal vs. Tumor). It invokes the R-based differential
engine to identify shifts in signaling strength, information flow, and
conserved/specific signaling pathways across biological conditions.

Author: Loai Eletr
"""

import subprocess
import config
import time
from config import CellChatConfig
from src.utils import save_adata, load_r_environment

# ==============================================================================
# MAIN EXECUTION BLOCK
# ==============================================================================

if __name__ == "__main__":
    print(f"\n{'=' * 60}")
    print(f"🧬 STARTING DIFFERENTIAL CELL-CELL COMMUNICATION ANALYSIS")
    print(f"{'=' * 60}")
    start_time = time.time()

    # --- 1. ENVIRONMENT SETUP ---
    # Re-initialize the R path and environment variables for subprocess execution.
    print(f"[*] Initializing R environment at: {config.R_HOME}")
    load_r_environment(config.R_HOME)

    print("[*] Preparing directories for differential results...")
    # Prepare environment handles directory creation and system path configurations
    custom_env = CellChatConfig.prepare_environment()

    # --- 2. R-SUBSYSTEM INVOCATION ---
    # Launching the comparative R script (09_differential_communication.R).
    # This script performs the 'Merge' operation and calculates differential weights.
    print(f"\n[!] Launching Rscript: 09_differential_communication.R")
    print(f"{'-' * 60}")

    try:
        # Subprocess call passes the necessary directories:
        # - individual_cellchat_dir: Where the separate RDS files are stored.
        # - diff_cellchat_dir: Where the merged comparative object will be saved.
        # - fig_dir: Output for comparison heatmaps, circle plots, and river plots.
        subprocess.run([
            "Rscript",
            f"{config.BASE_DIR}/scripts/08_differential_communication.R",
            str(config.BASE_DIR),
            str(CellChatConfig.CELLCHAT_IO["individual_cellchat_dir"]),
            str(CellChatConfig.CELLCHAT_IO["diff_cellchat_dir"]),
            str(CellChatConfig.CELLCHAT_IO["fig_dir"] / "differential_communication"),
        ], env=custom_env, check=True)

        # --- 3. PERFORMANCE MONITORING & LOGGING ---
        elapsed = (time.time() - start_time) / 60

        print(f"{'-' * 60}")
        print(f"✅ ANALYSIS COMPLETE in {elapsed:.2f} minutes")
        print(f"📂 Merged RDS:  {CellChatConfig.CELLCHAT_IO['diff_cellchat_dir']}")
        print(f"📊 Diff Figures: {CellChatConfig.CELLCHAT_IO['fig_dir'] / 'differential_communication'}")
        print(f"{'=' * 60}\n")

    except subprocess.CalledProcessError as e:
        # ERROR LOGGING: Captures R failures (e.g., non-overlapping cell types
        # or mismatched factor levels between objects).
        print(f"\n❌ ERROR: Rscript failed with exit code {e.returncode}")
        print(f"Ensure that individual CellChat objects exist and share common metadata.")
        print(f"{'=' * 60}\n")