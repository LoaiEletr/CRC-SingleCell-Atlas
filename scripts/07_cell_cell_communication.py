"""
Orchestration Script: 07_cell_cell_communication
Project: Colorectal Cancer (CRC) Single-Cell Atlas
Stage: Cell-Cell Communication (CCC) Analysis

This script serves as the bridge between the Python/Scanpy workflow and
the R/CellChat ecosystem. It handles data partitioning by biological
condition (Normal vs. Tumor), environment preparation, and the execution
of the R-based signaling inference engine.

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
    print(f"🧬 STARTING CELL-CELL COMMUNICATION ANALYSIS")
    print(f"{'=' * 60}")
    start_time = time.time()

    # --- 1. PREPARATION & DATA SPLITTING ---
    # Load R paths and ensure the environment is ready for subprocess execution.
    print(f"[*] Initializing R environment at: {config.R_HOME}")
    load_r_environment(config.R_HOME)

    print("[*] Preparing directories and K-selection configuration...")
    custom_env = CellChatConfig.prepare_environment()

    # Partitioning: CellChat performs best when comparing separate objects
    # (e.g., one for Normal, one for Primary Tumor).
    print("[*] Splitting integrated AnnData into condition-specific files...")
    CellChatConfig.split_data_by_condition(save_adata)
    print(f"    --> Inputs saved to: {CellChatConfig.CELLCHAT_IO['input_dir']}")

    # --- 2. R-SUBSYSTEM INVOCATION ---
    # Launching the R script provided in previous steps.
    # We pass all configuration parameters (Species, Min Cells, IO paths)
    # as command-line arguments to maintain a "Single Source of Truth".
    print(f"\n[!] Launching Rscript: scripts/08_cell_cell_communication.R")
    print(f"    Species: {CellChatConfig.CELLCHAT_PARAMS['species']}")
    print(f"    Min Cells: {CellChatConfig.CELLCHAT_PARAMS['min_cells']}")
    print(f"{'-' * 60}")

    try:
        # The core hand-off to the R CellChat implementation
        subprocess.run([
            "Rscript",
            f"{config.BASE_DIR}/scripts/07_cell_cell_communication.R",
            str(config.BASE_DIR),
            str(CellChatConfig.CELLCHAT_IO["input_dir"]),
            str(CellChatConfig.CELLCHAT_IO["individual_cellchat_dir"]),
            str(CellChatConfig.CELLCHAT_IO["fig_dir"]),
            CellChatConfig.CELLCHAT_PARAMS["species"],
            str(CellChatConfig.CELLCHAT_PARAMS["min_cells"]),
            str(CellChatConfig.CELLCHAT_IO["k_json"])
        ], env=custom_env, check=True)

        # --- 3. TELEMETRY & RESULTS ---
        elapsed = (time.time() - start_time) / 60

        print(f"{'-' * 60}")
        print(f"✅ ANALYSIS COMPLETE in {elapsed:.2f} minutes")
        print(f"📂 RDS Objects: {CellChatConfig.CELLCHAT_IO['individual_cellchat_dir']}")
        print(f"📊 Figures:     {CellChatConfig.CELLCHAT_IO['fig_dir']}")
        print(f"{'=' * 60}\n")

    except subprocess.CalledProcessError as e:
        # ERROR HANDLING: Captures issues in the R script (e.g., missing libraries,
        # lack of signaling in a specific subset) and reports it to the Python log.
        print(f"\n❌ ERROR: Rscript failed with exit code {e.returncode}")
        print(f"Check the R logs above for specific CellChat errors.")
        print(f"{'=' * 60}\n")