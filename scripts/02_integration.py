"""
Orchestration Script: 02_Integration
Project: Colorectal Cancer (CRC) Single-Cell Atlas
Stage: Individual h5ad Files to Unified Integrated Atlas

This script manages the second stage of the pipeline: Batch Correction.
It aggregates all 'cleaned' samples from Stage 01 and utilizes the
'run_full_integration_pipeline' to project them into a shared latent space,
ensuring that biological signals are preserved while technical batch
effects between patients/runs are minimized.

Author: Loai Eletr
"""

from src.integration import run_full_integration_pipeline
import config

# ==============================================================================
# MAIN EXECUTION BLOCK
# ==============================================================================

if __name__ == "__main__":
    """
    Entry point for the integration workflow.

    The workflow follows these steps:
    1. Loads hyperparameters from IntegrationConfig.
    2. Dynamically locates all cleaned h5ad files from the preprocessing output.
    3. Triggers the integration engine (e.g., scVI) using 
       parameter unpacking to maintain a clean execution logic.
    """

    # 1. INITIALIZATION
    # Loads specific integration constants (e.g., categorical_covariate_keys)
    integr_cfg = config.IntegrationConfig()

    # 2. FILE DISCOVERY
    # Uses glob pattern matching to find all cleaned AnnData objects.
    # Sorting ensures consistent sample ordering during the merge.
    SAMPLES = sorted(list(config.PREPROC_OUT.glob("**/*_cleaned.h5ad")))

    # 3. PIPELINE EXECUTION
    if SAMPLES:
        print("-" * 64)
        print(f"🚀 Found {len(SAMPLES)} samples. Initializing Atlas Integration...")
        print("-" * 64)

        # We pass the list of file paths to the pipeline.
        # Dictionary unpacking (**) allows for flexible hyperparameter tuning
        # without modifying the orchestration logic in this script.
        run_full_integration_pipeline(
            SAMPLES,
            **integr_cfg.MODEL_PARAMS,  # e.g., latent_dim, n_hidden
            **integr_cfg.EVAL_PARAMS,  # e.g., batch_mixing, bio_conservation
            **integr_cfg.IO_PARAMS,  # e.g., output_dir, final_filename
        )

        print("-" * 64)
        print(f"✅ Integration complete. Unified atlas saved to IO_PARAMS paths.")
        print("-" * 64)
    else:
        # FAIL-SAFE: Alert the user if Stage 01 hasn't been completed yet.
        print("❌ CRITICAL: No cleaned h5ad files found!")
        print(f"Please check your Stage 01 output at: {config.PREPROC_OUT}")