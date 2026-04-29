"""
Orchestration Script: 01_Preprocessing
Project: Colorectal Cancer (CRC) Single-Cell Atlas
Stage: Raw Data to Cleaned AnnData (h5ad)

This script manages the high-level workflow for the preprocessing phase.
It scans the raw data directory, initializes the R environment for specialized
tools, and iterates through samples—applying QC, normalization, and
initial clustering via the 'src.preprocessing' module.

Author: Loai Eletr
"""

import sys
import os
from datetime import datetime
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr

# 1. PATH RESOLUTION
# Ensures the project root is in the system path so 'config' and 'src'
# are discoverable regardless of where the script is invoked.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

# 2. PROJECT IMPORTS
import config
from config import PreprocessingConfig
import rpy2.rinterface_lib.callbacks as r_callbacks
from src.preprocessing import run_single_sample_preprocessing
from src.utils import load_r_environment, save_adata, r_write_console
from src.qc_utils import run_report_generation

# Initialize global configuration object
config = PreprocessingConfig()


def get_sample_map(folder_path):
    """
    Scans the data directory to create a mapping of Sample IDs to file paths.

    Expects folders to follow a naming convention where the suffix
    (e.g., 'N1', 'T2') represents the biological sample tag.

    Returns:
        dict: { 'SampleTag': '/absolute/path/to/data' }
    """
    sample_dict = {}
    base_path = Path(folder_path)

    if not base_path.exists():
        print(f"⚠️ Warning: Data directory not found at {folder_path}")
        return {}

    for item in base_path.iterdir():
        if item.is_dir():
            # Extract tag (e.g., from 'GSM12345_N1' -> 'N1')
            tag = item.name.split('_')[-1]
            sample_dict[tag] = str(item.absolute())
    return sample_dict


def run_preprocessing_pipeline(samples_to_run=None):
    """
    Executes the end-to-end preprocessing loop for a set of samples.

    This function handles:
    - R console callback redirection (for rpy2 logging)
    - Directory creation for results and figures
    - Per-sample logging (capturing all stdout/stderr to disk)
    - Exception handling to prevent single-sample failures from halting the loop.
    """

    # --- R INTEROP SETUP ---
    # Redirect R's console output to Python's logging system
    r_callbacks.consolewrite_print_update = r_write_console
    r_callbacks.consolewrite_warnerror = r_write_console
    load_r_environment(config.R_HOME)

    # --- SAMPLE SELECTION ---
    if samples_to_run is None:
        all_samples = get_sample_map(config.DATA_RAW)
        samples_to_run = all_samples

    # --- MAIN EXECUTION LOOP ---
    for tag, path in samples_to_run.items():

        # Define output sub-directories based on config paths
        SAMPLE_DIR = config.PREPROC_OUT / tag
        SAMPLE_DIR.mkdir(parents=True, exist_ok=True)

        sample_qc_dir = config.FIGURES_ROOT / "qc" / tag
        sample_cluster_dir = config.FIGURES_ROOT / "clustering" / tag
        log_file = SAMPLE_DIR / f"{tag}_processing_log.txt"

        print(f"Processing started for {tag} at {datetime.now()}")

        # ATOMIC LOGGING: Redirect all output (Scanpy, R, and print) to a log file
        with open(log_file, "w", encoding="utf-8") as f:
            with redirect_stdout(f), redirect_stderr(f):
                print(f"Processing started for {tag} at {datetime.now()}")
                print(f">>> Loading Sample {tag} from {path}...")

                try:
                    # Execute the worker function with parameter unpacking from config
                    adata = run_single_sample_preprocessing(
                        # --- Sample Context ---
                        adata_path=path,
                        sample_name=tag,
                        save_dir_qc=sample_qc_dir,
                        save_dir_clustering=sample_cluster_dir,

                        # --- Hyperparameters (Unpacked from config.py) ---
                        **config.GENERAL_PARAMS,
                        **config.QC_PARAMS,
                        **config.NORM_PARAMS,
                        **config.CLUSTERING_PARAMS
                    )

                    # PERSISTENCE: Save the cleaned AnnData object
                    save_adata(adata, SAMPLE_DIR / f"{tag}_cleaned.h5ad")
                    print(f"Status: SUCCESS for {tag}")

                except Exception as e:
                    # ERROR ISOLATION: Ensure the loop continues if one sample fails
                    print(f"ERROR in {tag}: {str(e)}")

                print(f"Processing ended for {tag} at {datetime.now()}")

        print(f"Processing ended for {tag} at {datetime.now()}")


# --- ENTRY POINT ---
if __name__ == "__main__":
    # Scans raw data and triggers the pipeline
    SAMPLES = get_sample_map(config.DATA_RAW)

    if SAMPLES:
        # Run the full preprocessing loop
        run_preprocessing_pipeline(SAMPLES)

        # AGGREGATION: Generate summary metrics (e.g., QC reports) across all samples
        run_report_generation(config.PREPROC_OUT, "Summary.csv")