"""
Log Parsing & Report Generation Module.

This module automates the extraction of pipeline metrics (cell counts,
rho values, runtimes) from text-based logs into structured CSV reports.

Author: Loai Eletr
"""

import re
import pandas as pd
from pathlib import Path
from datetime import datetime

# ==============================================================================
# 1. REGEX EXTRACTION ENGINE
# ==============================================================================

def extract_patterns(content: str, patterns: dict) -> dict:
    """
    Extract structured data from text using a dictionary of regex patterns.

    This utility iterates through a mapping of keys to regex strings,
    extracting the first capturing group for each. It includes specialized
    boolean logic for detecting pipeline 'skips' or specific flags.

    Parameters
    ----------
    content : str
        The raw text string (e.g., a log file or report) to be searched.
    patterns : dict
        A dictionary where keys are the desired output labels and values
        are regex strings containing at least one capturing group.

    Returns
    -------
    dict
        A dictionary of extracted results. If a pattern is not found,
        the value returns `None`.

    Notes
    -----
    The key 'Skip_Detected' is treated as a boolean flag. If the pattern
    is found anywhere in the content, it returns `True`, otherwise `False`.
    """
    # Ensure content is a string to prevent regex crashes
    if not isinstance(content, str):
        return {key: None for key in patterns}

    extracted = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, content)

        # Specialized handling for boolean pipeline flags
        if key == "Skip_Detected":
            extracted["Cleaning_Skipped"] = bool(match)
        else:
            # Extract first group and clean whitespace
            extracted[key] = match.group(1).strip() if match else None

    return extracted

# ==============================================================================
# 2. TIME & BENCHMARKING UTILITIES
# ==============================================================================

def calculate_duration(start_str: str, end_str: str) -> str | None:
    """
    Convert timestamp strings into human-readable duration strings.

    This utility calculates the delta between two timestamps and formats
    the output dynamically into seconds (s), minutes (m), or hours (h)
    depending on the magnitude of the duration.

    Parameters
    ----------
    start_str : str
        The starting timestamp in '%Y-%m-%d %H:%M:%S.%f' format.
    end_str : str
        The ending timestamp in '%Y-%m-%d %H:%M:%S.%f' format.

    Returns
    -------
    str or None
        A formatted string (e.g., '45.2s', '12.5m', '2.1h'). Returns
        None if inputs are missing or the format is invalid.
    """
    if not start_str or not end_str:
        return None

    try:
        # Standard format used by datetime.now() when including microseconds
        fmt = "%Y-%m-%d %H:%M:%S.%f"
        start = datetime.strptime(start_str, fmt)
        end = datetime.strptime(end_str, fmt)

        total_sec = (end - start).total_seconds()

        # Dynamic formatting logic
        if total_sec < 60:
            duration = f"{round(total_sec, 1)}s"
        elif total_sec < 3600:
            duration = f"{round(total_sec / 60, 2)}m"
        else:
            duration = f"{round(total_sec / 3600, 2)}h"

        return duration

    except (ValueError, TypeError, Exception):
        # Specific catch for formatting errors (ValueError) or non-string inputs
        return None

# ==============================================================================
# 3. SAMPLE-LEVEL COORDINATION
# ==============================================================================

def parse_sample_log(log_path: str | Path) -> dict:
    """
    Coordinate extraction and formatting of metrics from a single pipeline log.

    This function reads a specific log file, applies regex patterns to extract
    biological and technical metrics (cell counts, doublet rates, rho),
    and calculates the total processing time.

    Parameters
    ----------
    log_path : str or Path
        The path to the text-based log file.

    Returns
    -------
    dict
        A dictionary containing extracted metrics and metadata. Keys include:
        'Sample', 'Initial_Cells', 'Post_Filter_Cells', 'Final_Cells',
        'Ambient_Rho', 'Doublet_Rate', 'Duration', and 'Cleaning_Skipped'.
    """
    path = Path(log_path)

    # Defensive check: ensure the log exists before attempting to read
    if not path.exists():
        print(f"⚠️ Warning: Log file not found at {path}")
        return {"Sample": path.parent.name, "Error": "File not found"}

    with open(path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Define patterns specific to your pipeline's log output format
    patterns = {
        "Initial_Cells": r"Cells before filtering:\s+(\d+)",
        "Post_Filter_Cells": r"Cells after filtering:\s+(\d+)",
        "Final_Cells": r"Final cell count:\s+(\d+)",
        "Ambient_Rho": r"Estimated global rho of\s+([\d.]+)",
        "Doublet_Rate": r"\(([\d.]+)%\) doublets called",
        "Start_Time": r"Processing started for .* at ([\d\-\s:.]+)",
        "End_Time": r"Processing ended for .* at ([\d\-\s:.]+)",
        "Skip_Detected": r"⚠️ SKIPPING SoupX/Doublet removal"
    }

    # Initialize results with the sample name (assumed to be the parent directory)
    results = {"Sample": path.parent.name}

    # Use the extract_patterns utility (defined previously)
    data = extract_patterns(content, patterns)
    results.update(data)

    # Logic for handling skipped biological corrections
    if results.get("Cleaning_Skipped"):
        results["Ambient_Rho"] = "SKIPPED"
        results["Doublet_Rate"] = "SKIPPED"

    # Use the calculate_duration utility (defined previously)
    duration = calculate_duration(results.get("Start_Time"), results.get("End_Time"))
    results["Duration"] = duration

    return results


# ==============================================================================
# 4. GLOBAL REPORTING (MAIN ENTRY POINT)
# ==============================================================================

def run_report_generation(log_directory: Path, output_filename: str = "Summary.csv"):
    """
    Aggregate all sample logs in a directory and generate a consolidated report.

    This function recursively finds all processing log files, parses their
    content into a structured format, performs data type conversion for
    numeric analysis, and exports the final summary to a CSV file.

    Parameters
    ----------
    log_directory : Path
        The root directory where preprocessing results and logs are stored.
    output_filename : str, default "Summary.csv"
        The name of the resulting summary report.

    Returns
    -------
    pd.DataFrame or None
        Returns the summary DataFrame if successful, otherwise None.
    """
    # 1. Identify all log files recursively
    # Standardizing on a suffix like '_processing_log.txt' ensures we don't
    # accidentally parse irrelevant text files.
    log_files = list(log_directory.glob("**/*_processing_log.txt"))

    if not log_files:
        print("⚠️ No log files found. Skipping report generation.")
        return None

    print(f"📊 Found {len(log_files)} logs. Compiling report...")

    # 2. Process all logs into a list of dictionaries
    # Uses the 'parse_sample_log' coordinator function defined previously
    summary_data = [parse_sample_log(log) for log in log_files]
    df = pd.DataFrame(summary_data)

    if df.empty:
        print("⚠️ No data extracted from logs. Report not generated.")
        return None

    # 3. Data Cleaning & Type Casting
    # Converting these to numeric ensures that Excel/R recognizes them as numbers,
    # allowing for immediate plotting of filtering efficiency.
    numeric_cols = ['Initial_Cells', 'Post_Filter_Cells', 'Final_Cells']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # 4. Final Formatting
    # Select and order columns for a consistent biological report
    preferred_order = [
        'Sample', 'Initial_Cells', 'Post_Filter_Cells', 'Final_Cells',
        'Ambient_Rho', 'Doublet_Rate', 'Duration'
    ]

    # Only keep columns that are actually present in the data
    final_cols = [c for c in preferred_order if c in df.columns]
    df_final = df[final_cols].sort_values(by="Sample")

    # 5. Serialization
    report_path = log_directory / output_filename
    df_final.to_csv(report_path, index=False)

    print("-" * 64)
    print(f"✅ Summary Report generated at: {report_path}")
    print("-" * 64)

    return df_final