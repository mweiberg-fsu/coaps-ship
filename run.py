# File: run.py
# Description: Main script to run data retrieval and processing pipeline

import time
import logging
import os
from imports import *
from config import ships, years, logs_dir, output_csvs, proc_num
from data_get import run_parallel_downloads
from data_prc_optimized import process_ship_year
from utils import setup_logger, create_directory
from tqdm import tqdm


def main():
    # === Setup ===
    create_directory(logs_dir)
    create_directory(output_csvs)
    logger = setup_logger('run', f'{logs_dir}/run.log', level=logging.INFO)
    logger.info("=== STARTING DATA PIPELINE ===")
    print("\nStarting the data retrieval and processing pipeline.\n")

    start_time = time.time()

    # === Download ===
    logger.info(f"Downloading data for ships: {ships}, years: {years}")
    print(
        f"Downloading data for {len(ships)} ship(s) across {len(years)} year(s)...")

    try:
        run_parallel_downloads(max_threads=proc_num)
        logger.info("Data download completed.")
        print("Download complete.\n")
    except Exception as e:
        logger.error(f"Download failed: {e}")
        print(f"Download failed: {e}")
        return

    # === Process ===
    logger.info("Starting 96-core flux processing...")
    print("Processing data...\n")

    total_files = 0
    processed_files = 0

    # FIXED: Use fixed desc for outer loop, dynamic for inner
    for ship in tqdm(ships, desc="Ships", position=0):
        # Inner loop: use lambda or format inside
        for year in tqdm(years, desc="Years", position=1, leave=False):
            current_desc = f"Processing {ship} - Year {year}"
            # Update tqdm description dynamically
            tqdm.write(current_desc)
            logger.info(current_desc)

            input_dir = f"data/input/{ship}"
            if not os.path.exists(input_dir):
                logger.warning(f"Input directory not found: {input_dir}")
                continue

            csv_files = [f for f in os.listdir(
                input_dir) if f.endswith('.csv')]
            total_files += len(csv_files)

            try:
                process_ship_year(ship, str(year))
                processed_files += len(csv_files)
            except Exception as e:
                logger.error(f"Failed {ship}/{year}: {e}", exc_info=True)

    # === Summary ===
    elapsed = time.time() - start_time
    mins, secs = divmod(elapsed, 60)
    summary = f"""
=== PIPELINE COMPLETE ===
Files found:     {total_files}
Files processed: {processed_files}
Time:            {int(mins)} min {secs:.1f} sec
Output:          {output_csvs}
"""
    print(summary)
    logger.info(summary.strip())


if __name__ == "__main__":
    main()
