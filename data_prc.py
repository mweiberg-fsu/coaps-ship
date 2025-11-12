## File: data_prc_nonprl.py
## Description: Optimized data processing and flux calculations for ship data,
##              including Monte Carlo uncertainty estimation.
##              **Now respects config.years – only files for those years are processed.**

import numpy as np
import pandas as pd
import os
import time
import csv
import logging
from config import *
from constants import *
from utils import *
from MFT23 import *

logger = setup_logger('data_prc', f'{logs_dir}/data_prc.log', level=logging.INFO)

# ----------------------------------------------------------------------
# Helper functions (unchanged)
# ----------------------------------------------------------------------
def check_array(arr, min_val, max_val, replacement_val):
    """Clamp values to a reasonable range."""
    modified_arr = np.where(arr < min_val, replacement_val, arr)
    modified_arr = np.where(modified_arr > max_val, replacement_val, modified_arr)
    return modified_arr


def is_valid_input(*args):
    """Return True if none of the arguments are NaN."""
    return all(not np.isnan(x) for x in args)


# ----------------------------------------------------------------------
# Directories
# ----------------------------------------------------------------------
input_csvs = directory_destination
create_directory(output_csvs)

# ----------------------------------------------------------------------
# MAIN PROCESSING FUNCTION
# ----------------------------------------------------------------------
def process_data():
    logger.info("Starting the data processing and flux calculation")
    logger.info(f"Only years listed in config.years will be processed: {years}")

    start_time = time.time()
    file_count = 0

    for ship in ships:
        ship_start = time.time()
        directory_path = f'{input_csvs}/{ship}'

        if not os.path.isdir(directory_path):
            logger.warning(f"Ship directory not found: {directory_path}")
            continue

        logger.info(f"Scanning directory for ship {ship}: {directory_path}")

        # --------------------------------------------------------------
        # Iterate over EVERY file in the ship folder
        # --------------------------------------------------------------
        for filename in os.listdir(directory_path):
            if not filename.lower().endswith('.csv'):
                continue

            # ----- Extract year from filename (e.g. KAOU2026.csv → 2026) -----
            try:
                file_year_str = filename.replace(ship, '').replace('.csv', '')
                file_year = int(file_year_str)
            except ValueError:
                logger.warning(f"Cannot parse year from filename, skipping: {filename}")
                continue

            # ----- ONLY PROCESS YEARS THAT ARE IN config.years -----
            if file_year not in years:
                logger.info(f"Skipping {filename} (year {file_year} not in config.years)")
                continue

            file_count += 1
            loop_start = time.time()
            file_path = os.path.join(directory_path, filename)

            # ----- Load the raw CSV (skip the two metadata lines) -----
            try:
                data = pd.read_csv(file_path, skiprows=[1, 2])
                logger.info(f"Loaded {filename} – {len(data)} rows")
            except Exception as e:
                logger.error(f"Failed to read {file_path}: {e}")
                continue

            # ------------------------------------------------------------------
            # Chunking / overlapping window setup (unchanged)
            # ------------------------------------------------------------------
            overlap = window_size // 2
            if chunk_size <= 2 * overlap:
                logger.error(f"Invalid chunk_size ({chunk_size}) vs overlap ({overlap}) for {filename}")
                continue

            effective_chunk_size = chunk_size - 2 * overlap
            number_of_chunks = max(1,
                                   (len(data) // effective_chunk_size) +
                                   (1 if len(data) % effective_chunk_size else 0))
            logger.info(f"Will process {filename} in {number_of_chunks} chunk(s)")

            # Output CSV (one per input file)
            csv_filename = f'{output_csvs}/{filename[:-4]}_processed.csv'
            with open(csv_filename, "w", newline="") as f_out:
                writer = csv.writer(f_out)
                writer.writerow(headers)

            chunk_count = 0
            row_buffer = []                     # collect rows before bulk write
            window_radius = 5                   # hard-coded from original code

            # ------------------------------------------------------------------
            # CHUNK LOOP
            # ------------------------------------------------------------------
            for start_idx in range(0, len(data), effective_chunk_size):
                chunk_count += 1
                chunk_start = time.time()
                logger.debug(f'  → Chunk {chunk_count}/{number_of_chunks}')

                # ----- Define chunk boundaries (with overlap) -----
                end_idx = min(start_idx + chunk_size, len(data))
                chunk = data.iloc[start_idx:end_idx]

                # ----- Extract columns (once per chunk) -----
                time_list   = chunk['time'].values
                ship_id     = chunk['platform_call_sign'].values
                latitude    = chunk['latitude'].values
                longitude   = chunk['longitude'].values
                truew       = chunk['in_SPD'].values
                press       = chunk['in_P'].values
                hum         = chunk['in_RH'].values
                airT        = chunk['in_T'].values
                waterT      = chunk['in_TS'].values

                # ----- Rolling-window means & stds (window = 11) -----
                truew_mean = pd.Series(truew).rolling(window=window_size, center=True, min_periods=1).mean().values
                truew_std  = pd.Series(truew).rolling(window=window_size, center=True, min_periods=1).std().fillna(0).values
                press_mean = pd.Series(press).rolling(window=window_size, center=True, min_periods=1).mean().values
                press_std  = pd.Series(press).rolling(window=window_size, center=True, min_periods=1).std().fillna(0).values
                hum_mean   = pd.Series(hum).rolling(window=window_size, center=True, min_periods=1).mean().values
                hum_std    = pd.Series(hum).rolling(window=window_size, center=True, min_periods=1).std().fillna(0).values
                airT_mean  = pd.Series(airT).rolling(window=window_size, center=True, min_periods=1).mean().values
                airT_std   = pd.Series(airT).rolling(window=window_size, center=True, min_periods=1).std().fillna(0).values
                waterT_mean= pd.Series(waterT).rolling(window=window_size, center=True, min_periods=1).mean().values
                waterT_std = pd.Series(waterT).rolling(window=window_size, center=True, min_periods=1).std().fillna(0).values

                # ----- Monte-Carlo loop over every row in the chunk -----
                for k in range(len(chunk)):
                    shf, lhf, tau, m_o = [], [], [], []

                    # Generate n Gaussian realisations of the inputs
                    tspd_gauss  = np.random.normal(truew_mean[k], truew_std[k], n)
                    press_gauss = np.random.normal(press_mean[k], press_std[k], n)
                    hum_gauss   = np.random.normal(hum_mean[k],   hum_std[k],   n)
                    airT_gauss  = np.random.normal(airT_mean[k],  airT_std[k],  n)
                    waterT_gauss= np.random.normal(waterT_mean[k],waterT_std[k],n)

                    for j in range(n):
                        if not is_valid_input(tspd_gauss[j], press_gauss[j],
                                              hum_gauss[j], airT_gauss[j],
                                              waterT_gauss[j]):
                            continue

                        try:
                            flux = mft_fluxes(
                                dyn_in_prm, tspd_gauss[j], dyn_in_val2,
                                sfc_current1, sfc_current2,
                                convect, press_gauss[j], air_moist_prm,
                                hum_gauss[j], sfc_moist_prm, sfc_moist_val,
                                salinity, ss_prm, ss_val,
                                airT_gauss[j], sst_prm, waterT_gauss[j],
                                ref_ht_wind, ref_ht_tq, z_wanted, astab,
                                eqv_neut, net_heat_flux, warn, flux_model,
                                z0_mom_prm, z0_theta_q_prm, stable_prm,
                                oil_fract_area, dimensionless_m_o_length,
                                zo_m, missing
                            )
                            # flux ordering (from original code): [?, shf, lhf, tau, ?, ?, ?, dmo, ...]
                            m_o.append(flux[7])
                            lhf.append(flux[2])
                            shf.append(flux[1])
                            tau.append(flux[3][0])
                        except Exception as e:
                            logger.debug(f'Flux error at row {k} iter {j}: {e}')

                    # ----- Build output row -----
                    if shf:   # at least one valid MC realisation
                        row = [
                            time_list[k], ship_id[k], latitude[k], longitude[k],
                            np.std(shf), np.mean(shf),
                            np.std(lhf), np.mean(lhf),
                            np.std(tau), np.mean(tau),
                            np.mean(m_o),
                            airT_std[k], waterT_std[k],
                            hum_std[k], press_std[k], truew_std[k]
                        ]
                    else:
                        row = [time_list[k], ship_id[k], latitude[k], longitude[k],
                               -9999, -9999, -9999, -9999,
                               -9999, -9999, -9999,
                               -9999, -9999, -9999, -9999, -9999]

                    row_buffer.append(row)

                # ----- Write the whole chunk at once -----
                if row_buffer:
                    with open(csv_filename, "a", newline="") as f_out:
                        writer = csv.writer(f_out)
                        writer.writerows(row_buffer)
                    row_buffer.clear()

                chunk_time = time.time() - chunk_start
                logger.debug(f'Chunk {chunk_count} finished in {chunk_time:.2f}s')

            # ----- End of file -----
            file_time = time.time() - loop_start
            logger.info(f'Finished {filename} in {file_time:.2f}s')

        ship_time = time.time() - ship_start
        logger.info(f'Ship {ship} total time: {ship_time:.2f}s')

    total_time = time.time() - start_time
    mins, secs = divmod(total_time, 60)
    logger.info(f'=== ALL PROCESSING DONE in {int(mins)} min {secs:.2f} sec ===')