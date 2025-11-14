# Filename: data_prc.py

# Description: This script processes ship data to calculate latent and sensible heat fluxes
#              using Dr. Bourassa's MFT calculations. It reads input CSV files, computes
#              rolling statistics, performs flux calculations, and writes the results to output CSV files.

import os
import csv
import time
import logging
import numpy as np
import pandas as pd
from multiprocessing import Pool
import tempfile
import shutil
from config import *
from constants import *
from utils import *
from MFT23 import *

logger = setup_logger('data_prc', f'{logs_dir}/data_prc.log', level=logging.INFO)

# Constants for incremental writing
BATCH_SIZE = 500  # Number of rows to write at a time for chunk size

####################################################################################################################
# Step 1: Checks that each element in the data array is within reasonable bounds, validate it
####################################################################################################################
def check_array(arr, min_val, max_val, replacement_val):
    modified_arr = np.where(arr < min_val, replacement_val, arr)
    modified_arr = np.where(modified_arr > max_val, replacement_val, modified_arr)
    return modified_arr

# Validates that all inputs are not NaN
def is_valid_input(*args):
    return all(not np.isnan(x) for x in args)

# Set input and output directories
input_csvs = directory_destination
create_directory(output_csvs)

####################################################################################################################
## Step 2: Function to process a single chunk (to be run in parallel)
####################################################################################################################
def process_chunk(args):
    chunk, chunk_idx, number_of_chunks, window_radius, headers, output_filename = args
    # Define temporary file name with ship, filename, and chunk index
    temp_filename = f"{output_filename[:-4]}_chunk{chunk_idx}.csv"
    
    # Check if temporary file already exists (in case of crashes/terminations)
    if os.path.exists(temp_filename):
        logger.info(f'Skipping chunk {chunk_idx+1}/{number_of_chunks} as {temp_filename} already exists')
        return temp_filename

    logger.info(f'Processing chunk {chunk_idx+1}/{number_of_chunks} in process {os.getpid()}')
    chunk_start = time.time()
    
    # Extract variables from chunk
    airTemp = chunk['in_T'].to_numpy().astype(float)
    hum = chunk['in_RH'].to_numpy().astype(float) / 100
    pressure = chunk['in_P'].to_numpy().astype(float) * 100
    waterTemp = chunk['in_TS'].to_numpy().astype(float)
    truew = chunk['in_SPD'].to_numpy().astype(float)
    time_list = chunk['time'].to_numpy()
    latitude = chunk['latitude'].to_numpy()
    longitude = chunk['longitude'].to_numpy()
    ship_id = chunk['platform_call_sign'].to_numpy()

    # Validate and clean data arrays
    airTemp = check_array(airTemp, -60, 100, np.nan)
    waterTemp = check_array(waterTemp, -4, 100, np.nan)
    hum = check_array(hum, 0, 1, np.nan)
    pressure = check_array(pressure, 0, 120000, np.nan)
    truew = check_array(truew, 0, 100, np.nan)

    ####################################################################################################################
    ## Step 3: Compute rolling means and stds
    ####################################################################################################################
    rolling_window = 2 * window_radius + 1
    min_periods = 3

    # Air temperature rolling stats
    airT_series = pd.Series(airTemp)
    airT_mean = airT_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
    airT_std = airT_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

    # Water temperature rolling stats
    waterT_series = pd.Series(waterTemp)
    waterT_mean = waterT_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
    waterT_std = waterT_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

    # Humidity rolling stats
    hum_series = pd.Series(hum)
    hum_mean = hum_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
    hum_std = hum_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

    # Pressure rolling stats
    press_series = pd.Series(pressure)
    press_mean = press_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
    press_std = press_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

    # True wind speed rolling stats
    truew_series = pd.Series(truew)
    truew_mean = truew_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
    truew_std = truew_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

    ####################################################################################################################
    ## Step 4: Process each data point in the chunk and write incrementally — FIXED FOR NO DUPLICATES
    ####################################################################################################################
    row_buffer = []
    with open(temp_filename, 'w', newline='') as temp_file:
        writer = csv.writer(temp_file)

        # === CALCULATE VALID OUTPUT RANGE (NO OVERLAP DUPLICATES) ===
        total_rows = len(truew)
        overlap = window_size // 2  # From constants.py

        # Start after leading overlap + window_radius
        start_k = overlap + window_radius
        # End before trailing overlap (unless last chunk)
        if chunk_idx == number_of_chunks - 1:
            end_k = total_rows - window_radius  # Last chunk: go to end
        else:
            end_k = total_rows - overlap - window_radius

        if start_k >= end_k:
            logger.warning(f"Chunk {chunk_idx+1}: No output rows after overlap trim. Skipping.")
            return temp_filename

        for k in range(start_k, end_k):
            lhf, shf, m_o, tau = [], [], [], []

            mean_tspd = truew_mean[k]
            std_tspd = truew_std[k]

            if np.isnan(mean_tspd) or std_tspd < 0:
                row = [
                    time_list[k], ship_id[k], latitude[k], longitude[k],
                    -9999, -9999, -9999, -9999,
                    -9999, -9999, -9999,
                    -9999, -9999, -9999, -9999, -9999
                ]
            else:
                tspd_gauss = np.random.normal(mean_tspd, std_tspd, n)
                press_gauss = np.random.normal(press_mean[k], press_std[k], n)
                hum_gauss = np.random.normal(hum_mean[k], hum_std[k], n)
                airT_gauss = np.random.normal(airT_mean[k], airT_std[k], n)
                waterT_gauss = np.random.normal(waterT_mean[k], waterT_std[k], n)

                for j in range(n):
                    if is_valid_input(tspd_gauss[j], press_gauss[j], hum_gauss[j], airT_gauss[j], waterT_gauss[j]):
                        try:
                            flux = mft_fluxes(
                                dyn_in_prm, tspd_gauss[j], dyn_in_val2, sfc_current1, sfc_current2,
                                convect, press_gauss[j], air_moist_prm, hum_gauss[j], sfc_moist_prm,
                                sfc_moist_val, salinity, ss_prm, ss_val, airT_gauss[j], sst_prm,
                                waterT_gauss[j], ref_ht_wind, ref_ht_tq, z_wanted, astab, eqv_neut,
                                net_heat_flux, warn, flux_model, z0_mom_prm, z0_theta_q_prm, stable_prm,
                                oil_fract_area, dimensionless_m_o_length, zo_m, missing
                            )
                            m_o.append(flux[7])
                            lhf.append(flux[2])
                            shf.append(flux[1])
                            tau.append(flux[3][0])
                        except Exception as e:
                            logger.error(f'Error in flux calculation at index {k}: {e}')
                    else:
                        if j == 0:
                            break

            if len(shf) > 0:
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
                row = [
                    time_list[k], ship_id[k], latitude[k], longitude[k],
                    -9999, -9999, -9999, -9999,
                    -9999, -9999, -9999,
                    -9999, -9999, -9999, -9999, -9999
                ]
            row_buffer.append(row)

            if len(row_buffer) >= BATCH_SIZE:
                writer.writerows(row_buffer)
                temp_file.flush()
                os.fsync(temp_file.fileno())
                row_buffer = []
                logger.debug(f'Wrote {BATCH_SIZE} rows to {temp_filename}')

        if row_buffer:
            writer.writerows(row_buffer)
            temp_file.flush()
            os.fsync(temp_file.fileno())
            logger.debug(f'Wrote {len(row_buffer)} final rows to {temp_filename}')

    chunk_end = time.time()
    logger.info(f'Completed chunk {chunk_idx+1}/{number_of_chunks} in {chunk_end - chunk_start:.2f} sec')
    return temp_filename


####################################################################################################################
## Step 5: Main function to process data and calculate fluxes using Dr. Bourassa's MFT calculations
####################################################################################################################
def process_data_prl():
    logger.info("Starting the data processing and flux calculation")
    
    start_time = time.time()
    file_count = 0

    create_directory(output_csvs)

    for ship in ships:
        ship_start = time.time()
        directory_path = f'{input_csvs}/{ship}'
        ship_output_dir = f'{output_csvs}/{ship}'
        create_directory(ship_output_dir)
       
        for filename in os.listdir(directory_path):
            file_count += 1
            loop_start = time.time()
            file_path = os.path.join(directory_path, filename)
            
            try:
                data = pd.read_csv(file_path, skiprows=[1, 2])
                logger.info(f"Loaded file {file_path} with {len(data)} rows")
            except Exception as e:
                logger.error(f"Failed to read {file_path}: {e}")
                continue

            # === DEFINE window_radius BEFORE chunking ===
            window_radius = window_size // 2  # e.g., 5

            # Define chunk parameters
            overlap = window_size // 2
            if chunk_size <= 2 * overlap:
                logger.error(f"Invalid chunk_size ({chunk_size}) or overlap ({overlap}) for file {file_path}.")
                continue
            if chunk_size <= 0 or overlap < 0:
                logger.error(f"Invalid chunk_size ({chunk_size}) or overlap ({overlap}) for file {file_path}.")
                continue
            
            effective_chunk_size = chunk_size - 2 * overlap
            number_of_chunks = max(1, (len(data) // effective_chunk_size) + (1 if len(data) % effective_chunk_size > 0 else 0))
            logger.info(f"Processing file {file_path} with {number_of_chunks} chunks")

            output_filename = os.path.join(ship_output_dir, f"{filename[:-4]}_processed.csv")

            if os.path.exists(output_filename):
                try:
                    output_data = pd.read_csv(output_filename)
                    if len(output_data) >= len(data) - 2 * window_radius:
                        logger.info(f"Skipping {output_filename} as it appears complete with {len(output_data)} rows")
                        continue
                except Exception as e:
                    logger.warning(f"Could not verify {output_filename}: {e}. Proceeding with processing.")

            with open(output_filename, "w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(headers)

            # === Prepare chunks ===
            chunks = []
            for start_idx in range(0, len(data), effective_chunk_size):
                chunk_start_idx = max(0, start_idx - overlap)
                chunk_end_idx = min(len(data), start_idx + chunk_size + overlap)
                chunk = data.iloc[chunk_start_idx:chunk_end_idx].copy()
                chunks.append((chunk, len(chunks), number_of_chunks, window_radius, headers, output_filename))

            # === Parallel processing ===
            with Pool(processes=proc_num) as pool:
                temp_files = pool.map(process_chunk, chunks)

            # === Merge temp files in order ===
            with open(output_filename, "a", newline="") as file:
                writer = csv.writer(file)
                for temp_file_path in sorted(temp_files, key=lambda x: int(x.split('_chunk')[-1].split('.csv')[0])):
                    try:
                        with open(temp_file_path, "r", newline="") as temp_file:
                            reader = csv.reader(temp_file)
                            writer.writerows(reader)

                        file.flush()
                        os.fsync(file.fileno())
                        logger.debug(f'Appended {temp_file_path}')
                    except Exception as e:
                        logger.error(f"Failed to append {temp_file_path}: {e}")
                    finally:
                        if os.path.exists(temp_file_path):
                            os.unlink(temp_file_path)

            loop_end = time.time()
            logger.info(f'File {filename} processed in {loop_end - loop_start:.2f} sec')

        ship_end = time.time()
        logger.info(f'Ship {ship} completed in {ship_end - ship_start:.2f} sec')
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    logger.info(f"Data processing completed in {int(minutes)} min {seconds:.2f} sec")