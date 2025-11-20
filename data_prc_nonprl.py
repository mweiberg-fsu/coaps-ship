## File: data_prc.py
## Description: Optimized data processing and flux calculations for ship data, including uncertainty estimations and output generation.

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

logger = setup_logger('data_prc', f'{logs_dir}/data_prc.log', level=logging.INFO)  # Set to INFO for production

# Checks that each element in the data array is within reasonable bounds for the variable
def check_array(arr, min_val, max_val, replacement_val):
    modified_arr = np.where(arr < min_val, replacement_val, arr)
    modified_arr = np.where(modified_arr > max_val, replacement_val, modified_arr)
    return modified_arr

# Validates that all inputs are not NaN
def is_valid_input(*args):
    return all(not np.isnan(x) for x in args)  # Updated to handle np.nan

# Set input and output directories
input_csvs = directory_destination
create_directory(output_csvs)

# Main function to process data and calculate fluxes
def process_data():
    logger.info("Starting the data processing and flux calculation")
    
    start_time = time.time()  # Start timer
    file_count = 0  # Initialize file counter

    for ship in ships:  # Iterate through each ship
        ship_start = time.time()  # Start timer for ship
        directory_path = f'{input_csvs}/{ship}'  # Directory for current ship
       
        # Iterate through each file in ship directory
        for filename in os.listdir(directory_path):
            file_count += 1
            loop_start = time.time()
            file_path = os.path.join(directory_path, filename)
            
            # Read the CSV file once, skipping metadata rows
            try:
                data = pd.read_csv(file_path, skiprows=[1, 2])  # Adjust skiprows based on file structure
                logger.info(f"Loaded file {file_path} with {len(data)} rows")
            except Exception as e:
                logger.error(f"Failed to read {file_path}: {e}")
                continue

            # Define overlap as in original code
            overlap = window_size // 2  # Define overlap size
            # Validate chunk_size and overlap
            if chunk_size <= 2 * overlap:
                logger.error(f"Invalid chunk_size ({chunk_size}) or overlap ({overlap}) for file {file_path}. "
                             f"Ensure chunk_size > 2 * overlap.")
                continue
            if chunk_size <= 0 or overlap < 0:
                logger.error(f"Invalid chunk_size ({chunk_size}) or overlap ({overlap}) for file {file_path}. "
                             f"chunk_size must be positive and overlap non-negative.")
                continue

            # Calculate number of chunks
            effective_chunk_size = chunk_size - 2 * overlap
            number_of_chunks = max(1, (len(data) // effective_chunk_size) + (1 if len(data) % effective_chunk_size > 0 else 0))
            logger.info(f"Processing file {file_path} with {number_of_chunks} chunks")

            # Initialize output CSV with headers
            csv_filename = f'{output_csvs}/{filename[:-4]}_processed.csv'
            with open(csv_filename, "w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(headers)

            # Process file in chunks with overlap
            chunk_count = 0
            row_buffer = []  # Buffer to hold rows before writing to CSV

            window_radius = 5  # Hardcoded based on original range(5, len-5); adjust if window_size changes

            # Manually slice the DataFrame into chunks
            for start_idx in range(0, len(data), effective_chunk_size):
                chunk_count += 1
                chunk_start = time.time()
                logger.info(f'Processing chunk {chunk_count}/{number_of_chunks}')

                # Define chunk boundaries with overlap
                chunk_start_idx = max(0, start_idx - overlap)  # Include previous overlap
                chunk_end_idx = min(len(data), start_idx + chunk_size + overlap)  # Include next overlap
                chunk = data.iloc[chunk_start_idx:chunk_end_idx].copy()  # Use .copy() to avoid SettingWithCopyWarning

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

                # Validate and clean data arrays, replace invalid with np.nan
                airTemp = check_array(airTemp, -60, 100, np.nan)
                waterTemp = check_array(waterTemp, -4, 100, np.nan)
                hum = check_array(hum, 0, 1, np.nan)
                pressure = check_array(pressure, 0, 120000, np.nan)
                truew = check_array(truew, 0, 100, np.nan)

                # Compute rolling means and stds for each variable
                rolling_window = 2 * window_radius + 1
                min_periods = 3

                airT_series = pd.Series(airTemp)
                airT_mean = airT_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
                airT_std = airT_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

                waterT_series = pd.Series(waterTemp)
                waterT_mean = waterT_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
                waterT_std = waterT_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

                hum_series = pd.Series(hum)
                hum_mean = hum_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
                hum_std = hum_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

                press_series = pd.Series(pressure)
                press_mean = press_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
                press_std = press_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

                truew_series = pd.Series(truew)
                truew_mean = truew_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
                truew_std = truew_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

                # Process each data point in the chunk
                for k in range(window_radius, len(truew) - window_radius):
                    lhf, shf, m_o, tau = [], [], [], []  # Initialize lists for flux calculations

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
                        # Precompute Gaussian samples for all iterations
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
                                logger.warning(f'Skipped iteration due to invalid inputs at index {k}')
                                if j == 0:  # If first try invalid, assume all will be, skip
                                    break
                    
                    # Prepare row for output
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

                # Write buffered rows to CSV after processing chunk
                if row_buffer:
                    with open(csv_filename, "a", newline="") as file:
                        writer = csv.writer(file)
                        writer.writerows(row_buffer)
                    row_buffer = []  # Clear buffer after writing
                chunk_end = time.time()
                chunk_total = chunk_end - chunk_start
                logger.info(f'Time to complete chunk: {chunk_total}')

            # End of file processing
            loop_end = time.time()
            loop_total = loop_end - loop_start
            logger.info(f'Time to complete file: {loop_total}')

        # End of ship processing
        ship_end = time.time()
        ship_total = ship_end - ship_start
        logger.info(f'Time to complete ship {ship}: {ship_total}')
    
    # End of all processing
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    logger.info(f"Data processing completed in {int(minutes)} min {seconds:.2f} sec")