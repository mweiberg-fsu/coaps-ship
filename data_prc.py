# Filename: data_prc.py
#
# Description: This script processes ship data to calculate latent and sensible heat fluxes
#              using Dr. Bourassa's MFT calculations. It reads input CSV files, computes
#              rolling statistics, performs flux calculations, and writes the results to output CSV files.
#
#  *** CRASH-SAFE / RESUMABLE VERSION ***
#  - Temp files are kept until the final CSV is verified.
#  - A .progress JSON file tracks rows_written, chunks_done, total_rows_expected.
#  - On restart the driver skips every chunk that already has a matching temp file.

import os
import csv
import time
import json
import logging
import numpy as np
import pandas as pd
from multiprocessing import Pool
from pathlib import Path
from config import *
from constants import *
from utils import *
from MFT23 import *

logger = setup_logger(
    'data_prc', f'{logs_dir}/data_prc.log', level=logging.INFO)

# ----------------------------------------------------------------------
# Constants for incremental writing
# ----------------------------------------------------------------------
BATCH_SIZE = 500                     # rows per flush inside a chunk
PROGRESS_SUFFIX = ".progress"        # e.g.  myfile_processed.csv.progress
TEMP_PREFIX = "tmp_chunk_"           # not used – we embed the absolute start index

# ----------------------------------------------------------------------
# Helper functions for progress tracking
# ----------------------------------------------------------------------


def _progress_path(out_csv: str) -> Path:
    return Path(out_csv).with_suffix(PROGRESS_SUFFIX)


def _load_progress(out_csv: str) -> dict:
    p = _progress_path(out_csv)
    if p.exists():
        try:
            return json.loads(p.read_text())
        except Exception as e:
            logger.warning(f"Corrupt progress file {p}, starting fresh: {e}")
    return {}


def _save_progress(out_csv: str, prog: dict):
    p = _progress_path(out_csv)
    p.write_text(json.dumps(prog, indent=2))

# ----------------------------------------------------------------------
# Step 1: Validation utilities (unchanged)
# ----------------------------------------------------------------------


def check_array(arr, min_val, max_val, replacement_val):
    modified_arr = np.where(arr < min_val, replacement_val, arr)
    modified_arr = np.where(modified_arr > max_val,
                            replacement_val, modified_arr)
    return modified_arr


def is_valid_input(*args):
    return all(not np.isnan(x) for x in args)

# ----------------------------------------------------------------------
# Step 2: Process a single chunk (Monte-Carlo worker)
# ----------------------------------------------------------------------


def process_chunk(args):
    """
    args = (
        chunk_dict,
        chunk_start_idx,          # absolute start index inside the *valid* array
        total_chunks,
        window_radius,
        headers,
        output_csv_path
    )
    """
    chunk, chunk_start_idx, total_chunks, window_radius, headers, output_csv = args

    # Temp file lives next to the final CSV, name encodes the absolute start index
    temp_filename = f"{Path(output_csv).stem}_chunk{chunk_start_idx}.csv"
    temp_path = Path(temp_filename)

    if temp_path.exists():
        logger.info(
            f"Temp file {temp_path} already present – assuming chunk done")
        return str(temp_path)

    logger.info(
        f"Processing chunk starting at row {chunk_start_idx} (PID {os.getpid()})")
    start_t = time.time()

    # ------------------------------------------------------------------
    # Extract numpy arrays from the chunk dict
    # ------------------------------------------------------------------
    airTemp_mean = chunk['airT_mean']
    airTemp_std = chunk['airT_std']
    waterTemp_mean = chunk['waterT_mean']
    waterTemp_std = chunk['waterT_std']
    hum_mean = chunk['hum_mean']
    hum_std = chunk['hum_std']
    press_mean = chunk['press_mean']
    press_std = chunk['press_std']
    truew_mean = chunk['truew_mean']
    truew_std = chunk['truew_std']
    time_list = chunk['time']
    latitude = chunk['latitude']
    longitude = chunk['longitude']
    ship_id = chunk['ship_id']

    # ------------------------------------------------------------------
    # Monte-Carlo loop (identical to original logic)
    # ------------------------------------------------------------------
    row_buffer = []
    with open(temp_path, 'w', newline='') as temp_file:
        writer = csv.writer(temp_file)
        # header once per temp file
        writer.writerow(headers)

        # ---- safe output range inside the chunk (no overlap with neighbours) ----
        total_rows_in_chunk = len(truew_mean)
        overlap = window_size // 2
        start_k = overlap + window_radius
        end_k = total_rows_in_chunk - overlap - window_radius

        if start_k >= end_k:
            logger.warning(
                f"Chunk {chunk_start_idx}: no usable rows after overlap trim")
            return str(temp_path)

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
                hum_gauss = np.random.normal(hum_mean[k],   hum_std[k],   n)
                airT_gauss = np.random.normal(
                    airTemp_mean[k], airTemp_std[k], n)
                waterT_gauss = np.random.normal(
                    waterTemp_mean[k], waterTemp_std[k], n)

                for j in range(n):
                    if is_valid_input(tspd_gauss[j], press_gauss[j], hum_gauss[j],
                                      airT_gauss[j], waterT_gauss[j]):
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
                            logger.error(f'Flux error at index {k}: {e}')
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
                    airTemp_std[k], waterTemp_std[k],
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

        # flush any remaining rows
        if row_buffer:
            writer.writerows(row_buffer)
            temp_file.flush()
            os.fsync(temp_file.fileno())

    logger.info(
        f"Chunk {chunk_start_idx}/{total_chunks} finished in {time.time()-start_t:.1f}s")
    return str(temp_path)


# ----------------------------------------------------------------------
# Main driver – resumable, parallel processing
# ----------------------------------------------------------------------
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

            # ------------------------------------------------------------------
            # 1. Rolling statistics (full file)
            # ------------------------------------------------------------------
            window_radius = window_size // 2
            rolling_window = 2 * window_radius + 1
            min_periods = 3

            def rolling_stats(series):
                s = pd.Series(series)
                mean = s.rolling(rolling_window, center=True,
                                 min_periods=min_periods).mean()
                std = s.rolling(rolling_window, center=True,
                                min_periods=min_periods).std(ddof=0)
                return mean, std

            airTemp = check_array(data['in_T'].astype(float), -60, 100, np.nan)
            waterTemp = check_array(
                data['in_TS'].astype(float), -4, 100, np.nan)
            hum = check_array(data['in_RH'].astype(float)/100, 0, 1, np.nan)
            pressure = check_array(data['in_P'].astype(
                float)*100, 0, 120000, np.nan)
            truew = check_array(data['in_SPD'].astype(float), 0, 100, np.nan)

            airT_mean, airT_std = rolling_stats(airTemp)
            waterT_mean, waterT_std = rolling_stats(waterTemp)
            hum_mean, hum_std = rolling_stats(hum)
            press_mean, press_std = rolling_stats(pressure)
            truew_mean, truew_std = rolling_stats(truew)

            # ------------------------------------------------------------------
            # 2. Valid (edge-trimmed) slice
            # ------------------------------------------------------------------
            start_idx = rolling_window // 2 + window_radius
            end_idx = len(data) - window_radius
            if start_idx >= end_idx:
                logger.warning(
                    f"File {filename}: Not enough data after edge trim.")
                continue

            valid = slice(start_idx, end_idx)
            time_valid = data['time'].iloc[valid].to_numpy()
            lat_valid = data['latitude'].iloc[valid].to_numpy()
            lon_valid = data['longitude'].iloc[valid].to_numpy()
            ship_id_valid = data['platform_call_sign'].iloc[valid].to_numpy()

            def _extract(arr): return arr.iloc[valid].to_numpy()
            airT_mean_v,   airT_std_v = _extract(
                airT_mean),   _extract(airT_std)
            waterT_mean_v, waterT_std_v = _extract(
                waterT_mean), _extract(waterT_std)
            hum_mean_v,    hum_std_v = _extract(hum_mean),    _extract(hum_std)
            press_mean_v,  press_std_v = _extract(
                press_mean),  _extract(press_std)
            truew_mean_v,  truew_std_v = _extract(
                truew_mean),  _extract(truew_std)

            # ------------------------------------------------------------------
            # 3. Output file + progress init
            # ------------------------------------------------------------------
            output_filename = os.path.join(
                ship_output_dir, f"{filename[:-4]}_processed.csv")

            prog = _load_progress(output_filename)
            total_rows_expected = len(time_valid)

            # Skip if final file already complete
            if Path(output_filename).exists():
                try:
                    if len(pd.read_csv(output_filename)) >= total_rows_expected:
                        logger.info(
                            f"Skipping {output_filename} – already complete")
                        continue
                except Exception:
                    pass

            # Write header once
            if not Path(output_filename).exists():
                with open(output_filename, "w", newline="") as f:
                    csv.writer(f).writerow(headers)

            # Initialise progress the first time we see this file
            if "total_rows_expected" not in prog:
                prog = {
                    "total_rows_expected": total_rows_expected,
                    "rows_written": 0,
                    "chunks_done": []
                }
                _save_progress(output_filename, prog)

            # ------------------------------------------------------------------
            # 4. Build list of chunks that still need work
            # ------------------------------------------------------------------
            points_per_chunk = chunk_size
            n_points = len(time_valid)
            chunks_to_do = []

            for chunk_idx in range(0, n_points, points_per_chunk):
                if chunk_idx in prog["chunks_done"]:
                    temp_name = f"{Path(output_filename).stem}_chunk{chunk_idx}.csv"
                    if Path(temp_name).exists():
                        logger.info(
                            f"Chunk {chunk_idx} already done – skipping")
                        continue
                    else:
                        logger.warning(
                            f"Chunk {chunk_idx} marked done but temp missing – re-computing")

                idx = slice(chunk_idx, chunk_idx + points_per_chunk)
                chunk_data = {
                    'time': time_valid[idx],
                    'ship_id': ship_id_valid[idx],
                    'latitude': lat_valid[idx],
                    'longitude': lon_valid[idx],
                    'truew_mean': truew_mean_v[idx],
                    'truew_std': truew_std_v[idx],
                    'press_mean': press_mean_v[idx],
                    'press_std': press_std_v[idx],
                    'hum_mean': hum_mean_v[idx],
                    'hum_std': hum_std_v[idx],
                    'airT_mean': airT_mean_v[idx],
                    'airT_std': airT_std_v[idx],
                    'waterT_mean': waterT_mean_v[idx],
                    'waterT_std': waterT_std_v[idx],
                }
                args = (
                    chunk_data,
                    chunk_idx,
                    len(range(0, n_points, points_per_chunk)),
                    window_radius,
                    headers,
                    output_filename
                )
                chunks_to_do.append(args)

            if not chunks_to_do:
                logger.info(f"All chunks for {filename} already processed")
                continue

            # ------------------------------------------------------------------
            # 5. Parallel Monte-Carlo
            # ------------------------------------------------------------------
            with Pool(processes=proc_num) as pool:
                temp_files = pool.map(process_chunk, chunks_to_do)

            # ------------------------------------------------------------------
            # 6. Append new temp files to the final CSV
            # ------------------------------------------------------------------
            with open(output_filename, "a", newline="") as out_f:
                writer = csv.writer(out_f)
                for temp_file in temp_files:
                    with open(temp_file, "r") as tf:
                        reader = csv.reader(tf)
                        next(reader, None)               # skip header
                        writer.writerows(reader)

                    # update progress
                    chunk_idx = int(Path(temp_file).stem.split(
                        '_chunk')[-1].split('.')[0])
                    prog["chunks_done"].append(chunk_idx)
                    # -header
                    prog["rows_written"] += sum(1 for _ in open(temp_file)) - 1
                    _save_progress(output_filename, prog)

            # ------------------------------------------------------------------
            # 7. Final verification & clean-up
            # ------------------------------------------------------------------
            if prog["rows_written"] >= prog["total_rows_expected"]:
                actual = len(pd.read_csv(output_filename))
                if actual >= prog["total_rows_expected"]:
                    # remove all temp files for this output
                    stem = Path(output_filename).stem
                    for p in Path(".").glob(f"{stem}_chunk*.csv"):
                        p.unlink(missing_ok=True)
                    _progress_path(output_filename).unlink(missing_ok=True)
                    logger.info(
                        f"Finished {output_filename} – {actual} rows written")
                else:
                    logger.warning(
                        f"Row count mismatch ({actual} vs {prog['total_rows_expected']}) – keeping temps")
            else:
                logger.info(
                    f"Partial progress for {output_filename}: {prog['rows_written']}/{prog['total_rows_expected']} rows")

        logger.info(
            f"Ship {ship} done in {(time.time()-ship_start)/60:.1f} min")

# ----------------------------------------------------------------------
# (Optional) Clean stray temp/progress files from a previous run
# ----------------------------------------------------------------------
# for p in Path(".").glob("*_chunk*.csv"): p.unlink(missing_ok=True)
# for p in Path(".").glob("*.progress"):   p.unlink(missing_ok=True)


if __name__ == "__main__":
    process_data_prl()
