# data_prc_optimized.py
import os
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count
from constants import *
from config import *
from utils import *
from MFT23 import mft_fluxes
import logging

# Setup
logger = setup_logger('opt', f'{logs_dir}/opt.log', level=logging.INFO)
N_CORES = min(96, cpu_count())  # Use all 96
CHUNK_SIZE = 50_000             # Big chunks = less overhead
N_SAMPLES = n                   # 500 from constants

create_directory(output_csvs)


def process_ship_year(ship, year):
    ship_dir = f"{directory_destination}/{ship}"
    out_dir = f"{output_csvs}/{ship}"
    create_directory(out_dir)

    for filename in sorted(os.listdir(ship_dir)):
        if not filename.endswith('.csv'):
            continue
        filepath = os.path.join(ship_dir, filename)
        outpath = os.path.join(out_dir, filename.replace('.csv', '_flux.csv'))

        if os.path.exists(outpath):
            logger.info(f"Skipping {outpath} (exists)")
            continue

        logger.info(f"Processing {filename} → {outpath}")
        df = pd.read_csv(filepath, skiprows=[1, 2])

        # === 1. Clean & extract ===
        T = df['in_T'].astype(float).values
        TS = df['in_TS'].astype(float).values
        RH = df['in_RH'].astype(float).values / 100
        P = df['in_P'].astype(float).values * 100
        SPD = df['in_SPD'].astype(float).values

        # Bounds check
        T = np.clip(T, -60, 100)
        TS = np.clip(TS, -4, 100)
        RH = np.clip(RH, 0, 1)
        P = np.clip(P, 0, 120000)
        SPD = np.clip(SPD, 0, 100)

        # === 2. Rolling stats (ONCE) ===
        window = window_size
        T_mean = pd.Series(T).rolling(
            window, center=True, min_periods=3).mean().values
        T_std = pd.Series(T).rolling(window, center=True,
                                     min_periods=3).std(ddof=0).values
        TS_mean = pd.Series(TS).rolling(
            window, center=True, min_periods=3).mean().values
        TS_std = pd.Series(TS).rolling(window, center=True,
                                       min_periods=3).std(ddof=0).values
        RH_mean = pd.Series(RH).rolling(
            window, center=True, min_periods=3).mean().values
        RH_std = pd.Series(RH).rolling(window, center=True,
                                       min_periods=3).std(ddof=0).values
        P_mean = pd.Series(P).rolling(
            window, center=True, min_periods=3).mean().values
        P_std = pd.Series(P).rolling(window, center=True,
                                     min_periods=3).std(ddof=0).values
        SPD_mean = pd.Series(SPD).rolling(
            window, center=True, min_periods=3).mean().values
        SPD_std = pd.Series(SPD).rolling(
            window, center=True, min_periods=3).std(ddof=0).values

        # === 3. Trim edges ===
        edge = window // 2
        valid = slice(edge, len(df) - edge)
        idx = np.arange(len(df))[valid]

        # Extract valid
        time = df['time'].values[valid]
        lat = df['latitude'].values[valid]
        lon = df['longitude'].values[valid]
        ship_id = df['platform_call_sign'].values[valid]

        T_mean, T_std = T_mean[valid], T_std[valid]
        TS_mean, TS_std = TS_mean[valid], TS_std[valid]
        RH_mean, RH_std = RH_mean[valid], RH_std[valid]
        P_mean, P_std = P_mean[valid], P_std[valid]
        SPD_mean, SPD_std = SPD_mean[valid], SPD_std[valid]

        # === 4. Chunk for parallel Monte Carlo ===
        n_points = len(time)
        chunks = []
        for i in range(0, n_points, CHUNK_SIZE):
            sl = slice(i, i + CHUNK_SIZE)
            chunk = {
                'idx': idx[sl],
                'time': time[sl],
                'lat': lat[sl],
                'lon': lon[sl],
                'ship': ship_id[sl],
                'T_mean': T_mean[sl], 'T_std': T_std[sl],
                'TS_mean': TS_mean[sl], 'TS_std': TS_std[sl],
                'RH_mean': RH_mean[sl], 'RH_std': RH_std[sl],
                'P_mean': P_mean[sl], 'P_std': P_std[sl],
                'SPD_mean': SPD_mean[sl], 'SPD_std': SPD_std[sl],
            }
            chunks.append(chunk)

        # === 5. Parallel Monte Carlo ===
        with Pool(N_CORES) as pool:
            results = pool.map(monte_carlo_chunk, chunks)

        # === 6. Combine & write ===
        with open(outpath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            for chunk_result in results:
                writer.writerows(chunk_result)

        logger.info(f"Done {filename} → {len(time)} rows")


def monte_carlo_chunk(chunk):
    n = len(chunk['time'])
    rows = []

    # Pre-sample all Monte Carlo inputs
    SPD_s = np.random.normal(
        chunk['SPD_mean'][:, None], chunk['SPD_std'][:, None], (n, N_SAMPLES))
    P_s = np.random.normal(chunk['P_mean'][:, None],
                           chunk['P_std'][:, None],   (n, N_SAMPLES))
    RH_s = np.random.normal(
        chunk['RH_mean'][:, None],  chunk['RH_std'][:, None],  (n, N_SAMPLES))
    T_s = np.random.normal(chunk['T_mean'][:, None],
                           chunk['T_std'][:, None],   (n, N_SAMPLES))
    TS_s = np.random.normal(
        chunk['TS_mean'][:, None],  chunk['TS_std'][:, None],  (n, N_SAMPLES))

    for i in range(n):
        shf, lhf, tau, mo = [], [], [], []

        for j in range(N_SAMPLES):
            if not all(np.isfinite(x) for x in [SPD_s[i, j], P_s[i, j], RH_s[i, j], T_s[i, j], TS_s[i, j]]):
                continue
            try:
                flux = mft_fluxes(
                    dyn_in_prm, SPD_s[i,
                                      j], dyn_in_val2, sfc_current1, sfc_current2,
                    convect, P_s[i, j], air_moist_prm, RH_s[i,
                                                            j], sfc_moist_prm,
                    sfc_moist_val, salinity, ss_prm, ss_val, T_s[i,
                                                                 j], sst_prm,
                    TS_s[i, j], ref_ht_wind, ref_ht_tq, z_wanted, astab, eqv_neut,
                    net_heat_flux, warn, flux_model, z0_mom_prm, z0_theta_q_prm, stable_prm,
                    oil_fract_area, dimensionless_m_o_length, zo_m, missing
                )
                shf.append(flux[1])
                lhf.append(flux[2])
                tau.append(flux[3][0])
                mo.append(flux[7])
            except:
                pass

        if shf:
            rows.append([
                chunk['time'][i], chunk['ship'][i], chunk['lat'][i], chunk['lon'][i],
                np.std(shf), np.mean(shf),
                np.std(lhf), np.mean(lhf),
                np.std(tau), np.mean(tau),
                np.mean(mo),
                chunk['T_std'][i], chunk['TS_std'][i],
                chunk['RH_std'][i], chunk['P_std'][i], chunk['SPD_std'][i]
            ])
        else:
            rows.append([
                chunk['time'][i], chunk['ship'][i], chunk['lat'][i], chunk['lon'][i],
                -9999, -9999, -9999, -9999, -9999, -9999, -9999,
                -9999, -9999, -9999, -9999, -9999
            ])
    return rows
