# debug_one_file_parallel_80.py
import pandas as pd
import numpy as np
from MFT23 import mft_fluxes
import os
from config import *
from constants import *
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

# ------------------------------------------------------------------
#  USER SETTINGS
# ------------------------------------------------------------------
INPUT_FILE   = "data/input/KAOU/KAOU2011.csv"
DEBUG_CSV    = "DEBUG_SHF_ONLY.csv"
N_SAMPLES    = 20
NP_SEED      = 42
MAX_WORKERS  = 80
# ------------------------------------------------------------------

def compute_shf_for_row(args):
    """One row → Monte-Carlo → (mean, std, n_valid)"""
    i, t, ship, lat, lon, T, Ts, U, RH, Pr = args
    np.random.seed(NP_SEED + i)                 # reproducible per row

    if not all(np.isfinite(x) for x in [T, Ts, U, RH, Pr]):
        return [t, ship, lat, lon, T, Ts, U, RH, -9999, -9999, 0]

    shf_vals = []
    for _ in range(N_SAMPLES):
        T_g  = T  + np.random.normal(0, 0.1)
        Ts_g = Ts + np.random.normal(0, 0.1)
        U_g  = max(0.1, U + np.random.normal(0, 0.3))
        RH_g = np.clip(RH + np.random.normal(0, 0.01), 0.01, 0.99)
        Pr_g = Pr + np.random.normal(0, 50)

        try:
            flux = mft_fluxes(
                dyn_in_prm, U_g, dyn_in_val2, sfc_current1, sfc_current2,
                convect, Pr_g, air_moist_prm, RH_g, sfc_moist_prm,
                sfc_moist_val, salinity, ss_prm, ss_val, T_g, sst_prm,
                Ts_g, ref_ht_wind, ref_ht_tq, z_wanted, astab, eqv_neut,
                net_heat_flux, warn, flux_model, z0_mom_prm, z0_theta_q_prm,
                stable_prm, oil_fract_area, dimensionless_m_o_length, zo_m, missing
            )
            shf = flux[1]                     # SHF
            if np.isfinite(shf):
                shf_vals.append(shf)
        except Exception:
            continue

    if shf_vals:
        return [t, ship, lat, lon, T, Ts, U, RH,
                np.mean(shf_vals), np.std(shf_vals), len(shf_vals)]
    else:
        return [t, ship, lat, lon, T, Ts, U, RH, -9999, -9999, 0]


# ------------------------------------------------------------------
#  MAIN – **MUST** be inside `if __name__ == '__main__':`
# ------------------------------------------------------------------
if __name__ == '__main__':                     # <<< THIS LINE IS CRUCIAL
    assert os.path.exists(INPUT_FILE), "File not found!"

    print(f"Loading {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE, skiprows=[1, 2])
    print(f"Loaded {len(df)} rows")

    # ------------------------------------------------------------------
    #  Extract + clean columns
    # ------------------------------------------------------------------
    airT = df['in_T'].astype(float).values
    rh   = df['in_RH'].astype(float).values / 100
    P    = df['in_P'].astype(float).values * 100
    sst  = df['in_TS'].astype(float).values
    wind = df['in_SPD'].astype(float).values
    time = df['time'].values
    lat  = df['latitude'].values
    lon  = df['longitude'].values
    ship = df['platform_call_sign'].values

    airT = np.where((airT < -60) | (airT > 100), np.nan, airT)
    sst  = np.where((sst < -4)  | (sst > 100),  np.nan, sst)
    rh   = np.where((rh < 0)    | (rh > 1),     np.nan, rh)
    P    = np.where((P < 0)     | (P > 120000), np.nan, P)
    wind = np.where((wind < 0)  | (wind > 100), np.nan, wind)

    # ------------------------------------------------------------------
    #  Build list of inputs for the pool
    # ------------------------------------------------------------------
    inputs = [
        (i, time[i], ship[i], lat[i], lon[i],
         airT[i], sst[i], wind[i], rh[i], P[i])
        for i in range(len(df))
    ]

    # ------------------------------------------------------------------
    #  Write CSV header once
    # ------------------------------------------------------------------
    with open(DEBUG_CSV, 'w') as f:
        f.write("time,ship,lat,lon,airT,sst,wind,rh,SHF_mean,SHF_std,n_valid\n")

    # ------------------------------------------------------------------
    #  Parallel execution – 80 workers
    # ------------------------------------------------------------------
    print(f"\nComputing SHF with {MAX_WORKERS} workers ...")
    results = []

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_idx = {executor.submit(compute_shf_for_row, inp): inp[0] for inp in inputs}

        for future in as_completed(future_to_idx):
            idx    = future_to_idx[future]
            row    = future.result()
            results.append(row)

            # append line immediately (atomic on most filesystems)
            t, s, la, lo, T, Ts, U, RH, mean_shf, std_shf, n = row
            with open(DEBUG_CSV, 'a') as f:
                f.write(f"{t},{s},{la:.4f},{lo:.4f},{T:.2f},{Ts:.2f},{U:.2f},{RH:.3f},"
                        f"{mean_shf:.2f},{std_shf:.2f},{n}\n")

            # progress for the first 10, every 200th, and last 5 rows
            if idx < 10 or idx % 200 == 0 or idx >= len(df) - 5:
                print(f"Row {idx:5d} → SHF = {mean_shf:+6.2f} ± {std_shf:.2f} W/m²")

    # ------------------------------------------------------------------
    #  Final summary
    # ------------------------------------------------------------------
    df_out = pd.DataFrame(results,
        columns="time,ship,lat,lon,airT,sst,wind,rh,SHF_mean,SHF_std,n_valid".split(','))
    valid = df_out['SHF_mean'] > -999
    print("\n=== DONE ===")
    print(f"Valid rows : {valid.sum()} / {len(df)}")
    if valid.any():
        print(f"SHF range  : {df_out.loc[valid,'SHF_mean'].min():.1f} … "
              f"{df_out.loc[valid,'SHF_mean'].max():.1f} W/m²")
        print(f"Mean SHF   : {df_out.loc[valid,'SHF_mean'].mean():.2f} W/m²")
    print(f"\nOpen {DEBUG_CSV} in Excel → plot SHF vs (airT-sst)")