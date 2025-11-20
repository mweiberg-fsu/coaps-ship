import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# === CONFIGURATION ===
N = 500  # Window size used to compute stdv in processing

# --- Load flux runs ---
print("Loading flux runs...")
r1 = pd.read_csv("KAOU2025_flux_r1.csv")
r2 = pd.read_csv("KAOU2025_flux_r2.csv")

# Merge on timestamp
df = r1.merge(r2, on='time', suffixes=('_r1', '_r2'))
print(f"Merged rows: {len(df)}")

# Ensure 'time' is datetime
df['time'] = pd.to_datetime(df['time'])

# -------------------------------------------------
# --- Load RAW INPUT data (KAOU2025.csv) ---
# -------------------------------------------------
print("\nLoading KAOU2025.csv ...")

try:
    raw_df = pd.read_csv(
        "KAOU2025.csv",
        skiprows=[1],           # Skip units row (UTC, degrees_north, ...)
        engine="python",        # Handle malformed lines
        on_bad_lines="skip",    # Skip bad lines
        dtype=str               # Avoid mixed-type warnings
    )
    print(f"  → Raw CSV loaded: {len(raw_df)} rows (pre-cleaning)")

    # Parse ISO-8601 with 'Z' (UTC)
    raw_df['time'] = pd.to_datetime(
        raw_df['time'],
        utc=True,
        errors='coerce',
        format='%Y-%m-%dT%H:%M:%SZ'
    )

    invalid_times = raw_df['time'].isna().sum()
    raw_df = raw_df.dropna(subset=['time']).reset_index(drop=True)
    print(
        f"  → Valid timestamps: {len(raw_df)} (dropped {invalid_times} invalid)")

except Exception as e:
    print(f"ERROR loading KAOU2025.csv: {e}")
    raise

# === VARIABLES TO COMPARE ===
variables = ['hfss stdv', 'hfls stdv']

# Create output directory
os.makedirs("plots", exist_ok=True)

# === MAIN LOOP ===
for var in variables:
    print(f"\n{'='*70}")
    print(f"ANALYZING: {var.upper()}")
    print(f"{'='*70}")

    a = df[f'{var}_r1'].to_numpy()
    b = df[f'{var}_r2'].to_numpy()

    # --- REMOVE NaN PAIRS ---
    mask = ~np.isnan(a) & ~np.isnan(b)
    a_clean = a[mask]
    b_clean = b[mask]
    times_clean = df['time'].iloc[np.where(mask)[0]].to_numpy()

    print(f"  Valid pairs: {len(a_clean):,} / {len(df):,}")
    if len(a_clean) == 0:
        print("  No valid data → skipping")
        continue

    # --- BASIC STATISTICS ---
    mean_a = np.mean(a_clean)
    mean_b = np.mean(b_clean)
    mean_diff = mean_b - mean_a

    print(f"  Min Run1: {a_clean.min():.3f}")
    print(f"  Max Run1: {a_clean.max():.3f}")
    print(f"  Min Run2: {b_clean.min():.3f}")
    print(f"  Max Run2: {b_clean.max():.3f}")
    print(f"  Mean Run1: {mean_a:.3f}")
    print(f"  Mean Run2: {mean_b:.3f}")
    print(f"  Mean Diff (R2 - R1): {mean_diff:+.3f}")

    # --- CORRELATION & REGRESSION ---
    if len(a_clean) > 1 and np.var(a_clean) > 1e-12 and np.var(b_clean) > 1e-12:
        corr, _ = stats.pearsonr(a_clean, b_clean)
        slope, intercept, r, p, se = stats.linregress(a_clean, b_clean)
    else:
        corr = slope = intercept = r = np.nan

    print(f"  Correlation: {corr:.3f}" if not np.isnan(
        corr) else "  Correlation: nan")
    print(f"  Slope: {slope:.3f}, Intercept: {intercept:.3f}" if not np.isnan(
        slope) else "  Slope: nan")
    print(f"  R²: {r**2:.3f}" if not np.isnan(r) else "  R²: nan")

    # --- RELATIVE DIFFERENCE ---
    denom = a_clean + b_clean + 1e-8
    rel_diff = 200 * (b_clean - a_clean) / denom
    rel_diff = np.nan_to_num(rel_diff, nan=0.0)

    frac_50 = np.mean(np.abs(rel_diff) > 50)
    frac_20 = np.mean(np.abs(rel_diff) > 20)
    frac_10 = np.mean(np.abs(rel_diff) > 10)

    print(f"  Fraction >50% rel diff: {frac_50:.3f} ({frac_50*100:.1f}%)")
    print(f"  Fraction >20% rel diff: {frac_20:.3f} ({frac_20*100:.1f}%)")
    print(f"  Fraction >10% rel diff: {frac_10:.3f} ({frac_10*100:.1f}%)")

    # --- OUTLIERS (>50% rel diff) ---
    outlier_mask = np.abs(rel_diff) > 50
    outlier_times = times_clean[outlier_mask]
    outlier_a = a_clean[outlier_mask]
    outlier_b = b_clean[outlier_mask]
    print(f"  Outliers (>50% diff): {len(outlier_times)} points")

    # --- PLOT 1: Scatter ---
    plt.figure(figsize=(7, 6))
    plt.scatter(a_clean, b_clean, s=4, alpha=0.5, label='All points')
    if len(outlier_times):
        plt.scatter(outlier_a, outlier_b, s=8, color='red', label='>50% diff')

    max_val = max(a_clean.max(), b_clean.max()) * 1.05
    plt.plot([0, max_val], [0, max_val], 'r--', label='1:1 line')
    if not np.isnan(slope):
        x_fit = np.array([a_clean.min(), a_clean.max()])
        plt.plot(x_fit, slope * x_fit + intercept,
                 'k-', label=f'Fit (R²={r**2:.3f})')

    plt.xlabel(f'{var} — Run 1')
    plt.ylabel(f'{var} — Run 2')
    plt.title(f'{var} | Run Comparison (N={N})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"plots/{var.replace(' ', '_')}_scatter_N{N}.png", dpi=180)
    plt.close()

    # --- PLOT 2: Bland-Altman ---
    mean_ab = 0.5 * (a_clean + b_clean)
    diff_ab = b_clean - a_clean

    plt.figure(figsize=(7, 4.5))
    plt.scatter(mean_ab, diff_ab, s=4, alpha=0.5, label='All points')
    if len(outlier_times):
        plt.scatter(0.5*(outlier_a + outlier_b), outlier_b - outlier_a,
                    s=10, color='red', label='>50% diff')

    plt.axhline(0, color='k', linewidth=1)
    plt.axhline(mean_diff, color='blue', linestyle='--',
                label=f'Mean diff = {mean_diff:+.3f}')
    plt.axhline(mean_diff + 1.96*np.std(diff_ab), color='gray', linestyle=':',
                label=f'±1.96σ = {1.96*np.std(diff_ab):+.3f}')
    plt.axhline(mean_diff - 1.96*np.std(diff_ab), color='gray', linestyle=':')

    plt.xlabel('Mean of Run1 and Run2')
    plt.ylabel('Run2 − Run1')
    plt.title(f'Bland–Altman: {var}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(
        f"plots/{var.replace(' ', '_')}_bland_altman_N{N}.png", dpi=180)
    plt.close()

    # --- PLOT 3: Time Series ---
    plt.figure(figsize=(12, 5))
    plt.plot(times_clean, a_clean, label='Run 1', alpha=0.7, linewidth=0.8)
    plt.plot(times_clean, b_clean, label='Run 2', alpha=0.7, linewidth=0.8)
    if len(outlier_times):
        plt.scatter(outlier_times, outlier_a, color='red',
                    s=10, label='>50% diff (Run1)')
        plt.scatter(outlier_times, outlier_b, color='orange',
                    s=10, label='>50% diff (Run2)')
    plt.xlabel('Time')
    plt.ylabel(var)
    plt.title(f'Time Series: {var} | Outliers >50% Difference')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(
        f"plots/{var.replace(' ', '_')}_time_series_outliers_N{N}.png", dpi=180)
    plt.close()

    # --- SAVE OUTLIERS: FULL FLUX + RAW INPUT DATA (FIXED) ---
    if len(outlier_times) > 0:
        # 1. Full flux data for these times
        outlier_flux_df = df[df['time'].isin(outlier_times)].copy()

        # 2. Diagnostics for current variable
        diag_df = pd.DataFrame({
            'time': outlier_times,
            f'{var}_r1': outlier_a,
            f'{var}_r2': outlier_b,
            'diff': outlier_b - outlier_a,
            'rel_diff_%': rel_diff[outlier_mask]
        })

        # 3. Avoid merge conflict: drop current var from flux data
        flux_cols_to_keep = [
            c for c in outlier_flux_df.columns if not c.startswith(f'{var}_')]
        outlier_flux_df = outlier_flux_df[flux_cols_to_keep]

        # Merge
        outlier_full = diag_df.merge(outlier_flux_df, on='time', how='left')
        outlier_full = outlier_full.merge(raw_df, on='time', how='left')

        # 4. Reorder safely
        raw_cols = [c for c in raw_df.columns if c != 'time']
        r1_cols = [c for c in df.columns if c.endswith(
            '_r1') and c != 'time' and not c.startswith(f'{var}_')]
        r2_cols = [c for c in df.columns if c.endswith(
            '_r2') and c != 'time' and not c.startswith(f'{var}_')]

        order = (
            ['time'] +
            raw_cols +
            [f'{var}_r1', f'{var}_r2'] +
            r1_cols +
            r2_cols +
            ['diff', 'rel_diff_%']
        )
        order = [col for col in order if col in outlier_full.columns]
        outlier_full = outlier_full[order]

        # Save
        csv_path = f"plots/{var.replace(' ', '_')}_outliers_FULL_N{N}.csv"
        outlier_full.to_csv(csv_path, index=False)
        print(f"  FULL OUTLIER CSV saved → {csv_path}")
        print(
            f"    → {len(outlier_full):,} outliers, {len(outlier_full.columns)} columns")
        print(f"      - Raw inputs: {len(raw_cols)}")
        print(f"      - Run 1 flux: {len(r1_cols) + 1}")
        print(f"      - Run 2 flux: {len(r2_cols) + 1}")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("Plots and FULL outlier CSVs saved in ./plots/")
print("="*70)
