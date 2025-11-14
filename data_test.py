import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

N = 100

# --- load both runs ---
r1 = pd.read_csv(f"data/output/KAOU/KAOU2026_processed_{N}_r1.csv")
r2 = pd.read_csv(f"data/output/KAOU/KAOU2026_processed_{N}_r2.csv")

# make sure timestamps align
df = r1.merge(r2, on='time', suffixes=('_r1', '_r2'))
print("Merged rows:", len(df))

# --- variables to compare ---
for var in ['hfss stdv', 'hfls stdv']:
    a = df[f'{var}_r1'].to_numpy()
    b = df[f'{var}_r2'].to_numpy()
    diff = b - a
    mean_a = np.nanmean(a)
    mean_diff = np.nanmean(diff)
    corr = np.corrcoef(a, b)[0,1]
    slope, intercept, r, p, se = stats.linregress(a, b)
    print(f'\n{var}')
    print(f'  mean run1: {mean_a:.2f}')
    print(f'  mean run2: {np.nanmean(b):.2f}')
    print(f'  mean diff (r2 - r1): {mean_diff:.2f}')
    print(f'  correlation: {corr:.3f}')
    print(f'  slope: {slope:.3f}, intercept: {intercept:.3f}')
    print(f'  r^2: {r**2:.3f}')

    expected_se = (0.5*(a+b)) / np.sqrt(2*N)
    z = diff / expected_se
    frac_large = np.mean(np.abs(z) > 3)
    print(f'  fraction |z|>3: {frac_large:.3f}')

    # --- visuals ---
    plt.figure(figsize=(6,6))
    plt.scatter(a, b, s=5, alpha=0.3)
    plt.plot([0, max(a.max(), b.max())],
             [0, max(a.max(), b.max())],
             'r--', label='1:1')
    plt.plot(a, slope*a + intercept, 'k-', label=f'fit r²={r**2:.2f}')
    plt.xlabel(f'{var} Run 1')
    plt.ylabel(f'{var} Run 2')
    plt.title(f'{var} comparison (N={N})')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'plots/{var.replace(" ","_")}_run_comparison_{N}.png', dpi=150)
    # plt.show()

    # Bland–Altman
    m = 0.5*(a+b)
    plt.figure(figsize=(6,4))
    plt.scatter(m, diff, s=5, alpha=0.3)
    plt.axhline(0, color='k')
    plt.xlabel('Mean of runs')
    plt.ylabel('Run2 - Run1')
    plt.title(f'Difference between runs against mean: {var}')
    plt.savefig(f'plots/{var.replace(" ","_")}_diff_{N}.png', dpi=150)
    # plt.show()
