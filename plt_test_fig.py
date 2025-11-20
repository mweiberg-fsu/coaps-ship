# ------------------------------------------------------------
# plot_outliers_enhanced.py
# Enhanced diagnostic plots for outlier CSV (with MANY env vars)
# ------------------------------------------------------------
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import cm

# ------------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------------
CSV_FILE = "plots/hfls_stdv_outliers_FULL_N500.csv"   # Change if needed
OUT_DIR = "outlier_plots_enhanced"
os.makedirs(OUT_DIR, exist_ok=True)

# ------------------------------------------------------------------
# 2. LOAD DATA
# ------------------------------------------------------------------
print(f"Loading {CSV_FILE} ...")
df = pd.read_csv(CSV_FILE)
df["time"] = pd.to_datetime(df["time"], utc=True)
print(f"  → {len(df)} outliers loaded")

# ------------------------------------------------------------------
# 3. PLOT 1 – Scatter Run1 vs Run2 (hfls stdv)
# ------------------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(8, 7))
sc = ax1.scatter(df["hfls stdv_r1"], df["hfls stdv_r2"],
                 c=df["rel_diff_%"], cmap="RdBu_r", edgecolor="k", s=70, alpha=0.9)
max_val = max(df["hfls stdv_r1"].max(), df["hfls stdv_r2"].max()) * 1.1
ax1.plot([0, max_val], [0, max_val], "r--", lw=2, label="1:1 line")
ax1.set_xlabel("hfls stdv – Run 1")
ax1.set_ylabel("hfls stdv – Run 2")
ax1.set_title("Outliers: hfls stdv Run 1 vs Run 2")
ax1.legend()
ax1.grid(True, alpha=0.3)
cbar = plt.colorbar(sc, ax=ax1)
cbar.set_label("rel_diff_%")
plt.tight_layout()
fig1.savefig(os.path.join(
    OUT_DIR, "1_scatter_colored_by_reldiff.png"), dpi=220)
print("→ Plot 1 saved")

# ------------------------------------------------------------------
# 4. PLOT 2 – Time series of disagreement
# ------------------------------------------------------------------
fig2, (ax2a, ax2b) = plt.subplots(2, 1, figsize=(13, 7), sharex=True)
ax2a.plot(df["time"], df["diff"], marker="o", ls="-",
          color="#2ca02c", label="diff (R2−R1)")
ax2a.set_ylabel("diff")
ax2a.grid(True, alpha=0.3)
ax2a.legend()

ax2b.plot(df["time"], df["rel_diff_%"], marker="^",
          ls="--", color="#d62728", label="rel_diff_%")
ax2b.set_ylabel("rel_diff_%")
ax2b.set_xlabel("Time (UTC)")
ax2b.grid(True, alpha=0.3)
ax2b.legend()
ax2b.tick_params(axis='x', rotation=30)

plt.suptitle("Temporal Evolution of hfls stdv Disagreement")
plt.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, "2_time_series_disagreement.png"), dpi=220)
print("→ Plot 2 saved")

# ------------------------------------------------------------------
# 5. PLOT 3 – ENVIRONMENTAL DRIVERS vs RELATIVE DIFFERENCE (8 PANELS!)
# ------------------------------------------------------------------
# List of variables to plot
env_vars = [
    ("in_SPD", "Wind Speed (m/s)", "viridis"),
    ("in_T", "Air Temp (°C)", "plasma"),
    ("in_TS", "Sea Surface Temp (°C)", "inferno"),
    ("in_RH", "Relative Humidity (%)", "cividis"),
    ("in_P", "Pressure (mb)", "magma"),
    ("hfls mean_r1", "Latent Heat Flux Mean R1 (W/m²)", "coolwarm"),
    ("hfss mean_r1", "Sensible Heat Flux Mean R1 (W/m²)", "turbo"),
    ("tau mean_r1", "Wind Stress Mean R1 (N/m²)", "hot"),
]

fig3, axes = plt.subplots(4, 2, figsize=(16, 20))
axes = axes.flatten()

for idx, (var, label, cmap_name) in enumerate(env_vars):
    if var not in df.columns:
        axes[idx].text(0.5, 0.5, f"{var}\nnot found", ha='center',
                       va='center', transform=axes[idx].transAxes)
        axes[idx].set_title(label)
        continue

    sc = axes[idx].scatter(df[var], df["rel_diff_%"],
                           c=df["diff"], cmap=cmap_name, edgecolor="k", s=60, alpha=0.85)
    axes[idx].set_xlabel(label)
    axes[idx].set_ylabel("rel_diff_%")
    axes[idx].set_title(f"{label} vs Relative Difference")
    axes[idx].grid(True, alpha=0.3)
    cbar = plt.colorbar(sc, ax=axes[idx], shrink=0.8)
    cbar.set_label("diff (R2−R1)")

# Hide unused subplots
for idx in range(len(env_vars), len(axes)):
    fig3.delaxes(axes[idx])

plt.suptitle(
    "Environmental Drivers vs hfls stdv Relative Difference", fontsize=16, y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])
fig3.savefig(os.path.join(
    OUT_DIR, "3_env_drivers_vs_reldiff_8panel.png"), dpi=250)
print("→ Plot 3 (8-panel) saved")

# ------------------------------------------------------------------
# BONUS: One giant correlation heatmap (optional)
# ------------------------------------------------------------------
corr_vars = ["in_T", "in_RH", "in_TS", "in_SPD", "in_P",
             "hfls stdv_r1", "hfls stdv_r2", "hfss stdv_r1", "hfss stdv_r2",
             "hfls mean_r1", "hfss mean_r1", "tau mean_r1", "rel_diff_%"]

available = [v for v in corr_vars if v in df.columns]
corr_matrix = df[available].corr()

fig4, ax4 = plt.subplots(figsize=(10, 8))
im = ax4.imshow(corr_matrix, cmap="coolwarm", vmin=-1, vmax=1)
ax4.set_xticks(np.arange(len(available)))
ax4.set_yticks(np.arange(len(available)))
ax4.set_xticklabels(available, rotation=45, ha="right")
ax4.set_yticklabels(available)
plt.colorbar(im, ax=ax4, label="Correlation")
for i in range(len(available)):
    for j in range(len(available)):
        ax4.text(j, i, f"{corr_matrix.iloc[i, j]:.2f}",
                 ha="center", va="center", color="black" if abs(corr_matrix.iloc[i, j]) < 0.5 else "white")

ax4.set_title("Correlation Matrix of Key Variables")
plt.tight_layout()
fig4.savefig(os.path.join(OUT_DIR, "4_correlation_heatmap.png"), dpi=220)
print("→ Bonus: Correlation heatmap saved")

# ------------------------------------------------------------------
print("\nALL DONE!")
print(f"Plots saved in: {OUT_DIR}")
print("""
Files created:
├─ 1_scatter_colored_by_reldiff.png
├─ 2_time_series_disagreement.png
├─ 3_env_drivers_vs_reldiff_8panel.png   ← THE BIG ONE!
└─ 4_correlation_heatmap.png
""")
