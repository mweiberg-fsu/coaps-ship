#!/usr/bin/env python
import polars as pl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datashader as ds
import datashader.transfer_functions as tf
import colorcet
from scipy.stats import linregress
import re

# ————————————————————————  CONFIG  ————————————————————————
FILE_1 = "data/output/KAOU/KAOU2025_processed_v1.csv"
FILE_2 = "data/output/KAOU/KAOU2025_processed_v2.csv"

COLS = ["time", "hfss stdv", "hfls stdv"]

# sensible heat filter (your request)
HFSS_MIN, HFSS_MAX = -50.0, 50.0

# optional latent heat filter (uncomment if you want)
# HFLS_MIN, HFLS_MAX = 0.0, 300.0
# ——————————————————————————————————————————————————————————

match = re.search(r'/([A-Za-z]+)(\d{4})_', FILE_1)
if match:
    ship_name = match.group(1)
    year = match.group(2)
    print("ship_name:", ship_name)
    print("year:", year)
else:
    raise ValueError("Could not parse ship/year from filename")

def load(f):
    return (pl.read_csv(f, columns=COLS, dtypes={"time": pl.Utf8})
            .with_columns(pl.col("time").str.to_datetime("%Y-%m-%dT%H:%M:%SZ", time_zone="UTC"))
            .set_sorted("time"))

df1 = load(FILE_1).rename({"hfss stdv": "hfss_1", "hfls stdv": "hfls_1"})
df2 = load(FILE_2).rename({"hfss stdv": "hfss_2", "hfls stdv": "hfls_2"})

print(f"df1 rows: {len(df1):,}")
print(f"df2 rows: {len(df2):,}")

# Unique timestamps
t1 = set(df1["time"].dt.strftime("%Y-%m-%d %H:%M:%S"))
t2 = set(df2["time"].dt.strftime("%Y-%m-%d %H:%M:%S"))
print(f"Unique timestamps in v1: {len(t1):,}")
print(f"Unique timestamps in v2: {len(t2):,}")
print(f"Overlap: {len(t1 & t2):,}")

for name, df in [("Run-1", df1), ("Run-2", df2)]:
    n_mc = df.height // df["time"].n_unique()
    print(f"{name}:  {df.height:,} rows  →  {n_mc:,} Monte-Carlo samples per timestamp")

# ————————  JOIN + FILTER OUTLIERS  ————————
df_full = df1.join(df2, on="time", how="inner")
print(f"Matched (pre-filter): {len(df_full):,} rows")

df = (df_full
      .filter(pl.col("hfss_1").is_between(HFSS_MIN, HFSS_MAX))
      .filter(pl.col("hfss_2").is_between(HFSS_MIN, HFSS_MAX))
      # .filter(pl.col("hfls_1").is_between(HFLS_MIN, HFLS_MAX))   # <-- uncomment if needed
      # .filter(pl.col("hfls_2").is_between(HFLS_MIN, HFLS_MAX))
)

removed = len(df_full) - len(df)
print(f"Matched & filtered: {len(df):,} rows  (-{removed:,} outliers removed)")

# ————————  DATASHADER SCATTER  ————————
def ds_scatter(x, y):
    cvs = ds.Canvas(plot_width=1000, plot_height=1000)
    agg = cvs.points(pd.DataFrame({"x": x, "y": y}), "x", "y")
    img = tf.shade(agg, cmap=colorcet.fire, how="log")
    img = tf.dynspread(img, threshold=0.5, max_px=6)
    return img

vars = [
    ("hfls stdv", "hfls_1", "hfls_2", "Latent Heat Std-Dev"),
    ("hfss stdv", "hfss_1", "hfss_2", "Sensible Heat Std-Dev"),
]

fig, axs = plt.subplots(1, 2, figsize=(16, 8), constrained_layout=True)

for ax, (name, c1, c2, label) in zip(axs, vars):
    x = df[c1].to_numpy()
    y = df[c2].to_numpy()
    img = ds_scatter(x, y)
    img = tf.spread(img, px=2)
    lims = [min(x.min(), y.min()), max(x.max(), y.max())]

    ax.imshow(np.flipud(img.to_pil()),
              extent=[lims[0], lims[1], lims[0], lims[1]],
              origin="lower", aspect="equal")

    ax.plot(lims, lims, "r--", lw=2, label="1:1")
    mask = ~np.isnan(x) & ~np.isnan(y)
    slope, intercept, r, _, _ = linregress(x[mask], y[mask])
    ax.plot(x[mask], slope * x[mask] + intercept, "k-", lw=1.5,
            label=f"fit (r²={r**2:.4f})")

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel(f"{label} – Run 1")
    ax.set_ylabel(f"{label} – Run 2")
    ax.set_title(name)
    ax.legend()
    ax.grid(True, ls=":", alpha=0.5)

# ————————  TITLE & SAVE  ————————
plt.suptitle(f"Uncertainty Calculations for {ship_name} "
             f"({year}, {len(df):,} matched & filtered obs)", 
             fontsize=16, y=1.04)

out_png = f"flux_comparison_{ship_name}{year}_filtered.png"
plt.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"Saved: {out_png}")
# plt.show()