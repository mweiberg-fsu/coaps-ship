import polars as pl
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pytz

# Step 0: Load and Merge Data
file1 = 'data/output/KAOU/KAOU2025_processed_v1.csv'
file2 = 'data/output/KAOU/KAOU2025_processed_v2.csv'
title_text = "(KAOU 2025 HFSS stdv)"

df1 = pl.read_csv(file1, columns=['time', 'hfss stdv', 'latitude', 'SPD stdv'], 
                  schema_overrides={'hfss stdv': pl.Float32, 'latitude': pl.Float32, 'SPD stdv': pl.Float32}, 
                  try_parse_dates=True)
df2 = pl.read_csv(file2, columns=['time', 'hfss stdv', 'latitude', 'SPD stdv'], 
                  schema_overrides={'hfss stdv': pl.Float32, 'latitude': pl.Float32, 'SPD stdv': pl.Float32}, 
                  try_parse_dates=True)

merged = df1.rename({'hfss stdv': 'hfss_stdv_file1', 'latitude': 'latitude_file1', 'SPD stdv': 'SPD_stdv_file1'}).join(
    df2.rename({'hfss stdv': 'hfss_stdv_file2', 'latitude': 'latitude_file2', 'SPD stdv': 'SPD_stdv_file2'}), 
    on='time', how='inner'
).drop_nulls()

# Convert to pandas
merged_pd = merged.to_pandas()
merged_pd.set_index('time', inplace=True)

# Step 1: Validate Data
print("Data Validation - Descriptive Stats:")
print(merged_pd.describe())

# Optional: Normalize if scale seems off
if merged_pd['hfss_stdv_file1'].max() > 10:
    print("Warning: High values detected. Normalizing by dividing by 100.")
    merged_pd['hfss_stdv_file1'] /= 100
    merged_pd['hfss_stdv_file2'] /= 100

print("\nNumber of duplicate timestamps:", merged_pd.index.duplicated().sum())

# Step 2: Investigate Spikes
spike_threshold = 200  # Adjust based on data (e.g., 0.5 if normalized)
spikes = merged_pd[merged_pd['hfss_stdv_file1'] > spike_threshold]
print("\nSpike Periods (hfss_stdv_file1 > {}):".format(spike_threshold))
print(spikes[['hfss_stdv_file1', 'hfss_stdv_file2', 'latitude_file1', 'SPD_stdv_file1']])

if not spikes.empty:
    corr_spd = spikes['hfss_stdv_file1'].corr(spikes['SPD_stdv_file1'])
    print("Correlation between hfss stdv and SPD stdv during spikes:", corr_spd)

# 2b: Zoom-in plot for a specific period (handle timezone-aware index)
zoom_start = pd.Timestamp('2025-05-10', tz=pytz.UTC)  # Set to UTC
zoom_end = pd.Timestamp('2025-05-20', tz=pytz.UTC)    # Set to UTC
merged_pd_sorted = merged_pd.sort_index()
zoomed = merged_pd_sorted.loc[zoom_start:zoom_end]

plt.figure(figsize=(12, 6))
plt.plot(zoomed.index, zoomed['hfss_stdv_file1'], label='File 1 hfss stdv')
plt.plot(zoomed.index, zoomed['hfss_stdv_file2'], label='File 2 hfss stdv')
plt.xlabel('Time')
plt.ylabel('hfss stdv')
plt.title('Zoomed Comparison ({} to {})'.format(zoom_start, zoom_end))
plt.legend()
plt.xticks(rotation=45)
plt.grid(True)
plt.tight_layout()
plt.show()

# Step 3: Quantify Differences
spike_date_start = pd.Timestamp('2025-05-15', tz=pytz.UTC)
spike_date_end = pd.Timestamp('2025-05-15', tz=pytz.UTC)
spike_diff = merged_pd_sorted.loc[spike_date_start:spike_date_end]
if not spike_diff.empty:
    mean_diff = (spike_diff['hfss_stdv_file1'] - spike_diff['hfss_stdv_file2']).mean()
    print("\nMean difference during spike ({} to {}): {}".format(spike_date_start, spike_date_end, mean_diff))

overall_corr = merged_pd['hfss_stdv_file1'].corr(merged_pd['hfss_stdv_file2'])
print("Overall correlation between File 1 and File 2 hfss stdv:", overall_corr)

non_spikes = merged_pd[merged_pd['hfss_stdv_file1'] <= spike_threshold]
if not non_spikes.empty:
    corr_non_spikes = non_spikes['hfss_stdv_file1'].corr(non_spikes['hfss_stdv_file2'])
    print("Correlation during non-spike periods:", corr_non_spikes)
if not spikes.empty:
    corr_spikes = spikes['hfss_stdv_file1'].corr(spikes['hfss_stdv_file2'])
    print("Correlation during spike periods:", corr_spikes)

# Step 4: Refine Plot
daily = merged_pd.resample('D').mean()

plt.figure(figsize=(12, 6))
plt.plot(daily.index, daily['hfss_stdv_file1'], label='File 1 hfss stdv (Daily Mean)')
plt.plot(daily.index, daily['hfss_stdv_file2'], label='File 2 hfss stdv (Daily Mean)')
plt.xlabel('Time')
plt.ylabel('hfss stdv')
plt.title(f'Daily Downsampled Comparison {title_text}')
plt.legend()
plt.xticks(rotation=45)
plt.grid(True)
plt.tight_layout()
plt.show()

merged_pd['diff'] = merged_pd['hfss_stdv_file1'] - merged_pd['hfss_stdv_file2']

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(merged_pd.index, merged_pd['hfss_stdv_file1'], 'b-', label='File 1 hfss stdv')
ax1.plot(merged_pd.index, merged_pd['hfss_stdv_file2'], 'orange', label='File 2 hfss stdv')
ax1.set_xlabel('Time')
ax1.set_ylabel('hfss stdv', color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.legend(loc='upper left')

ax2 = ax1.twinx()
ax2.plot(merged_pd.index, merged_pd['diff'], 'g--', label='Difference (File1 - File2)')
ax2.set_ylabel('Difference', color='green')
ax2.tick_params(axis='y', labelcolor='green')
ax2.axhline(0, color='red', linestyle='--')
ax2.legend(loc='upper right')

plt.title('Refined Comparison with Difference Overlay')
plt.xticks(rotation=45)
plt.grid(True)
plt.tight_layout()
plt.show()