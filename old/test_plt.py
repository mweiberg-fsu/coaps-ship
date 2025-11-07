import polars as pl
import matplotlib.pyplot as plt
import numpy as np

# Load only needed columns selectively (fast, low memory)
df1 = pl.read_csv('data/output/KAOU/KAOU2025_processed_v1.csv', columns=['time', 'hfss stdv'], dtypes={'hfss stdv': pl.Float32}, try_parse_dates=True)
df2 = pl.read_csv('data/output/KAOU/KAOU2025_processed_v2.csv', columns=['time', 'hfss stdv'], dtypes={'hfss stdv': pl.Float32}, try_parse_dates=True)

# Rename for clarity and join on 'time' (inner join for common timestamps only; use 'outer' if needed)
merged = df1.rename({'hfss stdv': 'hfss_stdv_file1'}).join(
    df2.rename({'hfss stdv': 'hfss_stdv_file2'}), on='time', how='inner'
)

# Drop any NA (rare, but handles mismatches)
merged = merged.drop_nulls()

# Convert to pandas only for plotting (Polars plotting is limited; this step is quick on filtered data)
merged_pd = merged.to_pandas()

# Option 1: Line plot (time-series comparison) - Use datetimes directly
plt.figure(figsize=(12, 6))
plt.plot(merged_pd['time'], merged_pd['hfss_stdv_file1'], label='File 1 hfss stdv')
plt.plot(merged_pd['time'], merged_pd['hfss_stdv_file2'], label='File 2 hfss stdv')
plt.xlabel('Time')
plt.ylabel('hfss stdv')
plt.title('Comparison of hfss stdv Over Time')
plt.legend()
plt.xticks(rotation=45)
plt.grid(True)
plt.tight_layout()
plt.show()

# Option 2: Scatter plot (point-by-point)
plt.figure(figsize=(8, 8))
plt.scatter(merged_pd['hfss_stdv_file1'], merged_pd['hfss_stdv_file2'], alpha=0.5)  # Alpha for density visibility
min_val = min(merged_pd['hfss_stdv_file1'].min(), merged_pd['hfss_stdv_file2'].min())
max_val = max(merged_pd['hfss_stdv_file1'].max(), merged_pd['hfss_stdv_file2'].max())
plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='y=x')
plt.xlabel('File 1 hfss stdv')
plt.ylabel('File 2 hfss stdv')
plt.title('Point-by-Point Comparison')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Option 3: Difference plot
merged_pd['diff'] = merged_pd['hfss_stdv_file1'] - merged_pd['hfss_stdv_file2']
plt.figure(figsize=(12, 6))
plt.plot(merged_pd['time'], merged_pd['diff'], color='green')
plt.axhline(0, color='red', linestyle='--')
plt.xlabel('Time')
plt.ylabel('Difference')
plt.title('Difference Over Time')
plt.xticks(rotation=45)
plt.grid(True)
plt.tight_layout()
plt.show()

# Stats (mean, std, corr) - Fast in Polars
print(merged.describe())
corr = merged.select(pl.corr('hfss_stdv_file1', 'hfss_stdv_file2')).item()
print(f"Correlation: {corr}")