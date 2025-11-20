import polars as pl
import matplotlib.pyplot as plt

# Define the actual column names (from your first header line)
columns = [
    'time', 'platform_call_sign', 'latitude', 'longitude',
    'in_T', 'in_RH', 'in_TS', 'in_SPD', 'in_P'
]

# Lazy load: skip_rows=1 skips the header+units? No—skip the first two lines (header used manually, units skipped)
# We use has_header=False and rename to our columns
lf = (
    pl.scan_csv(
        'data/input/KAOU/KAOU2025.csv',  # Adjust path if needed
        skip_rows=2,  # Skip header and units lines; data starts at line 3
        has_header=False,
        try_parse_dates=False  # We'll parse time manually
    )
    .with_columns([
        pl.col('column_1').str.to_datetime('%Y-%m-%dT%H:%M:%SZ')  # Parse time from first column
    ])
    .rename(dict(zip([f'column_{i+1}' for i in range(len(columns))], columns)))  # Rename all to proper names
)

# Select, sort, and collect efficiently
data = (
    lf
    .select(['time', 'in_TS'])
    .sort('time')
    .collect(streaming=True)  # Chunked for low memory on large data
)

# Plot with Matplotlib (convert to numpy for speed)
plt.figure(figsize=(12, 6))
plt.plot(data['time'].to_numpy(), data['in_TS'].to_numpy(), label='Wind Speed (m/s)')
plt.xlabel('Time (UTC)')
plt.ylabel('Wind Speed (m/s)')
plt.title('Wind Speed Over Time')
plt.legend()
plt.grid(True)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# Optional: Resample for large spans (e.g., hourly mean to reduce plot clutter)
# resampled = data.group_by_dynamic('time', every='1h').agg(pl.col('in_SPD').mean())
# Then plot resampled['time'], resampled['in_SPD']