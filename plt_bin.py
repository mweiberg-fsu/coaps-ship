import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 1. Load CSV
df = pd.read_csv("processed/KAOU2011_processed_v1.csv")

# 2. Decide bin width
bin_width = 1000

# 3. Automatically create bins from min to max of 'mean dmo'
min_val = df['mean dmo'].min()
max_val = df['mean dmo'].max()
# +bin_width to include the last edge
bins = np.arange(min_val, max_val + bin_width, bin_width)

# 4. Bin the data
df['dmo_bin'] = pd.cut(df['mean dmo'], bins=bins, right=False)

# 5. Count per bin
bin_counts = df['dmo_bin'].value_counts().sort_index()

# 6. Plot
plt.figure(figsize=(10, 6))
bin_counts.plot(kind='bar', color='skyblue')
plt.xlabel('mean dmo bins')
plt.ylabel('Count')
plt.title('Distribution of mean dmo')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
