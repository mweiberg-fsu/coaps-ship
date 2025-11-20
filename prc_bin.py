import pandas as pd
import matplotlib.pyplot as plt

# Load CSV
df = pd.read_csv("data\input\KAOU\KAOU2011_eval.csv")

# -----------------------------
# 1. Bin SHF in intervals of 5
# -----------------------------
bins = list(range(-40, 45, 5))  # from -40 to 40, step 5
# intervals like [-40, -35)
df['shf_bin'] = pd.cut(df['shf'], bins=bins, right=False)

# Count the number of values in each bin
bin_counts = df['shf_bin'].value_counts().sort_index()

# -----------------------------
# 2. Save binned data to CSV
# -----------------------------
df.to_csv("KAOU2011_eval_binned.csv", index=False)

# -----------------------------
# 3. Filter extreme SHF values
# -----------------------------
high_shf = df[(df['shf'] < -40) | (df['shf'] > 40)]

# Print extreme rows
print(high_shf)

# Save extreme rows to a separate CSV
high_shf.to_csv("KAOU2011_eval_high_shf.csv", index=False)

# -----------------------------
# 4. Plot binned counts
# -----------------------------
plt.figure(figsize=(12, 6))
plt.bar(bin_counts.index.astype(str), bin_counts.values,
        color='skyblue', edgecolor='black')
plt.xlabel('SHF Bins')
plt.ylabel('Count')
plt.title('Bin Distribution of SHF values (intervals of 5)')
plt.xticks(rotation=45)
plt.show()
