import pandas as pd

# Load your data
df = pd.read_csv("data/output/KAOU/KAOU2025_processed_v1_50k.csv")

# Columns you want to summarize
cols = [
    "hfss stdv","hfss mean","hfls stdv","hfls mean",
    "tau stdv","tau mean","mean dmo","T stdv",
    "TS stdv","RH stdv","P stdv","SPD stdv"
]

# Compute min, max, and mean for each column
summary = df[cols].agg(['min', 'max', 'mean']).transpose()

# Display neatly
print(summary)
