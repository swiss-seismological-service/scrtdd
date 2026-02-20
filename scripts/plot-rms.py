#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) != 2:
    print("plot-rms.py reloc-event.csv")
    exit(0)

# Load the data
df = pd.read_csv(sys.argv[1])

# Define common bins for both datasets
min_rms = min(df['startRms'].min(), df['finalRms'].min())
max_rms = max(df['startRms'].max(), df['finalRms'].max())
bins = np.linspace(min_rms, max_rms, 20)  # 20 bins

# Calculate the histogram counts for both
start_counts, _ = np.histogram(df['startRms'], bins=bins)
final_counts, _ = np.histogram(df['finalRms'], bins=bins)

# Convert counts to cumulative percentages
start_cum_pct = np.cumsum(start_counts) / len(df) * 100
final_cum_pct = np.cumsum(final_counts) / len(df) * 100

# Bar chart setup
bin_centers = (bins[:-1] + bins[1:]) / 2
width = (bins[1] - bins[0]) * 0.35  # Adjusted width of the bars

plt.figure(figsize=(10, 6))

# Plot the bars side-by-side
# We shift startRms left and finalRms right by half the width
plt.bar(bin_centers - width/2, start_cum_pct, width=width, 
        label='startRms', color='skyblue', edgecolor='black', alpha=0.8)

plt.bar(bin_centers + width/2, final_cum_pct, width=width, 
        label='finalRms', color='salmon', edgecolor='black', alpha=0.8)

plt.title('Cumulative RMS Histogram: Starting Lcation vs DD Relocated')
plt.xlabel('RMS [s]')
plt.ylabel('Percentage of Events (%)')
plt.xticks(bins.round(3), rotation=45) # Show bin edges on x-axis
plt.legend(loc='upper left')
plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.ylim(0, 110)

plt.savefig('rms-comparison.png')
