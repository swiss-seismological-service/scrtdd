#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

# Set the seaborn theme for aesthetics
sns.set_theme(style="darkgrid")

num_quantiles = 100
ylim = None

if len(sys.argv) < 2:
    print(f"{sys.argv[0]} reloc-event.csv [num-quantiles] [y-limit]")
    exit(0)

if len(sys.argv) >= 3:
    num_quantiles = int(sys.argv[2])

if len(sys.argv) >= 4:
    ylim = float(sys.argv[3])

# Load the data
df = pd.read_csv(sys.argv[1])

# Define the quantiles we want to calculate (from 0 to 1)
quantiles = np.linspace(0, 1, num_quantiles)
percentiles = quantiles * 100

# Calculate the RMS values at each quantile
start_rms_quantiles = df['startRms'].quantile(quantiles).rename('rms')
final_rms_quantiles = df['finalRms'].quantile(quantiles).rename('rms')

# Create a tidy DataFrame suitable for seaborn
plot_df_start = pd.DataFrame({'percentile': percentiles, 'rms': start_rms_quantiles, 'RMS Type': 'startRms'})
plot_df_final = pd.DataFrame({'percentile': percentiles, 'rms': final_rms_quantiles, 'RMS Type': 'finalRms'})
plot_df = pd.concat([plot_df_start, plot_df_final])

# Line plot setup
plt.figure(figsize=(12, 7))

# Plot
sns.lineplot(data=plot_df, x='percentile', y='rms', hue='RMS Type', style='RMS Type', markers=True, dashes=False)

plt.title('RMS Quantile Plot: Starting Location vs DD Relocated')
plt.xlabel('Percentage of Events (%)')
plt.ylabel('RMS [s]')
plt.legend(title='RMS Type', loc='upper left')
plt.ylim(bottom=0)  # RMS can't be negative
if ylim is not None:
    plt.ylim(top=ylim)
plt.xlim(0, 100)

# Set ticks for better readability
plt.xticks(np.arange(0, 101, 10))

plt.savefig('rms-comparison.png')

