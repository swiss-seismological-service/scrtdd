#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys

sns.set_theme(style="darkgrid")

def visualize_split_boxplots(file1_path, file2_path, num_bins=5, isXcorr=False):

    # Load data
    df1 = pd.read_csv(file1_path)
    df2 = pd.read_csv(file2_path)

    # Sort by distance first so the bin labels are created in order
    df1 = df1.sort_values('interEventDistance').reset_index(drop=True)
    df2 = df2.sort_values('interEventDistance').reset_index(drop=True)

    # Apply Equal-Frequency Binning (Quantiles) independently
    # pd.qcut divides the data into bins containing the same number of points
    # Create bins - keep as categorical to preserve mathematical order
    df1['DistBin'] = pd.qcut(df1['interEventDistance'], q=num_bins, duplicates='drop')
    df2['DistBin'] = pd.qcut(df2['interEventDistance'], q=num_bins, duplicates='drop')

    # Get the unique bins in the order they appear in the sorted data
    bin_order1 = df1['DistBin'].unique()
    bin_order2 = df2['DistBin'].unique()

    # Setup Side-by-Side Figure
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
    sns.set_theme(style="whitegrid")

    # Plot File 1
    sns.boxplot(
        ax=axes[0],
        data=df1,
        x='DistBin',
        y=('xcorrCoefficient' if isXcorr else 'doubleDifferenceResidual'),
        order=bin_order1, # Explicitly force the sorted order
        showfliers=False,
        color="skyblue"
    )
    axes[0].set_title('Starting Lcation', fontsize=14)
    axes[0].set_xlabel('Inter-Event Distance Bin [km]')
    axes[0].set_ylabel('Correlation Coefficient' if isXcorr else 'Double Difference Residual [s]')
    axes[0].tick_params(axis='x', rotation=30)
    # Add a dashed line at zero for reference
    axes[0].axhline(0, color='red', linestyle='--', alpha=0.5)

    # Plot File 2
    sns.boxplot(
        ax=axes[1],
        data=df2,
        x='DistBin',
        y=('xcorrCoefficient' if isXcorr else 'doubleDifferenceResidual'),
        order=bin_order2, # Explicitly force the sorted order
        showfliers=False,
        color="salmon"
    )
    axes[1].set_title('DD Relocated', fontsize=14)
    axes[1].set_xlabel('Inter-Event Distance Bin [km]')
    axes[1].set_ylabel('')
    axes[1].tick_params(axis='x', rotation=30)
    axes[1].axhline(0, color='red', linestyle='--', alpha=0.5)

    if isXcorr:
        plt.suptitle('Comparison of Correlation Coefficient By Inter-Event distance', fontsize=16, y=1.05)
    else:
        plt.suptitle('Comparison of Double-Difference Residual By Inter-Event distance', fontsize=16, y=1.05)

    fig.savefig('xcorr-comparison.png'  if isXcorr else 'dd-comparison.png', bbox_inches='tight')

if len(sys.argv) != 3:
    print("plot-dd.py cluster-1-initial-double-difference.csv cluster-1-final-double-difference.csv")
    exit(0)

visualize_split_boxplots(sys.argv[1], sys.argv[2], num_bins=20, isXcorr=False)
visualize_split_boxplots(sys.argv[1], sys.argv[2], num_bins=20, isXcorr=True)
