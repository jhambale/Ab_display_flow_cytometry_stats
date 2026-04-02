#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import re
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import argparse


def analyze_flow_cytometry_expression(csv_file):
    """
    Analyze flow cytometry data to calculate:
    1. Expression fraction: (+/+ and +/-) / (all 4 quadrants)
    2. Percentage of expressing cells that are in R6
    
    Parameters:
    csv_file: path to the CSV file containing flow cytometry statistics
    
    Returns:
    DataFrame with results for each sample
    """
    
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Get unique samples
    unique_samples = df['Sample'].unique()
    
    # Initialize results list
    results = []
    
    for sample in unique_samples:
        # Filter data for current sample
        sample_data = df[df['Sample'] == sample]
        
        # Initialize values for quadrants (using Count column for absolute numbers)
        concentration = None
        pos_pos_count = None  # +/+(1)
        pos_neg_count = None  # +/-(1)
        neg_pos_count = None  # -/+(1)
        neg_neg_count = None  # -/-(1)
        r6_count = None
        
        # Extract counts from Count column for specific gates
        for _, row in sample_data.iterrows():
            gate = row['Gate']
            
            if gate == '+/+(1)':
                pos_pos_count = row['Count']
            elif gate == '+/-(1)':
                pos_neg_count = row['Count']
            elif gate == '-/+(1)':
                neg_pos_count = row['Count']
            elif gate == '-/-(1)':
                neg_neg_count = row['Count']
            elif gate == 'R6':
                r6_count = row['Count']
        
        # Calculate total cells in all 4 quadrants
        if all(x is not None for x in [pos_pos_count, pos_neg_count, neg_pos_count, neg_neg_count]):
            total_quadrant_cells = pos_pos_count + pos_neg_count + neg_pos_count + neg_neg_count
            expressing_cells = pos_pos_count + neg_pos_count
            
            # Calculate expression fraction
            if total_quadrant_cells > 0:
                expression_fraction = (expressing_cells / total_quadrant_cells) * 100
            else:
                expression_fraction = 0
            
            # Calculate R6 percentage of expressing cells
            if expressing_cells > 0 and r6_count is not None:
                r6_percent_of_expressing = (r6_count / expressing_cells) * 100
            else:
                r6_percent_of_expressing = None if r6_count is None else 0
        else:
            expression_fraction = None
            expressing_cells = None
            r6_percent_of_expressing = None
        
        # Extract Ag concentration from sample name
        if '120nMAg' in sample:
            concentration = 120
        elif '30nMAg' in sample:
            concentration = 30
        elif '0nMAg' in sample:
            concentration = 0
        # If still None, try regex as backup
        if concentration is None:
            pattern = r'_(\d+)nMAg'
            match = re.search(pattern, sample)
            if match:
                concentration = int(match.group(1))
        
        # Extract sample ID (e.g., y66, y67, etc.)
        sample_id = sample.split('_')[0] if '_' in sample else sample
        
        # Append to results
        results.append({
            'Sample_ID': sample_id,
            'Sample_Full': sample,
            'Ag_Concentration_nM': concentration,
            '+/+(1)_Count': pos_pos_count,
            '+/-(1)_Count': pos_neg_count,
            '-/+(1)_Count': neg_pos_count,
            '-/-(1)_Count': neg_neg_count,
            'Total_Quadrant_Cells': total_quadrant_cells if 'total_quadrant_cells' in locals() else None,
            'Expressing_Cells': expressing_cells,
            'Expression_Fraction_%': expression_fraction,
            'R6_Count': r6_count,
            'R6_Percent_of_Expressing': r6_percent_of_expressing
        })
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results)
    
    # Sort by Sample_ID and Ag concentration
    results_df = results_df.sort_values(['Sample_ID', 'Ag_Concentration_nM'])
    
    return results_df

def analyze_flow_cytometry_mfi(csv_file):
    """
    Complete MFI-based flow cytometry analysis function
    
    Analyzes Mean Fluorescence Intensity (MFI) data from flow cytometry,
    calculating MFI ratios, fold changes, and binding metrics.
    
    Parameters:
    -----------
    csv_file : str
        Path to the CSV file containing flow cytometry statistics
    
    Returns:
    --------
    tuple : (mfi_results_df, summary_stats_df, plots)
        - mfi_results_df: Detailed MFI analysis for each sample
        - summary_stats_df: Summary statistics by concentration
        - plots: Generated visualization figures
    """
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Get unique samples
    unique_samples = df['Sample'].unique()
    print(f"Analyzing {len(unique_samples)} samples for MFI...")
    
    # Initialize results list
    results = []
    
    for sample in unique_samples:
        # Filter data for current sample
        sample_data = df[df['Sample'] == sample]
        
        # Initialize MFI storage
        mfi_data = {}
        
        # Gates of interest for quadrant analysis
        quadrant_gates = ['+/+(1)', '+/-(1)', '-/+(1)', '-/-(1)']
        
        # Extract MFI values for each gate
        for _, row in sample_data.iterrows():
            gate = row['Gate']
            
            # Store X and Y MFI for relevant gates
            if gate in quadrant_gates or gate == 'R6':
                mfi_data[f'{gate}_X_MFI'] = row['X Mean']
                mfi_data[f'{gate}_Y_MFI'] = row['Y Mean']
                mfi_data[f'{gate}_Count'] = row['Count']  # Keep count for reference
        
        # Calculate expressing population MFI (geometric mean of +/+ and +/-)
        pos_pos_x = mfi_data.get('+/+(1)_X_MFI', 0)
        pos_pos_y = mfi_data.get('+/+(1)_Y_MFI', 0)
        pos_neg_x = mfi_data.get('+/-(1)_X_MFI', 0)
        pos_neg_y = mfi_data.get('+/-(1)_Y_MFI', 0)
        
        # Geometric mean for expressing cells
        if pos_pos_x > 0 and pos_neg_x > 0:
            expressing_x_mfi_geomean = np.sqrt(pos_pos_x * pos_neg_x)
        else:
            expressing_x_mfi_geomean = np.mean([pos_pos_x, pos_neg_x])
        
        if pos_pos_y > 0 and pos_neg_y > 0:
            expressing_y_mfi_geomean = np.sqrt(pos_pos_y * pos_neg_y)
        else:
            expressing_y_mfi_geomean = np.mean([pos_pos_y, pos_neg_y])
        
        # Calculate non-expressing population MFI (geometric mean of -/- and -/+)
        neg_neg_x = mfi_data.get('-/-(1)_X_MFI', 0)
        neg_neg_y = mfi_data.get('-/-(1)_Y_MFI', 0)
        neg_pos_x = mfi_data.get('-/+(1)_X_MFI', 0)
        neg_pos_y = mfi_data.get('-/+(1)_Y_MFI', 0)
        
        if neg_neg_x > 0 and neg_pos_x > 0:
            nonexpressing_x_mfi_geomean = np.sqrt(neg_neg_x * neg_pos_x)
        else:
            nonexpressing_x_mfi_geomean = np.mean([neg_neg_x, neg_pos_x])
        
        if neg_neg_y > 0 and neg_pos_y > 0:
            nonexpressing_y_mfi_geomean = np.sqrt(neg_neg_y * neg_pos_y)
        else:
            nonexpressing_y_mfi_geomean = np.mean([neg_neg_y, neg_pos_y])
        
        # Calculate MFI fold change (expressing/non-expressing)
        if nonexpressing_x_mfi_geomean > 0:
            x_mfi_fold_change = expressing_x_mfi_geomean / nonexpressing_x_mfi_geomean
        else:
            x_mfi_fold_change = None
        
        if nonexpressing_y_mfi_geomean > 0:
            y_mfi_fold_change = expressing_y_mfi_geomean / nonexpressing_y_mfi_geomean
        else:
            y_mfi_fold_change = None
        
        # R6 MFI analysis
        r6_x_mfi = mfi_data.get('R6_X_MFI', 0)
        r6_y_mfi = mfi_data.get('R6_Y_MFI', 0)
        
        # Calculate R6 MFI enrichment over expressing population
        if expressing_x_mfi_geomean > 0 and r6_x_mfi > 0:
            r6_x_enrichment = r6_x_mfi / expressing_x_mfi_geomean
        else:
            r6_x_enrichment = None
        
        if expressing_y_mfi_geomean > 0 and r6_y_mfi > 0:
            r6_y_enrichment = r6_y_mfi / expressing_y_mfi_geomean
        else:
            r6_y_enrichment = None
        
        # Calculate combined MFI score (product of X and Y MFI)
        expressing_combined_mfi = expressing_x_mfi_geomean * expressing_y_mfi_geomean
        r6_combined_mfi = r6_x_mfi * r6_y_mfi if r6_x_mfi and r6_y_mfi else 0
        
        # Extract Ag concentration
        concentration = None
        if '0nMAg' in sample:
            concentration = 0
        elif '30nMAg' in sample:
            concentration = 30
        elif '120nMAg' in sample:
            concentration = 120
        
        # Extract sample ID
        sample_id = sample.split('_')[0] if '_' in sample else sample
        
        # Build results dictionary
        result_dict = {
            'Sample_ID': sample_id,
            'Sample_Full': sample,
            'Ag_Concentration_nM': concentration,
            # Expressing population MFI
            'Expressing_X_MFI': expressing_x_mfi_geomean,
            'Expressing_Y_MFI': expressing_y_mfi_geomean,
            'Expressing_Combined_MFI': expressing_combined_mfi,
            # Non-expressing population MFI
            'NonExpressing_X_MFI': nonexpressing_x_mfi_geomean,
            'NonExpressing_Y_MFI': nonexpressing_y_mfi_geomean,
            # Fold changes
            'X_MFI_Fold_Change': x_mfi_fold_change,
            'Y_MFI_Fold_Change': y_mfi_fold_change,
            # R6 gate MFI
            'R6_X_MFI': r6_x_mfi,
            'R6_Y_MFI': r6_y_mfi,
            'R6_Combined_MFI': r6_combined_mfi,
            # R6 enrichment
            'R6_X_Enrichment': r6_x_enrichment,
            'R6_Y_Enrichment': r6_y_enrichment,
        }
        
        # Add individual gate MFI values
        result_dict.update(mfi_data)
        
        results.append(result_dict)
    
    # Create DataFrame
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(['Sample_ID', 'Ag_Concentration_nM'])
    
    return results_df

def calculate_mfi_statistics(results_df):
     """
    Calculate summary statistics for MFI data grouped by concentration
    """
    # Group by concentration
    grouped = results_df.groupby('IL17_Concentration_nM')
    
    # Calculate statistics
    stats_dict = {}
    
    metrics = [
        'Expressing_X_MFI', 'Expressing_Y_MFI', 'Expressing_Combined_MFI',
        'X_MFI_Fold_Change', 'Y_MFI_Fold_Change',
        'R6_X_MFI', 'R6_Y_MFI', 'R6_Combined_MFI',
        'R6_X_Enrichment', 'R6_Y_Enrichment'
    ]
    
    for metric in metrics:
        if metric in results_df.columns:
            stats_dict[f'{metric}_mean'] = grouped[metric].mean()
            stats_dict[f'{metric}_std'] = grouped[metric].std()
            stats_dict[f'{metric}_sem'] = grouped[metric].sem()
            stats_dict[f'{metric}_median'] = grouped[metric].median()
    
    summary_df = pd.DataFrame(stats_dict)
    
    return summary_df

def plot_mfi_analysis(results_df):
     """
    Create comprehensive MFI visualization plots
    """
    # Set style
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Create figure with subplots
    fig = plt.figure(figsize=(20, 12))
    
    # 1. Expression MFI by concentration
    ax1 = plt.subplot(2, 3, 1)
    plot_concentration_response(results_df, 'Expressing_Combined_MFI', 
                                ax1, 'Expression Combined MFI')
    
    # 2. MFI Fold Change by concentration
    ax2 = plt.subplot(2, 3, 2)
    plot_concentration_response(results_df, 'X_MFI_Fold_Change', 
                                ax2, 'X-Channel MFI Fold Change')
    
    # 3. R6 MFI by concentration
    ax3 = plt.subplot(2, 3, 3)
    plot_concentration_response(results_df, 'R6_Combined_MFI', 
                                ax3, 'R6 Combined MFI')
    
    # 4. R6 Enrichment
    ax4 = plt.subplot(2, 3, 4)
    plot_concentration_response(results_df, 'R6_X_Enrichment', 
                                ax4, 'R6 X-Channel Enrichment')
    
    # 5. Heatmap of MFI values
    ax5 = plt.subplot(2, 3, 5)
    plot_mfi_heatmap(results_df, ax5)
    
    # 6. Correlation plot
    ax6 = plt.subplot(2, 3, 6)
    plot_mfi_correlation(results_df, ax6)
    
    plt.suptitle('Flow Cytometry MFI Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    return fig

def plot_concentration_response(df, metric, ax, title):
    """
    Plot concentration-response curves for a specific MFI metric
    """
    # Get unique samples
    df['Sample_Base'] = df['Sample_ID'].str.extract(r'(y\d+)', expand=False)
    unique_samples = df['Sample_Base'].dropna().unique()
    
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_samples)))
    concentrations = [0, 30, 120]
    
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = df[df['Sample_Base'] == sample].copy()
        sample_data = sample_data.sort_values('IL17_Concentration_nM')
        
        x_vals = []
        y_vals = []
        for conc in concentrations:
            conc_data = sample_data[sample_data['IL17_Concentration_nM'] == conc]
            if not conc_data.empty and pd.notna(conc_data[metric].iloc[0]):
                x_vals.append(conc)
                y_vals.append(conc_data[metric].iloc[0])
        
        if len(x_vals) > 0:
            ax.plot(x_vals, y_vals, marker='o', label=sample, 
                   color=colors[i], linewidth=2, markersize=8, alpha=0.7)
    
    ax.set_xlabel('IL17 Concentration (nM)', fontsize=10)
    ax.set_ylabel(title, fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.set_xticks(concentrations)
    ax.grid(True, alpha=0.3)
    
    # Add mean line
    mean_vals = []
    for conc in concentrations:
        conc_data = df[df['IL17_Concentration_nM'] == conc][metric].dropna()
        if not conc_data.empty:
            mean_vals.append((conc, conc_data.mean()))
    
    if mean_vals:
        x_mean, y_mean = zip(*mean_vals)
        ax.plot(x_mean, y_mean, 'k--', linewidth=2, label='Mean', alpha=0.8)

def plot_mfi_heatmap(df, ax):
    """
    Create a heatmap of MFI values across samples and concentrations
    """
    # Pivot data for heatmap
    pivot_data = df.pivot_table(
        values='Expressing_Combined_MFI',
        index='Sample_ID',
        columns='IL17_Concentration_nM',
        aggfunc='first'
    )
    
    # Create heatmap
    sns.heatmap(pivot_data, annot=True, fmt='.0f', cmap='YlOrRd', 
                ax=ax, cbar_kws={'label': 'Combined MFI'})
    ax.set_title('Expression MFI Heatmap', fontsize=11, fontweight='bold')
    ax.set_xlabel('IL17 Concentration (nM)', fontsize=10)
    ax.set_ylabel('Sample ID', fontsize=10)

def plot_mfi_correlation(df, ax):
    """
    Plot correlation between expression MFI and R6 MFI
    """
    # Filter valid data
    valid_data = df.dropna(subset=['Expressing_Combined_MFI', 'R6_Combined_MFI'])
    
    if not valid_data.empty:
        x = valid_data['Expressing_Combined_MFI']
        y = valid_data['R6_Combined_MFI']
        
        # Color by concentration
        colors = valid_data['IL17_Concentration_nM'].map({0: 'blue', 30: 'green', 120: 'red'})
        
        ax.scatter(x, y, c=colors, alpha=0.6, s=50)
        
        # Add trend line
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        ax.plot(x, p(x), "k--", alpha=0.5)
        
        # Calculate correlation
        corr = x.corr(y)
        ax.text(0.05, 0.95, f'Correlation: {corr:.3f}', 
                transform=ax.transAxes, fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_xlabel('Expression Combined MFI', fontsize=10)
        ax.set_ylabel('R6 Combined MFI', fontsize=10)
        ax.set_title('Expression vs R6 Binding MFI', fontsize=11, fontweight='bold')
        
        # Add legend for concentrations
        for conc, color in [(0, 'blue'), (30, 'green'), (120, 'red')]:
            ax.scatter([], [], c=color, label=f'{conc} nM', s=50)
        ax.legend(title='IL17', loc='lower right')

# Main execution function
def run_mfi_analysis(csv_file, output_prefix='mfi_analysis'):
    """
    Complete MFI analysis pipeline
    """
    print("="*80)
    print("FLOW CYTOMETRY MFI ANALYSIS")
    print("="*80)
    
    # Analyze MFI data
    results_df = analyze_flow_cytometry_mfi(csv_file)
    
    # Calculate statistics
    summary_df = calculate_mfi_statistics(results_df)
    
    # Display results
    print("\nSample Results (first 10):")
    print("-"*80)
    display_cols = ['Sample_ID', 'IL17_Concentration_nM', 
                   'Expressing_Combined_MFI', 'X_MFI_Fold_Change', 
                   'R6_Combined_MFI', 'R6_X_Enrichment']
    print(results_df[display_cols].head(10).to_string(index=False))
    
    print("\n" + "-"*80)
    print("Summary Statistics by IL17 Concentration:")
    print("-"*80)
    print(summary_df.round(2))
    
    # Save results
    results_df.to_csv(f'{output_prefix}_results.csv', index=False)
    summary_df.to_csv(f'{output_prefix}_summary.csv')
    print(f"\n✓ Results saved to '{output_prefix}_results.csv'")
    print(f"✓ Summary saved to '{output_prefix}_summary.csv'")
    
    # Create plots
    fig = plot_mfi_analysis(results_df)
    plt.savefig(f'{output_prefix}_plots.png', dpi=300, bbox_inches='tight')
    print(f"✓ Plots saved to '{output_prefix}_plots.png'")
    
    plt.show()
    
    return results_df, summary_df


def summarize_by_concentration(results_df):
    """
    Summarize results by Ag concentration
    """
    # Filter out rows with None values for summary
    valid_df = results_df.dropna(subset=['Expression_Fraction_%', 'R6_Percent_of_Expressing'])
    
    if not valid_df.empty:
        summary = valid_df.groupby('Ag_Concentration_nM').agg({
            'Expression_Fraction_%': ['mean', 'std', 'count'],
            'R6_Percent_of_Expressing': ['mean', 'std', 'count']
        }).round(2)
    else:
        summary = pd.DataFrame()
    
    return summary


def plot_concentration_curves(results_df, dir_prefix, save_plots=True):
    """
    Plot concentration-response curves for each sample showing:
    1. Expression fraction vs Ag concentration
    2. R6% of expressing cells vs Ag concentration
    
    Parameters:
    results_df: DataFrame with analysis results
    save_plots: Whether to save plots to files
    """
    
    # Set up the style
    plt.style.use('seaborn-v0_8-darkgrid')
    
    # Get unique sample IDs (without concentration)
    results_df['Sample_Base'] = results_df['Sample_ID'].str.extract(r'(y\d+)', expand=False)
    unique_samples = results_df['Sample_Base'].dropna().unique()
    
    # Define concentration points
    concentrations = [0, 30, 120]
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    
    # Color palette for different samples
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_samples)))
    
    # Plot 1: Expression Fraction vs Concentration
    ax1 = axes[0]
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = results_df[results_df['Sample_Base'] == sample].copy()
        sample_data = sample_data.sort_values('Ag_Concentration_nM')
        
        # Get data points for this sample
        x_vals = []
        y_vals = []
        for conc in concentrations:
            conc_data = sample_data[sample_data['Ag_Concentration_nM'] == conc]
            if not conc_data.empty and pd.notna(conc_data['Expression_Fraction_%'].iloc[0]):
                x_vals.append(conc)
                y_vals.append(conc_data['Expression_Fraction_%'].iloc[0])
        
        if len(x_vals) > 0:
            ax1.plot(x_vals, y_vals, marker='o', label=sample, 
                    color=colors[i], linewidth=2, markersize=8)
    
    ax1.set_xlabel('Ag Concentration (nM)', fontsize=12)
    ax1.set_ylabel('Expression Fraction (%)', fontsize=12)
    ax1.set_title('Expression Fraction vs Ag Concentration', fontsize=14, fontweight='bold')
    ax1.set_xticks(concentrations)
    ax1.grid(True, alpha=0.3)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)
    
    # Plot 2: R6% of Expressing Cells vs Concentration
    ax2 = axes[1]
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = results_df[results_df['Sample_Base'] == sample].copy()
        sample_data = sample_data.sort_values('Ag_Concentration_nM')
        
        # Get data points for this sample
        x_vals = []
        y_vals = []
        for conc in concentrations:
            conc_data = sample_data[sample_data['Ag_Concentration_nM'] == conc]
            if not conc_data.empty and pd.notna(conc_data['R6_Percent_of_Expressing'].iloc[0]):
                x_vals.append(conc)
                y_vals.append(conc_data['R6_Percent_of_Expressing'].iloc[0])
        
        if len(x_vals) > 0:
            ax2.plot(x_vals, y_vals, marker='s', label=sample, 
                    color=colors[i], linewidth=2, markersize=8)
    
    ax2.set_xlabel('Ag Concentration (nM)', fontsize=12)
    ax2.set_ylabel('R6% of Expressing Cells', fontsize=12)
    ax2.set_title('R6 Binding (% of Expressing Cells) vs Ag Concentration', 
                  fontsize=14, fontweight='bold')
    ax2.set_xticks(concentrations)
    ax2.grid(True, alpha=0.3)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)
    
    plt.tight_layout()
    
    if save_plots:
        plt.savefig(dir_prefix +'/concentration_response_curves.png', dpi=300, bbox_inches='tight')
        print("\n✓ Plots saved to 'concentration_response_curves.png'")
    
    plt.show()
    
    # Create individual plots for each sample
    create_individual_sample_plots(results_df, unique_samples, dir_prefix, concentrations, save_plots)

def create_individual_sample_plots(results_df, unique_samples, dir_prefix, concentrations, save_plots=True):
    """
    Create individual plots for each sample showing both metrics
    """
    
    # Calculate grid dimensions
    n_samples = len(unique_samples)
    n_cols = 4
    n_rows = (n_samples + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
    axes = axes.flatten() if n_samples > 1 else [axes]
    
    for idx, sample in enumerate(sorted(unique_samples)):
        ax = axes[idx]
        
        sample_data = results_df[results_df['Sample_Base'] == sample].copy()
        sample_data = sample_data.sort_values('Ag_Concentration_nM')
        
        # Get expression fraction data
        x_expr = []
        y_expr = []
        for conc in concentrations:
            conc_data = sample_data[sample_data['Ag_Concentration_nM'] == conc]
            if not conc_data.empty and pd.notna(conc_data['Expression_Fraction_%'].iloc[0]):
                x_expr.append(conc)
                y_expr.append(conc_data['Expression_Fraction_%'].iloc[0])
        
        # Get R6 data
        x_r6 = []
        y_r6 = []
        for conc in concentrations:
            conc_data = sample_data[sample_data['Ag_Concentration_nM'] == conc]
            if not conc_data.empty and pd.notna(conc_data['R6_Percent_of_Expressing'].iloc[0]):
                x_r6.append(conc)
                y_r6.append(conc_data['R6_Percent_of_Expressing'].iloc[0])
        
        # Create twin axis
        ax2 = ax.twinx()
        
        # Plot data
        if len(x_expr) > 0:
            line1 = ax.plot(x_expr, y_expr, 'b-o', label='Expression %', 
                           linewidth=2, markersize=8)
        if len(x_r6) > 0:
            line2 = ax2.plot(x_r6, y_r6, 'r-s', label='Binding %', 
                            linewidth=2, markersize=8)
        
        ax.set_xlabel('Ag (nM)', fontsize=10)
        ax.set_ylabel('Expression Fraction (%)', color='b', fontsize=10)
        ax.set_ylim(0,100) # full scale for expression from no cells to all cells expressing
        ax2.set_ylabel('Binding Fraction (%)', color='r', fontsize=10)
        ax2.set_ylim(0,100) # scale binding from no expressing cells to all expressing cells binding
        ax.set_title(f'Sample {sample}', fontsize=12, fontweight='bold')
        ax.set_xticks(concentrations)
        ax.grid(True, alpha=0.3)
        
        # Color the y-axis labels
        ax.tick_params(axis='y', labelcolor='b')
        ax2.tick_params(axis='y', labelcolor='r')
        
        # Add legend
        if len(x_expr) > 0 and len(x_r6) > 0:
            lines = line1 + line2
            labels = [l.get_label() for l in lines]
            ax.legend(lines, labels, loc='upper left', fontsize=8)
    
    # Hide empty subplots
    for idx in range(n_samples, len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle('Individual Sample Concentration-Response Curves', 
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    if save_plots:
        plt.savefig(dir_prefix + '/individual_sample_curves.png', dpi=300, bbox_inches='tight')
        print("✓ Individual sample plots saved to 'individual_sample_curves.png'")
    
    plt.show()

# Main execution
if __name__ == "__main__":
    # Specify your CSV file path
    # csv_file = "2026-01-14_Agscfv_stats.csv"
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs='+', help='input mutation analysis csv files')
    parser.add_argument(-'m', help='type of metric for analysis (Count, Mean, or Median)')

    args = parser.parse_args()

    csv_file = args.i
    # print(csv_file[0])
    dir_prefix =  '/'.join(csv_file[0].split('/')[:-1])
    print(f'output directory: {dir_prefix}')

    # determine desired metric
    metric = str(args.m)
    
    # Analyze the data
    print("Analyzing flow cytometry data...")
    print("Calculating expression fractions and R6 binding percentages...")
    if metric = 'Count':
        results = analyze_flow_cytometry_expression(csv_file[0])
    else metric = 'Mean'
        results = analyze_flow_cytometry_mfi(csv_file[0])
    else metric = 'Median'
        results = analyze_flow_cytometry_median(csv_file[0])
    
    # Display results
    print("\n" + "="*100)
    print("FLOW CYTOMETRY EXPRESSION ANALYSIS RESULTS")
    print("="*100)
    
    # Show detailed results
    print("\nDetailed Results (first 15 samples):")
    print("-"*100)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.float_format', '{:.2f}'.format)
    
    display_cols = ['Sample_ID', 'Ag_Concentration_nM', 
                    'Expression_Fraction_%', 'R6_Percent_of_Expressing']
    print(results[display_cols].head(15).to_string(index=False))
    
    # Summary statistics by concentration
    print("\n" + "-"*100)
    print("Summary Statistics by Ag Concentration:")
    print("-"*100)
    summary = summarize_by_concentration(results)
    if not summary.empty:
        print(summary)
    
    # Save detailed results to CSV
    output_file = dir_prefix + "/flow_cytometry_expression_analysis.csv"
    results.to_csv(output_file, index=False)
    print(f"\n✓ Detailed results saved to '{output_file}'")
    
    # Analyze trends
    print("\n" + "-"*100)
    print("Key Findings:")
    print("-"*100)
    
    # Find samples with highest expression
    top_expressors = results.nlargest(5, 'Expression_Fraction_%', keep='all')
    print("\nTop 5 Samples by Expression Fraction:")
    print(top_expressors[display_cols].to_string(index=False))
    
    # Find samples with highest R6 binding among expressing cells
    valid_r6 = results.dropna(subset=['R6_Percent_of_Expressing'])
    if not valid_r6.empty:
        top_r6 = valid_r6.nlargest(5, 'R6_Percent_of_Expressing', keep='all')
        print("\nTop 5 Samples by R6% of Expressing Cells:")
        print(top_r6[display_cols].to_string(index=False))
    
    # Calculate average expression by concentration
    print("\n" + "-"*100)
    print("Average Expression Fraction by Ag Concentration:")
    print("-"*100)
    
    for conc in sorted(results['Ag_Concentration_nM'].dropna().unique()):
        conc_data = results[results['Ag_Concentration_nM'] == conc]
        valid_data = conc_data.dropna(subset=['Expression_Fraction_%'])
        if not valid_data.empty:
            mean_expr = valid_data['Expression_Fraction_%'].mean()
            std_expr = valid_data['Expression_Fraction_%'].std()
            n_samples = len(valid_data)
            print(f"  {conc:3.0f} nM: {mean_expr:6.2f}% ± {std_expr:5.2f}% (n={n_samples})")
    
    # Calculate average R6 binding by concentration
    print("\n" + "-"*100)
    print("Average R6% of Expressing Cells by Ag Concentration:")
    print("-"*100)
    
    for conc in sorted(results['Ag_Concentration_nM'].dropna().unique()):
        conc_data = results[results['Ag_Concentration_nM'] == conc]
        valid_data = conc_data.dropna(subset=['R6_Percent_of_Expressing'])
        if not valid_data.empty:
            mean_r6 = valid_data['R6_Percent_of_Expressing'].mean()
            std_r6 = valid_data['R6_Percent_of_Expressing'].std()
            n_samples = len(valid_data)
            print(f"  {conc:3.0f} nM: {mean_r6:6.2f}% ± {std_r6:5.2f}% (n={n_samples})")
    
    # Check for correlation
    print("\n" + "-"*100)
    print("Correlation Analysis:")
    print("-"*100)
    
    valid_corr = results.dropna(subset=['Expression_Fraction_%', 'R6_Percent_of_Expressing'])
    if len(valid_corr) > 1:
        correlation = valid_corr['Expression_Fraction_%'].corr(valid_corr['R6_Percent_of_Expressing'])
        print(f"Correlation between Expression Fraction and R6% of Expressing: {correlation:.3f}")
    
    print("\n" + "="*100)
    print("Analysis Complete!")
    print("="*100)
    
    # Optional: Create a summary CSV with just the key metrics
    summary_output = dir_prefix + "/flow_cytometry_summary.csv"
    summary_df = results[['Sample_ID', 'Ag_Concentration_nM', 
                          'Expression_Fraction_%', 'R6_Percent_of_Expressing']].copy()
    summary_df.to_csv(summary_output, index=False)
    print(f"\n✓ Summary results saved to '{summary_output}'")

    # After getting results
    results = analyze_flow_cytometry_expression(csv_file[0])
    
    # Add plotting
    print("\n" + "-"*80)
    print("Generating Concentration-Response Curves...")
    print("-"*80)
    
    # Create the plots
    plot_concentration_curves(results, dir_prefix, save_plots=True)
    
    # Optional: Create a summary table for easy viewing
    print("\n" + "-"*80)
    print("Sample Summary by Concentration:")
    print("-"*80)
    
    pivot_expr = results.pivot_table(
        values='Expression_Fraction_%',
        index='Sample_ID',
        columns='Ag_Concentration_nM',
        aggfunc='first'
    ).round(2)
    
    pivot_r6 = results.pivot_table(
        values='R6_Percent_of_Expressing',
        index='Sample_ID',
        columns='Ag_Concentration_nM',
        aggfunc='first'
    ).round(2)
    
    print("\nExpression Fraction (%) by Concentration:")
    print(pivot_expr)
    
    print("\nR6% of Expressing Cells by Concentration:")
    print(pivot_r6)