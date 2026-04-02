#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

def plot_concentration_curves(results_df, save_plots=True):
    """
    Plot concentration response curves for expression and R6 binding.
    
    Parameters:
    -----------
    results_df : pd.DataFrame
        DataFrame with analysis results from analyze_flow_cytometry_expression
    save_plots : bool
        Whether to save plots to files
    
    Returns:
    --------
    fig : matplotlib.figure.Figure
        The generated figure
    """
    
    # Set style
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Get unique samples - handle different sample ID formats
    results_df['Sample_Base'] = results_df['Sample_ID'].str.extract(r'(y\d+)', expand=False)
    unique_samples = results_df['Sample_Base'].dropna().unique()
    
    # Check if we have data
    if len(unique_samples) == 0:
        print("Warning: No valid samples found in data")
        return fig
    
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_samples)))
    
    concentrations = [0, 30, 120]
    
    # Plot 1: Expression Fraction vs Concentration
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = results_df[results_df['Sample_Base'] == sample]
        
        x_vals = []
        y_vals = []
        
        for conc in concentrations:
            conc_data = sample_data[sample_data['Ag_Concentration_nM'] == conc]
            if not conc_data.empty:
                x_vals.append(conc)
                y_vals.append(conc_data['Expression_Fraction_%'].iloc[0])
        
        if len(x_vals) > 0:  # Only plot if we have data
            ax1.plot(x_vals, y_vals, marker='o', label=sample, 
                    color=colors[i], linewidth=2, markersize=8, alpha=0.7)
    
    ax1.set_xlabel('Ag Concentration (nM)', fontsize=12)
    ax1.set_ylabel('Expression Fraction (%)', fontsize=12)
    ax1.set_title('Ag Dose Response - Expression', fontsize=14, fontweight='bold')
    ax1.set_xticks(concentrations)
    ax1.grid(True, alpha=0.3)
    
    # Only add legend if we have data
    if len(unique_samples) > 0:
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Plot 2: R6 Binding vs Concentration
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = results_df[results_df['Sample_Base'] == sample]
        
        x_vals = []
        y_vals = []
        
        for conc in concentrations:
            conc_data = sample_data[sample_data['Ag_Concentration_nM'] == conc]
            if not conc_data.empty and pd.notna(conc_data['R6_Percent_of_Expressing'].iloc[0]):
                x_vals.append(conc)
                y_vals.append(conc_data['R6_Percent_of_Expressing'].iloc[0])
        
        if len(x_vals) > 0:  # Only plot if we have data
            ax2.plot(x_vals, y_vals, marker='s', label=sample,
                    color=colors[i], linewidth=2, markersize=8, alpha=0.7)
    
    ax2.set_xlabel('Ag Concentration (nM)', fontsize=12)
    ax2.set_ylabel('R6 Binding (% of Expressing)', fontsize=12)
    ax2.set_title('Ag Dose Response - R6 Binding', fontsize=14, fontweight='bold')
    ax2.set_xticks(concentrations)
    ax2.grid(True, alpha=0.3)
    
    plt.suptitle('Flow Cytometry Cell Fraction Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    if save_plots:
        plt.savefig('cell_fraction_curves.png', dpi=300, bbox_inches='tight')
        print("Saved: cell_fraction_curves.png")
    
    return fig

def create_heatmap(results_df):
    """
    Create heatmaps for expression and R6 binding data.
    
    Parameters:
    -----------
    results_df : pd.DataFrame
        DataFrame with analysis results
    
    Returns:
    --------
    fig : matplotlib.figure.Figure
        The generated figure
    """
    
    # Check if we have data
    if results_df.empty:
        print("Warning: No data available for heatmap")
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        return fig
    
    # Create pivot tables for heatmaps
    # Handle case where pivot might fail
    try:
        pivot_expr = results_df.pivot_table(
            values='Expression_Fraction_%',
            index='Sample_ID',
            columns='Ag_Concentration_nM',
            aggfunc='first'
        )
        
        pivot_r6 = results_df.pivot_table(
            values='R6_Percent_of_Expressing',
            index='Sample_ID', 
            columns='Ag_Concentration_nM',
            aggfunc='first'
        )
    except Exception as e:
        print(f"Warning: Could not create pivot tables: {e}")
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        return fig
    
    # Check if pivots are empty
    if pivot_expr.empty or pivot_r6.empty:
        print("Warning: Pivot tables are empty")
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        return fig
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Heatmap 1: Expression Fraction
    if not pivot_expr.empty:
        sns.heatmap(pivot_expr, annot=True, fmt='.1f', cmap='YlOrRd', 
                    ax=axes[0], cbar_kws={'label': 'Expression %'})
        axes[0].set_title('Expression Fraction Heatmap', fontweight='bold')
        axes[0].set_xlabel('Ag Concentration (nM)')
        axes[0].set_ylabel('Sample ID')
    
    # Heatmap 2: R6 Binding
    if not pivot_r6.empty:
        sns.heatmap(pivot_r6, annot=True, fmt='.1f', cmap='YlGnBu',
                    ax=axes[1], cbar_kws={'label': 'R6 %'})
        axes[1].set_title('R6 Binding Heatmap', fontweight='bold')
        axes[1].set_xlabel('Ag Concentration (nM)')
        axes[1].set_ylabel('Sample ID')
    
    plt.tight_layout()
    plt.savefig('cell_fraction_heatmaps.png', dpi=300, bbox_inches='tight')
    print("Saved: cell_fraction_heatmaps.png")
    
    return fig

# Main execution
if __name__ == "__main__":
    # Check if results file exists
    import os
    
    if not os.path.exists('cell_fraction_results.csv'):
        print("Error: cell_fraction_results.csv not found!")
        print("Please run analyze_flow_cytometry_expression first.")
    else:
        # Load results
        try:
            results_df = pd.read_csv('cell_fraction_results.csv')
            
            if results_df.empty:
                print("Warning: Results file is empty!")
            else:
                print(f"Loaded {len(results_df)} samples for plotting")
                print(f"Unique concentrations: {sorted(results_df['Ag_Concentration_nM'].unique())}")
                print(f"Unique samples: {results_df['Sample_ID'].unique()[:5]}...")  # Show first 5
                
                # Create plots
                print("\nGenerating concentration curves...")
                fig1 = plot_concentration_curves(results_df)
                
                print("\nGenerating heatmaps...")
                fig2 = create_heatmap(results_df)
                
                plt.show()
                
                print("\nAll plots generated successfully!")
                
        except Exception as e:
            print(f"Error loading or processing data: {e}")
            import traceback
            traceback.print_exc()