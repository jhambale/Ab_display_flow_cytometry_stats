#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os

def analyze_flow_cytometry_expression(csv_file):
    """
    Analyzes flow cytometry data to calculate expression fractions and R6 percentages
    using cell counts.
    
    Parameters:
    -----------
    csv_file : str
        Path to the CSV file containing flow cytometry statistics
    
    Returns:
    --------
    pd.DataFrame : Results with expression fractions and R6 percentages
    """
    
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Get unique samples
    unique_samples = df['Sample'].unique()
    print(f"Found {len(unique_samples)} unique samples")
    print("-" * 80)
    
    # Initialize results list
    results = []
    
    for sample in unique_samples:
        # Filter data for current sample
        sample_data = df[df['Sample'] == sample]
        
        # Initialize all variables at the start
        pos_pos_count = None
        pos_neg_count = None
        neg_pos_count = None
        neg_neg_count = None
        r6_count = None
        total_quadrant_cells = None
        expressing_cells = None
        expression_fraction = None
        r6_percent_of_expressing = None
        concentration = None
        
        # Extract counts from relevant gates
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
        
        # Calculate total quadrant cells and expression fraction only if all counts exist
        if all(x is not None for x in [pos_pos_count, pos_neg_count, neg_pos_count, neg_neg_count]):
            total_quadrant_cells = pos_pos_count + pos_neg_count + neg_pos_count + neg_neg_count
            expressing_cells = pos_pos_count + neg_pos_count
            
            if total_quadrant_cells > 0:
                expression_fraction = (expressing_cells / total_quadrant_cells) * 100
            else:
                expression_fraction = 0
            
            # Calculate R6 as percentage of expressing cells
            if expressing_cells > 0 and r6_count is not None:
                r6_percent_of_expressing = (r6_count / expressing_cells) * 100
            else:
                r6_percent_of_expressing = None
        
        # Extract Ag concentration from sample name
        if '120nMAg' in sample:
            concentration = 120
        elif '30nMAg' in sample:
            concentration = 30
        elif '0nMAg' in sample:
            concentration = 0
        else:
            concentration = None
        
        # Extract sample ID
        sample_id = sample.split('_')[0] if '_' in sample else sample
        
        # Append results
        results.append({
            'Sample_ID': sample_id,
            'Sample_Full': sample,
            'Ag_Concentration_nM': concentration,
            '+/+(1)_Count': pos_pos_count,
            '+/-(1)_Count': pos_neg_count,
            '-/+(1)_Count': neg_pos_count,
            '-/-(1)_Count': neg_neg_count,
            'Total_Quadrant_Cells': total_quadrant_cells,
            'Expressing_Cells': expressing_cells,
            'Expression_Fraction_%': expression_fraction,
            'R6_Count': r6_count,
            'R6_Percent_of_Expressing': r6_percent_of_expressing
        })
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results)
    
    # Sort by sample ID and concentration
    results_df = results_df.sort_values(['Sample_ID', 'Ag_Concentration_nM'])
    
    return results_df

# Main execution
if __name__ == "__main__":
    # Specify your CSV file path
    csv_file = "2026-01-14_Agscfv_stats.csv"
    
    print("=" * 80)
    print("FLOW CYTOMETRY CELL FRACTION ANALYSIS")
    print("=" * 80)
    
    # Analyze data
    results = analyze_flow_cytometry_expression(csv_file)
    
    # Save results
    output_file = "cell_fraction_results.csv"
    results.to_csv(output_file, index=False)
    
    # Display summary
    print("\nAnalysis Complete!")
    print(f"Results saved to: {output_file}")
    print("\nSample Results (first 10 rows):")
    print("-" * 80)
    print(results[['Sample_ID', 'Ag_Concentration_nM', 'Expression_Fraction_%', 
                  'R6_Percent_of_Expressing']].head(10).to_string(index=False))
    
    # Display summary statistics by concentration
    print("\n" + "-" * 80)
    print("Summary by Ag Concentration:")
    print("-" * 80)
    
    for conc in [0, 30, 120]:
        conc_data = results[results['Ag_Concentration_nM'] == conc]
        if not conc_data.empty:
            expr_mean = conc_data['Expression_Fraction_%'].mean()
            expr_std = conc_data['Expression_Fraction_%'].std()
            r6_mean = conc_data['R6_Percent_of_Expressing'].mean()
            r6_std = conc_data['R6_Percent_of_Expressing'].std()
            
            print(f"\n{conc} nM Ag:")
            print(f"  Expression: {expr_mean:.2f} ± {expr_std:.2f}%")
            print(f"  R6 Binding: {r6_mean:.2f} ± {r6_std:.2f}%")