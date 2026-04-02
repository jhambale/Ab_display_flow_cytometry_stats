#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

# import os


def analyze_flow_cytometry_mfi(csv_file):
    """
    Complete MFI-based flow cytometry analysis function.
    Analyzes Mean Fluorescence Intensity (MFI) data from flow cytometry.

    Parameters:
    -----------
    csv_file : str
        Path to the CSV file containing flow cytometry statistics

    Returns:
    --------
    pd.DataFrame : MFI analysis results
    """

    # Read the CSV file
    df = pd.read_csv(csv_file)

    # Get unique samples
    unique_samples = df["Sample"].unique()
    print(f"Analyzing MFI for {len(unique_samples)} samples...")

    # Initialize results list
    results = []

    for sample in unique_samples:
        # Filter data for current sample
        sample_data = df[df["Sample"] == sample]

        # Initialize MFI storage
        mfi_data = {}

        # Gates of interest
        quadrant_gates = ["+/+(1)", "+/-(1)", "-/+(1)", "-/-(1)"]

        # Extract MFI values for each gate
        for _, row in sample_data.iterrows():
            gate = row["Gate"]

            if gate in quadrant_gates or gate == "R6":
                mfi_data[f"{gate}_X_MFI"] = row["X Mean"]
                mfi_data[f"{gate}_Y_MFI"] = row["Y Mean"]
                mfi_data[f"{gate}_Count"] = row["Count"]

        # Calculate expressing population MFI (geometric mean of +/+ and -/+)
        pos_pos_x = mfi_data.get("+/+(1)_X_MFI", 0)
        pos_pos_y = mfi_data.get("+/+(1)_Y_MFI", 0)
        neg_pos_x = mfi_data.get("-/+(1)_X_MFI", 0)
        neg_pos_y = mfi_data.get("-/+(1)_Y_MFI", 0)

        # Geometric mean for expressing cells
        if pos_pos_x > 0 and neg_pos_x > 0:
            expressing_x_mfi_geomean = np.sqrt(pos_pos_x * neg_pos_x)
        else:
            expressing_x_mfi_geomean = np.mean([pos_pos_x, neg_pos_x])

        if pos_pos_y > 0 and neg_pos_y > 0:
            expressing_y_mfi_geomean = np.sqrt(pos_pos_y * neg_pos_y)
        else:
            expressing_y_mfi_geomean = np.mean([pos_pos_y, neg_pos_y])

        # Calculate non-expressing population MFI
        neg_neg_x = mfi_data.get("-/-(1)_X_MFI", 0)
        neg_neg_y = mfi_data.get("-/-(1)_Y_MFI", 0)
        pos_neg_x = mfi_data.get("+/-(1)_X_MFI", 0)
        pos_neg_y = mfi_data.get("+/-(1)_Y_MFI", 0)

        if neg_neg_x > 0 and pos_neg_x > 0:
            nonexpressing_x_mfi_geomean = np.sqrt(neg_neg_x * pos_neg_x)
        else:
            nonexpressing_x_mfi_geomean = np.mean([neg_neg_x, pos_neg_x])

        if neg_neg_y > 0 and pos_neg_y > 0:
            nonexpressing_y_mfi_geomean = np.sqrt(neg_neg_y * pos_neg_y)
        else:
            nonexpressing_y_mfi_geomean = np.mean([neg_neg_y, pos_neg_y])

        # Calculate MFI fold change
        if nonexpressing_x_mfi_geomean > 0:
            x_mfi_fold_change = expressing_x_mfi_geomean / nonexpressing_x_mfi_geomean
        else:
            x_mfi_fold_change = None

        if nonexpressing_y_mfi_geomean > 0:
            y_mfi_fold_change = expressing_y_mfi_geomean / nonexpressing_y_mfi_geomean
        else:
            y_mfi_fold_change = None

        # R6 MFI analysis
        r6_x_mfi = mfi_data.get("R6_X_MFI", 0)
        r6_y_mfi = mfi_data.get("R6_Y_MFI", 0)

        # Convert 0 to a small non-zero value for log scale compatibility
        if r6_x_mfi == 0:
            r6_x_mfi = 0.1  # Use small value instead of 0
        if r6_y_mfi == 0:
            r6_y_mfi = 0.1

        # Calculate R6 MFI enrichment
        if expressing_x_mfi_geomean > 0 and r6_x_mfi > 0:
            r6_x_enrichment = r6_x_mfi / expressing_x_mfi_geomean
        else:
            r6_x_enrichment = None

        # if expressing_y_mfi_geomean > 0 and r6_y_mfi > 0:
        #     r6_y_enrichment = r6_y_mfi / expressing_y_mfi_geomean
        # else:
        #     r6_y_enrichment = None

        # # Combined MFI score
        # expressing_combined_mfi = expressing_x_mfi_geomean * expressing_y_mfi_geomean
        # r6_combined_mfi = r6_x_mfi * r6_y_mfi if r6_x_mfi and r6_y_mfi else 0

        # Extract Ag concentration
        concentration = None
        if "120nMAg" in sample:
            concentration = 120
        elif "30nMAg" in sample:
            concentration = 30
        elif "0nMAg" in sample:
            concentration = 0

        # Extract sample ID
        sample_id = sample.split("_")[0] if "_" in sample else sample

        # Build results dictionary
        result_dict = {
            "Sample_ID": sample_id,
            "Sample_Full": sample,
            "Ag_Concentration_nM": concentration,
            # Expression metrics (Y-axis MFI)
            "Expression_Y_MFI": expressing_y_mfi_geomean,  # PRIMARY EXPRESSION METRIC
            "NonExpressing_Y_MFI": nonexpressing_y_mfi_geomean,
            "Expression_Y_Fold_Change": y_mfi_fold_change,
            # Binding metrics (X-axis MFI)
            "R6_X_MFI": r6_x_mfi,  # PRIMARY BINDING METRIC
            "PosPos_X_MFI": mfi_data.get("+/+(1)_X_MFI", 0),  # +/+ quadrant X MFI
            "R6_X_Enrichment": r6_x_enrichment,
            # Keep other channel data for reference
            "Expressing_X_MFI": expressing_x_mfi_geomean,
            "NonExpressing_X_MFI": nonexpressing_x_mfi_geomean,
            "X_MFI_Fold_Change": x_mfi_fold_change,
        }

        # Add individual gate MFI values
        result_dict.update(mfi_data)

        results.append(result_dict)

    # Create DataFrame
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(["Sample_ID", "Ag_Concentration_nM"])

    return results_df


# Main execution
if __name__ == "__main__":
    csv_file = "2026-01-14_Agscfv_stats.csv"

    print("=" * 80)
    print("FLOW CYTOMETRY MFI ANALYSIS")
    print("=" * 80)

    # Analyze MFI data
    mfi_results = analyze_flow_cytometry_mfi(csv_file)

    # Save results
    output_file = "mfi_results.csv"
    mfi_results.to_csv(output_file, index=False)

    print(f"\nMFI analysis complete! Results saved to: {output_file}")

    # Display sample results
    display_cols = [
        "Sample_ID",
        "Ag_Concentration_nM",
        "Expression_Y_MFI",
        "Expression_Y_Fold_Change",
        "R6_X_MFI",
        "R6_X_Enrichment",
    ]

    print("\nSample MFI Results (first 10):")
    print("-" * 80)
    print(mfi_results[display_cols].head(10).to_string(index=False))

    # Summary by concentration
    print("\n" + "-" * 80)
    print("MFI Summary by Ag Concentration:")
    print("-" * 80)

    for conc in [0, 30, 120]:
        conc_data = mfi_results[mfi_results["Ag_Concentration_nM"] == conc]
        if not conc_data.empty:
            expr_mfi_mean = conc_data["Expression_Y_MFI"].mean()
            fold_change_mean = conc_data["X_MFI_Fold_Change"].mean()
            r6_enrichment_mean = conc_data["R6_X_Enrichment"].mean()

            print(f"\n{conc} nM Ag:")
            print(f"  Expression MFI: {expr_mfi_mean:.0f}")
            print(f"  Fold Change: {fold_change_mean:.2f}x")
            print(f"  R6 Enrichment: {r6_enrichment_mean:.2f}x")
