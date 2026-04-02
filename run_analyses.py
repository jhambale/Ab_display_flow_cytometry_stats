#!/usr/bin/env python3
"""
Master script to run complete flow cytometry analysis pipeline.
This script imports and executes both cell fraction and MFI analyses,
along with their respective plotting functions.
"""

import argparse
import os
import sys
import traceback
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import utils_cb

# Import your analysis modules
try:
    # Cell fraction analysis
    from cell_fraction_analysis import analyze_flow_cytometry_expression

    # MFI analysis
    from mfi_analysis import analyze_flow_cytometry_mfi

    # Plotting functions
    from plot_cell_fractions import create_heatmap, plot_concentration_curves
    from plot_mfi import plot_mfi_analysis

except ImportError as e:
    print(f"Error importing modules: {e}")
    print("Make sure all script files are in the same directory!")
    sys.exit(1)


def create_output_directory(dir_prefix):
    """Create timestamped output directory for results."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = dir_prefix + f"flow_analysis_results_{timestamp}"
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def analyze_binding_significance(mfi_results, sd_threshold=2, output_dir=None):
    """
    Analyze binding significance by comparing non-zero concentration samples
    to the baseline (0 nM) using standard deviation thresholds.

    Parameters:
    -----------
    mfi_results : pd.DataFrame
        MFI analysis results containing 'Binding_X_MFI_Ratio' or 'X_MFI_Fold_Change'
    sd_threshold : int or list
        Number of standard deviations for significance (default=2, can be [2,3])
    output_dir : str
        Directory to save results

    Returns:
    --------
    dict : Analysis results including baseline stats and significant samples
    """

    print("\n" + "=" * 80)
    print("BINDING SIGNIFICANCE ANALYSIS")
    print("=" * 80)

    # Handle different column names based on your setup
    fold_change_col = None
    if "Binding_X_MFI_Ratio" in mfi_results.columns:
        fold_change_col = "Binding_X_MFI_Ratio"
    elif "X_MFI_Fold_Change" in mfi_results.columns:
        fold_change_col = "X_MFI_Fold_Change"
    else:
        print("ERROR: Could not find X MFI fold change column!")
        return None

    # Filter for 0 nM samples (baseline)
    baseline_data = mfi_results[mfi_results["Ag_Concentration_nM"] == 0].copy()

    if baseline_data.empty:
        print("ERROR: No 0 nM samples found for baseline!")
        return None

    # Calculate baseline statistics
    baseline_values = baseline_data[fold_change_col].dropna()

    if len(baseline_values) == 0:
        print("ERROR: No valid fold change values in baseline!")
        return None

    baseline_mean = baseline_values.mean()
    baseline_std = baseline_values.std()
    baseline_sem = baseline_values.sem()  # Standard error of mean
    n_baseline = len(baseline_values)

    print("\nBaseline Statistics (0 nM Ag):")
    print("-" * 60)
    print(f"  Number of samples: {n_baseline}")
    print(f"  Mean {fold_change_col}: {baseline_mean:.3f}")
    print(f"  Standard Deviation: {baseline_std:.3f}")
    print(f"  Standard Error: {baseline_sem:.3f}")
    print(f"  Range: [{baseline_values.min():.3f}, {baseline_values.max():.3f}]")

    # Handle multiple thresholds
    if not isinstance(sd_threshold, list):
        sd_thresholds = [sd_threshold]
    else:
        sd_thresholds = sd_threshold

    results = {
        "baseline_mean": baseline_mean,
        "baseline_std": baseline_std,
        "baseline_sem": baseline_sem,
        "baseline_n": n_baseline,
        "baseline_samples": baseline_data["Sample_ID"].tolist(),
        "analysis": {},
    }

    # Analyze non-zero concentrations
    for threshold in sd_thresholds:
        print("\n" + "-" * 60)
        print(f"Analysis with {threshold} Standard Deviation Threshold:")
        print("-" * 60)

        upper_bound = baseline_mean + (threshold * baseline_std)
        lower_bound = baseline_mean - (threshold * baseline_std)

        print(f"  Threshold bounds: [{lower_bound:.3f}, {upper_bound:.3f}]")

        threshold_results = {
            "upper_bound": upper_bound,
            "lower_bound": lower_bound,
            "above": [],
            "below": [],
            "within": [],
        }

        # Check each non-zero concentration
        for conc in [30, 120]:
            conc_data = mfi_results[mfi_results["Ag_Concentration_nM"] == conc].copy()

            if conc_data.empty:
                continue

            print(f"\n  {conc} nM Ag:")

            above_samples = []
            below_samples = []
            within_samples = []

            for _, row in conc_data.iterrows():
                sample_id = row["Sample_ID"]
                value = row[fold_change_col]

                if pd.isna(value):
                    continue

                # Calculate z-score
                z_score = (
                    (value - baseline_mean) / baseline_std if baseline_std > 0 else 0
                )

                sample_info = {
                    "sample_id": sample_id,
                    "concentration": conc,
                    "value": value,
                    "z_score": z_score,
                    "deviation_from_mean": value - baseline_mean,
                }

                if value > upper_bound:
                    above_samples.append(sample_info)
                    threshold_results["above"].append(sample_info)
                elif value < lower_bound:
                    below_samples.append(sample_info)
                    threshold_results["below"].append(sample_info)
                else:
                    within_samples.append(sample_info)
                    threshold_results["within"].append(sample_info)

            # Print summary for this concentration
            total_samples = (
                len(above_samples) + len(below_samples) + len(within_samples)
            )

            if total_samples > 0:
                print(f"    Total samples: {total_samples}")
                print(
                    f"    Above threshold (>{upper_bound:.3f}): {len(above_samples)} samples"
                )
                if above_samples:
                    for s in above_samples[:3]:  # Show first 3
                        print(
                            f"      • {s['sample_id']}: {s['value']:.3f} (z={s['z_score']:.2f})"
                        )
                    if len(above_samples) > 3:
                        print(f"      ... and {len(above_samples) - 3} more")

                print(
                    f"    Below threshold (<{lower_bound:.3f}): {len(below_samples)} samples"
                )
                if below_samples:
                    for s in below_samples[:3]:  # Show first 3
                        print(
                            f"      • {s['sample_id']}: {s['value']:.3f} (z={s['z_score']:.2f})"
                        )
                    if len(below_samples) > 3:
                        print(f"      ... and {len(below_samples) - 3} more")

                print(f"    Within normal range: {len(within_samples)} samples")

        results["analysis"][f"{threshold}sd"] = threshold_results

    # Create detailed report
    if output_dir:
        report_path = os.path.join(output_dir, "binding_significance_analysis.csv")

        # Prepare data for CSV
        report_data = []

        # Add baseline samples
        for _, row in baseline_data.iterrows():
            report_data.append(
                {
                    "Sample_ID": row["Sample_ID"],
                    "Ag_Concentration_nM": 0,
                    "X_MFI_Fold_Change": row[fold_change_col],
                    "Category": "Baseline",
                    "Z_Score": 0,
                    "Deviation_from_Baseline": 0,
                    "2SD_Significance": "Baseline",
                    "3SD_Significance": "Baseline",
                }
            )

        # Add test samples
        test_data = mfi_results[mfi_results["Ag_Concentration_nM"].isin([30, 120])]
        for _, row in test_data.iterrows():
            value = row[fold_change_col]
            if pd.isna(value):
                continue

            z_score = (value - baseline_mean) / baseline_std if baseline_std > 0 else 0
            deviation = value - baseline_mean

            # Determine significance at different thresholds
            sig_2sd = "Within Range"
            if value > baseline_mean + sd_threshold[0] * baseline_std:  # type: ignore
                sig_2sd = "Above 2SD"
            elif value < baseline_mean - sd_threshold[0] * baseline_std:  # type: ignore
                sig_2sd = "Below 2SD"

            sig_3sd = "Within Range"
            if value > baseline_mean + sd_threshold[1] * baseline_std:  # type: ignore
                sig_3sd = "Above 3SD"
            elif value < baseline_mean - sd_threshold[1] * baseline_std:  # type: ignore
                sig_3sd = "Below 3SD"

            report_data.append(
                {
                    "Sample_ID": row["Sample_ID"],
                    "Ag_Concentration_nM": row["Ag_Concentration_nM"],
                    "X_MFI_Fold_Change": value,
                    "Category": "Test",
                    "Z_Score": z_score,
                    "Deviation_from_Baseline": deviation,
                    "2SD_Significance": sig_2sd,
                    "3SD_Significance": sig_3sd,
                }
            )

        # Save to CSV
        report_df = pd.DataFrame(report_data)
        report_df.to_csv(report_path, index=False)
        print(f"\n✓ Detailed report saved to: {report_path}")

        # Create visualization
        create_significance_plot(
            mfi_results,
            fold_change_col,
            baseline_mean,
            baseline_std,
            output_dir,
            sd_threshold,
        )

    return results


def create_significance_plot(
    mfi_results, fold_change_col, baseline_mean, baseline_std, output_dir, sd_threshold
):
    """
    Create a visualization of binding significance analysis.
    """
    # import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Prepare data
    concentrations = [0, 30, 120]
    colors = {0: "#3498db", 30: "#e74c3c", 120: "#f39c12"}

    # Plot 1: Scatter plot with significance bands
    for conc in concentrations:
        conc_data = mfi_results[mfi_results["Ag_Concentration_nM"] == conc]
        values = conc_data[fold_change_col].dropna()

        if len(values) > 0:
            # Add jitter for better visibility
            x_positions = np.random.normal(conc, 2, len(values))
            ax1.scatter(
                x_positions,
                values,
                alpha=0.6,
                s=50,
                color=colors[conc],
                label=f"{conc} nM",
            )

    # Add significance bands
    x_range = [-10, 130]
    ax1.fill_between(
        x_range,
        baseline_mean - sd_threshold[1] * baseline_std,
        baseline_mean + sd_threshold[1] * baseline_std,
        alpha=0.1,
        color="gray",
        label=f"±{sd_threshold[1]} SD",
    )
    ax1.fill_between(
        x_range,
        baseline_mean - sd_threshold[0] * baseline_std,
        baseline_mean + sd_threshold[0] * baseline_std,
        alpha=0.2,
        color="gray",
        label=f"±{sd_threshold[0]} SD",
    )
    ax1.axhline(
        baseline_mean, color="black", linestyle="--", linewidth=1, label="Baseline Mean"
    )

    ax1.set_xlim(x_range)
    ax1.set_xticks(concentrations)
    ax1.set_xlabel("Ag Concentration (nM)")
    ax1.set_ylabel(fold_change_col.replace("_", " "))
    ax1.set_title("Binding Significance Analysis")
    ax1.legend(loc="best")
    ax1.grid(True, alpha=0.3)

    # Plot 2: Box plot comparison
    data_for_box = []
    labels_for_box = []

    for conc in concentrations:
        conc_data = mfi_results[mfi_results["Ag_Concentration_nM"] == conc]
        values = conc_data[fold_change_col].dropna()
        if len(values) > 0:
            data_for_box.append(values)
            labels_for_box.append(f"{conc} nM")

    bp = ax2.boxplot(data_for_box, tick_labels=labels_for_box, patch_artist=True)

    # Color the boxes
    for patch, conc in zip(bp["boxes"], concentrations):
        patch.set_facecolor(colors[conc])
        patch.set_alpha(0.6)

    # Add significance lines
    ax2.axhline(
        baseline_mean + sd_threshold[0] * baseline_std,
        color="orange",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
    )
    ax2.axhline(
        baseline_mean - sd_threshold[0] * baseline_std,
        color="orange",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
    )
    ax2.axhline(baseline_mean, color="black", linestyle="--", linewidth=1, alpha=0.7)

    ax2.set_ylabel(fold_change_col.replace("_", " "))
    ax2.set_title("Distribution by Concentration")
    ax2.grid(True, alpha=0.3)

    plt.suptitle(
        "Binding Significance Analysis - X MFI Fold Change",
        fontsize=14,
        fontweight="bold",
    )
    plt.tight_layout()

    if output_dir:
        plt.savefig(
            os.path.join(output_dir, "binding_significance_plot.png"),
            dpi=300,
            bbox_inches="tight",
        )
        print("✓ Significance plot saved")

    plt.close()


def run_complete_analysis(csv_file, sd_threshold):
    """
    Run complete flow cytometry analysis pipeline.

    Parameters:
    -----------
    csv_file : str
        Path to the CSV file containing flow cytometry statistics
    """

    print("=" * 80)
    print("COMPLETE FLOW CYTOMETRY ANALYSIS PIPELINE")
    print("=" * 80)
    print(f"\nInput file: {csv_file}")

    # Verify file exists
    if not os.path.exists(csv_file):
        print(f"ERROR: File '{csv_file}' not found!")
        return

    # Create output directory
    output_dir = create_output_directory(dir_prefix)
    print(f"Output directory: {output_dir}")
    print("-" * 80)

    # =============================================
    # PART 1: CELL FRACTION ANALYSIS
    # =============================================
    print("\n" + "=" * 80)
    print("PART 1: CELL FRACTION ANALYSIS")
    print("=" * 80)

    try:
        # Run cell fraction analysis
        fraction_results = analyze_flow_cytometry_expression(csv_file)

        # Save results
        fraction_output = os.path.join(output_dir, "cell_fraction_results.csv")
        fraction_results.to_csv(fraction_output, index=False)
        print(f"✓ Cell fraction results saved to: {fraction_output}")

        # Display summary
        print("\nCell Fraction Summary by Ag Concentration:")
        print("-" * 60)
        for conc in [0, 30, 120]:
            conc_data = fraction_results[
                fraction_results["Ag_Concentration_nM"] == conc
            ]
            if not conc_data.empty:
                expr_mean = conc_data["Expression_Fraction_%"].mean()
                expr_std = conc_data["Expression_Fraction_%"].std()
                r6_mean = conc_data["R6_Percent_of_Expressing"].mean()
                r6_std = conc_data["R6_Percent_of_Expressing"].std()

                print(f"\n{conc} nM Ag:")
                print(f"  Expression: {expr_mean:.2f} ± {expr_std:.2f}%")
                print(f"  R6 Binding: {r6_mean:.2f} ± {r6_std:.2f}%")

        # Create plots for cell fractions
        print("\nGenerating cell fraction plots...")

        # Concentration curves
        fig1 = plot_concentration_curves(fraction_results, save_plots=False)
        fig1.savefig(
            os.path.join(output_dir, "cell_fraction_curves.png"),
            dpi=300,
            bbox_inches="tight",
        )
        print("✓ Concentration curves saved")

        # Heatmaps
        fig2 = create_heatmap(fraction_results)
        fig2.savefig(
            os.path.join(output_dir, "cell_fraction_heatmaps.png"),
            dpi=300,
            bbox_inches="tight",
        )
        print("✓ Heatmaps saved")

        plt.close("all")  # Close all figures to save memory

    except Exception as e:
        print(f"ERROR in cell fraction analysis: {e}")
        fraction_results = None

    # =============================================
    # PART 2: MFI ANALYSIS
    # =============================================
    print("\n" + "=" * 80)
    print("PART 2: MFI ANALYSIS")
    print("=" * 80)

    try:
        # Run MFI analysis
        mfi_results = analyze_flow_cytometry_mfi(csv_file)

        # Save results
        mfi_output = os.path.join(output_dir, "mfi_results.csv")
        mfi_results.to_csv(mfi_output, index=False)
        print(f"✓ MFI results saved to: {mfi_output}")

        # Display summary
        print("\nMFI Summary by Ag Concentration:")
        print("-" * 60)
        for conc in [0, 30, 120]:
            conc_data = mfi_results[mfi_results["Ag_Concentration_nM"] == conc]
            if not conc_data.empty:
                expr_mfi = conc_data["Expression_Y_MFI"].mean()
                fold_change = conc_data["X_MFI_Fold_Change"].mean()
                r6_enrichment = conc_data["R6_X_Enrichment"].mean()

                print(f"\n{conc} nM Ag:")
                print(f"  Expression MFI: {expr_mfi:.0f}")
                print(f"  Fold Change: {fold_change:.2f}x")
                print(f"  R6 Enrichment: {r6_enrichment:.2f}x")

        # Create plots for MFI
        print("\nGenerating MFI plots...")
        fig3 = plot_mfi_analysis(mfi_results)
        fig3.savefig(
            os.path.join(output_dir, "mfi_analysis_plots.png"),
            dpi=300,
            bbox_inches="tight",
        )
        print("✓ MFI plots saved")

        plt.close("all")  # Close all figures to save memory

    except Exception as e:
        print(f"ERROR in MFI analysis: {e}")
        mfi_results = None

    # =============================================
    # PART 3: COMBINED ANALYSIS
    # =============================================
    print("\n" + "=" * 80)
    print("PART 3: COMBINED ANALYSIS")
    print("=" * 80)

    if fraction_results is not None and mfi_results is not None:
        try:
            # Merge results
            combined_results = fraction_results.merge(
                mfi_results,
                on=["Sample_ID", "Sample_Full", "Ag_Concentration_nM"],
                suffixes=("_fraction", "_mfi"),
            )

            # Save combined results
            combined_output = os.path.join(output_dir, "combined_analysis.csv")
            combined_results.to_csv(combined_output, index=False)
            print(f"✓ Combined results saved to: {combined_output}")

            # Create comparison plot
            create_comparison_plot(fraction_results, mfi_results, output_dir)

        except Exception as e:
            print(f"ERROR in combined analysis: {e}")

    # =============================================
    # PART 4: STATISTICAL ANALYSIS
    # =============================================
    print("\n" + "=" * 80)
    print("PART 4: STATISTICAL SUMMARY")
    print("=" * 80)

    # Create summary statistics file
    summary_stats = create_summary_statistics(fraction_results, mfi_results)
    if summary_stats is not None:
        summary_output = os.path.join(output_dir, "summary_statistics.csv")
        summary_stats.to_csv(summary_output)
        print(f"✓ Summary statistics saved to: {summary_output}")

    # Create report
    create_analysis_report(output_dir, csv_file, fraction_results, mfi_results)

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"\n📁 All results saved in: {output_dir}")
    print("\nFiles generated:")
    print("  • cell_fraction_results.csv - Cell population percentages")
    print("  • mfi_results.csv - Mean fluorescence intensity analysis")
    print("  • combined_analysis.csv - Merged dataset")
    print("  • summary_statistics.csv - Statistical summaries")
    print("  • analysis_report.txt - Complete analysis report")
    print("  • Various .png files - Visualization plots")

    # =============================================
    # PART 5: BINDING SIGNIFICANCE ANALYSIS
    # =============================================
    print("\n" + "=" * 80)
    print("PART 5: BINDING SIGNIFICANCE ANALYSIS")
    print("=" * 80)

    if mfi_results is not None:
        try:
            # Analyze binding significance using 2 and 3 SD thresholds
            significance_results = analyze_binding_significance(
                mfi_results,
                sd_threshold=[int(sd_threshold), int(sd_threshold) + 1],
                output_dir=output_dir,
            )

            if significance_results:
                # Summary
                print("\n" + "-" * 60)
                print("SIGNIFICANCE SUMMARY:")
                print("-" * 60)

                for threshold in [int(sd_threshold), int(sd_threshold) + 1]:
                    analysis = significance_results["analysis"].get(
                        f"{threshold}sd", {}
                    )
                    n_above = len(analysis.get("above", []))
                    n_below = len(analysis.get("below", []))
                    n_within = len(analysis.get("within", []))
                    total = n_above + n_below + n_within

                    if total > 0:
                        print(f"\n{threshold} SD Threshold:")
                        print(
                            f"  Samples above: {n_above} ({100 * n_above / total:.1f}%)"
                        )
                        print(
                            f"  Samples below: {n_below} ({100 * n_below / total:.1f}%)"
                        )
                        print(
                            f"  Samples within: {n_within} ({100 * n_within / total:.1f}%)"
                        )

        except Exception as e:
            print(f"ERROR in binding significance analysis: {e}")
            traceback.print_exc()

    return output_dir, fraction_results, mfi_results


def create_comparison_plot(fraction_results, mfi_results, output_dir):
    """Create a comparison plot between fraction and MFI analyses."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    concentrations = [0, 30, 120]

    # Get unique samples
    fraction_results["Sample_Base"] = fraction_results["Sample_ID"].str.extract(
        r"(y\d+)", expand=False
    )
    unique_samples = fraction_results["Sample_Base"].dropna().unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_samples)))  # type: ignore

    # Store line objects for legend
    legend_lines = []
    legend_labels = []

    # Plot 1: Expression Fraction
    ax1 = axes[0, 0]
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = fraction_results[fraction_results["Sample_Base"] == sample]
        x_vals, y_vals = [], []
        for conc in concentrations:
            conc_data = sample_data[sample_data["Ag_Concentration_nM"] == conc]
            if not conc_data.empty:
                x_vals.append(conc)
                y_vals.append(conc_data["Expression_Fraction_%"].iloc[0])
        if x_vals:
            (line,) = ax1.plot(x_vals, y_vals, marker="o", color=colors[i], alpha=0.7)
            # Store first plot's lines for legend
            if len(legend_lines) < len(unique_samples):
                legend_lines.append(line)
                legend_labels.append(sample)

    ax1.set_title("Expression Fraction (%)")
    ax1.set_xlabel("Ag (nM)")
    ax1.set_ylabel("% Expressing Cells")
    ax1.grid(True, alpha=0.3)

    # Plot 2: Expression MFI
    ax2 = axes[0, 1]
    mfi_results["Sample_Base"] = mfi_results["Sample_ID"].str.extract(
        r"(y\d+)", expand=False
    )
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = mfi_results[mfi_results["Sample_Base"] == sample]
        x_vals, y_vals = [], []
        for conc in concentrations:
            conc_data = sample_data[sample_data["Ag_Concentration_nM"] == conc]
            if not conc_data.empty and pd.notna(conc_data["Expression_Y_MFI"].iloc[0]):
                x_vals.append(conc)
                y_vals.append(conc_data["Expression_Y_MFI"].iloc[0])
        if x_vals:
            ax2.plot(x_vals, y_vals, marker="s", color=colors[i], alpha=0.7)
    ax2.set_title("Expression MFI")
    ax2.set_xlabel("Ag (nM)")
    ax2.set_ylabel("MFI (AU)")
    ax2.grid(True, alpha=0.3)

    # Plot 3: R6 Binding %
    ax3 = axes[1, 0]
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = fraction_results[fraction_results["Sample_Base"] == sample]
        x_vals, y_vals = [], []
        for conc in concentrations:
            conc_data = sample_data[sample_data["Ag_Concentration_nM"] == conc]
            if not conc_data.empty and pd.notna(
                conc_data["R6_Percent_of_Expressing"].iloc[0]
            ):
                x_vals.append(conc)
                y_vals.append(conc_data["R6_Percent_of_Expressing"].iloc[0])
        if x_vals:
            ax3.plot(x_vals, y_vals, marker="^", color=colors[i], alpha=0.7)
    ax3.set_title("R6 Binding (% of Expressing)")
    ax3.set_xlabel("Ag (nM)")
    ax3.set_ylabel("R6 %")
    ax3.grid(True, alpha=0.3)

    # Plot 4: X MFI Enrichment
    ax4 = axes[1, 1]
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = mfi_results[mfi_results["Sample_Base"] == sample]
        x_vals, y_vals = [], []
        for conc in concentrations:
            conc_data = sample_data[sample_data["Ag_Concentration_nM"] == conc]
            if not conc_data.empty and pd.notna(conc_data["X_MFI_Fold_Change"].iloc[0]):
                x_vals.append(conc)
                y_vals.append(conc_data["X_MFI_Fold_Change"].iloc[0])
        if x_vals:
            ax4.plot(x_vals, y_vals, marker="d", color=colors[i], alpha=0.7)
    ax4.set_title("X MFI Fold_Change (Expressing/Non-Expressing)")
    ax4.set_xlabel("Ag (nM)")
    ax4.set_ylabel("Fold Change")
    ax4.grid(True, alpha=0.3)

    plt.suptitle(
        "Comparison: Cell Fractions vs MFI Analysis", fontsize=14, fontweight="bold"
    )

    # Add super legend at the bottom of the figure
    fig.legend(
        legend_lines,
        legend_labels,
        loc="center",  # Center the legend
        bbox_to_anchor=(0.5, -0.05),  # Position below the plots
        ncol=min(4, len(legend_labels)),  # Use up to 4 columns
        frameon=True,
        fancybox=True,
        shadow=True,
        title="Sample IDs",
        title_fontsize=10,
        fontsize=9,
    )

    plt.tight_layout(rect=[0, 0.05, 1, 0.98])

    plt.savefig(
        os.path.join(output_dir, "comparison_plot.png"), dpi=300, bbox_inches="tight"
    )
    plt.close()
    print("✓ Comparison plot saved")


def create_summary_statistics(fraction_results, mfi_results):
    """Create summary statistics combining both analyses."""

    if fraction_results is None or mfi_results is None:
        return None

    summary = []

    for conc in [0, 30, 120]:
        # Fraction statistics
        frac_data = fraction_results[fraction_results["Ag_Concentration_nM"] == conc]

        # MFI statistics
        mfi_data = mfi_results[mfi_results["Ag_Concentration_nM"] == conc]

        if not frac_data.empty and not mfi_data.empty:
            summary.append(
                {
                    "Ag_Concentration_nM": conc,
                    "n_samples": len(frac_data),
                    # Cell fractions
                    "Expression_Fraction_Mean": frac_data[
                        "Expression_Fraction_%"
                    ].mean(),
                    "Expression_Fraction_SD": frac_data["Expression_Fraction_%"].std(),
                    "R6_Percent_Mean": frac_data["R6_Percent_of_Expressing"].mean(),
                    "R6_Percent_SD": frac_data["R6_Percent_of_Expressing"].std(),
                    # MFI values
                    "Expression_MFI_Mean": mfi_data["Expression_Y_MFI"].mean(),
                    "Expression_MFI_SD": mfi_data["Expression_Y_MFI"].std(),
                    "MFI_Fold_Change_Mean": mfi_data["X_MFI_Fold_Change"].mean(),
                    "MFI_Fold_Change_SD": mfi_data["X_MFI_Fold_Change"].std(),
                    "R6_Enrichment_Mean": mfi_data["R6_X_Enrichment"].mean(),
                    "R6_Enrichment_SD": mfi_data["R6_X_Enrichment"].std(),
                }
            )

    return pd.DataFrame(summary)


def create_analysis_report(output_dir, csv_file, fraction_results, mfi_results):
    """Create a text report summarizing the analysis."""

    report_path = os.path.join(output_dir, "analysis_report.txt")

    with open(report_path, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("FLOW CYTOMETRY ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input File: {csv_file}\n")
        f.write(f"Output Directory: {output_dir}\n\n")

        f.write("-" * 80 + "\n")
        f.write("ANALYSIS SUMMARY\n")
        f.write("-" * 80 + "\n\n")

        if fraction_results is not None:
            f.write("Cell Fraction Analysis:\n")
            f.write(f"  • Total samples analyzed: {len(fraction_results)}\n")
            f.write(
                f"  • Unique sample IDs: {fraction_results['Sample_ID'].nunique()}\n"
            )
            f.write("  • Ag concentrations tested: 0, 30, 120 nM\n\n")

        if mfi_results is not None:
            f.write("MFI Analysis:\n")
            f.write(f"  • Total samples analyzed: {len(mfi_results)}\n")
            f.write("  • Metrics calculated: Fold change, enrichment, combined MFI\n\n")

        f.write("-" * 80 + "\n")
        f.write("KEY FINDINGS\n")
        f.write("-" * 80 + "\n\n")

        if fraction_results is not None:
            for conc in [0, 30, 120]:
                conc_data = fraction_results[
                    fraction_results["Ag_Concentration_nM"] == conc
                ]
                if not conc_data.empty:
                    f.write(f"At {conc} nM Ag:\n")
                    f.write(
                        f"  Expression fraction: {conc_data['Expression_Fraction_%'].mean():.2f}%\n"
                    )
                    f.write(
                        f"  R6 binding: {conc_data['R6_Percent_of_Expressing'].mean():.2f}%\n\n"
                    )

        f.write("-" * 80 + "\n")
        f.write("FILES GENERATED\n")
        f.write("-" * 80 + "\n\n")

        # List all files in output directory
        for file in os.listdir(output_dir):
            f.write(f"  • {file}\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")

    print("✓ Analysis report saved")


# Main execution
if __name__ == "__main__":
    # Specify your CSV file path
    # csv_file = "2026-01-14_Agscfv_stats.csv"
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", nargs="+", help="input mutation analysis csv files")
    parser.add_argument(
        "-t",
        default=2,
        help="threshold of standard devs for determining significant binding above bg",
    )
    # parser.add_argument(-'m', help='type of metric for analysis (Count, Mean, or Median)')

    args = parser.parse_args()

    csv_file = args.i
    # print(csv_file[0])
    dir_prefix = "/".join(csv_file[0].split("/")[:-1])
    print(f"output directory: {dir_prefix}")

    sd_threshold = args.t

    # Check if file exists
    if not os.path.exists(csv_file[0]):
        print(f"ERROR: Cannot find '{csv_file}'")
        # print("Please make sure the file is in the current directory.")
        sys.exit(1)

    # Run the complete analysis
    try:
        output_dir, fraction_results, mfi_results = run_complete_analysis(
            csv_file[0], sd_threshold
        )

        # Optional: Open the output directory
        if sys.platform == "win32":
            os.startfile(output_dir)
        elif sys.platform == "darwin":
            os.system(f"open {output_dir}")
        else:  # linux
            os.system(f"xdg-open {output_dir}")

    except Exception as e:
        print(f"\nERROR: Analysis failed with error: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
