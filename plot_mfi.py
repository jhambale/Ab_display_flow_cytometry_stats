#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_mfi_analysis(results_df):
    """
    Create comprehensive MFI visualization plots.
    """
    # Set style
    plt.style.use("seaborn-v0_8-whitegrid")

    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))

    # Get unique samples
    results_df["Sample_Base"] = results_df["Sample_ID"].str.extract(
        r"(y\d+)", expand=False
    )
    unique_samples = results_df["Sample_Base"].dropna().unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_samples)))  # type: ignore
    concentrations = [0, 30, 120]

    # 1. Expression Combined MFI by concentration
    ax1 = plt.subplot(2, 3, 1)
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = results_df[results_df["Sample_Base"] == sample]
        x_vals = []
        y_vals = []
        for conc in concentrations:
            conc_data = sample_data[sample_data["Ag_Concentration_nM"] == conc]
            if not conc_data.empty and pd.notna(conc_data["Expression_Y_MFI"].iloc[0]):
                x_vals.append(conc)
                y_vals.append(conc_data["Expression_Y_MFI"].iloc[0])
        if len(x_vals) > 0:
            ax1.plot(
                x_vals,
                y_vals,
                marker="o",
                label=sample,
                color=colors[i],
                linewidth=2,
                markersize=6,
                alpha=0.7,
            )
    ax1.set_xlabel("Ag (nM)")
    ax1.set_ylabel("Combined MFI")
    ax1.set_title("Expression Combined MFI")
    ax1.set_xticks(concentrations)
    ax1.grid(True, alpha=0.3)

    # 2. MFI Fold Change by concentration
    ax2 = plt.subplot(2, 3, 2)
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = results_df[results_df["Sample_Base"] == sample]
        x_vals = []
        y_vals = []
        for conc in concentrations:
            conc_data = sample_data[sample_data["Ag_Concentration_nM"] == conc]
            if not conc_data.empty and pd.notna(conc_data["X_MFI_Fold_Change"].iloc[0]):
                x_vals.append(conc)
                y_vals.append(conc_data["X_MFI_Fold_Change"].iloc[0])
        if len(x_vals) > 0:
            ax2.plot(
                x_vals,
                y_vals,
                marker="s",
                label=sample,
                color=colors[i],
                linewidth=2,
                markersize=6,
                alpha=0.7,
            )
    ax2.set_xlabel("Ag (nM)")
    ax2.set_ylabel("Fold Change")
    ax2.set_title("X-Channel MFI Fold Change")
    ax2.set_xticks(concentrations)
    ax2.grid(True, alpha=0.3)

    # 3. R6 Combined MFI
    ax3 = plt.subplot(2, 3, 3)
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = results_df[results_df["Sample_Base"] == sample]
        x_vals = []
        y_vals = []
        for conc in concentrations:
            conc_data = sample_data[sample_data["Ag_Concentration_nM"] == conc]
            if not conc_data.empty:
                # Handle None/NaN values
                r6_val = conc_data["R6_X_MFI"].iloc[0]
                if pd.notna(r6_val):
                    x_vals.append(conc)
                    y_vals.append(max(r6_val, 0.1))

        if len(x_vals) > 0:
            ax3.plot(
                x_vals,
                y_vals,
                marker="^",
                label=sample,
                color=colors[i],
                linewidth=2,
                markersize=6,
                alpha=0.7,
            )
    # Use symlog scale for y-axis to handle values from 0 to high values
    # ax3.set_yscale('symlog', linthresh=10)  # linthresh sets the linear threshold
    ax3.set_xlabel("Ag (nM)")
    ax3.set_ylabel("R6 X MFI")
    ax3.set_title("R6 X MFI")
    ax3.set_xticks(concentrations)
    ax3.grid(True, alpha=0.3)

    # 4. R6 Enrichment
    ax4 = plt.subplot(2, 3, 4)
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = results_df[results_df["Sample_Base"] == sample]
        x_vals = []
        y_vals = []
        for conc in concentrations:
            conc_data = sample_data[sample_data["Ag_Concentration_nM"] == conc]
            if not conc_data.empty and pd.notna(conc_data["R6_X_Enrichment"].iloc[0]):
                x_vals.append(conc)
                y_vals.append(conc_data["R6_X_Enrichment"].iloc[0])
        if len(x_vals) > 0:
            ax4.plot(
                x_vals,
                y_vals,
                marker="d",
                label=sample,
                color=colors[i],
                linewidth=2,
                markersize=6,
                alpha=0.7,
            )
    ax4.set_xlabel("Ag (nM)")
    ax4.set_ylabel("R6 Enrichment (fold)")
    ax4.set_title("R6 X-Channel Enrichment")
    ax4.set_xticks(concentrations)
    ax4.grid(True, alpha=0.3)
    ax4.axhline(y=1, color="k", linestyle="--", alpha=0.3)

    # 5. Y-Channel MFI Fold Change
    ax5 = plt.subplot(2, 3, 5)
    for i, sample in enumerate(sorted(unique_samples)):
        sample_data = results_df[results_df["Sample_Base"] == sample]
        x_vals = []
        y_vals = []
        for conc in concentrations:
            conc_data = sample_data[sample_data["Ag_Concentration_nM"] == conc]
            if not conc_data.empty and pd.notna(
                conc_data["Expression_Y_Fold_Change"].iloc[0]
            ):
                x_vals.append(conc)
                y_vals.append(conc_data["Expression_Y_Fold_Change"].iloc[0])
        if len(x_vals) > 0:
            ax5.plot(
                x_vals,
                y_vals,
                marker="p",
                label=sample,
                color=colors[i],
                linewidth=2,
                markersize=6,
                alpha=0.7,
            )
    ax5.set_xlabel("Ag (nM)")
    ax5.set_ylabel("Fold Change")
    ax5.set_title("Y-Channel MFI Fold Change")
    ax5.set_xticks(concentrations)
    ax5.grid(True, alpha=0.3)

    # # 6. R6 Y-Channel Enrichment
    # ax6 = plt.subplot(2, 3, 6)
    # for i, sample in enumerate(sorted(unique_samples)):
    #     sample_data = results_df[results_df['Sample_Base'] == sample]
    #     x_vals = []
    #     y_vals = []
    #     for conc in concentrations:
    #         conc_data = sample_data[sample_data['Ag_Concentration_nM'] == conc]
    #         if not conc_data.empty and pd.notna(conc_data['R6_Y_Enrichment'].iloc[0]):
    #             x_vals.append(conc)
    #             y_vals.append(conc_data['R6_Y_Enrichment'].iloc[0])
    #     if len(x_vals) > 0:
    #         ax6.plot(x_vals, y_vals, marker='h', label=sample,
    #                 color=colors[i], linewidth=2, markersize=6, alpha=0.7)
    # ax6.set_xlabel('Ag (nM)')
    # ax6.set_ylabel('R6 Enrichment (fold)')
    # ax6.set_title('R6 Y-Channel Enrichment')
    # ax6.set_xticks(concentrations)
    # ax6.grid(True, alpha=0.3)
    # ax6.axhline(y=1, color='k', linestyle='--', alpha=0.3)

    # Add legend to first plot
    ax1.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize="small")

    plt.suptitle("MFI Analysis Dashboard", fontsize=14, fontweight="bold")
    plt.tight_layout()

    return fig


def create_mfi_heatmaps(results_df):
    """
    Create heatmaps for MFI data visualization.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Create pivot tables for heatmaps
    # 1. Expressing Combined MFI
    pivot_combined = results_df.pivot_table(
        values="Expression_Y_MFI",
        index="Sample_ID",
        columns="Ag_Concentration_nM",
        aggfunc="first",
    )

    sns.heatmap(
        pivot_combined,
        annot=True,
        fmt=".0f",
        cmap="RdYlBu_r",
        ax=axes[0, 0],
        cbar_kws={"label": "MFI"},
    )
    axes[0, 0].set_title("Expressing MFI", fontweight="bold")
    axes[0, 0].set_xlabel("Ag Concentration (nM)")
    axes[0, 0].set_ylabel("Sample ID")

    # 2. X-Channel Fold Change
    pivot_fold = results_df.pivot_table(
        values="X_MFI_Fold_Change",
        index="Sample_ID",
        columns="Ag_Concentration_nM",
        aggfunc="first",
    )

    sns.heatmap(
        pivot_fold,
        annot=True,
        fmt=".2f",
        cmap="coolwarm",
        center=1,
        ax=axes[0, 1],
        cbar_kws={"label": "Fold Change"},
    )
    axes[0, 1].set_title("X-Channel MFI Fold Change", fontweight="bold")
    axes[0, 1].set_xlabel("Ag Concentration (nM)")
    axes[0, 1].set_ylabel("Sample ID")

    # 3. R6 Combined MFI
    pivot_r6 = results_df.pivot_table(
        values="R6_X_MFI",
        index="Sample_ID",
        columns="Ag_Concentration_nM",
        aggfunc="first",
    )

    sns.heatmap(
        pivot_r6,
        annot=True,
        fmt=".0f",
        cmap="viridis",
        ax=axes[1, 0],
        cbar_kws={"label": "MFI"},
    )
    axes[1, 0].set_title("R6 Combined MFI", fontweight="bold")
    axes[1, 0].set_xlabel("Ag Concentration (nM)")
    axes[1, 0].set_ylabel("Sample ID")

    # 4. R6 X-Channel Enrichment
    pivot_enrichment = results_df.pivot_table(
        values="R6_X_Enrichment",
        index="Sample_ID",
        columns="Ag_Concentration_nM",
        aggfunc="first",
    )

    sns.heatmap(
        pivot_enrichment,
        annot=True,
        fmt=".2f",
        cmap="PuOr",
        center=1,
        ax=axes[1, 1],
        cbar_kws={"label": "Enrichment"},
    )
    axes[1, 1].set_title("R6 X-Channel Enrichment", fontweight="bold")
    axes[1, 1].set_xlabel("Ag Concentration (nM)")
    axes[1, 1].set_ylabel("Sample ID")

    plt.suptitle("MFI Analysis Heatmaps", fontsize=14, fontweight="bold")
    plt.tight_layout()

    return fig


def create_mfi_bar_plots(results_df):
    """
    Create bar plots comparing MFI values across concentrations.
    """
    # Set up the plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    concentrations = [0, 30, 120]

    # 1. Average Expressing MFI by concentration
    ax1 = axes[0, 0]
    avg_expr_mfi = []
    std_expr_mfi = []
    for conc in concentrations:
        conc_data = results_df[results_df["Ag_Concentration_nM"] == conc]
        avg_expr_mfi.append(conc_data["Expression_Y_MFI"].mean())
        std_expr_mfi.append(conc_data["Expression_Y_MFI"].std())

    bars1 = ax1.bar(
        concentrations,
        avg_expr_mfi,
        yerr=std_expr_mfi,
        capsize=5,
        color=["#3498db", "#e74c3c", "#f39c12"],
    )
    ax1.set_xlabel("Ag Concentration (nM)")
    ax1.set_ylabel("Average Combined MFI")
    ax1.set_title("Average Expressing Combined MFI", fontweight="bold")
    ax1.set_xticks(concentrations)

    # Add value labels on bars
    for bar, val in zip(bars1, avg_expr_mfi):
        height = bar.get_height()
        ax1.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            f"{val:.0f}",
            ha="center",
            va="bottom",
        )

    # 2. Average Fold Change by concentration
    ax2 = axes[0, 1]
    avg_fold = []
    std_fold = []
    for conc in concentrations:
        conc_data = results_df[results_df["Ag_Concentration_nM"] == conc]
        avg_fold.append(conc_data["X_MFI_Fold_Change"].mean())
        std_fold.append(conc_data["X_MFI_Fold_Change"].std())

    bars2 = ax2.bar(
        concentrations,
        avg_fold,
        yerr=std_fold,
        capsize=5,
        color=["#2ecc71", "#9b59b6", "#e67e22"],
    )
    ax2.set_xlabel("Ag Concentration (nM)")
    ax2.set_ylabel("Average Fold Change")
    ax2.set_title("Average X-Channel Fold Change", fontweight="bold")
    ax2.set_xticks(concentrations)
    ax2.axhline(y=1, color="k", linestyle="--", alpha=0.3)

    # Add value labels
    for bar, val in zip(bars2, avg_fold):
        height = bar.get_height()
        ax2.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            f"{val:.2f}",
            ha="center",
            va="bottom",
        )

    # 3. Average R6 MFI by concentration
    ax3 = axes[1, 0]
    avg_r6_mfi = []
    std_r6_mfi = []
    for conc in concentrations:
        conc_data = results_df[results_df["Ag_Concentration_nM"] == conc]
        avg_r6_mfi.append(conc_data["R6_X_MFI"].mean())
        std_r6_mfi.append(conc_data["R6_X_MFI"].std())

    bars3 = ax3.bar(
        concentrations,
        avg_r6_mfi,
        yerr=std_r6_mfi,
        capsize=5,
        color=["#1abc9c", "#34495e", "#f1c40f"],
    )
    ax3.set_xlabel("Ag Concentration (nM)")
    ax3.set_ylabel("Average R6 Combined MFI")
    ax3.set_title("Average R6 Combined MFI", fontweight="bold")
    ax3.set_xticks(concentrations)

    # Add value labels
    for bar, val in zip(bars3, avg_r6_mfi):
        height = bar.get_height()
        ax3.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            f"{val:.0f}",
            ha="center",
            va="bottom",
        )

    # 4. Average R6 Enrichment by concentration
    ax4 = axes[1, 1]
    avg_enrichment = []
    std_enrichment = []
    for conc in concentrations:
        conc_data = results_df[results_df["Ag_Concentration_nM"] == conc]
        avg_enrichment.append(conc_data["R6_X_Enrichment"].mean())
        std_enrichment.append(conc_data["R6_X_Enrichment"].std())

    bars4 = ax4.bar(
        concentrations,
        avg_enrichment,
        yerr=std_enrichment,
        capsize=5,
        color=["#16a085", "#8e44ad", "#d35400"],
    )
    ax4.set_xlabel("Ag Concentration (nM)")
    ax4.set_ylabel("Average R6 Enrichment")
    ax4.set_title("Average R6 X-Channel Enrichment", fontweight="bold")
    ax4.set_xticks(concentrations)
    ax4.axhline(y=1, color="k", linestyle="--", alpha=0.3)

    # Add value labels
    for bar, val in zip(bars4, avg_enrichment):
        height = bar.get_height()
        ax4.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            f"{val:.2f}",
            ha="center",
            va="bottom",
        )

    plt.suptitle("MFI Average Values by Concentration", fontsize=14, fontweight="bold")
    plt.tight_layout()

    return fig


def create_mfi_scatter_plots(results_df):
    """
    Create scatter plots showing relationships between MFI metrics.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Define colors for concentrations
    color_map = {0: "#3498db", 30: "#e74c3c", 120: "#f39c12"}

    # 1. Expressing X MFI vs Y MFI
    ax1 = axes[0, 0]
    for conc in [0, 30, 120]:
        conc_data = results_df[results_df["Ag_Concentration_nM"] == conc]
        ax1.scatter(
            conc_data["Expressing_X_MFI"],
            conc_data["Expressing_Y_MFI"],
            label=f"{conc} nM",
            alpha=0.6,
            s=100,
            color=color_map[conc],
        )
    ax1.set_xlabel("Expressing X MFI")
    ax1.set_ylabel("Expressing Y MFI")
    ax1.set_title("Expressing X vs Y MFI", fontweight="bold")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. R6 X MFI vs R6 Y MFI
    ax2 = axes[0, 1]
    for conc in [0, 30, 120]:
        conc_data = results_df[results_df["Ag_Concentration_nM"] == conc]
        ax2.scatter(
            conc_data["R6_X_MFI"],
            conc_data["R6_Y_MFI"],
            label=f"{conc} nM",
            alpha=0.6,
            s=100,
            color=color_map[conc],
        )
    ax2.set_xlabel("R6 X MFI")
    ax2.set_ylabel("R6 Y MFI")
    ax2.set_title("R6 X vs Y MFI", fontweight="bold")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. Combined MFI vs Fold Change
    ax3 = axes[1, 0]
    for conc in [0, 30, 120]:
        conc_data = results_df[results_df["Ag_Concentration_nM"] == conc]
        ax3.scatter(
            conc_data["Expression_Y_MFI"],
            conc_data["X_MFI_Fold_Change"],
            label=f"{conc} nM",
            alpha=0.6,
            s=100,
            color=color_map[conc],
        )
    ax3.set_xlabel("Expressing Combined MFI")
    ax3.set_ylabel("X MFI Fold Change")
    ax3.set_title("Combined MFI vs Fold Change", fontweight="bold")
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=1, color="k", linestyle="--", alpha=0.3)

    # 4. R6 Combined MFI vs R6 Enrichment
    ax4 = axes[1, 1]
    for conc in [0, 30, 120]:
        conc_data = results_df[results_df["Ag_Concentration_nM"] == conc]
        ax4.scatter(
            conc_data["R6_X_MFI"],
            conc_data["R6_X_Enrichment"],
            label=f"{conc} nM",
            alpha=0.6,
            s=100,
            color=color_map[conc],
        )
    ax4.set_xlabel("R6 X MFI")
    ax4.set_ylabel("R6 X Enrichment")
    ax4.set_title("R6 Combined MFI vs Enrichment", fontweight="bold")
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.axhline(y=1, color="k", linestyle="--", alpha=0.3)

    plt.suptitle("MFI Relationship Scatter Plots", fontsize=14, fontweight="bold")
    plt.tight_layout()

    return fig


# Main execution
if __name__ == "__main__":
    # Load MFI results
    results_df = pd.read_csv("mfi_results.csv")

    print("Creating MFI visualizations...")
    print("-" * 60)

    # Create all plot types
    fig1 = plot_mfi_analysis(results_df)
    plt.savefig("mfi_analysis_curves.png", dpi=300, bbox_inches="tight")
    print("✓ MFI analysis curves saved")

    fig2 = create_mfi_heatmaps(results_df)
    plt.savefig("mfi_heatmaps.png", dpi=300, bbox_inches="tight")
    print("✓ MFI heatmaps saved")

    fig3 = create_mfi_bar_plots(results_df)
    plt.savefig("mfi_bar_plots.png", dpi=300, bbox_inches="tight")
    print("✓ MFI bar plots saved")

    fig4 = create_mfi_scatter_plots(results_df)
    plt.savefig("mfi_scatter_plots.png", dpi=300, bbox_inches="tight")
    print("✓ MFI scatter plots saved")

    plt.show()

    print("\nAll MFI plots generated successfully!")
