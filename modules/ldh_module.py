# Standard library imports
import os
from pathlib import Path

# Third-party library imports
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd


class LDHModule:
    def __init__(self):
        self.od = pd.read_excel(
            Path(os.getcwd()) / "data" / "ldh" / "ldh_data.xlsx", sheet_name="OD"
        )
        self.ldh = pd.read_excel(
            Path(os.getcwd()) / "data" / "ldh" / "ldh_data.xlsx", sheet_name="LDH"
        )
        self.od_adjust()
        self.calculate_cytotoxicity()

    def od_adjust(self):
        self.od = self.od
        self.od["Adjusted_OD"] = self.od["Abs 490"] - self.od["Abs 680"]
        self.od["replicate"] = (
            self.od["bio. Rep."].astype(str) + "_" + self.od["tech. Rep"].astype(str)
        )

    # def check_technical_replicates(self):

    def calculate_cytotoxicity(self):
        mergedf = self.od.merge(self.ldh, on=["Time (h)"])
        max_ldh = mergedf["Maximum LDH"].max()
        mergedf["cytotoxicity"] = (
            mergedf["Adjusted_OD"] - mergedf["Spontaneous LDH"]
        ) / (max_ldh - mergedf["Spontaneous LDH"])

        # mergedf['cytotoxicity'] = (mergedf['Adjusted_OD'] - mergedf['Spontaneous LDH']) / (mergedf['Maximum LDH'] - mergedf['Spontaneous LDH'])
        self.data = mergedf

    def plot_replicates(self, filepath: str = None):
        # Define colors for treatments
        default_colors = sns.color_palette()
        colors = {
            "ATP": default_colors[0],
            "MSU": default_colors[1],
            "Nigericin": default_colors[2],
        }
        # Plotting
        for treatment, color in colors.items():
            subset = self.data[self.data["Treatment"] == treatment]
            plt.scatter(
                subset["Time (h)"], subset["cytotoxicity"], color=color, label=treatment
            )

            # Connect data points across time for same bio. Rep. and tech. Rep
            for bio_rep in subset["bio. Rep."].unique():
                for tech_rep in subset["tech. Rep"].unique():
                    data_points = subset[
                        (subset["bio. Rep."] == bio_rep)
                        & (subset["tech. Rep"] == tech_rep)
                    ].sort_values(by="Time (h)")
                    plt.plot(
                        data_points["Time (h)"],
                        data_points["cytotoxicity"],
                        color=color,
                        linestyle="--",
                        alpha=0.5,
                    )

        plt.xlabel("Time (h)")
        plt.ylabel("Cytotoxicity")
        plt.title("Scatter plot of Cytotoxicity")
        plt.legend()
        if filepath is not None:
            plt.savefig(filepath, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def impute_value(self, treatment, time, replicate):
        # get all values for treatment, time, and replicate
        subset = self.data[
            (self.data["Treatment"] == treatment) & (self.data["Time (h)"] == time)
        ]
        # get subset where replicate is not the one we want to impute
        subset_norep = subset[subset["replicate"] != replicate]
        # get mean of cytotoxicity for subset
        mean = subset_norep["cytotoxicity"].mean()
        # set value of replicate to mean
        self.data.loc[
            (self.data["Treatment"] == treatment)
            & (self.data["Time (h)"] == time)
            & (self.data["replicate"] == replicate),
            "cytotoxicity",
        ] = mean

    def plot_bar_chart(self, filepath: str = None):
        # Group by Treatment and Time (h) and calculate mean and standard error
        grouped = self.data.groupby(["Treatment", "Time (h)"])["cytotoxicity"].agg(
            ["mean", "std", "count"]
        )
        grouped["SE"] = grouped["std"] / np.sqrt(grouped["count"])

        default_colors = sns.color_palette()
        colors = {
            "ATP": default_colors[0],
            "MSU": default_colors[1],
            "Nigericin": default_colors[2],
        }

        # Plotting
        width = 0.2  # width of the bars
        positions = list(range(len(self.data["Time (h)"].unique())))
        for idx, (treatment, color) in enumerate(colors.items()):
            means = grouped.loc[treatment]["mean"]
            errors = grouped.loc[treatment]["SE"]
            plt.bar(
                [p + idx * width for p in positions],
                means,
                yerr=errors,
                color=color,
                width=width,
                label=treatment,
                capsize=5,
                align="center",
            )

            # Format the Y-axis as percentages
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))

        plt.xlabel("Time (h)")
        plt.ylabel("Cytotoxicity (%)")  # Update ylabel to indicate percentages
        plt.title("Percentage Cytotoxicity")
        plt.xticks(positions, self.data["Time (h)"].unique())
        plt.legend()
        plt.tight_layout()
        if filepath is not None:
            plt.savefig(filepath, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def plot_separate_bar_charts(self, filepath: str = None):
        # Group by Treatment and Time (h) and calculate mean and standard error
        grouped = self.data.groupby(["Treatment", "Time (h)"])["cytotoxicity"].agg(
            ["mean", "std", "count"]
        )
        grouped["SE"] = grouped["std"] / np.sqrt(grouped["count"])

        default_colors = sns.color_palette()
        colors = {
            "ATP": default_colors[0],
            "MSU": default_colors[1],
            "Nigericin": default_colors[2],
        }

        # Plotting
        width = 0.2  # width of the bars
        positions = list(range(len(self.data["Time (h)"].unique())))

        for idx, (treatment, color) in enumerate(colors.items()):
            plt.figure(idx)  # create a new figure for each treatment
            means = grouped.loc[treatment]["mean"]
            errors = grouped.loc[treatment]["SE"]
            plt.bar(
                positions,
                means,
                yerr=errors,
                color=color,
                width=width,
                label=treatment,
                capsize=5,
                align="center",
            )

            # Format the Y-axis as percentages
            plt.gca().yaxis.set_major_formatter(PercentFormatter(1))

            plt.xlabel("Time (h)")
            plt.ylabel("Cytotoxicity (%)")  # Update ylabel to indicate percentages
            plt.title(f"Percentage Cytotoxicity - {treatment}")
            plt.xticks(positions, self.data["Time (h)"].unique())
            plt.legend()
            plt.tight_layout()
            if filepath is not None:
                file, ext = os.path.splitext(filepath)
                new_filepath = f"{file}_{treatment}{ext}"
                plt.savefig(new_filepath, dpi=300, bbox_inches="tight")
                plt.close()
            else:
                plt.show()

    def anova_tukey_analysis(self, output_filename="reports//anova_tukey_results.xlsx"):
        # Create a Pandas Excel writer using XlsxWriter as the engine
        with pd.ExcelWriter(output_filename, engine="xlsxwriter") as writer:
            # For each time point
            for time in self.data["Time (h)"].unique():
                subset = self.data[self.data["Time (h)"] == time]

                # Add a "zero" group
                zero_group = pd.DataFrame(
                    {
                        "cytotoxicity": [0 for _ in range(subset.shape[0])],
                        "Treatment": ["Null" for _ in range(subset.shape[0])],
                    }
                )
                subset = pd.concat([subset, zero_group], ignore_index=True)

                # ANOVA
                model = ols("cytotoxicity ~ Treatment", data=subset).fit()
                anova_table = sm.stats.anova_lm(model, typ=2)
                p_value = anova_table["PR(>F)"]["Treatment"]

                # Create a new DataFrame to store results for this time point
                results_df = pd.DataFrame()
                results_df = results_df.append(anova_table)

                # If p-value is significant (e.g., < 0.05), proceed with Tukey's HSD test
                if p_value < 0.05:
                    tukey = pairwise_tukeyhsd(
                        endog=subset["cytotoxicity"],
                        groups=subset["Treatment"],
                        alpha=0.05,
                    )
                    tukey_df = pd.DataFrame(
                        data=tukey._results_table.data[1:],
                        columns=tukey._results_table.data[0],
                    )
                    results_df = results_df.append(tukey_df)

                # Write results DataFrame to Excel
                results_df.to_excel(writer, sheet_name=f"Time {time}h")

        print(f"Results exported to {output_filename}")
