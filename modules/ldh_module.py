# Standard library imports
import os
from pathlib import Path

# Third-party library imports
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator, PercentFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
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


LDH = LDHModule()
plots_path = Path(os.getcwd()) / "plots"
for extension in [".png", ".tiff"]:
    LDH.plot_replicates(filepath=plots_path / f"ldh_replicates{extension}")
# Due to bad data point, impute value for Nig 1_1_21
LDH.impute_value(treatment="Nigericin", time=21, replicate="1_1")
for extension in [".png", ".tiff"]:
    LDH.plot_replicates(filepath=plots_path / f"ldh_replicates_imputed{extension}")
    LDH.plot_bar_chart(filepath=plots_path / f"bar_chart_ldh{extension}")

LDH.plot_separate_bar_charts(filepath=plots_path / f"bar_chart_ldh.png")
LDH.anova_tukey_analysis()


##New Playground  -  MUST RUN MAIN FIRST##


def plot_everything(
    analysis_module,
    ldh_module,
    ratio_name,
    treatments=None,
    filepath=None,
    padding: tuple = (0, 0),
    invert=False,
):
    treatment = treatments[0]
    subset = ldh_module.data[ldh_module.data["Treatment"] == treatment]
    ratio_name_x = ratio_name.split(":")[0]
    ratio_name_y = ratio_name.split(":")[1]

    def get_color(treatment, palette):
        colors = {"ATP": palette[0], "MSU": palette[1], "Nigericin": palette[2]}
        return colors[treatment]

    color = get_color(treatment, sns.color_palette())
    pastel_color = get_color(treatment, sns.color_palette("pastel"))
    darker_pastel = tuple(max(0, x - 0.2) for x in pastel_color)

    TAS = analysis_module
    speck_info = TAS.modules["TS_Speck"]
    ratio_info = TAS.modules["TS_Cyto"].ratio_data[ratio_name]
    speck_data = speck_info.data[speck_info.data["Treatment"].isin(treatments)]
    ratio_data = ratio_info[ratio_info["Treatment"].isin(treatments)]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_ylabel("Normalized Speck Formation")
    ax.yaxis.labelpad = padding[0]

    ax2 = ax.twinx()
    sns.lineplot(
        data=speck_data,
        x="Time (hrs)",
        y="Measurement",
        color=color,
        ax=ax,
        errorbar="se",
        legend=False,
    )

    if invert == True:
        y_name = f"Ratio_{ratio_name_y}_{ratio_name_x}"
        clean_name = f"Ratio {ratio_name_y}:{ratio_name_x}".replace("b", "β")
        ax2.set_ylabel(f"Ratio of Cytokine Measurements, {ratio_name_y}:{ratio_name_x}")
    else:
        y_name = f"Ratio_{ratio_name_x}_{ratio_name_y}"
        clean_name = f"Ratio {ratio_name_x}:{ratio_name_y}".replace("b", "β")
        ax2.set_ylabel(f"Ratio of Cytokine Measurements, {ratio_name_x}:{ratio_name_y}")
    sns.pointplot(
        data=ratio_data,
        x="Time (hrs)",
        y=y_name,
        color=pastel_color,
        ax=ax2,
        errorbar="se",
    )

    ax2.yaxis.labelpad = padding[1]
    ax2.yaxis.tick_left()
    ax2.yaxis.set_label_position("left")
    ax2.tick_params(axis="y", direction="in", pad=-27.5, colors=darker_pastel)
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax2.set_ylabel(clean_name, color=darker_pastel)

    ax3 = ax.twinx()
    ax3.spines["right"].set_position(("outward", 0))
    ax3.set_ylabel("Cytotoxicity", color="red")
    ax3.tick_params(axis="y", colors="red")
    ax3.scatter(
        subset["Time (h)"], subset["cytotoxicity"], color="red", label=treatment, s=12
    )

    for bio_rep in subset["bio. Rep."].unique():
        for tech_rep in subset["tech. Rep"].unique():
            data_points = subset[
                (subset["bio. Rep."] == bio_rep) & (subset["tech. Rep"] == tech_rep)
            ].sort_values(by="Time (h)")
            ax3.plot(
                data_points["Time (h)"],
                data_points["cytotoxicity"],
                color=color,
                linestyle="--",
                alpha=0.25,
            )
    ax.set_ylim(
        bottom=ax.get_ylim()[0], top=ax.get_ylim()[1] * 1.1
    )  # Increase the upper limit by 10%
    ax2.set_ylim(
        bottom=0, top=ax2.get_ylim()[1] * 1.1
    )  # Increase the upper limit by 10%
    ax3.set_ylim(
        bottom=ax3.get_ylim()[0], top=ax3.get_ylim()[1] * 1.1
    )  # Increase the upper limit by 10%
    # After all your plotting and after setting the limits
    ticks = ax2.get_yticks()
    ticks = ticks[ticks != 0]  # remove 0 from ticks
    ticks = ticks[ticks < max(ticks)]  # this removes the highest tick
    ax2.set_yticks(ticks)  # set the new ticks

    legend_elements = [
        Line2D([0], [0], color=color, lw=2, label="Normalized Speck Formation"),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=pastel_color,
            markersize=10,
            label=clean_name,
            linestyle="None",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="red",
            markersize=8,
            label="Cytotoxicity",
            linestyle="None",
        ),
    ]

    ax.legend(
        handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, 1), ncol=3
    )
    plt.xlim(-1.25, 21.2)
    plt.xlabel("Time (hrs)")
    if filepath is not None:
        plt.savefig(filepath, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()
        return fig, ax


plot_everything(
    TAS,
    LDH,
    ratio_name="IL18:IL1b",
    treatments=["ATP"],
    padding=(15, 30),
    invert=True,
    filepath=Path("plots", "ratio_speck_ldh_combined_atp.png"),
)
plot_everything(
    TAS,
    LDH,
    ratio_name="IL18:IL1b",
    treatments=["MSU"],
    padding=(20, 40),
    invert=True,
    filepath=Path("plots", "ratio_speck_ldh_combined_msu.png"),
)
plot_everything(
    TAS,
    LDH,
    ratio_name="IL18:IL1b",
    treatments=["Nigericin"],
    padding=(20, 40),
    invert=True,
    filepath=Path("plots", "ratio_speck_ldh_combined_nigericin.png"),
)
