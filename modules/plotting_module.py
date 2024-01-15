from typing import Type, List, Dict, Callable
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from modules.stats_functions import bootstrap_t_test
import textwrap

class plotting_module:
    
    treatment_indeces = {
        "ATP": 0,
        "MSU": 1,
        "Nigericin": 2,
    }
    linestyles = list(mlines.Line2D.lineStyles.keys())

    @staticmethod
    def show_colors(index):
        """
        Display the colors from the custom palette as a bar plot.

        :param index: Index of color in seaborn's default palette.
        :type index: int
        """
        palette = plotting_module.custom_palette(index)

        # Create a new figure
        fig, ax = plt.subplots(figsize=(3, 1))  # 3 bars, each with width of 1

        # Plot a bar for each color in the palette
        for i, color in enumerate(palette):
            ax.bar(i, 1, color=color)

        # Remove y-axis and x-axis ticks and labels
        ax.set_yticks([])
        ax.set_xticks([])

        # Remove the x and y axis spines
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        plt.show()

    @staticmethod
    def custom_palette(index):
        """
        Generate a custom color palette based on a color index from seaborn's default palette.

        :param index: Index of color in seaborn's default palette.
        :type index: int
        :return: List of three color codes.
        :rtype: list
        """
        # Get the default color
        default_color = sns.color_palette()[index]

        # Get a pastel version of the color
        pastel_color = sns.color_palette("pastel")[index]

        # Get a dark version of the color
        dark_color = sns.color_palette("dark")[index]

        return [default_color, pastel_color, dark_color]

    @staticmethod
    def change_plot(
        module, analytes=None, treatments=None, mode="change", axis=None, filepath=None, normalized=False,
    ):
        sns.set_context("notebook", font_scale=2)
        if axis == None:
            fig, ax = plt.subplots(figsize=(16, 9))
        else:
            ax = axis
        if mode.lower() == "delta":
            y = "Delta"
        elif mode.lower() == "acceleration":
            y = "Acceleration"
        elif mode.lower() == "smoothed acceleration":
            y = "Acceleration_Smooth"
        elif mode.lower() == "combo":
            y1 = "Delta"
            y2 = "Acceleration"
            y3 = "Acceleration_Smooth"
            ax2 = ax.twinx()
        else:
            raise ValueError("Invalid Mode Specified")
        if normalized == True:
            y = f"Normalized_{y}"       
        if treatments is None:
            treatments = module.data["Treatment"].unique()
        if module.name == "TS_Speck":
            label = "Speck Count"
        elif module.name == "TS_Cyto":
            label = ""
        legend_elements = []
        if analytes is None:
            analyte_list = module.data["Analyte"].unique()
        else:
            analyte_list = analytes
        for i, analyte in enumerate(analyte_list):
            # print(i, analyte)
            linestyle = plotting_module.linestyles[i]
            if analytes is not None or len(analyte_list) > 1:
                analyte_tag = f" - {analyte}"
            else:
                analyte_tag = ""
            data = module.data[module.data["Analyte"] == analyte]
            for treatment in treatments:
                treatment_data = data[data["Treatment"] == treatment]
                # Get the color for the treatment
                index = plotting_module.treatment_indeces[treatment]
                if mode.lower() != "combo":
                    palette = (
                        plotting_module.custom_palette(index)
                        if module.name == "TS_Speck"
                        else plotting_module.custom_palette(index)[1:]
                    )
                    sns.lineplot(
                        data=treatment_data,
                        x="Time (hrs)",
                        y=y,
                        ax=ax,
                        errorbar="se",
                        color=palette[i],
                        legend=False,
                        linestyle=linestyle,
                    )
                    if mode.lower() != "delta":
                        legend_elements.append(
                        Line2D(
                            [0],
                            [0],
                            color=palette[i],
                            lw=2,
                            label=f"{label} {mode.title()} {analyte_tag.replace('b', 'β')}",
                            linestyle=linestyle,
                        ),
                    )
                    else:
                        legend_elements.append(
                        Line2D(
                            [0],
                            [0],
                            color=palette[i],
                            lw=2,
                            label=f"{label} Change {analyte_tag.replace('b', 'β')}",
                            linestyle=linestyle,
                        ),
                    )
                    if mode.lower() != "delta":
                        ax.set_ylabel(f"{mode.title()}")
                    else:
                        ax.set_ylabel(f"Change")
                else:
                    palette = plotting_module.custom_palette(index)
                    # ax2 = ax.twinx()
                    sns.lineplot(
                        data=treatment_data,
                        x="Time (hrs)",
                        y=y1,
                        ax=ax,
                        errorbar="se",
                        color=palette[0],
                        legend=False,
                        linestyle=linestyle,
                    )
                    sns.lineplot(
                        data=treatment_data,
                        x="Time (hrs)",
                        y=y2,
                        ax=ax2,
                        errorbar="se",
                        color=palette[1],
                        legend=False,
                        linestyle=linestyle,
                    )
                    sns.lineplot(
                        data=treatment_data,
                        x="Time (hrs)",
                        y=y3,
                        ax=ax2,
                        color=palette[2],
                        legend=False,
                        linestyle=linestyle,
                    )
                    legend_elements.append(
                        Line2D(
                            [0],
                            [0],
                            color=palette[0],
                            lw=2,
                            label=f"{label} {y1.replace('_', ' ')} {analyte_tag}",
                            linestyle=linestyle,
                        )
                    )
                    legend_elements.append(
                        Line2D(
                            [0],
                            [0],
                            color=palette[1],
                            lw=2,
                            label=f"{label} {y2.replace('_', ' ')} {analyte_tag}",
                            linestyle=linestyle,
                        ),
                    )
                    legend_elements.append(
                        Line2D(
                            [0],
                            [0],
                            color=palette[2],
                            lw=2,
                            label=f"{label} {y3.replace('_', ' ')} {analyte_tag}",
                            linestyle=linestyle,
                        ),
                    )
                    ax.set_ylabel("Total Change")
                    ax2.set_ylabel("Acceleration")

        if axis == None:
            ax.legend(
                handles=legend_elements,
                loc="upper center",
                bbox_to_anchor=(0.5, 1),
                ncol=3,
            )
            plt.xlim(-1.25, 24)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
            if filepath is not None:
                ax.set_title(treatment)
                #ax.set_title(filepath.split("/")[-1].split(".")[0].replace("_", " ").replace("Delta", "Change"))
                plt.savefig(filepath, dpi=300, bbox_inches="tight")
                plt.close()
            else:
                plt.show()

        else:
            return legend_elements

    @staticmethod
    def plot_lineplot(
        dataframe: pd.DataFrame,
        y: str,
        x: str = "Time (hrs)",
        hue: str = "Treatment",
        hue_order: List = None,
        errorbar: str = "se",
        filepath: str = None,
        show_p=False,
        manual_ax_modification=None,
    ):
        sns.set_context("notebook", font_scale=1.75)

        fig, ax = plt.subplots(figsize=(16, 9))
        lineplot = sns.lineplot(
            data=dataframe,
            x=x,
            y=y,
            hue=hue,
            hue_order=hue_order,
            errorbar=errorbar,
            ax=ax,
        )
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_title(f"{y} Over Time")
        ax.legend()
        if manual_ax_modification is not None:
            manual_ax_modification(
                ax
            )  # Call the callback function with the ax parameter
        if show_p == True:
            ax = plotting_module.show_p(data=dataframe, y=y, separator=hue, ax=ax)
        if filepath is not None:
            plt.savefig(filepath, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()
            return fig, ax


    @staticmethod
    def plot_speck_count(
        analysis_class,
        treatments=None,
        filepath=None,
        axis=None,
    ):
        if treatments is None:
            treatments = analysis_class.modules["TS_Speck"].data["Treatment"].unique()

        speck_info = analysis_class.modules["TS_Speck"]
        if axis == None:
            fig, ax = plt.subplots(figsize=(10, 6))

        else:
            ax = axis
        legend_elements = []
        ax.set_ylabel("Speck Count")
        # ax.yaxis.labelpad = padding[0]
        for treatment in treatments:
            # Get the color for the treatment
            index = plotting_module.treatment_indeces[treatment]
            palette = plotting_module.custom_palette(index)

            # Filter data for the specific treatmen
            speck_data = speck_info.data[speck_info.data["Treatment"] == treatment]

            sns.lineplot(
                data=speck_data,
                x="Time (hrs)",
                y="Measurement",
                color=palette[0],
                ax=ax,
                errorbar="se",
                legend=False,
            )
            # Adds some white space so the legend doesn't overlap with the plot
            if len(treatments) > 1:
                specifier = f" - {treatment}"
            else:
                specifier = ""

            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    color=palette[0],
                    lw=2,
                    label=f"Speck Count{specifier}",
                ),
            )

        if axis == None:
            ax.set_ylim(bottom=ax.get_ylim()[0], top=ax.get_ylim()[1] * 1.1)
            plt.xlim(-1.25, 21)
            plt.xlabel("Time (hrs)")
            ax.legend(
                handles=legend_elements,
                loc="upper center",
                bbox_to_anchor=(0.5, 1),
                ncol=3,
            )
            if filepath is not None:
                plt.savefig(filepath, dpi=300, bbox_inches="tight")
                plt.close()
            else:
                plt.show()
                return fig, ax

        else:
            return legend_elements

    @staticmethod
    def plot_cytokines(
        analysis_class,
        treatments=None,
        filepath=None,
        padding: tuple = (0, 0),
        axis=None,
    ):
        if treatments is None:
            treatments = analysis_class.modules["TS_Cyto"].data["Treatment"].unique()

        legend_elements = []

        Cyto = analysis_class.modules["TS_Cyto"].data

        if axis == None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            ax = axis

        for treatment in treatments:
            # Get the color for the treatment
            index = plotting_module.treatment_indeces[treatment]
            palette = plotting_module.custom_palette(index)[1:]
            sub_cyto = Cyto[Cyto["Treatment"] == treatment]

            if "Analyte" not in sub_cyto:
                raise ValueError("'Analyte' column not found in the dataset.")

            sns.lineplot(
                data=sub_cyto,
                x="Time (hrs)",
                y="Measurement",
                hue="Analyte",
                ax=ax,
                errorbar="se",
                legend=False,
                palette=palette,
                style="Analyte",  # using style
                # dashes=[(2, 2)] * len(sub_cyto["Analyte"].unique()),
            )
            if len(treatments) > 1:
                specifier = f" - {treatment}"
            else:
                specifier = ""
            linestyles = ["dotted", "dashed"]
            for i, analyte in enumerate(sub_cyto["Analyte"].unique()):
                legend_elements.append(
                    Line2D(
                        [0],
                        [0],
                        color=palette[i],
                        lw=2,
                        label=f"{analyte.replace('b', 'β')}{specifier}",
                        linestyle=linestyles[i],
                    )
                )
        ax.lines[0].set_linestyle("dotted")

        # Axis labels and formatting
        clean_name = f"Cytokine Concentration (ng/ml)".replace("b", "β")
        ax.set_ylabel(clean_name, color=palette[1], labelpad=padding[1])
        ax.set_xlabel("Time (hrs)")
        ax.yaxis.tick_left()
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
        ax.yaxis.set_label_position("left")
        ax.set_ylim(bottom=0, top=ax.get_ylim()[1] * 1.1)

        # Adjust Axis tick marks
        ax.tick_params(axis="y", direction="in", pad=padding[0], colors=palette[1])
        ticks = ax.get_yticks()
        ticks = ticks[ticks != 0]  # remove 0 from ticks
        ticks = ticks[ticks < max(ticks)]  # this removes the highest tick
        ax.set_yticks(ticks)  # set the new ticks
        for label in ax.get_yticklabels():
            label.set_horizontalalignment("left")
        if axis == None:
            # ax.set_ylim(bottom=ax.get_ylim()[0], top=ax.get_ylim()[1] * 1.1)
            # plt.xlim(-1.25, 21.2)
            # plt.xlabel("Time (hrs)")
            ax.legend(
                handles=legend_elements,
                loc="upper center",
                bbox_to_anchor=(0.5, 1),
                ncol=3,
            )
            plt.xlim(-1.25, 24)
            if filepath is not None:
                plt.savefig(filepath, dpi=300, bbox_inches="tight")
                plt.close()
            else:
                plt.show()
                return fig, ax

        else:
            return legend_elements

    @staticmethod
    def plot_cytotoxicity(
        analysis_class,
        treatments=None,
        filepath=None,
        axis=None,
    ):
        data = analysis_class.modules["LDH"].data
        if treatments is None:
            treatments = data["Treatment"].unique()

        if axis == None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            ax = axis
        # ax = ax.twinx()
        legend_elements = []

        for treatment in treatments:
            cytotoxicity = data[data["Treatment"] == treatment]
            # change cytotoxicity to a percentage
            cytotoxicity["cytotoxicity"] = cytotoxicity["cytotoxicity"] * 100
            cytotox_means = (
                cytotoxicity.groupby(["Time (h)"])["cytotoxicity"].mean().reset_index()
            )
            sns.scatterplot(
                data=cytotoxicity,
                x="Time (h)",
                y="cytotoxicity",
                color="red",
                ax=ax,
                #    errorbar="se",
                legend=False,
            )
            sns.scatterplot(
                data=cytotox_means,
                x="Time (h)",
                y="cytotoxicity",
                color="black",
                marker="_",
                s=500,
                ax=ax,
            )
        ax.spines["right"].set_position(("outward", 0))
        ax.set_ylabel("Cell Death (% of Maximum)", color="red")
        ax.tick_params(axis="y", colors="red")
        # Add % symbol to tick values on ax3
        ax.yaxis.set_major_formatter(ticker.PercentFormatter())

        legend_elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="red",
                markersize=8,
                label="Cell Death",
                linestyle="None",
            )
        )

        ax.set_ylim(bottom=ax.get_ylim()[0], top=ax.get_ylim()[1] * 1.1)
        if axis == None:
            # ax.set_ylim(bottom=ax.get_ylim()[0], top=ax.get_ylim()[1] * 1.1)
            # plt.xlim(-1.25, 21.2)
            # plt.xlabel("Time (hrs)")
            ax.legend(
                handles=legend_elements,
                loc="upper center",
                bbox_to_anchor=(0.5, 1),
                ncol=3,
            )
            if filepath is not None:
                plt.savefig(filepath, dpi=300, bbox_inches="tight")
                plt.close()
            else:
                plt.show()
                return fig, ax

        else:
            return legend_elements


    @staticmethod
    def combined_plot(analysis_class, treatments=None, filepath=None, padding=(0, 0), font_scale=1):
        sns.set_context("notebook", font_scale=font_scale)

        fig, ax = plt.subplots(figsize=(16, 9))
        ax2 = ax.twinx()
        ax3 = ax.twinx()
        all_legend_elements = []
        legend_elements = plotting_module.plot_speck_count(
            analysis_class, treatments=treatments, axis=ax
        )
        all_legend_elements += legend_elements
        legend_elements = plotting_module.plot_cytokines(
            analysis_class,
            treatments=treatments,
            axis=ax2,
            padding=padding,
        )
        all_legend_elements += legend_elements
        legend_elements = plotting_module.plot_cytotoxicity(
            analysis_class,
            treatments=treatments,
            axis=ax3,
        )
        all_legend_elements += legend_elements
        ax.legend(
            handles=all_legend_elements,
            loc="upper center",
            bbox_to_anchor=(0.5, 1),
            ncol=3,
        )
        plt.xlim(-1.25, 21.2)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
        plt.xlabel("Time (hrs)")
        ax.set_title(treatments[0])

        if filepath is not None:
            plt.savefig(filepath, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()
