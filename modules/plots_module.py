import pandas as pd
import numpy as np
from typing import Type, List, Dict, Callable
import matplotlib.pyplot as plt
import seaborn as sns
from modules.stats_functions import bootstrap_t_test
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker


class plotting_module:
    treatment_indeces = {
        "ATP": 0,
        "MSU": 1,
        "Nigericin": 2,
    }

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
        fig, ax = plt.subplots(figsize=(10, 6))
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
    def plot_pointplot(
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
        fig, ax = plt.subplots(figsize=(10, 6))
        pointplot = sns.pointplot(
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
        # for line in pointplot.lines:
        #     line.set_linestyle("dashed")
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
    def plot_multimodel(
        modelA,
        modelB,
        measurement_type="Normalized_Measurement",
        output_directory=None,
    ):
        ##Confirm that ModelA and ModelB have the same treatments
        if modelA.comp.treatments != modelB.comp.treatments:
            raise ValueError("Models do not have the same treatments.")
        for treatment in modelA.comp.treatments:
            dfA = modelA.data[modelA.data["Treatment"] == treatment]
            dfB = modelB.data[modelB.data["Treatment"] == treatment]
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.set_ylabel("Normalized Speck Formation")
            ax2 = ax.twinx()
            ax2.set_ylabel("Normalized Cytokine Measurement")
            for model in [dfA, dfB]:
                ax_to_use = ax if model is dfA else ax2
                if len(model["Analyte"].unique()) > 1:
                    sns.lineplot(
                        data=model,
                        x="Time (hrs)",
                        y=measurement_type,
                        hue="Analyte",
                        ax=ax_to_use,
                        errorbar="se",
                    )
                else:
                    sns.lineplot(
                        data=model,
                        x="Time (hrs)",
                        y=measurement_type,
                        color="green",
                        ax=ax_to_use,
                        errorbar="se",
                    )
            ax.set_ylim(1, None)  # Set the lower limit of ax's y-axis to 1
            ax2.set_ylim(1, None)  # Set the lower limit of ax2's y-axis to 1
            plt.xlabel("Time (hrs)")
            # plt.ylabel("Normalized Value")
            plt.title(
                str(
                    "Normalized Measurements for IL1b, IL18, and Normalized Speck Formation - "
                    + treatment
                )
            )
            if output_directory is not None:
                filepath = output_directory / f"{treatment}_multimodel.png"
                plt.savefig(filepath, dpi=300, bbox_inches="tight")
                plt.close()
            else:
                plt.show()

    @staticmethod
    def plot_count_against_ratio(
        analysis_module,
        ratio_name,
        treatments=None,
        invert=False,
        manual_ax_modification=None,
        filepath=None,
    ):
        TAS = analysis_module
        speck_info = TAS.modules["TS_Speck"]
        ratio_info = TAS.modules["TS_Cyto"].ratio_data[ratio_name]
        ratio_name_x = ratio_name.split(":")[0]
        ratio_name_y = ratio_name.split(":")[1]
        if treatments is None:
            treatments = speck_info.comp.treatments
        if len(treatments) == 1:
            treatment = treatments[0]
        elif treatments == speck_info.comp.treatments:
            treatment = "All Treatments"
        speck_data = speck_info.data[speck_info.data["Treatment"].isin(treatments)]
        ratio_data = ratio_info[ratio_info["Treatment"].isin(treatments)]
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.set_ylabel("Normalized Speck Formation")
        ax2 = ax.twinx()
        lineplot_ax1 = sns.lineplot(
            data=speck_data,
            x="Time (hrs)",
            y="Measurement",
            hue="Treatment",
            ax=ax,
            errorbar="se",
        )
        if invert == True:
            y_name = f"Ratio_{ratio_name_y}_{ratio_name_x}"
            ax2.set_ylabel(
                f"Ratio of Cytokine Measurements, {ratio_name_y}:{ratio_name_x}"
            )
        else:
            y_name = f"Ratio_{ratio_name_x}_{ratio_name_y}"
            ax2.set_ylabel(
                f"Ratio of Cytokine Measurements, {ratio_name_x}:{ratio_name_y}"
            )
        pointplot = sns.pointplot(
            data=ratio_data,
            x="Time (hrs)",
            y=y_name,
            hue="Treatment",
            ax=ax2,
            palette="pastel",
            errorbar="se",
        )
        line_legend = ax.legend(loc="upper left", title="Speck Formation")
        point_legend = ax2.legend(loc="upper right", title="Cytokine Ratio")
        if manual_ax_modification is not None:
            manual_ax_modification(
                ax, ax2
            )  # Call the callback function with the ax parameter
        # ax.set_ylim(1, None)  # Set the lower limit of ax's y-axis to 1
        # ax2.set_ylim(1, None)
        # ax.get_legend().remove()  # Set the lower limit of ax2's y-axis to 1
        # for line in pointplot.lines:
        #     line.set_linestyle("dotted")
        plt.xlim(0, 21)
        plt.xlabel("Time (hrs)")

        # plt.title(
        #     str(
        #         "Normalized Measurements for IL1b, IL18, and Normalized Speck Formation - "
        #         + treatment
        #     )
        # )
        if filepath is not None:
            plt.savefig(filepath, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()
            return fig, ax

    @staticmethod
    def plot_ratio_old(
        module, measurement_type="Measurement", normalize_start=False, invert=False
    ):
        ##If the module's data does not have two analytes, return an error
        if len(module.data["Analyte"].unique()) < 2:
            raise ValueError("Model does not have two analytes.")
        split_data = plotting_module.split_analytes(module)
        df1, df2 = split_data.values()
        if invert == True:
            df1, df2 = df2, df1
        analyte1, analyte2 = split_data.keys()
        merged_data = pd.merge(
            df1, df2, on=["Treatment", "Time (hrs)", "Experimental_Replicate"]
        )
        merged_data[f"Ratio_{analyte1}_{analyte2}"] = (
            merged_data[f"{measurement_type}_x"] / merged_data[f"{measurement_type}_y"]
        )
        # If normalize_start is True, normalize the ratio to the first time point for each treatment and experimental replicate
        if normalize_start:
            merged_data[f"Ratio_{analyte1}_{analyte2}"] = merged_data.groupby(
                ["Treatment", "Experimental_Replicate"]
            )[f"Ratio_{analyte1}_{analyte2}"].transform(lambda x: x / x.iloc[0])

        plotting_module.plot_lineplot(
            merged_data, y=f"Ratio_{analyte1}_{analyte2}", errorbar="se"
        )
        return merged_data

    @staticmethod
    def plot_ratio(
        module,
        invert=False,
        analyte_pair=None,
        manual_ax_modification=None,
        filepath=None,
    ):
        ##If module does not have .ratio_data, return an error
        if module.ratio_data is None:
            raise ValueError("Module does not have ratio_data.")
        ##If analyte_pair is None, use the first key in the ratio_data dictionary
        if analyte_pair is None:
            analyte_pair = list(module.ratio_data.keys())[0]
        merged_data = module.ratio_data[analyte_pair]
        if invert == True:
            analyte1 = analyte_pair.split(":")[1]
            analyte2 = analyte_pair.split(":")[0]
        else:
            analyte1 = analyte_pair.split(":")[0]
            analyte2 = analyte_pair.split(":")[1]
        plotting_module.plot_pointplot(
            merged_data,
            y=f"Ratio_{analyte1}_{analyte2}",
            errorbar="se",
            manual_ax_modification=manual_ax_modification,
            filepath=str(filepath),
        )
        # return fig, ax

    @staticmethod
    def split_analytes(module):
        """
        Split the data into separate dataframes for each analyte.

        Args:
            model (time_sequence): The time_sequence object to split.

        Returns:
            dict: A dictionary of dataframes for each analyte.
        """
        analytes = module.data["Analyte"].unique()
        split_data = {}
        for analyte in analytes:
            split_data[analyte] = module.data[module.data["Analyte"] == analyte]
        return split_data

    @staticmethod
    def show_p(data, y, separator, ax):
        # Get the current ylim values
        ymin, ymax = ax.get_ylim()
        # Set the new ylim with an increased upper limit
        ax.set_ylim(ymin, ymax * 1.25)
        # Get the treatments in the data
        treatments = data[separator].unique()
        p_values_text = ""
        for treatment in treatments:
            treatment_data = data[data[separator] == treatment]
            max_time = treatment_data.loc[treatment_data[y].idxmax(), "Time (hrs)"]

            # max_time_per_replicate = (
            #     treatment_data.groupby("Experimental_Replicate")[y]
            #     .idxmax()
            #     .reset_index()
            # )
            # # Merge the original treatment_data with max_time_per_replicate to get the values for each replicate at their max times
            # max_values_per_replicate = treatment_data.loc[max_time_per_replicate[y]]
            # # Filter the max_values_per_replicate DataFrame for the current treatment
            # treatment_max_values = max_values_per_replicate[
            #     max_values_per_replicate[separator] == treatment
            # ]

            time_zero_data = np.array(
                treatment_data[treatment_data["Time (hrs)"] == 0][y]
            )
            max_time_data = np.array(
                treatment_data[treatment_data["Time (hrs)"] == max_time][y]
            )
            # print(time_zero_data, max_time_data)
            p, ci = bootstrap_t_test(max_time_data, time_zero_data)
            # print(p)
            p_values_text += f"{treatment} (Max = {max_time} hrs):\n p = {p:.4f}\n"

        ax.text(
            0.02,
            0.98,
            p_values_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            horizontalalignment="left",
        )
        return ax

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
        ax.set_ylabel("Cytotoxicity", color="red")
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
                label="Cytotoxicity",
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
    def combined_plot(analysis_class, treatments=None, filepath=None, padding=(0, 0)):
        fig, ax = plt.subplots(figsize=(10, 6))
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
        ax.legend(
            handles=all_legend_elements,
            loc="upper center",
            bbox_to_anchor=(0.5, 1),
            ncol=3,
        )
        plt.xlim(-1.25, 21.2)
        plt.xlabel("Time (hrs)")

        if filepath is not None:
            plt.savefig(filepath, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()
