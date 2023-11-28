import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# Assuming custom_palette function is already defined


def show_colors(index):
    """
    Display the colors from the custom palette as a bar plot.

    :param index: Index of color in seaborn's default palette.
    :type index: int
    """
    palette = custom_palette(index)

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


def plot_thing_1(
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
    # change cytotoxicity to a percentage
    subset["cytotoxicity"] = subset["cytotoxicity"] * 100
    cyto_means = subset.groupby(["Time (h)"])["cytotoxicity"].mean().reset_index()

    ratio_name_x = ratio_name.split(":")[0]
    ratio_name_y = ratio_name.split(":")[1]

    def get_color(treatment, palette):
        colors = {"ATP": palette[0], "MSU": palette[1], "Nigericin": palette[2]}
        return colors[treatment]

    treatment_indeces = {"ATP": 0, "MSU": 1, "Nigericin": 2}

    color = get_color(treatment, sns.color_palette())
    pastel_color = get_color(treatment, sns.color_palette("pastel"))
    darker_pastel = tuple(max(0, x - 0.2) for x in pastel_color)

    TAS = analysis_module
    Cyto = TAS.modules["TS_Cyto"].data
    sub_cyto = Cyto[Cyto["Treatment"] == treatment]
    speck_info = TAS.modules["TS_Speck"]
    ratio_info = TAS.modules["TS_Cyto"].ratio_data[ratio_name]

    speck_data = speck_info.data[speck_info.data["Treatment"].isin(treatments)]
    ratio_data = ratio_info[ratio_info["Treatment"].isin(treatments)]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_ylabel("Speck Count")
    ax.yaxis.labelpad = padding[0]

    # ax2 = ax.twinx()
    sns.lineplot(
        data=speck_data,
        x="Time (hrs)",
        y="Measurement",
        color=color,
        ax=ax,
        errorbar="se",
        legend=False,
    )

    # if invert == True:
    #     y_name = f"Ratio_{ratio_name_y}_{ratio_name_x}"
    #     clean_name = f"Cytokine Concentration (ng/ml)".replace("b", "β")
    #     ax2.set_ylabel(f"Ratio of Cytokine Measurements, {ratio_name_y}:{ratio_name_x}")
    # else:
    #     y_name = f"Ratio_{ratio_name_x}_{ratio_name_y}"
    #     clean_name = f"Cytokine Concentration (ng/ml)".replace("b", "β")
    #     ax2.set_ylabel(f"Ratio of Cytokine Measurements, {ratio_name_x}:{ratio_name_y}")
    palette = custom_palette(treatment_indeces[treatment])[1:]
    # sns.lineplot(
    #     data=sub_cyto,
    #     x="Time (hrs)",
    #     y="Measurement",
    #     hue="Analyte",
    #     ax=ax2,
    #     errorbar="se",
    #     legend=False,
    #     palette=palette,
    #     style="Analyte",  # using style
    #     # dashes=[(2, 2)] * len(sub_cyto["Analyte"].unique()),
    # )

    # ax2.yaxis.labelpad = padding[1]
    # ax2.yaxis.tick_left()
    # ax2.yaxis.set_label_position("left")
    # ax2.lines[0].set_linestyle("dotted")
    # ax2.tick_params(axis="y", direction="in", pad=-27.5, colors=darker_pastel)
    # ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
    # ax2.yaxis.set_major_locator(MaxNLocator(nbins=6))
    # ax2.set_ylabel(clean_name, color=darker_pastel)

    # ax3 = ax.twinx()
    # ax3.spines["right"].set_position(("outward", 0))
    # ax3.set_ylabel("Cytotoxicity", color="red")
    # ax3.tick_params(axis="y", colors="red")
    # Add % symbol to tick values on ax3

    # ax3.yaxis.set_major_formatter(ticker.PercentFormatter())
    # ax3.scatter(
    #     subset["Time (h)"],
    #     subset["cytotoxicity"], color="red", label=treatment, s=12
    # )

    # sns.scatterplot(
    #     data=subset,
    #     x="Time (h)",
    #     y="cytotoxicity",
    #     color="red",
    #     ax=ax3,
    #     #    errorbar="se",
    #     legend=False,
    # )
    # sns.scatterplot(
    #     data=cyto_means,
    #     x="Time (h)",
    #     y="cytotoxicity",
    #     color="black",
    #     marker="_",
    #     s=500,
    #     ax=ax3,
    # )
    # for bio_rep in subset["bio. Rep."].unique():
    #     for tech_rep in subset["tech. Rep"].unique():
    #         data_points = subset[
    #             (subset["bio. Rep."] == bio_rep) & (subset["tech. Rep"] == tech_rep)
    #         ].sort_values(by="Time (h)")
    #         ax3.plot(
    #             data_points["Time (h)"],
    #             data_points["cytotoxicity"],
    #             color=color,
    #             linestyle="--",
    #             alpha=0.25,
    #         )
    ax.set_ylim(
        bottom=ax.get_ylim()[0], top=ax.get_ylim()[1] * 1.1
    )  # Increase the upper limit by 10%
    # ax2.set_ylim(
    #     bottom=0, top=ax2.get_ylim()[1] * 1.1
    # )  # Increase the upper limit by 10%
    # ax3.set_ylim(
    #     bottom=ax3.get_ylim()[0], top=ax3.get_ylim()[1] * 1.1
    # )  # Increase the upper limit by 10%
    # # After all your plotting and after setting the limits
    # ticks = ax2.get_yticks()
    # ticks = ticks[ticks != 0]  # remove 0 from ticks
    # ticks = ticks[ticks < max(ticks)]  # this removes the highest tick
    # ax2.set_yticks(ticks)  # set the new ticks
    legend_palette = custom_palette(treatment_indeces[treatment])
    # legend_elements = [
    #     Line2D(
    #         [0], [0], color=legend_palette[0], lw=2, label="Normalized Speck Formation"
    #     ),
    #     Line2D(
    #         [0],
    #         [0],
    #         color=legend_palette[1],
    #         lw=2,
    #         label="IL-18 Measurement",
    #         linestyle="dotted",
    #     ),
    #     Line2D(
    #         [0],
    #         [0],
    #         color=legend_palette[2],
    #         lw=2,
    #         label="IL-1β Measurement",
    #         linestyle="--",
    #     ),
    #     Line2D(
    #         [0],
    #         [0],
    #         marker="o",
    #         color="w",
    #         markerfacecolor="red",
    #         markersize=8,
    #         label="Cytotoxicity",
    #         linestyle="None",
    #     ),
    # ]

    # ax.legend(
    #     handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, 1), ncol=3
    # )
    plt.xlim(-1.25, 21.2)
    plt.xlabel("Time (hrs)")
    if filepath is not None:
        plt.savefig(filepath, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()
        return fig, ax


def plot_thing_2(
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
    # change cytotoxicity to a percentage
    subset["cytotoxicity"] = subset["cytotoxicity"] * 100
    cyto_means = subset.groupby(["Time (h)"])["cytotoxicity"].mean().reset_index()

    ratio_name_x = ratio_name.split(":")[0]
    ratio_name_y = ratio_name.split(":")[1]

    def get_color(treatment, palette):
        colors = {"ATP": palette[0], "MSU": palette[1], "Nigericin": palette[2]}
        return colors[treatment]

    treatment_indeces = {"ATP": 0, "MSU": 1, "Nigericin": 2}

    color = get_color(treatment, sns.color_palette())
    pastel_color = get_color(treatment, sns.color_palette("pastel"))
    darker_pastel = tuple(max(0, x - 0.2) for x in pastel_color)

    TAS = analysis_module
    Cyto = TAS.modules["TS_Cyto"].data
    sub_cyto = Cyto[Cyto["Treatment"] == treatment]
    speck_info = TAS.modules["TS_Speck"]
    ratio_info = TAS.modules["TS_Cyto"].ratio_data[ratio_name]

    speck_data = speck_info.data[speck_info.data["Treatment"].isin(treatments)]
    ratio_data = ratio_info[ratio_info["Treatment"].isin(treatments)]

    fig, ax2 = plt.subplots(figsize=(10, 6))
    # ax.set_ylabel("Speck Count")
    # ax.yaxis.labelpad = padding[0]

    # ax2 = ax.twinx()
    # sns.lineplot(
    #     data=speck_data,
    #     x="Time (hrs)",
    #     y="Measurement",
    #     color=color,
    #     ax=ax,
    #     errorbar="se",
    #     legend=False,
    # )

    if invert == True:
        y_name = f"Ratio_{ratio_name_y}_{ratio_name_x}"
        clean_name = f"Cytokine Concentration (ng/ml)".replace("b", "β")
        ax2.set_ylabel(f"Ratio of Cytokine Measurements, {ratio_name_y}:{ratio_name_x}")
    else:
        y_name = f"Ratio_{ratio_name_x}_{ratio_name_y}"
        clean_name = f"Cytokine Concentration (ng/ml)".replace("b", "β")
        ax2.set_ylabel(f"Ratio of Cytokine Measurements, {ratio_name_x}:{ratio_name_y}")
    palette = custom_palette(treatment_indeces[treatment])[1:]
    sns.lineplot(
        data=sub_cyto,
        x="Time (hrs)",
        y="Measurement",
        hue="Analyte",
        ax=ax2,
        errorbar="se",
        legend=False,
        palette=palette,
        style="Analyte",  # using style
        # dashes=[(2, 2)] * len(sub_cyto["Analyte"].unique()),
    )

    ax2.yaxis.labelpad = padding[1]
    ax2.yaxis.tick_left()
    ax2.yaxis.set_label_position("left")
    ax2.lines[0].set_linestyle("dotted")
    ax2.tick_params(axis="y", direction="in", pad=-27.5, colors=darker_pastel)
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax2.set_ylabel(clean_name, color=darker_pastel)

    # ax3 = ax.twinx()
    # ax3.spines["right"].set_position(("outward", 0))
    # ax3.set_ylabel("Cytotoxicity", color="red")
    # ax3.tick_params(axis="y", colors="red")
    # Add % symbol to tick values on ax3

    # ax3.yaxis.set_major_formatter(ticker.PercentFormatter())
    # ax3.scatter(
    #     subset["Time (h)"],
    #     subset["cytotoxicity"], color="red", label=treatment, s=12
    # )

    # sns.scatterplot(
    #     data=subset,
    #     x="Time (h)",
    #     y="cytotoxicity",
    #     color="red",
    #     ax=ax3,
    #     #    errorbar="se",
    #     legend=False,
    # )
    # sns.scatterplot(
    #     data=cyto_means,
    #     x="Time (h)",
    #     y="cytotoxicity",
    #     color="black",
    #     marker="_",
    #     s=500,
    #     ax=ax3,
    # )
    # for bio_rep in subset["bio. Rep."].unique():
    #     for tech_rep in subset["tech. Rep"].unique():
    #         data_points = subset[
    #             (subset["bio. Rep."] == bio_rep) & (subset["tech. Rep"] == tech_rep)
    #         ].sort_values(by="Time (h)")
    #         ax3.plot(
    #             data_points["Time (h)"],
    #             data_points["cytotoxicity"],
    #             color=color,
    #             linestyle="--",
    #             alpha=0.25,
    #         )
    # ax.set_ylim(
    #     bottom=ax.get_ylim()[0], top=ax.get_ylim()[1] * 1.1
    # )  # Increase the upper limit by 10%
    ax2.set_ylim(
        bottom=0, top=ax2.get_ylim()[1] * 1.1
    )  # Increase the upper limit by 10%
    # ax3.set_ylim(
    #     bottom=ax3.get_ylim()[0], top=ax3.get_ylim()[1] * 1.1
    # )  # Increase the upper limit by 10%
    # After all your plotting and after setting the limits
    ticks = ax2.get_yticks()
    ticks = ticks[ticks != 0]  # remove 0 from ticks
    ticks = ticks[ticks < max(ticks)]  # this removes the highest tick
    ax2.set_yticks(ticks)  # set the new ticks
    legend_palette = custom_palette(treatment_indeces[treatment])
    # legend_elements = [
    #     Line2D(
    #         [0], [0], color=legend_palette[0], lw=2, label="Normalized Speck Formation"
    #     ),
    #     Line2D(
    #         [0],
    #         [0],
    #         color=legend_palette[1],
    #         lw=2,
    #         label="IL-18 Measurement",
    #         linestyle="dotted",
    #     ),
    #     Line2D(
    #         [0],
    #         [0],
    #         color=legend_palette[2],
    #         lw=2,
    #         label="IL-1β Measurement",
    #         linestyle="--",
    #     ),
    #     Line2D(
    #         [0],
    #         [0],
    #         marker="o",
    #         color="w",
    #         markerfacecolor="red",
    #         markersize=8,
    #         label="Cytotoxicity",
    #         linestyle="None",
    #     ),
    # ]

    # ax.legend(
    #     handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, 1), ncol=3
    # )
    plt.xlim(-1.25, 21.2)
    plt.xlabel("Time (hrs)")
    if filepath is not None:
        plt.savefig(filepath, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()
        return fig, ax2


def plot_thing_3(
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
    # change cytotoxicity to a percentage
    subset["cytotoxicity"] = subset["cytotoxicity"] * 100
    cyto_means = subset.groupby(["Time (h)"])["cytotoxicity"].mean().reset_index()

    ratio_name_x = ratio_name.split(":")[0]
    ratio_name_y = ratio_name.split(":")[1]

    def get_color(treatment, palette):
        colors = {"ATP": palette[0], "MSU": palette[1], "Nigericin": palette[2]}
        return colors[treatment]

    treatment_indeces = {"ATP": 0, "MSU": 1, "Nigericin": 2}

    color = get_color(treatment, sns.color_palette())
    pastel_color = get_color(treatment, sns.color_palette("pastel"))
    darker_pastel = tuple(max(0, x - 0.2) for x in pastel_color)

    TAS = analysis_module
    Cyto = TAS.modules["TS_Cyto"].data
    sub_cyto = Cyto[Cyto["Treatment"] == treatment]
    speck_info = TAS.modules["TS_Speck"]
    ratio_info = TAS.modules["TS_Cyto"].ratio_data[ratio_name]

    speck_data = speck_info.data[speck_info.data["Treatment"].isin(treatments)]
    ratio_data = ratio_info[ratio_info["Treatment"].isin(treatments)]

    fig, ax3 = plt.subplots(figsize=(10, 6))
    # ax.set_ylabel("Speck Count")
    # ax.yaxis.labelpad = padding[0]

    # ax2 = ax.twinx()
    # sns.lineplot(
    #     data=speck_data,
    #     x="Time (hrs)",
    #     y="Measurement",
    #     color=color,
    #     ax=ax,
    #     errorbar="se",
    #     legend=False,
    # )

    # if invert == True:
    #     y_name = f"Ratio_{ratio_name_y}_{ratio_name_x}"
    #     clean_name = f"Cytokine Concentration (ng/ml)".replace("b", "β")
    #     ax2.set_ylabel(f"Ratio of Cytokine Measurements, {ratio_name_y}:{ratio_name_x}")
    # else:
    #     y_name = f"Ratio_{ratio_name_x}_{ratio_name_y}"
    #     clean_name = f"Cytokine Concentration (ng/ml)".replace("b", "β")
    #     ax2.set_ylabel(f"Ratio of Cytokine Measurements, {ratio_name_x}:{ratio_name_y}")
    palette = custom_palette(treatment_indeces[treatment])[1:]
    # sns.lineplot(
    #     data=sub_cyto,
    #     x="Time (hrs)",
    #     y="Measurement",
    #     hue="Analyte",
    #     ax=ax2,
    #     errorbar="se",
    #     legend=False,
    #     palette=palette,
    #     style="Analyte",  # using style
    #     # dashes=[(2, 2)] * len(sub_cyto["Analyte"].unique()),
    # )

    # ax2.yaxis.labelpad = padding[1]
    # ax2.yaxis.tick_left()
    # ax2.yaxis.set_label_position("left")
    # ax2.lines[0].set_linestyle("dotted")
    # ax2.tick_params(axis="y", direction="in", pad=-27.5, colors=darker_pastel)
    # ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
    # ax2.yaxis.set_major_locator(MaxNLocator(nbins=6))
    # ax2.set_ylabel(clean_name, color=darker_pastel)

    # ax3 = ax.twinx()
    ax3.spines["right"].set_position(("outward", 0))
    ax3.set_ylabel("Cytotoxicity", color="red")
    ax3.tick_params(axis="y", colors="red")
    # Add % symbol to tick values on ax3

    ax3.yaxis.set_major_formatter(ticker.PercentFormatter())
    # ax3.scatter(
    #     subset["Time (h)"],
    #     subset["cytotoxicity"], color="red", label=treatment, s=12
    # )

    sns.scatterplot(
        data=subset,
        x="Time (h)",
        y="cytotoxicity",
        color="red",
        ax=ax3,
        #    errorbar="se",
        legend=False,
    )
    sns.scatterplot(
        data=cyto_means,
        x="Time (h)",
        y="cytotoxicity",
        color="black",
        marker="_",
        s=500,
        ax=ax3,
    )
    # for bio_rep in subset["bio. Rep."].unique():
    #     for tech_rep in subset["tech. Rep"].unique():
    #         data_points = subset[
    #             (subset["bio. Rep."] == bio_rep) & (subset["tech. Rep"] == tech_rep)
    #         ].sort_values(by="Time (h)")
    #         ax3.plot(
    #             data_points["Time (h)"],
    #             data_points["cytotoxicity"],
    #             color=color,
    #             linestyle="--",
    #             alpha=0.25,
    #         )
    # ax.set_ylim(
    #     bottom=ax.get_ylim()[0], top=ax.get_ylim()[1] * 1.1
    # )  # Increase the upper limit by 10%
    # ax2.set_ylim(
    #     bottom=0, top=ax2.get_ylim()[1] * 1.1
    # )  # Increase the upper limit by 10%
    ax3.set_ylim(
        bottom=ax3.get_ylim()[0], top=ax3.get_ylim()[1] * 1.1
    )  # Increase the upper limit by 10%
    # After all your plotting and after setting the limits
    # ticks = ax2.get_yticks()
    # ticks = ticks[ticks != 0]  # remove 0 from ticks
    # ticks = ticks[ticks < max(ticks)]  # this removes the highest tick
    # ax2.set_yticks(ticks)  # set the new ticks
    # legend_palette = custom_palette(treatment_indeces[treatment])
    # legend_elements = [
    #     Line2D(
    #         [0], [0], color=legend_palette[0], lw=2, label="Normalized Speck Formation"
    #     ),
    #     Line2D(
    #         [0],
    #         [0],
    #         color=legend_palette[1],
    #         lw=2,
    #         label="IL-18 Measurement",
    #         linestyle="dotted",
    #     ),
    #     Line2D(
    #         [0],
    #         [0],
    #         color=legend_palette[2],
    #         lw=2,
    #         label="IL-1β Measurement",
    #         linestyle="--",
    #     ),
    #     Line2D(
    #         [0],
    #         [0],
    #         marker="o",
    #         color="w",
    #         markerfacecolor="red",
    #         markersize=8,
    #         label="Cytotoxicity",
    #         linestyle="None",
    #     ),
    # ]

    # ax.legend(
    #     handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, 1), ncol=3
    # )
    plt.xlim(-1.25, 21.2)
    plt.xlabel("Time (hrs)")
    if filepath is not None:
        plt.savefig(filepath, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()
        return fig, ax3


#
for treatment in ["ATP", "MSU", "Nigericin"]:
    for function in [plot_thing_1, plot_thing_2, plot_thing_3]:
        print(function.__name__)
        function(
            TAS,
            LDH,
            ratio_name="IL18:IL1b",
            treatments=[treatment],
            padding=(20, 30),
            invert=True,
            filepath=Path(
                "plots", "broken_apart", f"alt1_{treatment}_{function.__name__}.png"
            ),
        )
plot_everything(
    analysis_module=TAS,
    ldh_module=LDH,
    ratio_name="IL18:IL1b",
    treatments=["MSU"],
    padding=(20, 40),
    invert=True,
    filepath=Path("plots", "alt1_msu_ratio_speck_ldh_combined.png"),
)
plot_everything(
    TAS,
    LDH,
    ratio_name="IL18:IL1b",
    treatments=["Nigericin"],
    padding=(20, 40),
    invert=True,
    filepath=Path("plots", "alt1_nigericin_ratio_speck_ldh_combined.png"),
)


sns.scatterplot
