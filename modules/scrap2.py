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
    Cyto = TAS.modules["TS_Cyto"].data
    sub_cyto = Cyto[Cyto["Treatment"] == treatment]
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
    analysis_module=TAS,
    ldh_module=LDH,
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
