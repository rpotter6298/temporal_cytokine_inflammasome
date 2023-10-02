
treatment = 'ATP'
color = colors['ATP']
def plot_everything(analysis_module,ldh_module,ratio_name,treatments=None,invert=False,manual_ax_modification=None,filepath=None,):
    treatment = treatments[0]
    def get_color(treatment, palette):
        colors = {
            'ATP': palette[0],
            'MSU': palette[1],
            'Nigericin': palette[2]
        }
        return colors[treatment]
    color = get_color(treatment, sns.color_palette(None))
    pastel_color = get_color(treatment, sns.color_palette("pastel"))

    
    TAS = analysis_module
    LDH = ldh_module
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
    ax.yaxis.label.set_color("blue")
    # ax.yaxis.label.set_position((0, 0))
    ax2 = ax.twinx()
    lineplot_ax1 = sns.lineplot(
        data=speck_data,
        x="Time (hrs)",
        y="Measurement",
        hue="Treatment",
        ax=ax,
        errorbar="se",
        legend=False,
    )
    if invert == True:
        y_name = f"Ratio_{ratio_name_y}_{ratio_name_x}"
        ax2.set_ylabel(
            f"Ratio of Cytokine Measurements, {ratio_name_y}:{ratio_name_x}", color=pastel_color)
    else:
        y_name = f"Ratio_{ratio_name_x}_{ratio_name_y}"
        ax2.set_ylabel(
            f"Ratio of Cytokine Measurements, {ratio_name_x}:{ratio_name_y}", color=pastel_color)

    darker_pastel = tuple(max(0, x- 0.2) for x in pastel_color)
    ax2.tick_params(axis="y", colors=darker_pastel)
    pointplot = sns.pointplot(
        data=ratio_data,
        x="Time (hrs)",
        y=y_name,
        hue="Treatment",
        ax=ax2,
        palette="pastel",
        errorbar="se",
    )
    pointplot.legend_.remove()
    ax3 = ax.twinx()
    ax3.spines['right'].set_position(('outward', 60))  # moving ax3 to the right
    ax3.set_ylabel("Cytotoxicity", color='red')  # setting y-axis label color to red
    ax3.tick_params(axis='y', colors='red')  # setting y-axis tick color to red
    ax3.scatter(subset['Time (h)'], subset['cytotoxicity'], color='red', label=treatment, s=5)
    for bio_rep in subset['bio. Rep.'].unique():
        for tech_rep in subset['tech. Rep'].unique():
            data_points = subset[(subset['bio. Rep.'] == bio_rep) & (subset['tech. Rep'] == tech_rep)].sort_values(by='Time (h)')
            ax3.plot(data_points['Time (h)'], data_points['cytotoxicity'], color=color, linestyle='--', alpha=0.75)    
    
    from matplotlib.lines import Line2D
    def pastel_color(color):
        sns.set_palette(sns.color_palette("pastel"))
        return sns.color_palette()[color]
    legend_elements = [Line2D([0], [0], color='b', lw=2, label='Normalized Speck Formation'),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor=get_color(treatment, sns.color_palette("pastel")), markersize=10, label=f'Ratio {ratio_name}', linestyle='None'),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=5, label='Cytotoxicity', linestyle='None')]

    ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, .95), ncol=3)
    plt.xlim(0, 21.2)
    plt.xlabel("Time (hrs)")

# Call your function with the required arguments
plot_everything(TAS, LDH, ratio_name="IL18:IL1b", treatments=["ATP"])
plt.show()