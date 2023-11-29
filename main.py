from modules import *
from pathlib import Path

## Basic Parameters
base_treatments = ["ATP", "MSU", "Nigericin"]
modifier = "MCC950"
short_name_dict = {"nig": "Nigericin"}
time_limits = {"ALL": 21}


def prepare_analysis_class():
    ## Build the separate class objects for the different tests
    cyto = time_sequence_module(
        Cytokine_Adapter("Cyto", base_treatments, short_name_dict=short_name_dict)
    )
    cyto.structure_data("Concentration")

    speck = time_sequence_module(
        Count_Aggregator(
            "Speck", base_treatments, modifier=modifier, short_name_dict=short_name_dict
        )
    )
    speck.structure_data("SpeckFormation")
    speck.data = speck.limit_data(speck.data, time_limits=time_limits)

    # Removal of outlier data
    speck.data = speck.data[
        ~(
            (speck.data["Treatment"].str.contains("MCC950"))
            & (speck.data["Experimental_Replicate"] == "1")
        )
    ]
    speck.data = speck.data[
        ~(
            (speck.data["Treatment"] == "MSU")
            & (speck.data["Experimental_Replicate"] == "1")
        )
    ]
    speck.data = speck.data[
        ~(
            (speck.data["Treatment"] == "MSU")
            & (speck.data["Experimental_Replicate"] == "2")
        )
    ]

    ##Construct main analysis class object
    TAS = analysis_module([cyto, speck])
    TAS.time_compare()
    TAS.compute_ratio("TS_Cyto")
    TAS.time_compare_dict()

    return TAS


def generate_reports(TAS):
    ##REPORTS SECTION
    # check if reports directory exists, if not create it
    Path("reports").mkdir(parents=True, exist_ok=True)
    ## Create Summaries
    for module in TAS.modules:
        (TAS.create_summary_table(TAS.modules[module])).to_excel(
            Path("reports", "".join(["summary_", module, ".xlsx"])), index=False
        )
    ## Create Modifier Impact Reports
    for module in [
        module for module in TAS.modules if TAS.modules[module].comp.modifier
    ]:
        TAS.check_modifier_impact(TAS.modules[module]).T.to_excel(
            Path("reports", "".join(["modifier_impact_", module, ".xlsx"])),
            header=False,
        )
    ## Create Time Point Comparisons
    for module in TAS.modules:
        TAS.aggregate_time_comparisons(TAS.modules[module])


def generate_plots(TAS, file_extension=".tiff"):
    cyto = TAS.modules["TS_Cyto"]
    module = TAS.modules["TS_Speck"]
    # check that plots directory exists and create it if not
    Path("plots").mkdir(parents=True, exist_ok=True)

    def normalized_speck_counts_with_MCC950():
        ## Normalized Speck Counts with MCC950
        for treatment in module.comp.treatments:
            file_title = str(f"Normalized ASC-Speck Count with MCC950 - {treatment}")

            def ax_modification(ax):
                ax.set_title(file_title)
                ax.set_ylabel("Normalized ASC-Speck Count")
                # ymin, ymax = ax.get_ylim()
                # ax.set_ylim(ymin, ymax * 1.2)

            plotting_module.plot_lineplot(
                **TAS.prepare_lineplot(
                    TAS.modules["TS_Speck"],
                    treatments=[treatment, f"MCC950_{treatment}"],
                    measurement_type="Normalized_Measurement",
                ),
                manual_ax_modification=ax_modification,
                filepath=Path(
                    "plots",
                    "".join(
                        [file_title.replace("-", "").replace(" ", "_"), file_extension]
                    ),
                ),
            )

    def speck_curves_raw_counts():
        ##Speck Curves with Raw Counts
        file_title = str("Inflammasome Speck Count over Time")
        treatment_sets = ["all", ["ATP", "MSU", "Nigericin"], ["ATP"]]

        def ax_modification(ax):
            ax.set_title(file_title)
            ax.set_ylabel("ASC-Speck Count")

        for treatment_set in treatment_sets:
            plotting_module.plot_lineplot(
                **TAS.prepare_lineplot(
                    TAS.modules["TS_Speck"],
                    treatments=treatment_set,
                    measurement_type="Measurement",
                ),
                manual_ax_modification=ax_modification,
                filepath=Path(
                    "plots",
                    "".join(
                        [
                            file_title.replace("-", "").replace(" ", "_"),
                            "_",
                            str(treatment_set),
                            file_extension,
                        ]
                    ),
                ),
            )

    def cytokine_plots_all_treatments_single_analyte():
        ## CYTOKINES PLOTS, ALL TREATMENTS, SINGLE ANALYTE
        for analyte in cyto.data["Analyte"].unique():
            if analyte == "IL1b":
                aname = "IL-1β"
            elif analyte == "IL18":
                aname = "IL-18"

            def ax_modification(ax):
                ax.set_title(str(aname + " Concentration over Time"))
                ax.set_ylabel("Concentration (pg/mL)")
                ymin, ymax = ax.get_ylim()
                ax.set_ylim(ymin, ymax * 1.2)

            plotting_module.plot_lineplot(
                **TAS.prepare_lineplot(
                    TAS.modules["TS_Cyto"],
                    treatments="all",
                    measurement_type="Measurement",
                    analyte_selection=[analyte],
                ),
                manual_ax_modification=ax_modification,
                filepath=Path(
                    "plots",
                    "".join(
                        [
                            aname.replace("-", "").replace(" ", "_"),
                            "_concentration_over_time",
                            file_extension,
                        ]
                    ),
                ),
            )

    def ratio_IL18_IL1b():
        file_title = str("Ratio of IL-1β:IL-18 Concentration over Time")

        ## SPECK COUNT, CYTO RATIO and ALL TREATMENTS
        def ax_modification(
            ax,
        ):
            ax.set_title(file_title)
            ax.set_ylabel("IL-1β:IL-18 Concentration Ratio")

        plotting_module.plot_ratio(
            module=TAS.modules["TS_Cyto"],
            invert=True,
            manual_ax_modification=ax_modification,
            filepath=Path(
                "plots",
                "".join(
                    [
                        file_title.replace("-", "").replace(" ", "_").replace(":", "-"),
                        file_extension,
                    ]
                ),
            ),
        )

    def ratio_IL18_IL1b_with_specks():
        file_title = str(
            "Ratio of IL-1β:IL-18 Concentration over Time with Speck Count"
        )

        def ax_modification(ax, ax2):
            ax.set_title("Ratio of IL-1β:IL-18 Concentration over Time")
            ax.set_ylabel("Speck Formation")
            ymin, ymax = ax.get_ylim()
            ax.set_ylim(ymin, ymax * 1.25)
            ymin, ymax = ax2.get_ylim()
            ax2.set_ylim(ymin, ymax * 1.25)
            ax2.set_ylabel("IL-1β:IL-18 Concentration Ratio")

        plotting_module.plot_count_against_ratio(
            TAS,
            ratio_name="IL18:IL1b",
            invert=True,
            manual_ax_modification=ax_modification,
            filepath=Path(
                "plots",
                "".join(
                    [
                        file_title.replace("-", "").replace(" ", "_").replace(":", "-"),
                        file_extension,
                    ]
                ),
            ),
        )

    def ratio_IL18_IL1b_with_specks_individual():
        # SAME, BUT FOR EACH TREATMENT
        for treatment in cyto.data["Treatment"].unique():
            file_title = f"Ratio of IL-1β:IL-18 Concentration over Time - {treatment}"

            def ax_modification(ax, ax2):
                ax.set_title(file_title)
                ax.set_ylabel("Speck Formation")
                ymin, ymax = ax.get_ylim()
                ax.set_ylim(ymin, ymax * 1.15)
                ymin, ymax = ax2.get_ylim()
                ax2.set_ylim(ymin, ymax * 1.15)
                ax2.set_ylabel("IL-1β:IL-18 Concentration Ratio")

            plotting_module.plot_count_against_ratio(
                TAS,
                ratio_name="IL18:IL1b",
                treatments=[treatment],
                invert=True,
                manual_ax_modification=ax_modification,
                filepath=Path(
                    "plots",
                    "".join(
                        [
                            file_title.replace("-", "")
                            .replace(" ", "_")
                            .replace(":", "-"),
                            file_extension,
                        ]
                    ),
                ),
            )

    def speck_plots_individual():
        ## Line Plots of Speck data for only each MCC950 treatment
        for treatment in [
            "MCC950_ATP",
            "MCC950_MSU",
            "MCC950_Nigericin",
            "ATP",
            "MSU",
            "Nigericin",
        ]:
            treatment_name = treatment.replace("_", " + ")
            file_title = str(f"Inflammasome Speck Count over Time - {treatment_name}")

            def ax_modification(ax):
                ax.set_title(file_title)
                ymin, ymax = ax.get_ylim()
                ax.set_ylim(ymin, ymax * 1.1)
                # set y axis label
                ax.set_ylabel("ASC-Speck Formation")
                handles, labels = ax.get_legend_handles_labels()
                new_labels = [label.replace("_", " + ") for label in labels]
                ax.legend(handles, new_labels)

            plotting_module.plot_lineplot(
                **TAS.prepare_lineplot(
                    TAS.modules["TS_Speck"],
                    treatments=[treatment],
                    measurement_type="Measurement",
                ),
                manual_ax_modification=ax_modification,
                filepath=Path(
                    "plots",
                    "".join(
                        [file_title.replace("-", "").replace(" ", "_"), file_extension]
                    ),
                ),
            )

    ### RUN ALL FUNCTIONS
    ratio_IL18_IL1b()
    ratio_IL18_IL1b_with_specks()
    ratio_IL18_IL1b_with_specks_individual()
    cytokine_plots_all_treatments_single_analyte()
    speck_curves_raw_counts()
    speck_plots_individual()
    normalized_speck_counts_with_MCC950()


def generate_plots_ldh(TAS, file_extension=".tiff"):
    plot_folder = Path("plots", "ldh")
    # if plot folder doesn't exist, create it
    plot_folder.mkdir(parents=True, exist_ok=True)
    for treatment in TAS.modules["LDH"].data["Treatment"].unique():
        plotting_module.plot_speck_count(
            analysis_class=TAS,
            treatments=[treatment],
            filepath=f"{plot_folder}/{treatment}_speck_count{file_extension}",
        )
        plotting_module.plot_cytokines(
            analysis_class=TAS,
            treatments=[treatment],
            filepath=f"{plot_folder}/{treatment}_cytokines{file_extension}",
            padding=(-5, 10),
        )
        plotting_module.plot_cytotoxicity(
            analysis_class=TAS,
            treatments=[treatment],
            filepath=f"{plot_folder}/{treatment}_cytotoxicity{file_extension}",
        )
        plotting_module.combined_plot(
            analysis_class=TAS,
            treatments=[treatment],
            filepath=f"{plot_folder}/{treatment}_combined{file_extension}",
            padding=(-5, 50),
        )


def main():
    # check that all requirements in requirements.txt are installed
    TAS = prepare_analysis_class()
    generate_reports(TAS)
    generate_plots(TAS)
    TAS.modules["LDH"] = LDHModule()
    generate_plots_ldh(TAS)
    generate_plots_ldh(TAS, file_extension=".png")
    # generate_plots(TAS, file_extension=".png")
