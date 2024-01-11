from modules import *
from pathlib import Path
import pandas as pd

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
                ax.set_title(f"Extracellular {aname} Concentration")
                ax.set_ylabel("Concentration (pg/mL)")
                ymin, ymax = ax.get_ylim()
                ax.set_ylim(ymin, ymax * 1.2)
                ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
                ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
                ax.set_xlim(0, 24)

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
                        replace("β", "b"),
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
    #ratio_IL18_IL1b()
    #ratio_IL18_IL1b_with_specks()
    #ratio_IL18_IL1b_with_specks_individual()
    cytokine_plots_all_treatments_single_analyte()
    #speck_curves_raw_counts()
    #speck_plots_individual()
    #normalized_speck_counts_with_MCC950()


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


def build_correlation_df(
    analysis_module,
    parameter="Normalized_Measurement",
    ratio="IL18:IL1b",
    invert_ratio=True,
):
    """
    Builds a DataFrame containing the correlation between all cytokines and speck counts for each treatment.

    :param analysis_module: analysis_module object containing the data to be analyzed.
    :param parameter: parameter to be used for the correlation. Default is "Normalized_Measurement".
    :return: DataFrame containing the correlation between all cytokines and speck counts for each treatment.
    """
    cytokines = analysis_module.modules["TS_Cyto"].data["Analyte"].unique()
    if invert_ratio == True:
        ratio_name = ratio.split(":")[::-1]
        ratio_name = "_".join(ratio_name)
    else:
        ratio_name = ratio.split(":")
        ratio_name = "_".join(ratio_name)
    ratios = analysis_module.modules["TS_Cyto"].ratio_data[ratio][
        [
            "Analyte_x",
            "Time (hrs)",
            "Treatment",
            "Experimental_Replicate",
            f"Ratio_{ratio_name}",
        ]
    ]
    ratios.columns = [
        "Analyte",
        "Time (hrs)",
        "Treatment",
        "Experimental_Replicate",
        parameter,
    ]
    treatments = analysis_module.modules["TS_Cyto"].comp.treatments
    overall_corr_df = pd.DataFrame()
    time_corr_df = pd.DataFrame()

    def build_corr_df(corr, cytokine, treatment):
        temp_df = pd.DataFrame(
            {
                "Treatment": [treatment],
                "Cytokine": [cytokine],
                "Correlation": [corr["Correlation - Overall"]],
                "p-value": [corr["p-value - Overall"]],
                "Cosine Similarity": [corr["Cosine Similarity"]],
            }
        )
        return temp_df

    def build_time_corr_df(corr, cytokine, treatment):
        if type(cytokine) == tuple:
            cytokine = "_".join(cytokine)
        else:
            cytokine = "Speck Count - " + cytokine
        tcorr = pd.DataFrame(corr["Time Correlation"]).T

        temp_df = pd.DataFrame(
            {
                "Time (hrs)": tcorr.index,
                "Treatment": [treatment] * len(tcorr),
                "Cytokine": [cytokine] * len(tcorr),
                "Correlation": tcorr["Correlation"],
                "p-value": tcorr["p-value"],
            }
        )
        return temp_df

    ## Speck: Cytokine Correlation for each treatment
    for cytokine in cytokines:
        for treatment in treatments:
            subset_1, subset_2, corr = analysis_module.correlation(
                dataset_1=("TS_Cyto", cytokine, [treatment]),
                dataset_2=("TS_Speck", "NA", [treatment]),
                parameter=parameter,
            )
            overall_corr_df = pd.concat(
                [overall_corr_df, build_corr_df(corr, cytokine, treatment)]
            )
            time_corr_df = pd.concat(
                [time_corr_df, build_time_corr_df(corr, cytokine, treatment)]
            )
            # temp_df = pd.DataFrame(
            #     {
            #         "Treatment": [treatment],
            #         "Cytokine": [cytokine],
            #         "Correlation": [corr["Correlation - Overall"]],
            #         "p-value": [corr["p-value - Overall"]],
            #         "Cosine Similarity": [corr["Cosine Similarity"]],
            #     }
            # )
            # overall_corr_df = pd.concat([overall_corr_df, temp_df])
            # tcorr = pd.DataFrame(corr["Time Correlation"]).T
            # temp_df = pd.DataFrame(
            #     {
            #         "Time (hrs)": tcorr.index,
            #         "Treatment": [treatment] * len(tcorr),
            #         "Cytokine": [cytokine] * len(tcorr),
            #         "Correlation": tcorr["Correlation"],
            #         "p-value": tcorr["p-value"],
            #     }
            # )
            # time_corr_df = pd.concat([time_corr_df, temp_df])

    # for treatment in treatments:
    #     ## Speck: Cytokine Ratio Correlation for each treatment
    #     subset_1, subset_2, corr = analysis_module.correlation(
    #         dataset_1=("TS_Cyto", ratios["Analyte"].unique()[0], [treatment]),
    #         dataset_2=("TS_Speck", "NA", [treatment]),
    #         parameter=parameter,
    #     )
    #     overall_corr_df = pd.concat([overall_corr_df, build_corr_df(corr, ratio, treatment)])
    #     time_corr_df = pd.concat([time_corr_df, build_time_corr_df(corr, ratio, treatment)])
    ## Cytokine: Cytokine Correlation for each treatment
    # # create pairs of all cytokines to create correlation for
    # import itertools

    # cytokine_pairs = list(itertools.combinations(cytokines, 2))
    # for pair in cytokine_pairs:
    #     subset_1, subset_2, corr = analysis_module.correlation(
    #         dataset_1=("TS_Cyto", pair[0], [treatment]),
    #         dataset_2=("TS_Cyto", pair[1], [treatment]),
    #         parameter=parameter,
    #     )
    #     overall_corr_df = pd.concat([overall_corr_df, build_corr_df(corr, pair, treatment)])
    #     time_corr_df = pd.concat([time_corr_df, build_time_corr_df(corr, pair, treatment)])

    overall_corr_df = overall_corr_df.reset_index(drop=True)
    time_corr_df = time_corr_df.reset_index(drop=True)

    return overall_corr_df, time_corr_df


def generate_plots_changes(TAS, file_extension=".tiff"):
    file_dir = Path("plots", "change_plots")
    file_dir.mkdir(parents=True, exist_ok=True)
    # This part makes the change plots for the cytokines
    for treatment in ["ATP", "MSU", "Nigericin"]:
        for normalization in [False]:
            for method in ["Delta", "Acceleration", "Smoothed Acceleration"]:
                plotting_module.change_plot(
                    TAS.modules["TS_Cyto"],
                    treatments=[treatment],
                    mode=method,
                    filepath=f"{file_dir}/{treatment}_{method}_Cytokines.png",
                    normalized=normalization,
                )
                plotting_module.change_plot(
                    TAS.modules["TS_Speck"],
                    treatments=[treatment],
                    mode=method,
                    filepath=f"{file_dir}/{treatment}_{method}_Speck.png",
                    normalized=normalization,
                )
                for analyte in ["IL1b", "IL18"]:
                    plotting_module.change_plot(
                        TAS.modules["TS_Cyto"],
                        treatments=[treatment],
                        analytes=[analyte],
                        mode=method,
                        filepath=f"{file_dir}/{treatment}_{method}_{analyte}.png",
                        normalized=normalization,
                    )
    for method in ["Delta", "Acceleration", "Smoothed Acceleration"]:
        plotting_module.change_plot(
            TAS.modules["TS_Speck"],
            treatments=["ATP", "MSU", "Nigericin"],
            mode=method,
            filepath=f"{file_dir}/Speck_{method}.png",
        )


def check_time_point(treatment, analyte, time, method):
    import scipy

    example_1 = TAS.modules["TS_Cyto"].data[
        (TAS.modules["TS_Cyto"].data["Treatment"] == treatment)
        & (TAS.modules["TS_Cyto"].data["Analyte"] == analyte)
    ]
    example_2 = TAS.modules["TS_Speck"].data[
        (TAS.modules["TS_Speck"].data["Treatment"] == treatment)
    ]
    parameter = method
    # drop all columns except for time, experimental replicate, and the parameter
    subset_1 = example_1[["Time (hrs)", "Experimental_Replicate", parameter]]
    subset_2 = example_2[["Time (hrs)", "Experimental_Replicate", parameter]]

    def balance_times(subset_1, subset_2):
        time_points_1 = subset_1["Time (hrs)"].unique()
        time_points_2 = subset_2["Time (hrs)"].unique()

        # identify any time points outside the range of the other set
        # if there are any, remove them from the set
        # also remove them from the time points list
        if time_points_1.max() > time_points_2.max():
            subset_1 = subset_1[subset_1["Time (hrs)"] <= time_points_2.max()]
            time_points_1 = subset_1["Time (hrs)"].unique()
        elif time_points_2.max() > time_points_1.max():
            subset_2 = subset_2[subset_2["Time (hrs)"] <= time_points_1.max()]
            time_points_2 = subset_2["Time (hrs)"].unique()

        # if there are any time points that are not in both sets, make a dictionary of the missing times and which set they are missing from
        missing_times = {}
        for time in time_points_1:
            if time not in time_points_2:
                missing_times[time] = "subset_2"
        for time in time_points_2:
            if time not in time_points_1:
                missing_times[time] = "subset_1"
        import pandas as pd

        # for each time that is missing, find the nearest two times from the respective subset and average them to impute these rows
        for time, subset_name in missing_times.items():
            current_subset = subset_1 if subset_name == "subset_1" else subset_2

            # average the nearest two times to the missing time
            imputed_values = current_subset.groupby(["Experimental_Replicate"]).apply(
                lambda x:
                # get the two nearest times
                x.loc[
                    (x["Time (hrs)"] == x["Time (hrs)"].max())
                    | (x["Time (hrs)"] == x["Time (hrs)"].min())
                ][parameter].mean()
            )
            imputed_df = pd.DataFrame(
                {
                    "Time (hrs)": [time] * len(imputed_values),
                    "Experimental_Replicate": imputed_values.index,
                    parameter: imputed_values,
                }
            )
            if subset_name == "subset_1":
                subset_1 = pd.concat([subset_1, imputed_df])
            else:
                subset_2 = pd.concat([subset_2, imputed_df])

        # Sorting the subsets by time after imputation
        subset_1 = subset_1.sort_values(
            by=["Time (hrs)", "Experimental_Replicate"]
        ).reset_index(drop=True)
        subset_2 = subset_2.sort_values(
            by=["Time (hrs)", "Experimental_Replicate"]
        ).reset_index(drop=True)

        return subset_1, subset_2

    subset_A, subset_B = balance_times(subset_1, subset_2)

    time_correlation = {}
    # time = 7

    subset_1_time = subset_A[subset_A["Time (hrs)"] == time]
    subset_2_time = subset_B[subset_B["Time (hrs)"] == time]
    # run shapiro wilks test to determine normality
    # if normal, run pearson, if not, run spearman
    normality_test_1 = scipy.stats.shapiro(subset_1_time[parameter])[1]
    normality_test_2 = scipy.stats.shapiro(subset_2_time[parameter])[1]
    if normality_test_1 > 0.05 and normality_test_2 > 0.05:
        corr_type = "pearson"
        corr = scipy.stats.pearsonr(subset_1_time[parameter], subset_2_time[parameter])
    else:
        corr_type = "spearman"
        corr = scipy.stats.spearmanr(subset_1_time[parameter], subset_2_time[parameter])
    time_correlation[time] = {
        "Normality Test - Subset 1": normality_test_1,
        "Normality Test - Subset 2": normality_test_2,
        "Correlation Type": corr_type,
        "Correlation": corr[0],
        "p-value": corr[1],
    }
    print(subset_1_time)
    print(subset_2_time)
    print(time_correlation[time])


# TODO Add the change plots that are missing here


def main():
    # check that all requirements in requirements.txt are installed
    TAS = prepare_analysis_class()
    generate_reports(TAS)
    generate_plots(TAS)
    generate_plots(TAS, file_extension=".png")
    TAS.modules["LDH"] = LDHModule()
    generate_plots_ldh(TAS)
    generate_plots_ldh(TAS, file_extension=".png")
    # generate_plots(TAS, file_extension=".png")
    generate_plots_changes(TAS, file_extension=".png")
    overall, time = build_correlation_df(TAS, parameter="Delta")
    plotting_module.time_heatmap(time, filepath=Path("plots", "Delta_time_heatmap.png"))
    overall, time = build_correlation_df(TAS, parameter="Acceleration")
    plotting_module.time_heatmap(
        time, filepath=Path("plots", "Acceleration_time_heatmap.png")
    )


# module = TAS.modules["TS_Cyto"]
# speck = TAS.modules["TS_Speck"].data

cyto = TAS.modules["TS_Cyto"]
cyto.data
cyto_df = cyto.data[cyto.data["Treatment"] == "MSU"][cyto.data["Time (hrs)"] == 7]

# ###PLAYGROUND
# import matplotlib.pyplot as plt
# import seaborn as sns

check_time_point("ATP", "IL1b", 16.5, "Delta")

TAS.correlation()
# file_dir = Path("plots", "change_plots", "combos")

# for treatment in ["ATP", "MSU", "Nigericin"]:
#     for method in ["Change", "Acceleration"]:
#         plotting_module.combo_change_plot(
#             TAS,
#             treatments=[treatment],
#             mode=method,
#             filepath=f"{file_dir}/{treatment}_{method}_Speck_and_Cytokines.png",
#         )


# plotting_module.combo_change_plot(
#     TAS, treatments=["ATP"], mode="smoothed acceleration", filepath="test.png"
# )
# plotting_module.change_plot(
#     TAS.modules["TS_Cyto"], treatments=["ATP"], mode="Combo", filepath="test.png"
# )


# TAS = prepare_analysis_class()
# subset_1, subset_2, corr = TAS.correlation(
#     dataset_1=("TS_Cyto", "IL1b", ["ATP"]), dataset_2=("TS_Speck", "NA", ["MSU"])
# )
# time_corr_df = pd.DataFrame(corr["Time Correlation"]).T
# time_corr_df["Time ("] = time_corr_df.index

# import pandas as pd
# import matplotlib.patches as patches


# def time_heatmap(time_correlation_df):
#     from matplotlib.colors import LinearSegmentedColormap
#     import numpy as np
#     import textwrap

#     time_correlation_df["Correlation"] = pd.to_numeric(
#         time_correlation_df["Correlation"], errors="coerce"
#     )
#     pivot_df = time_correlation_df.pivot_table(
#         index="Time (hrs)", columns=["Treatment", "Cytokine"], values="Correlation"
#     )
#     pivot_p_values = time_correlation_df.pivot_table(
#         index="Time (hrs)", columns=["Treatment", "Cytokine"], values="p-value"
#     )
#     plt.figure(figsize=(12, 8))

#     yellow = plt.get_cmap("Spectral")(0.5)
#     purple = plt.get_cmap("Spectral")(1.0)
#     # mod_colors = spectral(np.linspace(0.5, 1.0, 256))
#     ytp = LinearSegmentedColormap.from_list("ytp", [yellow, purple])

#     ax = sns.heatmap(pivot_df, cmap=ytp, annot=True, fmt=".2f")
#     column_width = pivot_df.shape[1]
#     for i, (idx, row) in enumerate(pivot_p_values.iterrows()):
#         for j, p_value in enumerate(row):
#             text = ""
#             if p_value < 0.001:
#                 text = "***"
#             elif p_value < 0.01:
#                 text = "**"
#             elif p_value < 0.05:
#                 text = "*"

#             if text:
#                 # Get current text
#                 current_text = ax.texts[i * column_width + j].get_text()
#                 # Update text with asterisks
#                 ax.texts[i * column_width + j].set_text(current_text + text)

#     treatments = list(pivot_df.columns.get_level_values(0).unique())
#     treatment_colors = {
#         treatment: sns.color_palette()[plotting_module.treatment_indeces[treatment]]
#         for treatment in treatments
#     }
#     # Calculate the linewidth in terms of axis fraction
#     linewidth = 4
#     adjustment_factor = 0.0225
#     # Draw lines at the start and end of each treatment group
#     for i, treatment in enumerate(treatments):
#         treatment_columns = [col for col in pivot_df.columns if col[0] == treatment]
#         left = pivot_df.columns.get_loc(treatment_columns[0]) + adjustment_factor
#         right = pivot_df.columns.get_loc(treatment_columns[-1]) + 1 - adjustment_factor

#         # Draw a line at the left edge of the treatment group, shifted slightly to the left
#         ax.axvline(
#             x=left,
#             color=treatment_colors[treatment],
#             linewidth=linewidth,
#             zorder=1,
#         )

#         # Draw a line at the right edge of the treatment group, shifted slightly to the left
#         ax.axvline(
#             x=right,
#             color=treatment_colors[treatment],
#             linewidth=linewidth,
#             zorder=1,
#         )

#     def format_x_ticks(ticklabels):
#         wrap_width = max(
#             len(word) for label in ticklabels for word in label.get_text().split()
#         )
#         # Get the tick labels
#         format_labels = [label.get_text() for label in ticklabels]
#         # Split the tick labels into treatment and cytokine
#         format_labels = [label.split("-") for label in format_labels]
#         # Wrap each part of the label and join with newline
#         wrapped_labels = [
#             "\n".join(textwrap.wrap(part, width=wrap_width))
#             for label in format_labels
#             for part in label
#         ]
#         # Group the wrapped parts back into labels
#         stacked_labels = [
#             "\n".join(wrapped_labels[i : i + len(format_labels[0])])
#             for i in range(0, len(wrapped_labels), len(format_labels[0]))
#         ]
#         return stacked_labels

#     ticklabels = ax.get_xticklabels()
#     ax.set_xticklabels(format_x_ticks(ax.get_xticklabels()))  # Rotate by 90 degrees
#     plt.xticks(rotation=0)
#     plt.title("Correlation Heatmap with Significance Levels")
#     plt.savefig(
#         Path("plots", "change_plots", "heatmap.png"), dpi=300, bbox_inches="tight"
#     )


# time_heatmap(time)
# heatmap_df = time

# import pandas as pd
# import textwrap

# # Assuming heatmap_df is your DataFrame
# heatmap_df["Correlation"] = pd.to_numeric(heatmap_df["Correlation"], errors="coerce")


# # Pivot to get Treatment and Cytokine combinations as columns and Time (hrs) as rows
# pivot_df = heatmap_df.pivot_table(
#     index="Time (hrs)", columns=["Treatment", "Cytokine"], values="Correlation"
# )

# # Pivot p-values
# pivot_p_values = heatmap_df.pivot_table(
#     index="Time (hrs)", columns=["Treatment", "Cytokine"], values="p-value"
# )


# # Assuming pivot_df contains your correlation data and pivot_p_values contains the corresponding p-values

# plt.figure(figsize=(12, 8))
# # Use fmt='.2f' to format the numbers to two decimal places
# ax = sns.heatmap(
#     pivot_df, cmap="coolwarm", annot=True, fmt=".2f"
# )  # Format the annotations to two decimal places

# # Get the width of the columns in order to place text correctly
# column_width = pivot_df.shape[1]

# # Add asterisks based on p-value thresholds
# for i, (idx, row) in enumerate(pivot_p_values.iterrows()):
#     for j, p_value in enumerate(row):
#         text = ""
#         if p_value < 0.001:
#             text = "***"
#         elif p_value < 0.01:
#             text = "**"
#         elif p_value < 0.05:
#             text = "*"

#         if text:
#             # Get current text
#             current_text = ax.texts[i * column_width + j].get_text()
#             # Update text with asterisks
#             ax.texts[i * column_width + j].set_text(current_text + text)

# plt.title("Correlation Heatmap with Significance Levels")
# plt.show()


# plotting_module.change_plot(
#     TAS.modules["TS_Speck"],
#     treatments=["ATP", "MSU", "Nigericin"],
#     mode="Acceleration",
#     filepath="test.png",
# )
# plotting_module.change_plot(
#     TAS.modules["TS_Speck"],
#     treatments=["ATP", "MSU", "Nigericin"],
#     mode="Change",
#     filepath="Combined_Change_Plot.png",
# )

# plotting_module.change_plot(
#     module=TAS.modules["TS_Cyto"], treatments=["ATP"], mode="Acceleration"
# )
# plotting_module.change_plot(
#     module=TAS.modules["TS_Cyto"], treatments=["ATP"], analytes=["IL1b"], mode="change"
# )


# file_dir = Path("plots", "change_plots", "cytokines")
# for treatment in ["ATP", "MSU", "Nigericin"]:
#     for method in ["Change", "Acceleration", "Smoothed Acceleration", "combo"]:
#         plotting_module.change_plot(
#             TAS.modules["TS_Cyto"],
#             treatments=[treatment],
#             mode=method,
#             filepath=f"{file_dir}/{treatment}_{method}_Cytokines.png",
#         )
#         for analyte in ["IL1b", "IL18"]:
#             plotting_module.change_plot(
#                 TAS.modules["TS_Cyto"],
#                 treatments=[treatment],
#                 analytes=[analyte],
#                 mode=method,
#                 filepath=f"{file_dir}/{treatment}_{method}_{analyte}.png",
#             )
# for method in ["Change", "Acceleration", "Smoothed Acceleration"]:
#     plotting_module.change_plot(
#         TAS.modules["TS_Speck"],
#         treatments=["ATP", "MSU", "Nigericin"],
#         mode=method,
#         filepath=f"{file_dir}/Speck_{method}.png",
#     )


# def correlation(dataset_1, dataset_2, parameter="Normalized_Measurement"):
#     """
#     Runs pearson/spearman (programmatically decided) correlation between two specified datasets.
#     Also runs cosine similarity between two specified datasets.
#     Returns a dictionary containing the results of the correlation and cosine similarity.

#     Parameters:
#     dataset_1: tuple containing module, analyte, and treatment list (optional)
#     dataset_2: tuple containing module, analyte, and treatment list (optional)

#     Returns:
#     Dictionary containing the results of the correlation and cosine similarity.

#     Example:
#     correlation(("TS_Cyto", "IL1b", "ATP"), ("TS_Speck", "None", "ATP"))
#     """

#     def balance_times(subset_1, subset_2):
#         time_points_1 = subset_1["Time (hrs)"].unique()
#         time_points_2 = subset_2["Time (hrs)"].unique()

#         # identify any time points outside the range of the other set
#         # if there are any, remove them from the set
#         # also remove them from the time points list
#         if time_points_1.max() > time_points_2.max():
#             subset_1 = subset_1[subset_1["Time (hrs)"] <= time_points_2.max()]
#             time_points_1 = subset_1["Time (hrs)"].unique()
#         elif time_points_2.max() > time_points_1.max():
#             subset_2 = subset_2[subset_2["Time (hrs)"] <= time_points_1.max()]
#             time_points_2 = subset_2["Time (hrs)"].unique()

#         # if there are any time points that are not in both sets, make a dictionary of the missing times and which set they are missing from
#         missing_times = {}
#         for time in time_points_1:
#             if time not in time_points_2:
#                 missing_times[time] = "subset_2"
#         for time in time_points_2:
#             if time not in time_points_1:
#                 missing_times[time] = "subset_1"
#         import pandas as pd

#         # for each time that is missing, find the nearest two times from the respective subset and average them to impute these rows
#         for time, subset_name in missing_times.items():
#             current_subset = subset_1 if subset_name == "subset_1" else subset_2

#             # average the nearest two times to the missing time
#             imputed_values = current_subset.groupby(["Experimental_Replicate"]).apply(
#                 lambda x:
#                 # get the two nearest times
#                 x.loc[
#                     (x["Time (hrs)"] == x["Time (hrs)"].max())
#                     | (x["Time (hrs)"] == x["Time (hrs)"].min())
#                 ][parameter].mean()
#             )
#             imputed_df = pd.DataFrame(
#                 {
#                     "Time (hrs)": [time] * len(imputed_values),
#                     "Experimental_Replicate": imputed_values.index,
#                     parameter: imputed_values,
#                 }
#             )
#             if subset_name == "subset_1":
#                 subset_1 = pd.concat([subset_1, imputed_df])
#             else:
#                 subset_2 = pd.concat([subset_2, imputed_df])

#         # Sorting the subsets by time after imputation
#         subset_1 = subset_1.sort_values(
#             by=["Time (hrs)", "Experimental_Replicate"]
#         ).reset_index(drop=True)
#         subset_2 = subset_2.sort_values(
#             by=["Time (hrs)", "Experimental_Replicate"]
#         ).reset_index(drop=True)

#         return subset_1, subset_2

#     def balance_replicates(subset_1, subset_2):
#         ##TODO: add later if needed.
#         return subset_1, subset_2

#     def cosine_similarity(vector_a, vector_b):
#         # Compute the dot product
#         dot_product = np.dot(vector_a, vector_b)
#         # Compute the magnitude (norm) of each vector
#         norm_a = np.linalg.norm(vector_a)
#         norm_b = np.linalg.norm(vector_b)
#         # Compute cosine similarity
#         similarity = dot_product / (norm_a * norm_b)
#         return similarity

#     subset_1 = self.modules[dataset_1[0]].data[
#         self.modules[dataset_1[0]].data["Analyte"] == dataset_1[1]
#     ]
#     subset_2 = self.modules[dataset_2[0]].data[
#         self.modules[dataset_2[0]].data["Analyte"] == dataset_2[1]
#     ]

#     if dataset_1[2] != "None":
#         subset_1 = subset_1[subset_1["Treatment"].isin(dataset_1[2])]
#     if dataset_2[2] != "None":
#         subset_2 = subset_2[subset_2["Treatment"].isin(dataset_2[2])]

#     # drop all columns except for time, experimental replicate, and the parameter
#     subset_1 = subset_1[["Time (hrs)", "Experimental_Replicate", parameter]]
#     subset_2 = subset_2[["Time (hrs)", "Experimental_Replicate", parameter]]

#     subset_1, subset_2 = balance_replicates(subset_1, subset_2)
#     subset_1, subset_2 = balance_times(subset_1, subset_2)

#     # run shapiro wilks test to determine normality
#     # if normal, run pearson, if not, run spearman
#     normality_test_1 = scipy.stats.shapiro(subset_1[parameter])[1]
#     normality_test_2 = scipy.stats.shapiro(subset_2[parameter])[1]
#     if normality_test_1 > 0.05 and normality_test_2 > 0.05:
#         corr_type = "pearson"
#         corr = scipy.stats.pearsonr(subset_1[parameter], subset_2[parameter])
#     else:
#         corr_type = "spearman"
#         corr = scipy.stats.spearmanr(subset_1[parameter], subset_2[parameter])
#     cos = cosine_similarity(subset_1[parameter], subset_2[parameter])
#     corr_dict = {
#         "Normality Test - Subset 1": normality_test_1,
#         "Normality Test - Subset 2": normality_test_2,
#         "Correlation Type": corr_type,
#         "Correlation": corr[0],
#         "p-value": corr[1],
#         "Cosine Similarity": cos,
#     }

#     return subset_1, subset_2, corr_dict


# parameter = "Normalized_Measurement"


# # count number of experimental replicates for each treatment
# cyto.groupby(["Treatment"])["Experimental_Replicate"].nunique()
# speck.groupby(["Treatment"])["Experimental_Replicate"].nunique()
