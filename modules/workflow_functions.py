import pandas as pd
from pathlib import Path
## This Module contains all the functions used in the workflow just to keep the main script clean
## This Module contains all the functions used in the workflow just to keep the main script clean
def prepare_analysis_class(base_treatments, time_limits, short_name_dict, modifier:str=None, exclusion_list=[]):
    from modules.time_sequence import Time_sequence
    from modules.analysis_module import Analysis_module
    from modules.cytokine_adapter import Cytokine_Adapter
    from modules.speck_adapter import Speck_adapter
    from modules.ldh_module import LDH_module
    ## Build the separate class objects for the different tests
    cytokine = Time_sequence(
        Cytokine_Adapter("Cyto", base_treatments, short_name_dict=short_name_dict), variable="Concentration"
    )
    speck = Time_sequence(
        Speck_adapter(
            "Speck", base_treatments, modifier=modifier, short_name_dict=short_name_dict
        )
    )
    speck.structure_data("SpeckFormation")
    speck.data = speck.limit_data(speck.data, time_limits=time_limits)

    def remove_exclusions(data, exclusion_list, speck=speck, cytokine=cytokine):
        for exclusion in exclusion_list:
            dataset_name = exclusion[0]
            dataset_obj = eval(dataset_name)
            dataset = dataset_obj.data
            for exclusion_tuple in exclusion[1:]:
                #print(exclusion_tuple)
                treatment = exclusion_tuple[0]
                replicate = exclusion_tuple[1]
                if "*" in treatment:
                    dataset = dataset[
                        ~(
                            (dataset["Treatment"].str.contains(treatment.replace("*", "")
                            )) & (dataset["Experimental_Replicate"] == str(replicate))
                        )
                    ]
                else:
                    dataset = dataset[
                        ~(
                            (dataset["Treatment"] == treatment)
                            & (dataset["Experimental_Replicate"] == str(replicate))
                        )
                    ]
        return dataset

    speck.data = remove_exclusions(speck.data, exclusion_list)
    ##Construct main analysis class object
    TAS = Analysis_module([cytokine, speck])
    TAS.time_compare()
    TAS.compute_ratio("TS_Cyto")
    TAS.time_compare_dict()
    TAS.modules["LDH"] = LDH_module()
    return TAS

def generate_plots(TAS, file_extension=".png"):
    from modules.plotting_module import plotting_module
    import warnings
    warnings.filterwarnings('ignore')
    def generate_change_plots(TAS, file_extension=".png"):
        file_dir = Path("plots", "change_plots")
        file_dir.mkdir(parents=True, exist_ok=True)
        for treatment in ["ATP", "MSU", "Nigericin"]:
            for method in ["Delta", "Acceleration", "Smoothed Acceleration"]:
                plotting_module.change_plot(
                    TAS.modules["TS_Cyto"],
                    treatments=[treatment],
                    mode=method,
                    filepath=f"{file_dir}/{treatment}_{method}_Cytokines.png",
                    normalized=False,
                )
                plotting_module.change_plot(
                    TAS.modules["TS_Speck"],
                    treatments=[treatment],
                    mode=method,
                    filepath=f"{file_dir}/{treatment}_{method}_Speck.png",
                    normalized=False,
                )
                for analyte in ["IL1b", "IL18"]:
                    plotting_module.change_plot(
                        TAS.modules["TS_Cyto"],
                        treatments=[treatment],
                        analytes=[analyte],
                        mode=method,
                        filepath=f"{file_dir}/{treatment}_{method}_{analyte}.png",
                        normalized=False,
                    )
    
    def generate_cytokine_plots(TAS, file_extension=".png"):
            import matplotlib.ticker as ticker
            plot_dir = Path("plots", "cytokine_plots")
            plot_dir.mkdir(parents=True, exist_ok=True)
            ## CYTOKINES PLOTS, ALL TREATMENTS, SINGLE ANALYTE
            for analyte in TAS.modules["TS_Cyto"].data["Analyte"].unique():
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
                        plot_dir,
                        "".join(
                            [
                                aname.replace("-", "").replace(" ", "_"),
                                "_concentration_over_time",
                                file_extension,
                            ]
                        ),
                    ),
                )

    def generate_combined_speck_cyto_toxicity_plot(TAS, file_extension=".png"):
        plot_dir = Path("plots")
        plot_dir.mkdir(parents=True, exist_ok=True)
        for treatment in TAS.modules["LDH"].data["Treatment"].unique():
            plotting_module.combined_plot(
                analysis_class=TAS,
                treatments=[treatment],
                filepath=f"{plot_dir}/{treatment}_combined{file_extension}",
                padding=(-5, 75),
                font_scale=1.5,
            )
    generate_combined_speck_cyto_toxicity_plot(TAS, file_extension)
    generate_cytokine_plots(TAS, file_extension)
    generate_change_plots(TAS, file_extension)

def generate_reports(TAS):
    import warnings
    # To ignore all warnings
    warnings.filterwarnings('ignore')
    # check if reports directory exists, if not create it
    Path("reports").mkdir(parents=True, exist_ok=True)
    
    ## Create Summaries
    for module in ["TS_Cyto", "TS_Speck"]:
        print(TAS.modules[module].data.columns)
        (TAS.create_summary_table(TAS.modules[module])).to_excel(
            Path("reports", "".join(["summary_", module, ".xlsx"])), index=False
        )
    ## Create Modifier Impact Reports
    for module in [
        module for module in ["TS_Cyto", "TS_Speck"] if TAS.modules[module].comp.modifier
    ]:

        TAS.check_modifier_impact(TAS.modules[module]).T.to_excel(
            Path("reports", "".join(["modifier_impact_", module, ".xlsx"])),
            header=False,
        )
    ## Create Time Point Comparisons
    for module in ["TS_Cyto", "TS_Speck"]:
        TAS.aggregate_time_comparisons(TAS.modules[module])

