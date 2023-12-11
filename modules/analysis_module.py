from typing import Type, List, Dict, Any
from modules.stats_functions import *
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.stats
from itertools import combinations
from itertools import product
from scipy.stats import f_oneway
from scipy.stats import shapiro
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
import scikit_posthocs as sp
from statsmodels.stats.libqsturng import psturng
import itertools


class arbitrary_max:
    def __init__(self):
        self.real_max = None
        self.real_time = None
        self.window = None
        self.time = None
        self.max = None


class placeholder:
    def __init__(self):
        pass


class analysis_module:
    def __init__(self, component_list: List = []):
        self.modules = {module.name: module for module in component_list}
        self.time_compare()

    def time_compare_dict(self, reverse_lagging_change_rate=False):
        for module in self.modules:
            # Loop over both factors to compare
            data_set = self.modules[module].data
            data_set_diclist = []
            for analyte in data_set["Analyte"].unique():
                analyte_subset = data_set[data_set["Analyte"] == analyte]
                for treatment in analyte_subset["Treatment"].unique():
                    treatment_subset = analyte_subset[
                        analyte_subset["Treatment"] == treatment
                    ]
                    average_factor = (
                        treatment_subset.groupby(["Time (hrs)"])["Change_Rate_Lagging"]
                        .mean(numeric_only=True)
                        .reset_index()
                    )
                    # Find the index of the maximum value for each group of (Analyte, Treatment)
                    idx = average_factor["Change_Rate_Lagging"].idxmax()
                    # Select the corresponding rows from average_factor DataFrame
                    max_from_average_with_time = average_factor.loc[idx].reset_index(
                        drop=True
                    )
                    window_time = max_from_average_with_time[0]
                    window_subset = treatment_subset[
                        (treatment_subset["Time (hrs)"] <= window_time)
                        & (treatment_subset["Time (hrs)"] >= window_time - 2)
                    ]
                    window_averages = (
                        window_subset.groupby(["Time (hrs)"])
                        .mean(numeric_only=True)
                        .reset_index()
                    )
                    # get max change_rate from window_averages
                    max_change_rate_idx = window_averages["Change_Rate"].idxmax()
                    max_change_rate = window_averages.loc[max_change_rate_idx][
                        "Change_Rate"
                    ]
                    max_change_rate_time = window_averages.loc[max_change_rate_idx][
                        "Time (hrs)"
                    ]
                    # repeat the whole sliding window process for change_rate_diff
                    average_factor = (
                        treatment_subset.groupby(["Time (hrs)"])[
                            "Change_Rate_Diff_Lagging"
                        ]
                        .mean()
                        .reset_index()
                    )
                    # Find the index of the maximum value for each group of (Analyte, Treatment)
                    idx = average_factor["Change_Rate_Diff_Lagging"].idxmax()
                    # Select the corresponding rows from average_factor DataFrame
                    max_from_average_with_time = average_factor.loc[idx].reset_index(
                        drop=True
                    )
                    window_time = max_from_average_with_time[0]
                    window_subset = treatment_subset[
                        (treatment_subset["Time (hrs)"] <= window_time)
                        & (treatment_subset["Time (hrs)"] >= window_time - 2)
                    ]
                    window_averages = (
                        window_subset.groupby(["Time (hrs)"])
                        .mean(numeric_only=True)
                        .reset_index()
                    )

                    max_change_rate_diff_idx = window_averages[
                        "Change_Rate_Diff"
                    ].idxmax()
                    max_change_rate_diff = window_averages.loc[
                        max_change_rate_diff_idx
                    ]["Change_Rate_Diff"]
                    max_change_rate_diff_time = window_averages.loc[
                        max_change_rate_diff_idx
                    ]["Time (hrs)"]

                    # get the maximum measurement for each expirimental replicate in the treatment_subset and the time of each max
                    max_measurement = treatment_subset.groupby(
                        ["Experimental_Replicate"]
                    )["Measurement"].idxmax()
                    max_measurement_times = treatment_subset.loc[max_measurement][
                        "Time (hrs)"
                    ].mean()
                    max_measurement_values = treatment_subset.loc[max_measurement][
                        "Measurement"
                    ].mean()

                    max_normalized_measurement = treatment_subset.groupby(
                        ["Experimental_Replicate"]
                    )["Normalized_Measurement"].idxmax()
                    max_normalized_measurement_times = treatment_subset.loc[
                        max_normalized_measurement
                    ]["Time (hrs)"].mean(numeric_only=True)
                    max_normalized_measurement_values = treatment_subset.loc[
                        max_normalized_measurement
                    ]["Normalized_Measurement"].mean(numeric_only=True)

                    subset_dict = {
                        "Analyte": analyte,
                        "Treatment": treatment,
                        "Max_Change_Rate_Time (Proportion)": max_change_rate_time,
                        "Max_Change_Rate (Proportion)": max_change_rate,
                        "Max_Change_Rate (Absolute)": max_change_rate_diff,
                        "Max_Change_Rate_Time (Absolute)": max_change_rate_diff_time,
                        "Max_Measurement_Time": max_measurement_times,
                        "Max_Measurement": max_measurement_values,
                        "Max_Normalized_Measurement_Time": max_normalized_measurement_times,
                        "Max_Normalized_Measurement": max_normalized_measurement_values,
                    }
                    data_set_diclist.append(subset_dict)
                setattr(self.modules[module], "tdic", data_set_diclist)
            # window_time = max_from_average_with_time["Time (hrs)"].values[0]
            # window_subset = data_set[
            #     (data_set["Time (hrs)"] <= window_time)
            #     & (data_set["Time (hrs)"] >= window_time - 2)
            # ]

    def time_compare(self):
        """
        Compare time-dependent factors across modules.

        This method compares the "Normalized_Measurement" and "Change_Rate"
        factors across all modules in the instance. For each factor, it
        calculates the maximum value for each group of (Analyte, Treatment,
        Experimental_Replicate) in each module and adds a new DataFrame
        attribute to the corresponding module object with the maximum
        value for that factor.

        Args:
            self (TimeComparison): The instance to compare.

        Returns:
            None.
        """
        # Loop over all modules in the instance
        for module in self.modules:
            # Loop over both factors to compare
            for factor in [
                "Normalized_Measurement",
                "Change_Rate",
                "Change_Rate_Lagging",
            ]:
                # Group the module data by (Analyte, Treatment, Experimental_Replicate)
                # and calculate the maximum value of the factor for each group
                max_factor = (
                    self.modules[module]
                    .data.groupby(["Analyte", "Treatment", "Experimental_Replicate"])
                    .agg(Max_Factor=(factor, "max"))
                    .reset_index()
                )
                # Merge the maximum factor values back into the module data to identify
                # the rows with the maximum factor values for each group
                max_factor = self.modules[module].data.merge(
                    max_factor,
                    left_on=["Analyte", "Treatment", "Experimental_Replicate", factor],
                    right_on=[
                        "Analyte",
                        "Treatment",
                        "Experimental_Replicate",
                        "Max_Factor",
                    ],
                )
                max_factor = max_factor.drop_duplicates(
                    subset=["Analyte", "Treatment", "Experimental_Replicate"],
                    keep="first",
                )

                # Sort the resulting DataFrame by (Analyte, Treatment, Experimental_Replicate)
                sorted_data = max_factor.sort_values(
                    by=["Analyte", "Treatment", "Experimental_Replicate"]
                )
                # Set the attribute of the module object to the sorted DataFrame
                setattr(self.modules[module], f"Max_{factor}", sorted_data)

    def check_modifier_impact(self, module=None, reverse_lag_change_rate=False):
        if module is None:
            module = self.modules[next(iter(self.modules))]
        # Check if a modifier exists, if not, return None
        if module.comp.modifier is None:
            return None
        # Check that the modifier is a string, if not, return None
        if not isinstance(module.comp.modifier, str):
            return None
        results = []
        reference = pd.DataFrame(module.tdic)
        reference["Max_Change_Rate"] = reference["Max_Change_Rate (Proportion)"]
        reference["Max_Change_Rate_Time"] = reference[
            "Max_Change_Rate_Time (Proportion)"
        ]
        for treatment in module.comp.treatments:
            print(treatment)
            result = {"Treatment": treatment, "Modifier": module.comp.modifier}
            for factor in [
                "Normalized_Measurement",
                "Change_Rate",
            ]:
                # Get the data for the current treatment and factor
                base_time = reference[reference["Treatment"] == treatment][
                    f"Max_{factor}_Time"
                ].values[0]
                mbase_time = reference[
                    reference["Treatment"] == f"{module.comp.modifier}_{treatment}"
                ][f"Max_{factor}_Time"].values[0]

                nearest_timepoint = round(base_time * 2) / 2
                mnearest_timepoint = round(mbase_time * 2) / 2
                value = module.data[
                    (module.data["Treatment"] == treatment)
                    & (module.data["Time (hrs)"] == nearest_timepoint)
                ][factor]
                print(value)
                mvalue = module.data[
                    (
                        module.data["Treatment"]
                        == str(module.comp.modifier + "_" + treatment)
                    )
                    & (module.data["Time (hrs)"] == mnearest_timepoint)
                ][factor]
                print(mvalue)

                data = getattr(module, f"Max_{factor}")
                if factor == "Normalized_Measurement":
                    value = data[data["Treatment"] == treatment][factor]
                    mvalue = data[
                        data["Treatment"] == str(result["Modifier"] + "_" + treatment)
                    ][factor]
                time = reference[reference["Treatment"] == treatment][
                    f"Max_{factor}_Time"
                ].values[0]
                mtime = reference[
                    reference["Treatment"] == f"{module.comp.modifier}_{treatment}"
                ][f"Max_{factor}_Time"].values[0]
                p, ci = bootstrap_t_test(value, mvalue)
                # tp, tci = bootstrap_t_test(time, mtime)
                result[f"base_peak_{factor}"] = value.mean()
                result[f"mod_peak_{factor}"] = mvalue.mean()
                result[f"peak_{factor}_p"] = p
                result[f"peak_{factor}_ci"] = ci
                result[f"base_peak_{factor}_time"] = time.mean()
                result[f"mod_peak_{factor}_time"] = mtime.mean()
                # result[f"peak_{factor}_time_p"] = tp
                # result[f"peak_{factor}_time_ci"] = tci
            results.append(result)
        results_df = pd.DataFrame(results)
        return results_df

    def time_point_comparison(
        self,
        Treatment: str = None,
        Analyte: str = None,
        module=None,
        time1: float = 0,
    ):
        """
        Compare the mean speck formation of a treatment at time1 with the mean speck formation at all other time points.

        Args:
            treatment (str): The treatment to analyze.
            time1 (float, optional): The initial time point for comparison. Defaults to 0.

        Returns:
            pd.DataFrame: A DataFrame containing time points, mean speck formation, differences, and p-values.
        """
        if Treatment is None:
            Treatment = module.data["Treatment"].unique()[0]
        if Analyte is None:
            Analyte = module.data["Analyte"].unique()[0]
        data = module.data[
            (module.data["Treatment"] == Treatment)
            & (module.data["Analyte"] == Analyte)
        ]
        unique_timepoints = data["Time (hrs)"].unique()
        unique_timepoints = np.delete(
            unique_timepoints, np.where(unique_timepoints == time1)
        )
        time1_data = data[data["Time (hrs)"] == time1]

        # Store raw p-values
        raw_p_values = []
        results = pd.DataFrame(columns=["Time", "Mean", "Difference", "CI"])

        for time2 in unique_timepoints:
            time2_data = data[(data["Time (hrs)"] == time2)]
            p, ci = bootstrap_t_test(
                time1_data["Measurement"], time2_data["Measurement"]
            )
            raw_p_values.append(p)
            new_row = pd.DataFrame(
                {
                    "Time": float(time2),
                    "Mean": time2_data["Measurement"].mean(),
                    "Difference": -(
                        time1_data["Measurement"].mean()
                        - time2_data["Measurement"].mean()
                    ),
                    "P-Value": p,
                    "CI": str(str(ci[0]) + ":" + str(ci[1])),
                },
                index=[0],
            )

            results = pd.concat([results, new_row], ignore_index=True)
        adj_p_values = multipletests(raw_p_values, method="fdr_bh")[1]
        results["Adjusted P-Value"] = adj_p_values
        # Dunnett's correction
        k = len(raw_p_values)
        q_values = [np.abs(psturng(p * k, k, np.inf)) for p in raw_p_values]

        results["Dunnett Adjusted P-Value"] = q_values
        # results.set_index("Time", inplace=True)
        return results

    def aggregate_time_comparisons(
        self, module, output_directory: Path = Path("reports"), time1: float = 0
    ):
        # Get a list of all treatments, including modified treatments
        treatments = [treatment for treatment in module.data["Treatment"].unique()]
        analytes = [analyte for analyte in module.data["Analyte"].unique()]
        # Initialize an empty dictionary to store the significant times for each treatment
        significant_times = {}
        filename = output_directory / f"{module.name}_significant_times.xlsx"

        for analyte in analytes:
            data = module.data[module.data["Analyte"] == analyte]
            # Iterate over each treatment
            for treatment in treatments:
                # Calculate the time range with a significant difference from zero count
                table = self.time_point_comparison(treatment, analyte, module)
                significant_range = (
                    table[table["P-Value"] < 0.05].Time.min(),
                    table[table["P-Value"] < 0.05].Time.max(),
                )

                # Store the significant times in the dictionary
                significant_times[treatment] = significant_range

                # Write the significant times to the report
                print(
                    f"{treatment} is significant between {significant_range[0]} and {significant_range[1]}."
                )

                sheet_name = treatment
                if len(analytes) > 1:
                    sheet_name = str(analyte + "_" + treatment)
                # Write the significant times to an Excel file
                # If filename doesn't exist, create it
                if not filename.exists():
                    table.to_excel(filename, sheet_name=sheet_name, index=True)
                with pd.ExcelWriter(
                    filename, mode="a", if_sheet_exists="replace"
                ) as writer:
                    table.to_excel(writer, sheet_name=sheet_name, index=True)

    def create_summary_table(self, module=None):
        """
        Creates a summary table with treatment statistics including max measurement, max change rate, and time to max values.

        Returns:
            pd.DataFrame: A summary table with treatment statistics.
        """
        if module is None:
            module = self.modules[next(iter(self.modules))]

        columns = [
            "Treatment",
            "Max Measurement",
            "Max Measurement (Normalized)",
            "Max Change Rate",
            "Time to Max Measurement",
            "Time to Max Change Rate",
        ]
        # If there are more than 1 analytes, add analyte column to start of columns list
        if len(module.data["Analyte"].unique()) > 1:
            columns.insert(0, "Analyte")
        summary = pd.DataFrame(columns=columns)
        all_treatments = module.data["Treatment"].unique()
        for analyte in module.data["Analyte"].unique():
            for treatment in all_treatments:
                set = module.Max_Normalized_Measurement[
                    (module.Max_Normalized_Measurement["Treatment"] == treatment)
                    & (module.Max_Normalized_Measurement["Analyte"] == analyte)
                ]
                time_0_set = module.data[
                    (module.data["Treatment"] == treatment)
                    & (module.data["Analyte"] == analyte)
                    & (module.data["Time (hrs)"] == 0)
                ]
                # Test for normality
                stat, p = shapiro(time_0_set["Measurement"])
                if p < 0.05:
                    print(
                        f"Data for Treatment {treatment} and Analyte {analyte} is not normally distributed"
                    )
                    print(p)
                    # You can decide to continue or break based on this result
                    continue
                else:
                    print(
                        f"Data for Treatment {treatment} and Analyte {analyte} is normally distributed"
                    )
                    print(p)
                max_measurement = set["Measurement"].mean()

                max_measuremnt_p, ci = bootstrap_t_test(
                    time_0_set["Measurement"], set["Measurement"]
                )
                max_measurement_normalized = set["Normalized_Measurement"].mean()
                time_to_max_measurement = set["Time (hrs)"].mean()
                tset = module.Max_Change_Rate[
                    (module.Max_Change_Rate["Treatment"] == treatment)
                    & (module.Max_Change_Rate["Analyte"] == analyte)
                ]
                max_change_rate = tset["Change_Rate"].mean()
                max_change_rate_normalized = tset["Normalized_Change_Rate"].mean()
                time_to_max_change_rate = tset["Time (hrs)"].mean()
                new_row = pd.DataFrame(
                    {
                        "Analyte": analyte,
                        "Treatment": treatment,
                        "Max Measurement": max_measurement,
                        "Max Measurement P-Value": max_measuremnt_p,
                        "Max Measurement (Normalized)": max_measurement_normalized,
                        "Max Change Rate": max_change_rate,
                        "Max Change Rate (Normalized)": max_change_rate_normalized,
                        "Time to Max Measurement": time_to_max_measurement,
                        "Time to Max Change Rate": time_to_max_change_rate,
                    },
                    index=[0],
                )
                if len(module.data["Analyte"].unique()) == 1:
                    new_row.drop(columns="Analyte", inplace=True)
                summary = pd.concat([summary, new_row], ignore_index=True)
        return summary

    def split_analytes(self, module):
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

    def prepare_lineplot(
        self,
        module,
        treatments: List = None,
        output_path: Path = None,
        measurement_type=None,
        analyte_selection=None,
    ):
        if treatments is None:
            treatments = module.comp.treatments
        if treatments == "all":
            treatments = module.data["Treatment"].unique()
        data = module.data[
            module.data["Treatment"].isin(treatments)
        ]  # Filter the data to only include the specified treatments
        if analyte_selection is not None:
            data = data[data["Analyte"].isin(analyte_selection)]
        # Check if there is only one treatment and multiple analytes
        if len(treatments) == 1 and len(data.Analyte.unique()) > 1:
            hue = "Analyte"
        # Check if there is only one analyte and multiple treatments
        elif len(data.Analyte.unique()) == 1 and len(treatments) >= 1:
            hue = "Treatment"
        # Otherwise, raise an error
        else:
            raise ValueError(
                "Must specify either one treatment and multiple analytes, or one analyte and multiple treatments."
            )
        if measurement_type:
            y = measurement_type
        elif len(data[hue].unique()) == 1:
            y = "Measurement"
        else:
            y = "Normalized_Measurement"

        # Specify the order of treatments to control the colors
        hue_order = None  # Initialize hue_order as None

        if hue == "Treatment":
            hue_order = treatments  # Set hue_order only when hue is "Treatment"

        if output_path is not None:
            filename = output_path / f"{module.name}_{'_'.join(treatments)}.png"
            lineplot_dict = {
                "dataframe": data,
                "y": y,
                "hue": hue,
                "hue_order": hue_order,
                "filepath": filename,
            }
        else:
            lineplot_dict = {
                "dataframe": data,
                "y": y,
                "hue": hue,
                "hue_order": hue_order,
            }
        return lineplot_dict

    def find_arbitrary_max(self, module):
        self.time_compare()
        max_values = module.Max_Normalized_Measurement
        max_dict = {}
        module.arbitrary_maximums = placeholder()
        for analyte in max_values["Analyte"].unique():
            analyte_data = max_values[max_values["Analyte"] == analyte]
            mean_max_times = analyte_data.groupby(["Treatment"]).mean()["Time (hrs)"]
            # Round the mean max times to the nearest measurement time
            mean_max_times = mean_max_times.apply(
                lambda x: analyte_data["Time (hrs)"].tolist()[
                    np.abs(np.array(analyte_data["Time (hrs)"].tolist()) - x).argmin()
                ]
            )
            for treatment, time in mean_max_times.items():
                print(treatment)
                time_table = self.time_point_comparison(
                    treatment, module=module, time1=time
                )
                max_window = time_table[time_table["P-Value"] > 0.05]
                arbitrary_max_time = time_table[time_table["P-Value"] > 0.05][
                    "Time"
                ].min()
                arb_obj = arbitrary_max()
                arb_obj.real_time = time
                arb_obj.real_max = analyte_data.groupby(["Treatment"]).mean()[
                    "Normalized_Measurement"
                ][treatment]
                arb_obj.time = arbitrary_max_time
                arb_obj.window = max_window
                arb_obj.max = module.data[
                    (module.data["Time (hrs)"] == arbitrary_max_time)
                    & (module.data["Analyte"] == analyte)
                    & (module.data["Treatment"] == treatment)
                ]["Normalized_Measurement"].mean()
                setattr(module.arbitrary_maximums, treatment, arb_obj)

    def compare_max_time_distances(self, moduleA, moduleB):
        # Get tlist of all treatments which appear in both moduleA and moduleB data
        def get_times_by_treatment_analyte(df):
            result = {}
            analytes = df["Analyte"].unique()

            for analyte in analytes:
                analyte_dict = {}
                for treatment in df["Treatment"].unique():
                    times = df[
                        (df["Treatment"] == treatment) & (df["Analyte"] == analyte)
                    ]["Time (hrs)"].tolist()
                    analyte_dict[treatment] = times
                result[analyte] = analyte_dict

            return result

        def compare_differences(results, treatment1, treatment2):
            mean_diff1 = results[treatment1]["mean_diff"]
            mean_diff2 = results[treatment2]["mean_diff"]
            ci_diff1 = results[treatment1]["ci_diff"]
            ci_diff2 = results[treatment2]["ci_diff"]

            mean_diff_diff = mean_diff1 - mean_diff2

            se_diff1 = (ci_diff1[1] - ci_diff1[0]) / (
                2
                * scipy.stats.t.ppf(
                    (1 + 0.95) / 2,
                    len(times_A[treatment]) + len(times_B[treatment]) - 2,
                )
            )
            se_diff2 = (ci_diff2[1] - ci_diff2[0]) / (
                2
                * scipy.stats.t.ppf(
                    (1 + 0.95) / 2,
                    len(times_A[treatment]) + len(times_B[treatment]) - 2,
                )
            )

            se_diff_diff = np.sqrt(se_diff1**2 + se_diff2**2)
            ci_diff_diff = (ci_diff1[0] - ci_diff2[1], ci_diff1[1] - ci_diff2[0])
            return mean_diff_diff, ci_diff_diff

        max_A = get_times_by_treatment_analyte(moduleA.Max_Normalized_Measurement)
        max_B = get_times_by_treatment_analyte(moduleB.Max_Normalized_Measurement)
        treatments = list(
            set(moduleA.data["Treatment"].unique()).intersection(
                set(moduleB.data["Treatment"].unique())
            )
        )
        results = {}
        if len(max_A) == 1:
            max_A = {moduleA.name: max_A[next(iter(max_A.keys()))]}
        if len(max_B) == 1:
            max_B = {moduleB.name: max_B[next(iter(max_B.keys()))]}
        for treatment in treatments:
            for analyte_A, times_A in max_A.items():
                name_A = analyte_A
                for analyte_B, times_B in max_B.items():
                    name_B = analyte_B
                    key = f"{name_A}-{name_B}_{treatment}"
                    mean_A, se_A = mean_standard_error(times_A[treatment])
                    mean_B, se_B = mean_standard_error(times_B[treatment])

                    mean_diff = mean_A - mean_B
                    se_diff = np.sqrt(se_A**2 + se_B**2)
                    ci_diff = scipy.stats.t.interval(
                        0.95,
                        len(times_A[treatment]) + len(times_B[treatment]) - 2,
                        loc=mean_diff,
                        scale=se_diff,
                    )
                    results[key] = {"mean_diff": mean_diff, "ci_diff": ci_diff}

        comparison_results = {}
        pairwise_combinations = list(product(list(max_A.keys()), list(max_B.keys())))
        for pair in pairwise_combinations:
            prefix = f"{pair[0]}-{pair[1]}"
            for treatment_combination in combinations(treatments, 2):
                option1 = f"{prefix}_{treatment_combination[0]}"
                option2 = f"{prefix}_{treatment_combination[1]}"
                mean_diff_diff, ci_diff_diff = compare_differences(
                    results, option1, option2
                )
                comparison_results[f"{option1}:{option2}"] = {
                    "mean_diff_diff": mean_diff_diff,
                    "ci_diff_diff": ci_diff_diff,
                    "significance": not (ci_diff_diff[0] <= 0 <= ci_diff_diff[1]),
                }
        return comparison_results

    def compare_max_normalized_measurement_anova(self, df):
        anova_results = {}

        for analyte in df["Analyte"].unique():
            df_analyte = df[df["Analyte"] == analyte]
            groups = [
                df_analyte["Normalized_Measurement"][
                    df_analyte["Treatment"] == treatment
                ]
                for treatment in df_analyte["Treatment"].unique()
            ]
            f_val, p_val = f_oneway(*groups)
            anova_results[analyte] = {"f_val": f_val, "p_val": p_val}
            if p_val < 0.05:
                post_hoc = pairwise_tukeyhsd(
                    df_analyte["Normalized_Measurement"], df_analyte["Treatment"]
                )
                anova_results[analyte]["post_hoc"] = post_hoc.summary()

        return anova_results

    def compute_ratio(
        self,
        module_name,
        measurement_type="Measurement",
        normalize_start=False,
    ):
        module = self.modules[module_name]
        ##If the module's data does not have two analytes, return an error
        if len(module.data["Analyte"].unique()) < 2:
            raise ValueError("Model does not have two analytes.")
        self.modules[module_name].ratio_data = {}
        for analyte_pair in list(combinations(module.data["Analyte"].unique(), 2)):
            split_data = self.split_analytes(module_name, analyte_pair)
            df1, df2 = split_data.values()
            analyte1, analyte2 = split_data.keys()
            merged_data = pd.merge(
                df1, df2, on=["Treatment", "Time (hrs)", "Experimental_Replicate"]
            )
            merged_data[f"Ratio_{analyte1}_{analyte2}"] = (
                merged_data[f"{measurement_type}_x"]
                / merged_data[f"{measurement_type}_y"]
            )
            merged_data[f"Ratio_{analyte2}_{analyte1}"] = (
                merged_data[f"{measurement_type}_y"]
                / merged_data[f"{measurement_type}_x"]
            )

            # If normalize_start is True, normalize the ratio to the first time point for each treatment and experimental replicate
            if normalize_start:
                merged_data[f"Ratio_{analyte1}_{analyte2}"] = merged_data.groupby(
                    ["Treatment", "Experimental_Replicate"]
                )[f"Ratio_{analyte1}_{analyte2}"].transform(lambda x: x / x.iloc[0])
            self.modules[module_name].ratio_data[
                analyte1 + ":" + analyte2
            ] = merged_data

    def split_analytes(self, module_name, analyte_pair):
        """
        Split the data into separate dataframes for each analyte.

        Args:
            model (time_sequence): The time_sequence object to split.

        Returns:
            dict: A dictionary of dataframes for each analyte.
        """
        module = self.modules[module_name]
        split_data = {}
        for analyte in analyte_pair:
            split_data[analyte] = module.data[module.data["Analyte"] == analyte]
        return split_data

    def correlation(self, dataset_1, dataset_2, parameter="Normalized_Measurement"):
        """
        Runs pearson/spearman (programmatically decided) correlation between two specified datasets.
        Also runs cosine similarity between two specified datasets.
        Returns a dictionary containing the results of the correlation and cosine similarity.

        Parameters:
        dataset_1: tuple containing module, analyte, and treatment list (optional)
        dataset_2: tuple containing module, analyte, and treatment list (optional)

        Returns:
        Dictionary containing the results of the correlation and cosine similarity.

        Example:
        correlation(("TS_Cyto", "IL1b", "ATP"), ("TS_Speck", "None", "ATP"))
        """

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
                imputed_values = current_subset.groupby(
                    ["Experimental_Replicate"]
                ).apply(
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

        def balance_replicates(subset_1, subset_2):
            ##TODO: add later if needed.
            return subset_1, subset_2

        def cosine_similarity(vector_a, vector_b):
            # Compute the dot product
            dot_product = np.dot(vector_a, vector_b)
            # Compute the magnitude (norm) of each vector
            norm_a = np.linalg.norm(vector_a)
            norm_b = np.linalg.norm(vector_b)
            # Compute cosine similarity
            similarity = dot_product / (norm_a * norm_b)
            return similarity

        subset_1 = self.modules[dataset_1[0]].data[
            self.modules[dataset_1[0]].data["Analyte"] == dataset_1[1]
        ]
        subset_2 = self.modules[dataset_2[0]].data[
            self.modules[dataset_2[0]].data["Analyte"] == dataset_2[1]
        ]

        if dataset_1[2] != "None":
            subset_1 = subset_1[subset_1["Treatment"].isin(dataset_1[2])]
        if dataset_2[2] != "None":
            subset_2 = subset_2[subset_2["Treatment"].isin(dataset_2[2])]

        # drop all columns except for time, experimental replicate, and the parameter
        subset_1 = subset_1[["Time (hrs)", "Experimental_Replicate", parameter]]
        subset_2 = subset_2[["Time (hrs)", "Experimental_Replicate", parameter]]

        subset_1, subset_2 = balance_replicates(subset_1, subset_2)
        subset_1, subset_2 = balance_times(subset_1, subset_2)

        time_correlation = {}
        for time in subset_1["Time (hrs)"].unique():
            subset_1_time = subset_1[subset_1["Time (hrs)"] == time]
            subset_2_time = subset_2[subset_2["Time (hrs)"] == time]
            # run shapiro wilks test to determine normality
            # if normal, run pearson, if not, run spearman
            normality_test_1 = scipy.stats.shapiro(subset_1_time[parameter])[1]
            normality_test_2 = scipy.stats.shapiro(subset_2_time[parameter])[1]
            if normality_test_1 > 0.05 and normality_test_2 > 0.05:
                corr_type = "pearson"
                corr = scipy.stats.pearsonr(
                    subset_1_time[parameter], subset_2_time[parameter]
                )
            else:
                corr_type = "spearman"
                corr = scipy.stats.spearmanr(
                    subset_1_time[parameter], subset_2_time[parameter]
                )
            time_correlation[time] = {
                "Normality Test - Subset 1": normality_test_1,
                "Normality Test - Subset 2": normality_test_2,
                "Correlation Type": corr_type,
                "Correlation": corr[0],
                "p-value": corr[1],
            }

        # run shapiro wilks test to determine normality
        # if normal, run pearson, if not, run spearman
        normality_test_1 = scipy.stats.shapiro(subset_1[parameter])[1]
        normality_test_2 = scipy.stats.shapiro(subset_2[parameter])[1]
        if normality_test_1 > 0.05 and normality_test_2 > 0.05:
            corr_type = "pearson"
            corr = scipy.stats.pearsonr(subset_1[parameter], subset_2[parameter])
        else:
            corr_type = "spearman"
            corr = scipy.stats.spearmanr(subset_1[parameter], subset_2[parameter])
        cos = cosine_similarity(subset_1[parameter], subset_2[parameter])
        corr_dict = {
            "Normality Test - Subset 1 Overall": normality_test_1,
            "Normality Test - Subset 2 Overall": normality_test_2,
            "Correlation Type - Overall": corr_type,
            "Correlation - Overall": corr[0],
            "p-value - Overall": corr[1],
            "Cosine Similarity": cos,
            "Time Correlation": time_correlation,
        }

        return subset_1, subset_2, corr_dict
