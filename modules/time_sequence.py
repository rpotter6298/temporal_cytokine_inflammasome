import pandas as pd


class Time_sequence:
    """
    A class for processing and analyzing time sequence data from a data_class object.
    """

    def __init__(self, data_class, time_limits: dict = None, variable: str = None):
        """
        Initializes the time_sequence_module object.

        :param data_class: A data_class object containing the raw data to be analyzed.
        :param time_limits: An optional dictionary specifying time limits for the data. Defaults to None.
        """
        self.comp = data_class
        self.name = "TS_" + data_class.name
        self.raw_data = data_class.data
        self.time_limits = time_limits
        

        if variable is not None:
            self.variable = variable
            self.data = self.structure_data(variable)
        else:
            possible_variables = self.raw_data.columns.difference(
                [
                    "Analyte",
                    "Treatment",
                    "Time (hrs)",
                    "Experimental_Replicate",
                    "Technical_Replicate",
                ]
            )
            if len(possible_variables) == 1:
                self.variable = possible_variables[0]
                self.data = self.structure_data(self.variable)
            else:
                print(
                    "Please specify the variable to be analyzed and run .structure_data(). Possible variables are:"
                )
                print(possible_variables)
                return

    def limit_data(self, data, time_limits):
        """
        Limits the data based on the provided time limits.

        :param data: DataFrame containing the data to be limited.
        :param time_limits: Dictionary specifying the time limits for the data.
        :return: DataFrame with limited data.
        """
        limited_data = data
        if time_limits is None:
            print("No time limits specified. Returning all data.")
            return limited_data
        else:
            for treatment in time_limits.keys():
                if treatment == "ALL":
                    limited_data = limited_data[
                        limited_data["Time (hrs)"] <= time_limits[treatment]
                    ]
                else:
                    limited_data = limited_data[
                        (limited_data["Treatment"].str.contains(treatment, case=False))
                        & (limited_data["Time (hrs)"] <= time_limits[treatment])
                        | (
                            ~limited_data["Treatment"].str.contains(
                                treatment, case=False
                            )
                        )
                    ]
            return limited_data

    def structure_data(
        self, variable: str, robust=True, merge_technical_replicates=True
    ):
        """
        Normalizes the data by dividing the selected variable by the initial values at time 0.
        Renames columns to standardize the names.

        :param variable: The variable to be normalized. Should be a column name in the data.
        :param robust: Optional boolean to calculate growth rate. Defaults to True.
        :return: DataFrame with normalized data.
        """

        ##Check that variable is in the columns
        if variable not in self.raw_data.columns:
            print("Variable not in columns. Please check input.")
            return
        ##Check that "Analyte" is in the columns and add it with a default value if not
        if "Analyte" not in self.raw_data.columns:
            self.raw_data["Analyte"] = "NA"
        if "Technical_Replicate" not in self.raw_data.columns:
            self.raw_data["Technical_Replicate"] = "NA"
        if "Time (hrs)" not in self.raw_data.columns:
            self.raw_data["Time (hrs)"] = 0
            print("Time (hrs) not in columns. Results may be corrupted.")
        ##Make sure that time column is numeric
        self.raw_data["Time (hrs)"] = pd.to_numeric(self.raw_data["Time (hrs)"])
        ##Rename the variable column to "Measurement"
        data = self.raw_data.rename(columns={variable: "Measurement"})
        if self.time_limits is not None:
            data = self.limit_data(data, self.time_limits)
        if merge_technical_replicates == True:
            data = (
                data.groupby(
                    ["Analyte", "Treatment", "Time (hrs)", "Experimental_Replicate"]
                )
                .mean(numeric_only=True)
                .reset_index()
            )

        ##Define function to normalize data
        def normalize_data(data):
            """
            Normalizes the data by dividing the selected variable by the initial values at time 0.
            :return: DataFrame with normalized data.
            """
            initial_values = data[data["Time (hrs)"] == 0].set_index(
                ["Analyte", "Treatment", "Experimental_Replicate"]
            )["Measurement"]
            data["Normalized_Measurement"] = data.apply(
                lambda row: row["Measurement"]
                / initial_values[
                    (row["Analyte"], row["Treatment"], row["Experimental_Replicate"])
                ],
                axis=1,
            )
            return data

        data = normalize_data(data)

        if robust == True:
            def calc_change_rate(data):
                
                change_rate = data.copy()
                time_increment = (
                    change_rate.groupby(
                        ["Analyte", "Treatment", "Experimental_Replicate"]
                    )["Time (hrs)"]
                    .diff()
                    .mode()[0]
                )
                change_rate.sort_values(
                    by=["Analyte", "Treatment", "Experimental_Replicate", "Time (hrs)"],
                    inplace=True,
                )
                change_rate["Normalized_Change"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Normalized_Measurement"].diff()
                # # set change_rate["Change_Rate"] to the difference between sequential measurements
                change_rate["Change_Rate"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Measurement"].diff()

                change_rate["Absolute_Change"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Measurement"].diff()

                change_rate["Change_Rate"] = (
                    change_rate.groupby(
                        ["Analyte", "Treatment", "Experimental_Replicate"]
                    )["Measurement"]
                    .apply(lambda x: x / x.shift(1))
                    .reset_index(
                        level=["Analyte", "Treatment", "Experimental_Replicate"],
                        drop=True,
                    )
                )
                
                change_rate["Change_Rate"] = change_rate["Change_Rate"].fillna(1)
                # # Function to calculate the window size based on time increments
                # # Change Rate Lagging - average change rate for the last two hours
                change_rate["Change_Rate_Diff"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Measurement"].diff()

                change_rate["Change_Rate_Lagging"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Change_Rate"].transform(
                    lambda x: x.rolling(window=(int(2 / time_increment))).mean(
                        numeric_only=True
                    )
                )
                
                change_rate["Change_Rate_Diff"] = change_rate[
                    "Change_Rate_Diff"
                ].fillna(0)
                change_rate["Change_Rate_Diff_Lagging"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Change_Rate_Diff"].transform(
                    lambda x: x.rolling(window=(int(2 / time_increment))).mean(
                        numeric_only=True
                    )
                )

                ##New Change Rate Variables, duplicated code with different names for easier comprehension.
                ##Can't delete the old ones because they are used elsewhere.
                
                change_rate["Delta"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Measurement"].diff()
                change_rate["Delta"] = change_rate["Delta"].fillna(0)
                change_rate["Normalized_Delta"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Normalized_Measurement"].diff()
                change_rate["Normalized_Delta"] = change_rate[
                    "Normalized_Delta"
                ].fillna(0)
                change_rate["Delta_Smooth"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Delta"].transform(
                    lambda x: x.rolling(window=(int(2 / time_increment))).mean(
                        numeric_only=True
                    )
                )
                change_rate["Normalized_Delta_Smooth"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Normalized_Delta"].transform(
                    lambda x: x.rolling(window=(int(2 / time_increment))).mean(
                        numeric_only=True
                    )
                )
                # Acceleration == Delta of Delta
                change_rate["Acceleration"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Delta"].diff()
                change_rate["Acceleration"] = change_rate["Acceleration"].fillna(0)
                change_rate["Normalized_Acceleration"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Normalized_Delta"].diff()
                change_rate["Normalized_Acceleration"] = change_rate[
                    "Normalized_Acceleration"
                ].fillna(0)

                change_rate["Acceleration_Smooth"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Acceleration"].transform(
                    lambda x: x.rolling(window=(int(2 / time_increment))).mean(
                        numeric_only=True
                    )
                )
                change_rate["Normalized_Acceleration_Smooth"] = change_rate.groupby(
                    ["Analyte", "Treatment", "Experimental_Replicate"]
                )["Normalized_Acceleration"].transform(
                    lambda x: x.rolling(window=(int(2 / time_increment))).mean(
                        numeric_only=True
                    )
                )

                return change_rate

            data = calc_change_rate(data)
        if self.time_limits is not None:
            data = self.limit_data(data, self.time_limits)
        return data

