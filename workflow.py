## Basic Parameters
base_treatments = ["ATP", "MSU", "Nigericin"]
modifier = "MCC950"
short_name_dict = {"nig": "Nigericin"}
time_limits = {"ALL": 21}
exclusion_list = [("speck", ("MSU", 1), ("MSU", 2), ("MCC950*", 1))]

##Build the analysis module
from modules.workflow_functions import prepare_analysis_class
TAS = prepare_analysis_class(base_treatments=base_treatments, modifier=modifier, short_name_dict=short_name_dict, time_limits=time_limits, exclusion_list=exclusion_list)

##Generate plots and reports
from modules.workflow_functions import generate_plots
generate_plots(TAS)

from modules.workflow_functions import generate_reports
generate_reports(TAS)
