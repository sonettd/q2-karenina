#/usr/bin/env python
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2016, The Karenina Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.0.1-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

import karenina.fit_timeseries as k_fit_timeseries
import pkg_resources
import qiime2
import q2templates
from q2_types.ordination import PCoAResults
import pandas as pd
import os

def fit_timeseries(output_dir: str, pcoa : str, metadata:str, method : str, 
                individual_col: str, timepoint_col: str, treatment_col: str) -> None:
    #pcoa = PCoAResults.read(pcoa).to_dataframe()
    #metadata = metadata.to_dataframe()
    if 'None' in metadata:
	    metadata = None
    if 'None' in treatment_col:
	    treatment_col = None
    site, metadata = k_fit_timeseries.parse_pcoa(pcoa, individual_col, timepoint_col, treatment_col, metadata=None)
    input = k_fit_timeseries.parse_metadata(metadata, individual_col, timepoint_col, treatment_col, site)
    #site = _parse_pcoa(pcoa)
    #input = _parse_metadata(metadata, individual_col, timepoint_col, treatment_col, site)
    if treatment_col is not None:
        output, cohort_output = k_fit_timeseries.fit_input(input, individual_col, timepoint_col, treatment_col, method)
        cohort_output.to_csv(os.path.join(output_dir,"cohort_fit_timeseries.csv"), index=False)
    else:
        output = k_fit_timeseries.fit_input(input, individual_col, timepoint_col, treatment_col, method)
    output.to_csv(os.path.join(output_dir,"individual_fit_timeseries.csv"), index=False)
