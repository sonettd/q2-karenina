#/usr/bin/env python
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2016, The Karenina Project"
__credits__ = ["Jesse Zaneveld","Samuel L. Peoples"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from q2templates import render
from os.path import join
from pkg_resources import resource_filename
from karenina.fit_timeseries import parse_pcoa,parse_metadata,fit_input

def fit_timeseries(output_dir: str, pcoa : str, metadata:str, method : str, 
                individual_col: str, timepoint_col: str, treatment_col: str) -> None:
    
    #Handle missing parameters from user interface
    if 'None' in metadata:
	    metadata = None
    if 'None' in treatment_col:
	    treatment_col = None


    #readin and format      
    site, metadata = parse_pcoa(pcoa, individual_col,\
      timepoint_col, treatment_col, metadata=None)
    model_input = parse_metadata(metadata, individual_col,\
      timepoint_col, treatment_col, site)
   
     
    if treatment_col is not None:
        #If there is a non-empty treatment column, then fit a cohort model,
        # in which all individuals with the same treatment
        #have the same model parameters.
        
        output, cohort_output = fit_input(model_input,\
          individual_col, timepoint_col, treatment_col, method)
        
        result = cohort_output
        plot_name = 'Ornstein-Uhlenbeck cohort model fit results (treatment = {},method={})'.format(treatment_col,method)
    else:
        #IF no treatment is provided, fit a model to each individual
 
       result = fit_input(model_input, individual_col,\
          timepoint_col, treatment_col, method)
 
       plot_name = 'Ornstein-Uhlenbeck individual model fit results (method={})'.format(method)
   
    #Output the results (individual or cohort based model fit) to a CSV filte 
    result.to_csv(join(output_dir,"fit_timeseries_results.csv"), index=False)

    #Generate the index.html file required by QIIME2 by filling in the template
    render_index_html(output_dir,plot_name)


def render_index_html(output_dir,plot_name):

    #Find the filepath for the q2_emperor folder 'assets'.
    #NOTE: to understand what's happening here see e.g. stackoverflow example here:
    #https://stackoverflow.com/questions/39104/finding-a-file-in-a-python-module-distribution

    template_dir = resource_filename('q2_karenina', 'assets')

    #get the path to our basic, unfilled index.html file (in the assets folder of q2_emperor)
    index = join(template_dir, 'index.html')

    #Use q2_templates.render to fill in data specific to our output in this visualization.
    #Documentation for q2_templates.render is available here:
    # https://github.com/qiime2/q2templates/blob/master/q2templates/_templates.py

    render(index, output_dir, context={'plot_name': plot_name})


