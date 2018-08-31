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
	

def _parse_pcoa(pcoa):
    # Parse PCoA Contents
    pcoa = open(pcoa,"rb")
    lines = pcoa.readlines()
    site = []
    i = 0
    for line in lines:
        line = line.decode("utf-8")
        if line.startswith("Eigvals"):
            eigs = lines[i+1].decode("utf-8").strip('\n').strip('\r').split("\t")
        elif line.startswith("Proportion explained"):
            propEx = lines[i+1].decode("utf-8").strip('\n').strip('\r').split("\t")
        elif line.startswith("Species"):
            species = lines[i + 1].decode("utf-8").strip('\n').strip('\r').split("\t")

        elif line.startswith("Site"):
            # We don't want Site constraints.
            if ("constraints") in line:
                break
            max = int(line.split("\t")[1])+i
            j = i + 1
            while j <= max:
                t_line = lines[j].decode('utf-8')
                site.append(t_line.strip('\n').strip('\r').split("\t"))
                j += 1
        i += 1
    t_site = []
    for item in site:
        t_site.append([item[0],item[1],item[2],item[3]])
    site = t_site

    #We now have variables:
            #First three are stored for now, will be utilized later when standardization of axes can be supported.
        # eigs: every eigenvalue
        # propEx: every proportion explained
        # species: every species

        # site: Every site, with the highest three proportion-explained PCs, in the form:
            # subjectID, x1, x2, x3
    return site

def _parse_metadata(metadata, individual_col, timepoint_col, treatment_col, site):
    df = pd.read_csv(metadata,sep="\t")
    # Drop any rows that are informational
    while df.iloc[0][0].startswith("#"):
	    df.drop(df.index[:1], inplace=True)
	
    # Combine Individual columns if multiple subject identifiers defined (Such as individual and site)
    if "," in individual_col:
        individual_temp = individual_col.split(",")
        df_ind = df[individual_temp[0]]

        # Iterate over the columns and combine on each pass
        for i in range(len(individual_temp)):
            if i == 0:
                pass
            else:
                df_ind = pd.concat([df_ind,df[individual_temp[i]]],axis=1)

        # Hyphenate the white space between identifiers
        for i in range(len(individual_temp)):
            if i == 0:
                pass
            else:
                df_ind = df_ind[individual_temp[0]].astype(str)+"_"+df_ind[individual_temp[1]]

        # Remove any user-generated white space
        inds = []
        for row in df_ind.iteritems():
            inds.append(row[1].replace(" ","-"))
        df_ind = pd.DataFrame({individual_col:inds})
        # Reindex to match dataframes to be merged
        df_ind.index += 1
    else:
        df_ind = df[individual_col]

    df_tp = df[timepoint_col]
    
    #Copy individual column into treatment if none assigned
    if treatment_col is not None:
        df_tx = df[treatment_col]
    else:
        df_tx = df_ind.copy()
        df_tx.columns=[str(df_tx.columns.values[0])+"_tx"]

    # Force timepoint to numerics
    df_tp = df_tp.replace('[^0-9]','', regex=True)


    # Build the return dataframe
    df_ret = pd.concat([df[df.columns[0]], df_ind, df_tp, df_tx], axis=1)

    df_site = pd.DataFrame.from_records(site, columns=["#SampleID","pc1","pc2","pc3"])

    df_ret = pd.merge(df_ret, df_site, on="#SampleID", how='inner')

    if "," in individual_col:
        df_ret.rename(columns = {0:str(individual)}, inplace=True)

    #Now have a full dataframe with sampleIDs, Subject Identifiers, timepoints, treatment, and initial x coordinates.
        #Preserves only values which have initial positions, drops null coordinates.
        #In example files, L3S242 was in metadata File, but not in the PCOA file, so there was one less in the df_ret.
    return df_ret

