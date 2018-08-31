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

from karenina.individual import Individual
import karenina.visualization as k_visualization
import pkg_resources
import qiime2
import q2templates
from q2_types.ordination import PCoAResults
import pandas as pd

def visualization(output_dir: str, pcoa : str, metadata : str,
				 individual_col: str, timepoint_col: str, treatment_col: str):
    """
	Connect to karenina.visualization and save output
	
	:param output_dir: output directory
	:param pcoa: pcoa qza file
	:param metadata: metadata file location
	:param individual_col: individual column identifiers
	:param timepoint_col: timepoint column identifier
	:param treatment_col: treatment column identifier
	"""
	
	# Parse in pcoa and metadata as dataframes and inject to k_visualization
    ind = []
    site = _parse_pcoa(pcoa)
    df = _parse_metadata(metadata, individual_col, timepoint_col, treatment_col, site)
    i = 0
    colors = ['fuchsia', 'cyan', 'darkorange', 'blue', 'yellow']
    tx = treatment_col
    treatments = df[tx].unique()
    while len(colors) < len(treatments):
	    colors.append('lightgray')

    for row in df.iterrows():
        curr_subject_id = "%s_%i" % (df[individual_col], i)
        j = 0
        while row[1][3] != treatments[j]:
            j += 1
        color = colors[j]
        params = {'lambda': 0.2, 'delta': 0.25, 'interindividual_variation': 0.01}
        params['color'] = color
        curr_subject = Individual(subject_id = curr_subject_id,
                                  params = params, \
                                  metadata = {treatment_col: df[treatment_col]}, \
                                  interindividual_variation = .01)
        ind.append(curr_subject)
        i += 1

    k_visualization.save_simulation_figure(individuals = ind, 
	                output_folder = output_dir, n_timepoints = 50, 
	                perturbation_timepoint = 25, n_individuals = 50)
    k_visualization.save_simulation_movie(individuals = ind,
	                output_folder = output_dir, n_timepoints = 50, n_individuals = 50)
	
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
        df_tx.columns = [str(df_tx.columns.values[0])+"_tx"]

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