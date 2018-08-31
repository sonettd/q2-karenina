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

import karenina.spatial_ornstein_uhlenbeck as k_OU
from karenina.experiment import Experiment
import pkg_resources
import qiime2
import q2templates
from scipy.spatial import distance
from skbio import OrdinationResults
from q2_types.ordination import PCoAResults
from q2_types.ordination import PCoAResults
from q2_types.distance_matrix import DistanceMatrix
import pandas as pd
import tempfile
import os
        
def spatial_ornstein_uhlenbeck(perturbation_fp:str, treatment_names:str, n_individuals:str, 
                            n_timepoints:int, perturbation_timepoint:int,
                            perturbation_duration:int,interindividual_variation:float,
                            delta:float,lam:float,fixed_start_pos:str) -> (OrdinationResults, DistanceMatrix):
    #pass items to k_OU, 
    #ensure that metadata is saved with both PcOAResults, and DistanceMatrix!
    k_OU.check_perturbation_timepoint(perturbation_timepoint,n_timepoints)
    individual_base_params = {"lambda":lam,"delta":delta,
        "interindividual_variation":interindividual_variation}
    if "None" in fixed_start_pos:
        fixed_start_pos = None
    if fixed_start_pos:
        try:
            x,y,z = map(float,fixed_start_pos.split(","))
            individual_base_params['x']=x
            individual_base_params['y']=y
            individual_base_params['z']=z

        except:
            print ("Supplied value for fixed start position after parsing:",fixed_start_pos)
            raise ValueError('Problem with --fixed_start_pos. Got %s Please supply x,y,z values in the range (-1,1) separated by commas and '+
                'enclosed in quotes. Example: "0.1,-0.2,0.3"'% fixed_start_pos)
    perturbations = k_OU.parse_perturbation_file(perturbation_fp, perturbation_timepoint, perturbation_duration)
    treatments = [[], perturbations]
    treatment_names = treatment_names.split(",")
    n_individuals = list(map(int,n_individuals.split(",")))
    experiment = Experiment(treatment_names,n_individuals,n_timepoints,
        individual_base_params,treatments,interindividual_variation, verbose=False)
    experiment.simulate_timesteps(0,n_timepoints, verbose=False)
    data, ids = experiment.q2_data()
    return _simulation_data(data, ids)
    
def _simulation_data(data, ids):
    with open("ordination.txt","w", encoding='utf8') as ordination:
        ordination.write('Eigvals\t0'+'\n\n')
        ordination.write('Proportion explained\t0'+'\n\n')
        ordination.write('Species\t0\t0\n\n')
        ordination.write('Site\t'+str(len(data)*len(data[0][0]))+'\t3\n')
        dm = {}
        j=0
        for row in data:
            identifier = ids[j]
            for i in range(len(row[0])):
                ordination.write(str(identifier)+"_t"+str(i)+"\t"+str(row[0][i])+"\t"+str(row[1][i])+"\t"+str(row[2][i])+"\n")
                dm.update({str(identifier)+"."+str(i):[row[0][i],row[1][i],row[2][i]]})
            j+=1
        ordination.write("\n")
        ordination.write("Biplot\t0\t0\n\n")
        ordination.write("Site constraints\t0\t0\n")
        ordination_results = OrdinationResults.read("ordination.txt")
    ordination.close
    os.remove("ordination.txt")
    
    # Distance matrix (euclidean)
    dm_0 = []
    dm_0.append("")
    distance_matrix = []
    for key in dm.keys():
        dm_0.append(key)
    distance_matrix.append(dm_0)
    for key in dm.keys():
        dm_1 = []
        dm_1.append(key)
        for key1 in dm.keys():
            dm_1.append(str(distance.euclidean(dm[key],dm[key1])))
        distance_matrix.append(dm_1)

    #Mapping file
    md_0 = ["#SampleID","Subject","Treatment","Timepoint"]
    md_1 = ["#q2:types","categorical","categorical","numeric"]
    md = []
    for id in ids:
        for i in range(len(data[0][0])):
            md.append([id+"_t"+str(i),id,''.join([k for k in id if not k.isdigit()])[:-1],i])
    metadata = [md_0,md_1]
    for row in md:
        metadata.append(row)
    #ADD FUNCTIONALITY TO RETURN MAPPING FILE
    return ordination_results, distance_matrix