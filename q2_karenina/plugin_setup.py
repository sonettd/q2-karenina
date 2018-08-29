import qiime2.plugin

import q2_karenina
from q2_karenina._spatial_ornstein_uhlenbeck import spatial_ornstein_uhlenbeck
from q2_karenina._fit_timeseries import fit_timeseries
from q2_karenina._visualization import visualization
from q2_types.ordination import PCoAResults
from q2_types.distance_matrix import DistanceMatrix
from qiime2.plugin import Metadata, Str, Choices, Int, Float


plugin = qiime2.plugin.Plugin(
    name='karenina',
    version=q2_karenina.__version__,
    website='https://github.com/zaneveld/karenina',
    package='q2_karenina',
    user_support_text=None,
    description="This script simulates microbiome " +
    "change over time using Ornstein-Uhlenbeck (OU) models.  These are " +
    "similar to Brownian motion models, with the exception that they " +
    "include reversion to a mean. Output is a tab-delimited data table " +
    "and figures.",
    citation_text=None
)
""""""
# Method omitted temporarily until:
    # Output can be in qiime2 Metadata, PcOAResults, and DistanceMatrix format

plugin.methods.register_function(
    function=spatial_ornstein_uhlenbeck,
    inputs = {},
    parameters = {
        'perturbation_fp': Str,
        'treatment_names': Str,
        'n_individuals': Str,
        'n_timepoints': Int,
        'perturbation_timepoint': Int,
        'perturbation_duration': Int,
        'interindividual_variation': Float,
        'delta': Float,
        'lam': Float,
        'fixed_start_pos': Str
    },
    outputs = [
        ('ordination', PCoAResults),
        ('distance_matrix', DistanceMatrix)
    ],
    input_descriptions = {},
    parameter_descriptions = {
        'perturbation_fp': 'filepath for perturbation parameters for simulation results',
        'treatment_names': '[\'control,destabilizing_treatment\'] Names for simulation treatments',
        'n_individuals': '[\'35,35\'] Number of individuals per treatment',
        'n_timepoints': '[\'10\'] Number of simulation timepoints',
        'perturbation_timepoint': '[\'5\'] Timepoint at which to apply treatment (<n_timepoints)',
        'perturbation_duration': '[\'100\'] Duration of perturbation.',
        'interindividual_variation': '[\'0.01\']Starting variability between individuals',
        'delta': '[\'0.25\'] Starting Delta parameter for Brownian Motion/ OU models. '
            +'Higher values indicate greater variability over time',
        'lam': '[\'0.20\'] Starting Lambda value for OU process. '
            +'Higher values indicate a greater tendancy to revert to the mean value.',
        'fixed_start_pos': 'Starting x,y,z position for each point. '
            +'If not defined, starting positions will be randomized based on '
            +'interindividual_variation; type: string, eg: [\'0.0,0.1,0.2\'].'
    },
    output_descriptions = {
        'ordination': 'Sample PCoA file containing simulation data',
        'distance_matrix': 'Sample Distance Matrix containing simulation data'
    },
    name = 'Spatial Ornstein Uhlenbeck microbial community simulation',
    description = ('This method simulates microbial behavior over time using'
                +'Ornstein Uhlenbeck models. This are similar to Brownian Motion'
                +'with the exception that they include reversion to a mean.')
)
""""""
    
# Modify to allow for PCoAResults as input    
plugin.visualizers.register_function(
    function = fit_timeseries,
    inputs = {
        #'pcoa' : PCoAResults
    },
    parameters = {
        'pcoa': Str,
        'method': Str % Choices({'basinhopping'}),
        'metadata': Str,
        'individual_col': Str,
        'timepoint_col': Str,
        'treatment_col': Str
    },
    parameter_descriptions = {
        'pcoa': 'filepath to PCoA results',
        'method': 'global optimization method',
        'metadata': 'filepath to Sample metadata',
        'individual_col': 'individual column identifier',
        'timepoint_col':'timepoint column identifier',
        'treatment_col': 'treatment column identifier'
    },
    name = 'Fit OU Models to PCoA Ordination output',
    description = 'This visualizer generates OU model parameters for PCoA output'
                'data, for each individual and each defined treatment cohort.'
)

plugin.visualizers.register_function(
    function = visualization,
    inputs = {
        #'pcoa' : PCoAResults
    },
    parameters = {
        'pcoa': Str,
        'metadata': Str,
        'individual_col': Str,
        'timepoint_col': Str,
        'treatment_col': Str
    },
    parameter_descriptions = {
        'pcoa': 'filepath to PCoA results',
        'metadata': 'filepath to Sample metadata',
        'individual_col': 'individual column identifier',
        'timepoint_col': 'timepoint column identifier',
        'treatment_col': 'treatment column identifier'
    },
    name = 'Generates 3D animations of PCoA Timeseries',
    description = 'This visualizer generates 3D animations of PCoA Timeseries.'
)
