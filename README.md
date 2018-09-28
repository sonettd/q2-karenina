# q2-karenina
[![Build Status](https://travis-ci.org/zaneveld/q2-karenina.svg?branch=master)](https://travis-ci.org/zaneveld/q2-karenina)

**Tutorial**
* [Tutorial](https://github.com/zaneveld/q2-karenina/blob/master/doc/q2-karenina_tutorial.md)

**Other links**
* [karenina](https://github.com/zaneveld/karenina)
* [Documentation](https://zaneveld.github.io/karenina/html/index.html)
* [QIIME 2](https://qiime2.org)

### Approach

q2-karenina in a QIIME2 plug-in for the karenina python package, which provides simulation and model-fitting tools for studying stochastic variability in microbial communities. The package is named for Anna Karenina Principle effects, in which certain diseases and environmental stressors can cause increased stochastic variability in animal microbiomes. A recent review of this phenomenon can be found in [Zaneveld et al.,Nature Microbiology, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28836573). This pattern is in contrast to stressors that shift the mean *position* of microbial communities in PCoA space.

The approach of the q2-karenina package to modeling microbiome dynamics is top down: microbial communities in PCoA space are modeled as if they were physical particle in Euclidean space. For example, a microbiome that shifts in composition by a small amount at each time-step can be modeled like the random walk of a particle undergoing Brownian Motion. Most healthy animal microbiomes do not appear to wander randomly like this. Instead, they tend to revert over time to a composition that is average for a particular individual. q2-karenina models this phenomenon using Ornstein-Uhlenbeck models, which include both Brownian motion and reversion to a particular home position. 

### Model

This model is fairly simple, having only 3 tunable parameters per PC axes. (A separate model is fit to the  motion of a microbial community on each PC axis separately, as many papers observe increased variance along some axes but not others).

Mathematically, the model is as follows:

dx/dt = dW * sigma + (x - theta)* lambda

Where:

  dx/dt = change in the position of a microbial community on a particular axis of the PCoA plot 

  dW = the Weiner process which generates Brownian motion (effectively drawing the extent of displacement from a normal distribution)

  sigma = a constant that scales the Brownian motion component

  theta = the 'home position' or attractor for the community. Over time (in the absence of disturbance or e.g. seasonal effects that alter the mean composition) it is expected that this position will be the mean composition of microbiome samples from an individual. It is also expected that in a stable community this will be the centroid of samples taken from an individual longitudinally.

### Release

The initial q2-karenina release allows for fitting of Ornstein-Uhlenbeck models to QIIME2 ordination results using the fit_spatial_ornstein_uhlenbeck visualizer, and returns the parameters for the model fit. Future development will wrap other functionality in karenina such as simulation of QIIME2-compatible OrdinationResult and Metadata files based on Ornstein-Uhlenbeck model parameters, benchmarks of the model fitting procedure, and visualization.

Models can either be fit to timeseries from each subject separately (individual fit), or can be fit to all individuals within a treatment collectively (cohort fit) using the spatial_ornstein_uhlenbeck.py script

**Inputs**

  QIIME2 Ordination Results

  QIIME2 Metadata (if metadata are not embedded in Ordination Results)

In practice, the sample should have multiple time-points per individual.
