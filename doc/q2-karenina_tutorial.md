# QIIME2 Tutorial: [q2-karenina](https://github.com/zaneveld/q2-karenina/)

[Documentation](https://zaneveld.github.io/karenina/html/index.html)

This is a Community Tutorial for [q2-karenina](https://github.com/zaneveld/q2-karenina/), which is an external plugin for use in QIIME2. This plugin fits microbial time-series data using stochastic Ornstein-Uhlenbeck models with the motivation to treat Anna Karenina effects (reviewed in [Zaneveld et al. 2017](https://www.nature.com/articles/nmicrobiol2017121)) in a more mathematical way than before. The specific goal is to fit a simple model that captures the extent to which microbiome changes are due to predictable shifts in community composition versus changes in intrinsic microbiome variability versus the extent to which intrinsic microbiome variability tends to revert to a mean position; for example, due to host factors.

The q2-karenina plugin is in active development and will include additional features upon later releases. This tutorial will cover **fit_timeseries**, which accepts PCoA output and Metadata, and fits the timeseries to an Ornstein-Uhlenbeck model of microbiome change over time.

## Setup

### Installation

q2-karenina requires **Qiime2**, **python 3.5+**, **scipy**, **pandas**, **matplotlib**, and **seaborn**. Since these are fairly ubiquitous packages, let's begin by installing **ffmpeg** (a visualization dependency):

    sudo apt-get install ffmpeg -y

Next, let's install **karenina** and **q2-karenina**:

    pip install git+https://github.com/zaneveld/karenina.git
    pip install git+https://github.com/zaneveld/q2-karenina.git

Now let's update the Qiime2 cache to enable the newly installed plugin:

    qiime dev refresh-cache

Lastly, we ensure that q2-karenina has been properly installed:

    qiime karenina fit-timeseries --help

### Download Data

Next we'll manually download ordination results that will be used in model fitting:

Data for the Moving Pictures of the Human Microbiome dataset [weighted_unifrac_pcoa_results.qza can be found here.](https://github.com/zaneveld/karenina/blob/master/data/fit_timeseries/weighted_unifrac_pcoa_results.qza?raw=true). These data arise from a time-series study of microbiome change across different human body-sites [Caporaso et al., Genome Biology, 2011](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-5-r50)

Synthetic data (simulated under an Ornstein-Uhlenbeck model with known parameters) is used in the second tutorial. [pcoa_results.qza available here](https://github.com/zaneveld/karenina/blob/master/data/fit_timeseries/simulation.qza?raw=true)

## Fitting an Ornstein-Uhlenbeck model to timeseries microbiome data

### What is an Ornstein-Uhlenbeck model?

q2-karenina models microbial communities in PCoA space as if they were physical particles. Microbial community changes reflected by shifts in position in the first three axes of a PCoA plot (PC1,PC2,PC3) are *separately* fit using OU models.


If microbial communities were changing purely randomly that could be described using Brownian motion dx/dt = W*sigma, where W is the Weiner process (effectively a draw from a normal distribution), and sigma scales the velocity of the random displacement (higher sigma = more rapid movement). Ornstein-Uhlenbeck models are a simple extension of Brownian motion that introduce an additional concept:  each particle has a home position (represented by the variable theta) to which it is attracted at each timestep. The strength of that attraction is controlled by lambda.

Mathematically this is represented as follows:

dx/dt = W*sigma + (x - theta)*lambda

Where
dx/dt = change in x with time
W = the random change in position from the Weiner process
sigma = a scaling factor representing intrinsic volatility
theta = the 'home' position to which a particle is attracted
lambda = how strongly the particle is attracted to its home position.

These models have three parameters:
sigma - this describes the intrinsic volatility of a community at each timepoint (higher numbers = more volatile)
theta - this describes an attractor or 'home location' that the process will tend to return to over time.
lambda - the strength with which the particle is attracted back to its home location.

### Example 1: Fitting an Ornstein Uhlenbeck model to the Moving Pictures dataset


#### Motivation

The Moving Pictures dataset captured microbiome from 4 body-sites sampled from 2 individuals at 396 timepoints. Fitting Ornstein-Uhlenbeck models to this dataset can potentially help quantify the extent of stochastic temporal variation across body-sites and between individuals. If one wished to extend the tutorial, one could also segment the data (e.g. by season) to test for whether dynamics differ across body-sites/individuals. 


The dataset is also a good candidate as it has extended timeseries data. Initial benchmarking suggests that the accuracy of the inferred single-individual OU model continues to improve up to about 50 timepoints per individual (fewer timepoints are likely required for fitting cohort data with many individuals, but this has not yet been adeqeutly benchmarked). The q2-longitudinal package recommends at least 5-6 timepoints for a meaningful analysis, and that is probably a reasonable minimal baseline here until benchmarking results for cohort samples are available. 


#### Getting ready

To fit our model, we will use the q2-karenina plugin fit-timeseries. Given a PCoA and metadata for individuals sampled over time, karenina will output a .qzv file at the specified location (NOTE: the .qzv extension will be appended).


You can get a list of the parameters taken by fit-timeseries by running:

```
qiime karenina fit-timeseries --help
``` 

The parameters required for fit_timeseries include:

* **--p-pcoa**: Filepath to PCoA Results (.qza)
* **--p-metadata**: Filepath to metadata file ("None" if within .qza)
* **--p-method**: Global optimization method (only "Basinhopping" supported at this time)
* **--p-individual-col**: Individual column identifier within metadata. Multiple columns can be supplied, comma-separated, and will be combined before splitting data into individuals.
* **--p-timepoint-col**: Timepoint column identifier within metadata
* **--p-treatment-col**: Treatment column identifier within metadata
* **--output-dir**: Output filepath for _visualization.qza_


#### Running fit-timeseries

The Moving Pictures metadata separates their individuals baed on antibiotic usage and bodysite. By combining the individual and bodysite columns into a single identifier, we can track each bodysite, in each individual over time.  This approach is preferred because we strongly expect the various _bodysites_ will have different microbial compositions -- and temporal dynamics --  even if they are part of the same  _individual_. So we define our _individual-col_ as "**Subject,BodySite**", which will combine the two columns from the metadata into a single identifier.

> **Note:** fit_timeseries supports concatenation of individual identifiers with input as a comma-separated list.
> **Note:** in this case, the metadata are already embedded in the OrdinationResults .qza file. If your metadata are instead in an externally supplied metadata source, this can be supplied using the --p-metadata option.  

The cohorts from Moving Pictures contain single-point observations, where the antibiotic treatment is applied at timepoint zero for some of the individuals. Because of this, we cannot estimate cohort models, but can proceed with estimation of individual Ornstein Uhlenbeck models. We control this by defining our _treatment-col_ as "**None**". We also define the _timepoint-col_ as "**DaysSinceExperimentStart**" from the metadata file.

> **Note:** fit_timeseries will generate _only_ individual output if _treatment-col_ is "**None**". If defined, the plugin will generate model parameters for both individuals and treatment cohorts.

The _output-dir_ is defined as **./moving_pictures_fit_ts**, and will contain _visualization.qza_ after a successful run.

To generate a result .qsv visualization folder, we execute the following command from the same directory as the PCoA Results (.qza) file:

```
qiime karenina fit-timeseries --p-timepoint-col DaysSinceExperimentStart  --p-pcoa weighted_unifrac_pcoa_results.qza --p-metadata None --o-visualization ./moving_pictures_fit_timeseries_results --p-individual-col 'Subject,BodySite'  --p-treatment-col None --verbose --p-method basinhopping
```

(Note that the analysis may take a few minutes to run)

A successful visualization provides the following output:

    Saved Visualization to: ./moving_pictures_fit_ts/visualization.qzv

If you'd like to proceed with the tutorial while the results are generating, or to compare your results against what we get, a pre-calculated verison of our result for fitting the model to the Moving Pictures dataset can be found [here](https://github.com/zaneveld/q2-karenina/blob/master/data/moving_pictures_fit_ts), as a .csv text file. i

#### Accessing your results.

All QIIME2 visualizations output to the .qzv format. Under the hood, this is just a zipped folder with the .qsv extension. The folder is structued in a particular way to allow QIIME2 and tools like view.qiime2.org to understand and display the contents in a consistent way.

There are three main ways of viewing your results.

##### Method 1: Using the QIIME2 commandline

The resulting .qza file can be viewed on the commandline as with any other QIIME2 visualization.

One way to do this on the commandline is to use the 'qiime tools view' command. The usage is:

```
qiime tools view [path to your file]
```

So if you set the output folder as ./moving_pictures_timeseries_results above, then to view results for the tutorial file, type:

```
qiime tools view ./moving_pictures_timeseries_results.qzv
```


A webpage should open displaying a  single link to a text file containing the results of the model fit.

##### Method 2: Using view.qiime2.org

Alternatively, you can navigate a web browser to [view.qiime2.org](http://www.view.qiime2.org) and drag the .qzv file where indicated to view the results.


##### Method 3: Manually opening the text file

Finally, if you're having any trouble, it's worth noting that .qzv files are simply zipped directories. It is also possible to unzip them directly and navigate their internal file structure to find what you need (although this isn't necessary in most cases). The results will be in the /data/ folder.


#### Interpreting the Ornstein-Uhlenbeck model results

If you successfully ran the fit-timeseries command and opened the .qzv file, you should now have a result .csv (comma-separated values) text file that can be opened in Microsoft Excel, Google Sheets, etc.

![image|690x463](https://i.imgur.com/RDbiRQu.png) 

Let's take a look.

The results .csv file has columns for describing the metadata parameter(s) on which the model was based. Each row in the results represents an Ornstein-Uhlenbeck model fit to one PC axis, for a particular set of samples. Because the model is fit once per PC axis, a pc column says which PC axis each row. The sigma, theta, and lambda parameters are described above. The n parameters column reflects the number of parameters that were fit. The nLogLik column gives the negative log likelihood of the data given the Ornstein-Uhlenbeck model fit. The AIC column describes the Akake Information Criterion score for the model, which accounts for both the negative log likelihood, and the model complexity (more complex models are penalized because they will tend to obtain better nLogLik scores even in random data).

 

### Example 2: Fitting an Ornstein-Uhlenbeck model to Simulated Data


#### Motivation:
In the Moving Pictures example above, we don't know for sure what the biologically correct answer is, althought the strong differences in PCoA clusters between body-sites suggests we should see different theta values for each body site.

Data simulated under an OU model, or other models provide a chance to try out model fitting when we know what the correct answer should be.

The karenina package provides capabilities for simulating microbiome ordination results for an experiment under various assumptions. In this section of the tutorial we'll show how you could use q2-karenina to fit a model to the results of those simulations. That package also allows export of movies to display simulated datasets (using ffmpeg) 


#### The simulated dataset

The data we will be using simulated an 'Anna Karenina Principle' effect within a group of six individuals, sampled over 50 timepoints. Specifically, it simulates a case in which some perturbation is applied to a treatment group, and reduces their ability to maintain a consistent microbiome. This is implemented in practice by reducing the lambda parameter of the Ornstein-Uhlenbeck model. 


>**Note:** If you'd like to learn more about microbiome community simulation in karenina, a usage guide for the relevant script is available [here](https://zaneveld.github.io/karenina/html/index.html#spatial-ornstein-uhlenbeck) tutorial.


Within the Ornstein-Uhlenbeck model, both control and treatment populations are initially simulated with consistent model parameters.

These are applied to all 3 PC axes for each individual: _sigma_: **0.25**, _lambda_: **0.20**, and _theta_: **0.00**.

 At timepoint 15, the simulation reduces **lambda** to  **0.00** for 3 of the individuals (destabilizing_treatment). 

The outcome produces the following model:

![Animation of AKP effects](https://media.giphy.com/media/ce1T0un8hq8W7v8qcm/giphy.gif)

#### Fitting the simulated datset using a cohort model

> **Note:** If you didn't download both datasets at the beginning, The data used in this tutorial can be found [here](https://github.com/zaneveld/karenina/blob/master/data/fit_timeseries/simulation.qza?raw=true)

Using fit_timeseries, we will attempt to estimate these input parameters by fitting to the treatment cohorts.

After inspection of the metadata file, we define our individual, timepoint, and treatment column identifiers as "**Subject**", "**Timepoint**", and "**Treatment**" respectively, and run the following command from within the same directory as _simulation.qza_:

```
qiime karenina fit-timeseries --p-pcoa ./simulation.qza --p-metadata None --p-method basinhopping --p-individual-col Subject --p-timepoint-col Timepoint --p-treatment-col Treatment --output-dir ./simulation_ou_fit_ts/
```

A successful visualization provides the following output:

    Saved Visualization to: ./simulation_ou_fit_ts/visualization.qzv


This process takes some time, so we've provided the output data [here](https://github.com/zaneveld/q2-karenina/blob/master/data/simulation_ou_fit_ts).

#### View the results

Use one of the methods discribed above to open.

_visualization.qzv_ contains output files with the  cohort model parameters. The cohort output can be viewed below, displaying estimations approaching the input parameters for the simulation described above.

![image|690x141](https://i.imgur.com/liy98gt.png) 



## Coming up!

* **[Method] qiime karenina spatial-ornstein-uhlenbeck**: We are currently developing a Qiime2 method for generation of Ornstein Uhlenbeck simulations which output Q2 PCoA Results, Metadata, Euclidean Distance Matrix, and an output visualization (see [spatial_ornstein_uhlenbeck](https://zaneveld.github.io/karenina/html/index.html#spatial-ornstein-uhlenbeck)).
* **[Visualizer] qiime karenina visualization**: We are currently working on separating the visualization generated within the Ornstein Uhlenbeck simulation to allow for visualization of user PCoA timeseries.
* **Benchmarking of individual time-series fit accuracy***: We are currently benchmarking accuracy across a range of scenarios.
Preliminary results fitting an OU model to simulated data within a single individual/PC axis

>![image|241x250](https://github.com/SLPeoples/karenina/blob/master/data/outputs/benchmarking/benchmark_sigma_err.png?raw=true) 

>![image|241x250](https://github.com/SLPeoples/karenina/blob/master/data/outputs/benchmarking/benchmark_lambda_err.png?raw=true) 

>![image|241x250](https://github.com/SLPeoples/karenina/blob/master/data/outputs/benchmarking/benchmark_theta_err.png?raw=true) 

Please feel free to submit any issues or suggestions as we continue to develop this plugin!
