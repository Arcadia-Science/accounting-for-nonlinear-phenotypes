# Acccounting for nonlinear phenotypes

This repository contains analysis code to accompany the publication 'Harnessing genotype-phenotype nonlinearity to accelerate biological prediction'.<br>

Notebooks containing expanded methods and the code + analyses for generating all figures in the publication can be found in [01_code/notebooks/](01_code/notebooks/).<br>

## Directory structure

[01_code/](01_code/) Contains code for simulating phenotypes (Python), analyzing phenotypes (R/Jupyter Notebooks), and predicting phenotypes with an autoencoder (Python).<br>
[02_outputs/](02_outputs/) Outputs associated with the pub.<br>

## Data

All data needed to replicate these analyses are available via [Zenodo]([https://docs.conda.io/en/latest/miniconda.html](https://zenodo.org/record/8298808)).<br>

Data used to study empirical variation across sets of phenotypes are publicly available and associated with the following publications:

*Arabidopsis thaliana*: [Exposito-Alonso et al. 2019](https://www.nature.com/articles/s41586-019-1520-9)<br>
*Saccharomyces cerevisiae*: [Bloom et al. 2019](https://elifesciences.org/articles/49212)<br>
*Caenorhabditis elegans*: [Snoek et al. 2019](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-019-0642-8)<br>
*Mus musculus 1*: [Gonzales et al. 2018](https://www.nature.com/articles/s41467-018-07642-8)<br>
*Mus musculus 2*: [Bogue et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4602074/)<br>
*Drosophila melanogaster*: [Mackay et al. 2012](https://www.nature.com/articles/nature10811)<br>

More details on these datasets are available in the associated pub.

## Packages used in this repo
`ArcadiaColorBrewer v0.0.0.9000` <br/>
`entropy v1.3.1` <br/>
`gplots v3.1.3` <br/>
`lmtest v0.9-40` <br/>
`mgcv v1.8-41` <br/>
`missMDA v1.18` <br/>
`RColorBrewer v1.1-3` <br/>
`reticulate v1.27` <br/>
`scales v1.2.1` <br/>
`vioplot v0.4.0` <br/>
