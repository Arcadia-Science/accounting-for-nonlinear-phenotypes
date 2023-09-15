# Acccounting for nonlinear phenotypes

This repository contains analysis code to accompany the publication 'Harnessing genotype-phenotype nonlinearity to accelerate biological prediction'.<br>

This directory contains code for simulating phenotypes (tools_for_phen_gen_creation.py) and predicting phenotypes with an autoencoder (autoencoder_denoise_nohup.py and batch_runner.py).

##Usage
#tools_for_phen_gen_creation.py
This is a script containing a series of functions for simulating phenotypic and genetic data. It is intended to be loaded as a module into an interpreter or called from another script.

#autoencoder_denoise_nohup.py
This is a script containing a denoising autoencoder network for predicting phenotypes from phenotypes. It can be run at the command line as follows:
python3 autoencoder_denoise_nohup.py


## packages used by these scripts
### Python
`packages that are part of the python standard library are not listed` <br>
`argparse 1.1` <br>
`CUDA 11.7` <br>
`numpy 1.21.5` <br>
`Python 3.10.12` <br>
`scipy 1.8.0` <br>
`sklearn 1.3.0` <br>
`torch 1.13.0` <br>
`torchvision 0.14.0` <br>

