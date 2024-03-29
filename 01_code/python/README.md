# Acccounting for nonlinear phenotypes

This directory contains code for simulating phenotypes (tools_for_phen_gen_creation.py) and predicting phenotypes with an autoencoder (autoencoder_denoise_nohup.py and batch_runner.py).

## Usage <br>
### Installing required packages
The version of python3 used is 3.10.12.  Be sure python3 is installed. <br>
The version of PIP3 used is 22.0.2.  Be sure PIP3 is installed. <br>
Then, at the command line type:

```
pip install -r requirements.txt
```

### [tools_for_phen_gen_creation.py](tools_for_phen_gen_creation.py)
This is a script containing a series of functions for simulating phenotypic and genetic data. It is intended to be loaded as a module into an interpreter or called from another script.

For example, you might write a script to simulate a population with 10,000 individuals, 100 phenotypes, 2,000 total genetic loci, 100 loci influencing any individual phenotype and the rest of the parameters set to default.
You could write a script called "my_simulation.py" containing the following text: <br>

```
import tools_for_phen_gen_creation

train_data, test_data = make_genotype(n_animals = 10000, n_phens = 100, n_loci = 2000, n_loci_ip = 100)

pk.dump(train_data, open('training_data.pk', 'rb'))
pk.dump(test_data, open('test_data.pk', 'rb'))
```

### [autoencoder_denoise_nohup.py](autoencoder_denoise_nohup.py)

This is a script containing a denoising autoencoder for predicting phenotypes from phenotypes. It expects to be pointed to a folder containing 2 files formatted in the following way test_[SUFFIX].pk train_[SUFFIX].pk. These files should contain genetic and phenotypic data organized in the format that is created by the tools_for_phen_gen_creation.py functions.
It can be run at the command line as follows:<br>

```
python3 autoencoder_denoise_nohup.py --dataset_path [path to data folder]
```

The output is a single tab separated line organized as follows:
[probability of pleiotropy] [probability of interaction] [number of phenotypes predicted] [number of phenotypes analyzed] [mean absolute percentage error averaged across phenotypes] [pearsons r correlation between real phenotypes and predicted phenotypes averaged across phenotypes] [individual pearsons r for all phenotypes] [individual MAPE for all phenotypes]

### [batch_runner.py](batch_runner.py)
This is a script for parallelizing and running the autoencoder_denoise_nohup.py script. It is the script that was used to create the 30->5, 20->5, 10->5, and 5->5 phenotype predictions that are presented in the pub. <br>
It can be run at the command line as follows:<br>

```
python3 batch_runner.py [path_to_directory_containing_test_and_train_data]
```

## Packages used by these scripts
### Python
`packages that are part of the python standard library are not listed` <br>
`CUDA is only required if a GPU is used` <br>
<br>
`CUDA 11.7 NOTE: unnecessary if the system does not have a GPU` <br>
`numpy 1.21.5` <br>
`Python 3.10.12` <br>
`pip3 22.0.2` <br>
`scipy 1.8.0` <br>
`sklearn 1.3.0` <br>
`torch 1.13.0+cpu` <br>
`torchvision 0.14.0+cpu` <br>

## System
This code has been tested on 2 systems: <br>
`Operating system: Ubuntu 22.04.3 LTS`
`CPU: i7-1260P`
`Memory: 32 Gb` <br>
`Operating system: Ubuntu 22.04.3 LTS`
`CPU: AMD Ryzen 9 5950X`
`Memory: 128 Gb`
`GPU: GeForce RTX 3070` <br>

NOTE: A GPU is unncecessary to run this code.  A GPU provides an acceleration that is the difference between minutes and hours.
