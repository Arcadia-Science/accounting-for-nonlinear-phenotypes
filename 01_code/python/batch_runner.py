import os as os
import sys as sys
from multiprocessing import Pool

'''This is a script to run the denoising autoencoder on a set of inputs and outputs.
 This script can be run by providing a 'taret folder' containing a 'test' and a 'train' data set.
 By default that target folder is:
 /home/dmets/git/accounting-for-nonlinear-phenotypes/02_output/ppleio_pint_sweep_no_noise_no_downsample/
 For the pub, the ppleio_pint_sweep_no_noise_no_downsample/ folder contains the primary data. This data is at https://zenodo.org/record/8298808'''

if len(sys.argv)<2:
 target_folder = "/home/dmets/git/accounting-for-nonlinear-phenotypes/02_output/ppleio_pint_sweep_no_noise_no_downsample/"
else: 
 target_folder = sys.argv[1]

roots = os.listdir(
    target_folder
)

roots = [x.strip("train") for x in roots]
roots = [x.strip("est") for x in roots]
roots = list(set(roots))


def main_5_5(x):
    outpath = target_folder+"batch_out.txt"
    program = "python3 autoencoder_denoise_nohup.py"
    suffix = x
    os.system(
        program
        + " --n_phens_to_predict 5 --n_phens_to_analyze 5 --train_suffix train"
        + x
        + " --test_suffix test"
        + x
        + " >> "
        + outpath
    )


def main_5_30(x):
    outpath = (
        target_folder+"batch_out.txt"
    )
    program = "python3 autoencoder_denoise_nohup.py"
    suffix = x
    os.system(
        program
        + " --n_phens_to_predict 5 --n_phens_to_analyze 30 --train_suffix train"
        + x
        + " --test_suffix test"
        + x
        + " >> "
        + outpath
    )


def main_5_10(x):
    outpath = (
        target_folder+"batch_out.txt"
    )
    program = "python3 autoencoder_denoise_nohup.py"
    suffix = x
    os.system(
        program
        + " --n_phens_to_predict 5 --n_phens_to_analyze 10 --train_suffix train"
        + x
        + " --test_suffix test"
        + x
        + " >> "
        + outpath
    )


def main_5_20(x):
    outpath = (
        target_folder+"batch_out.txt" 
    )
    program = "python3 autoencoder_denoise_nohup.py"
    suffix = x
    os.system(
        program
        + " --n_phens_to_predict 5 --n_phens_to_analyze 20 --train_suffix train"
        + x
        + " --test_suffix test"
        + x
        + " >> "
        + outpath
    )


if __name__ == "__main__":
    pool = Pool(processes=5)
    pool.map(main_5_5, roots)
    pool.map(main_5_30, roots)
    pool.map(main_5_10, roots)
    pool.map(main_5_20, roots)
