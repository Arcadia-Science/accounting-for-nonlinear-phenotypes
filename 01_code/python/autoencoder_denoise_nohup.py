import torch
import math
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torchvision import datasets, transforms
import numpy as np
from torch.autograd import Variable
from argparse import ArgumentParser
from torch.utils.data.dataset import Dataset
import pickle as pk
import time as tm
import scipy as sc
from scipy import stats
from sklearn.metrics import mean_squared_error


# parse commandline arguments
args = ArgumentParser()
args.add_argument(
    "--n_alleles",
    type=int,
    default=3,
    help="number of segregating causal alleles at any given causal locus",
)
args.add_argument(
    "--n_locs",
    type=int,
    default=900,
    help="number of causal loci to model.  This is the same as the number of genes",
)
args.add_argument(
    "--n_env", type=int, default=3, help="number of influential continuous components"
)
args.add_argument("--n_phens", type=int, default=30, help="number of phenotypes")
args.add_argument(
    "--gen_lw",
    type=float,
    default=1,
    help="weight for the loss attributed to genetic features",
)
args.add_argument(
    "--eng_lw",
    type=float,
    default=0.1,
    help="weight for the loss attributed to env features",
)
args.add_argument(
    "--n_epochs", type=int, default=100, help="number of epochs of training"
)
args.add_argument("--batch_size", type=int, default=16, help="batch size")
args.add_argument(
    "--lr_r", type=float, default=0.01, help="reconstruction learning rate"
)
args.add_argument(
    "--b1", type=float, default=0.9, help="adam: gradient decay variables"
)
args.add_argument(
    "--b2", type=float, default=0.999, help="adam: gradient decay variables"
)
args.add_argument("--n_cpu", type=int, default=14, help="number of cpus")
args.add_argument(
    "--e_hidden_dim",
    type=int,
    default=4096,
    help="number of neurons in the hidden layers of encoder",
)
args.add_argument(
    "--d_hidden_dim",
    type=int,
    default=4096,
    help="number of neurons in the hidden layers of decoder",
)
args.add_argument(
    "--batchnorm_momentum",
    type=float,
    default=0.8,
    help="momentum for the batchnormalization layers",
)
args.add_argument(
    "--latent_dim", type=int, default=2048, help="number of neurons in the latent space"
)
args.add_argument(
    "--n_phens_to_analyze", type=int, default=30, help="number of phenotypes to analyze"
)
args.add_argument("--sd_noise", type=float, default=20, help="noise added to phens")
args.add_argument(
    "--n_phens_to_predict", type=int, default=30, help="number of phenotypes to predict"
)
args.add_argument(
    "--dataset_path",
    type=str,
    default=None,
    help="where the train and test files are located",
)
args.add_argument(
    "--train_suffix",
    type=str,
    default="train.pk",
    help="name of the training data file",
)
args.add_argument(
    "--test_suffix", type=str, default="test.pk", help="name of the test data file"
)


vabs = args.parse_args()


# define a torch dataset object
class PhenDataset(Dataset):
    """a class for importing simulated phenotype data.
    It expects a pickled object that is organized as a list of tensors:
    genotypes[n_animals, n_loci, n_alleles] (one hot at allelic state)
    gen_locs[n_animals, n_loci] (index of allelic state)
    weights[n_phens, n_loci, n_alleles] float weight for allelic contribution to phen
    phens[n_animals,n_phens] float value for phenotype
    indexes_of_loci_influencing_phen[n_phens,n_loci_ip] integer indicies of loci that influence a phenotype
    interaction_matrix[FILL THIS IN]
    pleiotropy_matrix[n_phens, n_phens, gen_index]"""

    def __init__(self, data_file, n_phens):
        self.datset = pk.load(open(data_file, "rb"))
        self.phens = self.datset["noisy_phens"]
        self.genotypes = self.datset["genotypes"]
        self.weights = self.datset["weights"]
        self.data_file = data_file
        self.n_phens = n_phens

    def __len__(self):
        return len(self.phens)

    def __getitem__(self, idx):
        phenotypes = torch.tensor(self.phens[idx][: self.n_phens], dtype=torch.float32)
        genotype = torch.tensor(self.genotypes[idx], dtype=torch.float32)
        return phenotypes, genotype


# error definition
def mean_absolute_percentage_error(y_true, y_pred):
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    return np.mean(np.abs((y_true - y_pred) / y_true)) * 100


# convert data to torch.FloatTensor
transform = transforms.ToTensor()

# load the training and test datasets
if vabs.dataset_path == None:
    dataset_path = "/home/dmets/git/accounting-for-nonlinear-phenotypes/02_output/ppleio_pint_sweep_no_noise_no_downsample/"

else:
    dataset_path = vabs.dataset_path

train_dat = vabs.train_suffix
test_dat = vabs.test_suffix

train_data = phen_dataset(dataset_path + train_dat, n_phens=vabs.n_phens_to_analyze)
test_data = phen_dataset(dataset_path + test_dat, n_phens=vabs.n_phens_to_analyze)

# setting device on GPU if available, else CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# how many samples per batch to load
batch_size = vabs.batch_size
num_workers = vabs.n_cpu

# prepare data loaders
train_loader = torch.utils.data.DataLoader(
    dataset=train_data, batch_size=batch_size, num_workers=num_workers, shuffle=True
)
test_loader = torch.utils.data.DataLoader(
    dataset=test_data, batch_size=batch_size, num_workers=num_workers, shuffle=True
)
data_iter = iter(train_loader)


# encoder
class QNet(nn.Module):
    def __init__(self, phen_dim=None, N=None):
        super(Q_net, self).__init__()
        if N == None:
            N = vabs.e_hidden_dim
        if phen_dim == None:
            phen_dim = vabs.n_phens_to_analyze

        batchnorm_momentum = vabs.batchnorm_momentum
        latent_dim = vabs.latent_dim
        self.encoder = nn.Sequential(
            nn.Linear(in_features=phen_dim, out_features=N),
            nn.BatchNorm1d(N, momentum=batchnorm_momentum),
            nn.LeakyReLU(0.01, inplace=True),
            nn.Linear(in_features=N, out_features=latent_dim),
            nn.LeakyReLU(0.01, inplace=True),
        )

    def forward(self, x):
        x = self.encoder(x)
        return x


# decoder
class PNet(nn.Module):
    def __init__(self, phen_dim=None, N=None):
        if N == None:
            N = vabs.d_hidden_dim
        if phen_dim == None:
            phen_dim = vabs.n_phens_to_analyze

        out_phen_dim = vabs.n_phens_to_predict
        n_env = vabs.n_env
        n_allelic_states = vabs.n_locs * vabs.n_alleles
        latent_dim = vabs.latent_dim

        batchnorm_momentum = vabs.batchnorm_momentum

        super(P_net, self).__init__()
        self.decoder = nn.Sequential(
            nn.Linear(in_features=latent_dim, out_features=N),
            nn.BatchNorm1d(N, momentum=batchnorm_momentum),
            nn.LeakyReLU(0.01),
            nn.Linear(in_features=N, out_features=out_phen_dim),
            nn.LeakyReLU(0.01),
        )

    def forward(self, x):
        x = self.decoder(x)
        return x


# main
EPS = 1e-15
Q = QNet()
P = PNet()

# put everything on the GPU if it is there
Q.to(device)
P.to(device)

# Set learning rates
reg_lr = vabs.lr_r

# adam betas
adam_b = (vabs.b1, vabs.b2)

# encode/decode optimizers
optim_P = torch.optim.Adam(P.parameters(), lr=reg_lr, betas=adam_b)
optim_Q_enc = torch.optim.Adam(Q.parameters(), lr=reg_lr, betas=adam_b)

num_epochs = vabs.n_epochs

torch.manual_seed(10)

# train
n_phens = vabs.n_phens_to_analyze
n_phens_pred = vabs.n_phens_to_predict
rcon_loss = []

start_time = tm.time()

for n in range(num_epochs):
    for i, (phens, gen) in enumerate(train_loader):
        phens = phens.to(device)  # move data to GPU if it is there
        batch_size = phens.shape[
            0
        ]  # redefine batch size here to allow for incomplete batches

        # reconstruction loss
        P.zero_grad()
        Q.zero_grad()

        noise_phens = phens + (vabs.sd_noise**0.5) * torch.randn(phens.shape).to(
            device
        )

        z_sample = Q(noise_phens)
        X_sample = P(z_sample)
        recon_loss = F.mse_loss(X_sample + EPS, phens[:, :n_phens_pred] + EPS)
        rcon_loss.append(float(recon_loss.detach()))
        recon_loss.backward()
        optim_P.step()
        optim_Q_enc.step()

    cur_time = tm.time() - start_time
    start_time = tm.time()

#test
phen_encodings = []
phens = []
phen_latent = []

for dat in test_loader:
    ph, gt = dat
    ph = ph.to(device)
    batch_size = ph.shape[0]
    z_sample = Q(ph)
    X_sample = P(z_sample)
    phens += list(ph.detach().cpu().numpy())
    phen_encodings += list(X_sample.detach().cpu().numpy())
    phen_latent += list(z_sample.detach().cpu().numpy())


phens = np.array(phens).T
phen_encodings = np.array(phen_encodings).T

#calculate stats
cors = [
    sc.stats.pearsonr(phens[n], phen_encodings[n])[0]
    for n in range(len(phens[:n_phens_pred]))
]

errs = [
    mean_squared_error(phens[n], phen_encodings[n])
    for n in range(len(phens[:n_phens_pred]))
]

errs = [
    mean_absolute_percentage_error(phens[n], phen_encodings[n])
    for n in range(len(phens[:n_phens_pred]))
]

# output
dats = vabs.train_suffix.strip(".pk")
dats = dats.split("_")
p_pleio = dats[2]
p_int = dats[4]
n_pred = str(vabs.n_phens_to_predict)
n_ana = str(vabs.n_phens_to_analyze)

merrs = str(np.mean(errs))
mcor = str(np.mean(cors))

sep = "\t"
output = (
    p_pleio
    + sep
    + p_int
    + sep
    + n_pred
    + sep
    + n_ana
    + sep
    + merrs
    + sep
    + mcor
    + sep
    + str(cors)
    + sep
    + str(errs)
)
print(output)
