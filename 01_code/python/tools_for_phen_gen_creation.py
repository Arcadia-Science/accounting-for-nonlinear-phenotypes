import numpy as np


def make_genotype(n_as=None,n_loci=None, n_loci_ip=None, n_env=None, n_animals=None, n_phens=None, env_weight=None, p_interact=None, p_pleio=None, noise=None):

 '''a simple tool for generating phenotypic and genetic data.  Currently, this 
 allows for the addition of non-linear gene-gene interactions, but the model is limited.
 An allele of a gene can only interact with one other allele of another gene. This now
 allows for pleiotropy.  For any individual gene known to influence a phenotype, there is a 
 p_pleio probability that it will also influence any other individual phenotype.  The degree
 of influence is set by p_pleio.
 n_as is the number of allelic states at all segregating loci
 n_loci is the number of segregating loci in the analysis
 n_loci_ip is the number of loci influencing the phenotype'''

 if n_as==None: n_as=3
 if n_loci==None: n_loci=3000
 if n_loci_ip==None: n_loci_ip=10
 if n_animals==None: n_animals=500 
 if n_phens==None: n_phens=30
 if p_interact==None: p_interact=0.1
 if p_pleio==None: p_pleio=0.1
 if n_env==None: n_env=2
 if env_weight==None: env_weight=0.2
 if noise==None: noise=0.1

 out_dct={}
 out_dct['n_as']=n_as
 out_dct['n_loci']=n_loci
 out_dct['n_loci_ip']=n_loci_ip
 out_dct['n_animals']=n_animals
 out_dct['n_phens']=n_phens
 out_dct['p_interact']=p_interact
 out_dct['p_pleio']=p_pleio
 out_dct['n_env']=n_env
 out_dct['env_weight']=env_weight
 out_dct['noise']=noise
 
 #set up n_animals for test and train data sets
 n_animals_train=n_animals  
 n_animals=int(n_animals*1.2)
 n_animals_test=n_animals-n_animals_train

 #set random seed
 np.random.seed(47)

 #make weights for genotypes
 inds=list(np.zeros((n_phens,n_loci_ip),dtype=int))
 weights=np.zeros((n_phens,n_loci,n_as))
 for n in range(n_loci_ip):
  for m in range(n_phens):
   ind=int(np.random.randint(n_loci))
   weights[m][ind]=np.random.rand(n_as)*10
   inds[m][n]=int(ind)
   #weights[m][np.random.randint(n_loci)]=np.random.rand(n_as)*10

 #make indices for locations in the genome that influence a phenotype this is a slopy way of 
 #getting info from the previous simulator...
 #inds=[[n for n in range(len(weights[0])) if sum(weights[y][n])>0] for y in range(n_phens)]

 #add pleiotropy.  For each locus where there is influence on any one phenotype, this 
 #will add an influence on each other phenotype with probability p_pleio.
 pleiotropy_mat=[]
 for m in range(n_phens):
  gens_mat=[]
  for n in list(range(0,m))+list(range(m+1,n_phens)):
   for z in inds[m]:
    if np.random.binomial(1,p_pleio):
     #weights[n][z]=np.random.rand(n_as)*10
     weights[n][z]=weights[m][z]
     gens_mat.append([m,n,z])
  pleiotropy_mat.append(gens_mat)

 #make genotypes
 genotypes,gen_locs=zip(*[make_genotype_ind(n_as,n_loci) for x in range(n_animals)])

 #make weighted genotypes
 weighted_genes=[]
 for n in range(n_phens):
  weighted_genes.append(genotypes*weights[n])
 
 #make gene-gene interactions
 interact=[]
 interacting_loci=[]
 for n in range(n_phens):
  phen_interact=[]
  for m in range(n_loci_ip):
   for z in range(m+1,n_loci_ip): 
    if np.random.binomial(1,p_interact):
     ind_1=inds[n][m]
     ind_2=inds[n][z]
     if ind_1 not in interacting_loci and ind_2 not in interacting_loci:
      interacting_loci.append(ind_1)
      interacting_loci.append(ind_2)
      allele=np.random.randint(3)
      phen_interact.append([ind_1,ind_2,allele])
  interact.append(phen_interact)

 #impose gene-gene interactions
 for n in range(n_phens):
  for m in range(n_animals):
   for intr in interact[n]:
    loc_1=weighted_genes[n][m][intr[0]][intr[2]]
    loc_2=weighted_genes[n][m][intr[1]][intr[2]]
    if loc_1 !=0 and loc_2 !=0:
     new_weight=loc_1*loc_2
     weighted_genes[n][m][intr[0]][intr[2]]=new_weight
     weighted_genes[n][m][intr[1]][intr[2]]=0

 #calculate individual phenotypes
 phens=[]
 for n in range(n_phens):
  phens.append(np.sum(weighted_genes[n],axis=(2,1)))

 #add environmental effects
 env_phens=phens
 env_vects=np.random.rand(n_phens,n_env,n_animals)
 for n in range(n_phens):
  for m in range(n_env):
   phen_vect=phens[n]
   env_vect=env_vects[n][m]
   out_phen=env_vect*(np.mean(phen_vect)/np.mean(env_vect))*env_weight #because this is done iteratively, this results in the first environmental variable having a lower contribution to the total variance than the last environmental variable
   env_phens[n]=out_phen

 #add noise
 noise_vects=np.random.rand(n_phens,n_animals)
 noisy_phens=[]
 for n in range(n_phens):
  phen_vect=env_phens[n]
  noise_vect=noise_vects[n]
  out_phens=noise_vect*(np.mean(phen_vect)/np.mean(noise_vect))*noise
  noisy_phens.append(out_phens)
 
 #format data for output
 out_dct['weights']=weights
 out_dct['inds_of_loci_influencing_phen']=inds
 out_dct['interact_matrix']=interact
 out_dct['pleiotropy_matrix']=pleiotropy_mat

 test_out_dct=out_dct.copy()

 #append genetic data
 out_dct['genotypes']=genotypes[:n_animals_train]
 out_dct['gen_locs']=gen_locs[:n_animals_train]

 test_out_dct['n_animals']=n_animals_test
 test_out_dct['genotypes']=genotypes[n_animals_train:]
 test_out_dct['gen_locs']=gen_locs[n_animals_train:]

 #re-structure phenotype data and append
 phens=list(np.array(phens).T)
 noisy_phens=list(np.array(noisy_phens).T)
 env_phens=list(np.array(env_phens).T)

 out_dct['phens']=phens[:n_animals_train]
 out_dct['noisy_phens']=noisy_phens[:n_animals_train]
 out_dct['env_phens']=env_phens[:n_animals_train]

 test_out_dct['phens']=phens[n_animals_train:]
 test_out_dct['noisy_phens']=noisy_phens[n_animals_train:]
 test_out_dct['env_phens']=env_phens[n_animals_train:]

 

 #return genotypes,gen_locs,weights,phens,inds,interact,pleiotropy_mat
 return out_dct, test_out_dct


def make_genotype_ind(n_as,n_loci):
 '''utlity function for make_genotype.  Creates a random set of genotypes
 for an individual.'''
 gen=np.zeros((n_loci,n_as))
 locs=[]
 for n in range(n_loci):
  loc=np.random.randint(n_as)
  gen[n][loc]=1
  locs.append(loc)
 return gen,locs


