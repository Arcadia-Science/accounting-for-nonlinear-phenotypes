import numpy as np


def make_genotype(n_as=None,n_loci=None, n_loci_ip=None, n_animals=None,n_phens=None,p_interact=None,p_pleio=None):

 '''a simple tool for generating phenotypic and genetic data.  Currently, this 
 allows for the addition of non-linear gene-gene interactions, but the model is limited.
 An allele of a gene can only interact with one other allele of another gene. This now
 allows for pleiotropy.  For any individual gene known to influence a phenotype, there is a 
 p_pleio probability that it will also influence any other individual phenotype.  The degree
 of influence is random.
 n_as is the number of allelic states at all segregating loci
 n_loci is the number of segregating loci in the analysis
 n_loci_ip is the number of loci influencing the phenotype'''

 if n_as==None: n_as=3
 if n_loci==None: n_loci=3000
 if n_loci_ip==None: n_loci_ip=10
 if n_animals==None: n_animals=500 
 if n_phens==None: n_phens=1
 if p_interact==None: p_interact=0.1
 if p_pleio==None: p_pleio=0.1

 #make weights for genotypes
 weights=np.zeros((n_phens,n_loci,n_as))
 for n in range(n_loci_ip):
  for m in range(n_phens):
   weights[m][np.random.randint(n_loci)]=np.random.rand(n_as)*10

 #make indices for locations in the genome that influence a phenotype this is a slopy way of 
 #getting info from the previous simulator...
 inds=[[n for n in range(len(weights[0])) if sum(weights[y][n])>0] for y in range(n_phens)]

 #add pleiotropy.  For each locus where there is influence on any one phenotype, this 
 #will add an influence on each other phenotype with probability p_pleio
 pleiotropy_mat=[]
 for m in range(n_phens):
  gens_mat=[]
  for n in list(range(0,m))+list(range(m+1,n_phens)):
   for z in inds[m]:
    if np.random.binomial(1,p_pleio):
     weights[n][z]=np.random.rand(n_as)*10
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
   for int in interact[n]:
    loc_1=weighted_genes[n][m][int[0]][int[2]]
    loc_2=weighted_genes[n][m][int[1]][int[2]]
    if loc_1 !=0 and loc_2 !=0:
     new_weight=loc_1*loc_2
     weighted_genes[n][m][int[0]][int[2]]=new_weight
     weighted_genes[n][m][int[1]][int[2]]=0


 #calculate individual phenotypes
 phens=[]
 for n in range(n_phens):
  phens.append(np.sum(weighted_genes[n],axis=(2,1)))
 return genotypes,gen_locs,weights,phens,inds,interact,pleiotropy_mat

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
