import numpy as np


def make_genotype(n_as=None,n_loci=None, n_loci_ip=None, n_animals=None,n_phens=None,p_interact=None,p_pleio=None):

 '''a simple tool for generating phenotypic and genetic data.  Currently, this 
 allows for the addition of non-linear gene-gene interactions, but the model is limited.
 An allele of a gene can only interact with one other allele of another gene.
 n_as is the number of allelic states at all segregating loci
 n_loci is the number of segregating loci in the analysis
 n_loci_ip is the number of loci influencing the phenotype'''

 if n_as==None: n_as=3
 if n_loci==None: n_loci=3000
 if n_loci_ip==None: n_loci_ip=10
 if n_animals==None: n_animals=500 
 if n_phens==None: n_phens=1
 if p_interact==None: p_interact=0.0
 if p_pleio==None: p_pleio=0.2

 #make weights for genotypes
 weights=np.zeros((n_phens,n_loci,n_as))
 for n in range(n_loci_ip):
  for m in range(n_phens):
   weights[m][np.random.randint(n_loci)]=np.random.rand(n_as)*10

 #make genotypes
 genotypes,gen_locs=zip(*[make_genotype_ind(n_as,n_loci) for x in range(n_animals)])

 #make weighted genotypes
 weighted_genes=[]
 for n in range(n_phens):
  weighted_genes.append(genotypes*weights[n])

 #make indices for locations in the genome that influence a phenotype
 inds=[[n for n in range(len(weights[0])) if sum(weights[y][n])>0] for y in range(n_phens)]
 
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
      #new_weights=(weights[n][ind_1]*weights[n][ind_2])
      #weights[n][ind_1][allele]=new_weights[allele]
      #weights[n][ind_2][allele]=new_weights[allele]
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


 phens=[]
 for n in range(n_phens):
  #phens.append(np.sum(genotypes*weights[n],axis=(2,1))/n_loci)
  phens.append(np.sum(weighted_genes[n],axis=(2,1)))
 return genotypes,gen_locs,weights,phens,inds,interact

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
