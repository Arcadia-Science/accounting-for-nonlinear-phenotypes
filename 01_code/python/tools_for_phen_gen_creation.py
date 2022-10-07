import numpy as np


def make_genotype(n_as=None,n_loci=None, n_loci_ip=None, n_animals=None,n_phens=None,p_interact=None):

 '''a simple tool for generating phenotypic and genetic data.  Currently,
 the genetic model is purely additive and independent.
 n_as is the number of allelic states at all segregating loci
 n_loci is the number of segregating loci in the analysis
 n_loci_ip is the number of loci influencing the phenotype'''

 if n_as==None: n_as=3
 if n_loci==None: n_loci=3000
 if n_loci_ip==None: n_loci_ip=10
 if n_animals==None: n_animals=500 
 if n_phens==None: n_phens=1
 if p_interact==None: p_interact=0.3

 weights=np.zeros((n_phens,n_loci,n_as))
 for n in range(n_loci_ip):
  for m in range(n_phens):
   weights[m][np.random.randint(n_loci)]=np.random.rand(n_as)

 genotypes,gen_locs=zip(*[make_genotype_ind(n_as,n_loci) for x in range(n_animals)])

 weighted_gens=[]
 for n in range(n_phens):
  weighted_gens.append(genotypes*weights[n])


 inds=[[n for n in range(len(weights[0])) if sum(weights[y][n])>0] for y in range(n_phens)]
 
 
 '''for phn in weighted_gens:
  p_ind=[]
  #for ind in phn:
  for n in range(n_loci):
   loc=ind[n]
   if sum(loc)!=0:
    p_ind.append(n)
  inds.append(p_ind)'''
 

 ''' phens=[]
 for n in range(n_phens):
  phens.append(np.sum(genotypes*weights[n],axis=(2,1))/n_loci)'''

 return genotypes,gen_locs,weights,phens

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
