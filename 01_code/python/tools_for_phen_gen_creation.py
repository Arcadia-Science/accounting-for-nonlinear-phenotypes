import numpy as np

def make_genotype(n_as=None,n_loci=None, n_loci_ip=None, n_animals=None):
 '''a simple tool for generating phenotypic and genetic data.  Currently,
 the genetic model is purely additive. 
 n_as is the number of allelic states at all segregating loci
 n_loci is the number of segregating loci in the analysis
 n_loci_ip is the number of loci influencing the phenotype'''

 if n_as==None: n_as=3
 if n_loci==None: n_loci=3000
 if n_loci_ip==None: n_loci_ip=10
 if n_animals==None: n_animals=500 

 loci_that_matter=np.zeros(n_loci)
 for n in range(n_alleles_that_matter):
  loci_that_matter[np.random.randint(n_loci)]=1

 weights=np.zeros((n_loci,n_allelic_states))
 for n in range(n_alleles_that_matter):
  weights[np.random.randint(n_loci)]=np.random.rand(n_allelic_states)

 genotypes,gen_locs=zip(*[make_genotype(n_allelic_states,n_loci) for x in range(n_animals)])

 phens=np.sum(genotypes*weights,axis=(2,1))/n_loci
 return genotypes,gen_locs,phens

def make_genotype_ind(n_as,n_loci):
 '''utlity function for make_genotype'''
 gen=np.zeros((n_loci,n_as))
 locs=[]
 for n in range(n_loci):
  loc=np.random.randint(n_as)
  gen[n][loc]=1
  locs.append(loc)
 return gen,locs
