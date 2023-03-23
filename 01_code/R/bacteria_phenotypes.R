setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')
source('01_code/R/nonlinear_phenotype_utils.R')

#######################################
#####Load and clean phenotype data#####
#######################################
#Set working directory
setwd('~/Desktop/nonlinear_phenotype_review/00_data/bacteria_phenotypes/fitness_browser_data/')

#List directories
files = list.files()

#Loop through and load
phenos = list()
for(i in 1:length(files)){
  print(paste(i, 'out of', length(files)))
  x = read.delim(files[i])
  x = x[x$used == TRUE,]
  image(cov(x[,6:ncol(x)]))
  phenos[[as.character(files[i])]] = as.data.frame(x[,7:ncol(x)])
}

###########################
#####Calculate entropy#####
###########################
#Calculate for each study
entropies = list()
for(i in 1:length(phenos)){
  print(paste(i, 'out of', length(phenos)))
  x = cov(phenos[[i]])
  e = eigen(x)
  
  l = list(entropy::entropy.empirical(e$values, unit = 'log2'),
           entropy::entropy.empirical(e$values, unit = 'log2')/log2(length(e$values)),
           det(x),
           ncol(x))
  names(l) = c('entropy', 'normalized_entropy', 'determinant', 'n_phenotypes')
  entropies[[i]] = l
}

plot(unlist(lapply(entropies, function(x) x$normalized_entropy)),
     unlist(lapply(entropies, function(x) x$n_phenotypes)))

#Calculate entropy over subsamples
all_entropies = list()
for(i in 1:length(phenos)){
  print(names(phenos)[i])
  all_entropies[[names(phenos)[i]]] = subsample.entropy(phenos[[i]],
                                                        numBins = 100,
                                                        permutations = 100,
                                                        normalize_entropy = TRUE,
                                                        sample_sizes = round(seq(0.1, 0.9, 0.05)*ncol(phenos[[i]])),
                                                        plot = TRUE)}

#############################################
#####Calculate phenotypic non-linearity######
#############################################
#Run function on all experiments
res = list()
for(i in 1:length(phenos)){
  print(paste(i, 'out of', length(phenos)))
  res[[as.character(i)]] = compare_nonlinear(phenos[[i]], verbose = TRUE)
}