library(RColorBrewer)
library(gplots)
library(missMDA)
library(mgcv)
library(lmtest)
library(scales)
library(reticulate)
library(ArcadiaColorBrewer)
library(entropy)
library(vioplot)
library(RNOmni)
np <- import("numpy")
pd <- import("pandas")

###################
##### Functions#####
###################
# Function to compare linear vs. non-linear models for a phenotype matrix
compare_nonlinear <- function(dat,
                              verbose = FALSE,
                              return_models = FALSE,
                              k = 5) {
  # Get all comparisons
  all <- expand.grid(1:ncol(dat), 1:ncol(dat))

  # Get just unique combos
  all <- all[!all[, 1] == all[, 2], ]

  # Remove reciprocal pairings
  for (i in 1:nrow(all)) {
    all[i, ] <- sort(unlist(all[i, ]))
  }
  all <- all[-duplicated(all), ]

  # Convert to list
  all <- split(all, rownames(all))

  # Convert data to numeric matrix
  z <- data.matrix(dat)

  # Test linear vs. nonlinear models via AIC for all traits
  if (verbose == TRUE) {
    
    # Set progress bar
    pb <- txtProgressBar(
      min = 1,
      max = length(all),
      style = 3,
      width = 50,
      char = "."
    )
    counter <- 0

    res <- lapply(all, function(x) {
      counter <<- counter + 1
      setTxtProgressBar(pb, counter)

      mod1 <- lm(z[, x[, 1]] ~ z[, x[, 2]])
      if (length(unique(z[, x[, 2]])) <= 10) {
        mod2 <- gam(z[, x[, 1]] ~ s(z[, x[, 2]],
          k = k
        ))
      } else {
        mod2 <- gam(z[, x[, 1]] ~ s(z[, x[, 2]]))
      }

      out <- lrtest(mod1, mod2)
      a <- AIC(mod1, mod2)

      l <- list(out, a)
      names(l) <- c("lrtest", "AIC")
      l
    })
  } else {
    res <- lapply(all, function(x) {
      mod1 <- lm(z[, x[, 1]] ~ z[, x[, 2]])
      if (length(unique(z[, x[, 2]])) <= 10) {
        mod2 <- gam(z[, x[, 1]] ~ s(z[, x[, 2]],
          k = k
        ))
      } else {
        mod2 <- gam(z[, x[, 1]] ~ s(z[, x[, 2]]))
      }

      out <- lrtest(mod1, mod2)
      a <- AIC(mod1, mod2)

      l <- list(out, a)
      names(l) <- c("lrtest", "AIC")
      l
    })
  }

  # Calculate proportions
  x <- unlist(lapply(res, function(x) abs(signif(x$lrtest$LogLik[1], 3)) - abs(signif(x$lrtest$LogLik[2], 3))))
  linear <- sum(x > 0) / length(x)
  nonlinear <- sum(x < 0) / length(x)

  # Return
  if (return_models == TRUE) {
    r <- list(res, x, linear, nonlinear)
    names(r) <- c("models", "loglik_diffs", "linear", "nonlinear")
    return(r)
  } else {
    r <- list(linear, nonlinear)
    names(r) <- c("linear", "nonlinear")
    return(r)
  }
}

# Function calculate nonlinearity ratio using AIC
aic.ratios <- function(data) {
  a <- c()
  for (i in 1:length(data)) {
    tmp <- unlist(lapply(data[[i]]$models, function(x) abs(x$AIC[2, 2]) - abs(x$AIC[1, 2])))
    a <- c(a, sum(tmp < 0) / length(tmp))
  }
  return(a)
}

# Function to calculate nonlinearity via permutations
nonlinear.permutation <- function(dat,
                                  permutation_number = 10,
                                  data_proportion = 0.2,
                                  counter = FALSE,
                                  ...) {
  # Generate list to save results into
  out <- list()

  # Print counter if verbose
  if (counter == TRUE) {
    pb <- txtProgressBar(
      min = 1,
      max = permutation_number,
      style = 3,
      width = 50,
      char = "."
    )
  }

  # Loop through and calculate
  for (i in 1:permutation_number) {
    if (counter == TRUE) {
      setTxtProgressBar(pb, i)
    }

    # Run compare_nonlinear
    tmp <- compare_nonlinear(dat[, sample(1:ncol(dat), round(ncol(dat) * data_proportion))], ...)
    out[[i]] <- tmp
  }

  # Calculate aic ratios
  aic <- aic.ratios(out)

  # Return
  return(aic)
}

# Function to calculate entropy for a phenotype matrix
overall.entropy <- function(data) {
  entropies <- list()
  for (i in 1:length(data)) {
    # print(paste(i, 'out of', length(data)))
    x <- cov(data[[i]] / max(data[[i]]))
    e <- eigen(x)

    l <- list(
      entropy::entropy.empirical(e$values, unit = "log2"),
      entropy::entropy.empirical(e$values, unit = "log2") / log2(length(e$values)),
      det(x),
      ncol(x)
    )
    names(l) <- c("entropy", "normalized_entropy", "determinant", "n_phenotypes")
    entropies[[i]] <- l
  }

  return(entropies)
}

# Function to calculate mutual information for a phenotype matrix
overall.mi <- function(data) {
  mi <- lapply(data, function(x) {
    m <- infotheo::mutinformation(infotheo::discretize(x, nbins = 100))
    mean(m[lower.tri(m)])
  })
  return(mi)
}

# Function to calculate entropy for subsampled portions of a phenotype matrix
subsample.entropy <- function(data,
                              sample_sizes = seq(0.1, 0.9, 0.1),
                              permutations = 10,
                              numBins = 100,
                              normalize_entropy = FALSE,
                              plot = FALSE) {
  out <- list()
  distributions <- list()

  sample_sizes <- round(sample_sizes * ncol(data))
  sample_sizes <- sample_sizes[sample_sizes > 1]

  for (a in 1:length(sample_sizes)) {

    entropies <- list()
    seqs <- list()
    for (j in 1:permutations) {
      dat <- data[, sample(ncol(data), sample_sizes[a])]

      x <- cov(dat)
      e <- eigen(x)
      if (normalize_entropy == TRUE) {
        ent <- entropy::entropy.empirical(e$values, unit = "log2") / log2(length(e$values))
      } else {
        ent <- entropy::entropy.empirical(e$values, unit = "log2")
      }

      entropies[[as.character(j)]] <- ent
      seqs[[as.character(j)]] <- e$values
    }
    out[[as.character(sample_sizes[a])]] <- entropies
    distributions[[as.character(sample_sizes[a])]] <- seqs
  }

  out <- lapply(out, function(x) unlist(x))

  # Plot if desired
  if (plot == TRUE) {
    vioplot::vioplot(out,
      ylab = "Entropy",
      cex.axis = 1.5,
      cex.lab = 1.5,
      cex.sub = 1.5,
      xlab = "n phenotypes"
    )
  }

  # Return
  l <- list(out, distributions)
  names(l) <- c("entropies", "distributions")
  return(l)
}

# Function to calculate slopes of entropy curves from 'subsample.entropy'
entropy.slopes <- function(data) {
  slopes <-
    lapply(data, function(z) {
      y <- colMeans(do.call(cbind, z$entropies))
      coef(lm(y ~ seq(1, length(y), 1)))[2]
    })
  return(slopes)
}

# Function to calculate all measures
compute.phenotype.stats <- function(data,
                                    run_nonlinear = FALSE,
                                    sample_sizes = seq(0.1, 0.9, 0.1),
                                    entropy_permutations = 10,
                                    verbose = FALSE,
                                    normalize_entropy = FALSE,
                                    ...) {
  # Calculate overall entropy
  if (verbose == TRUE) {
    print("Overall entropy")
  }
  overall_entropy <- overall.entropy(data)

  # Calculate mutual information'
  if (verbose == TRUE) {
    print("Mutual information")
  }
  mutual_information <- overall.mi(data)

  # Calculate subsampled entropy
  if (verbose == TRUE) {
    print("Subsampled entropy")
  }
  subsampled_entropy <- lapply(data, function(x) {
    subsample.entropy(x,
      permutations = entropy_permutations,
      sample_sizes = sample_sizes,
      normalize_entropy = normalize_entropy
    )
  })

  # Calculate slopes
  if (verbose == TRUE) {
    print("Subsampled entropy slopes")
  }
  subsampled_entropy_slopes <- entropy.slopes(subsampled_entropy)

  if (run_nonlinear == TRUE) {
    # Compare non-linearity
    if (verbose == TRUE) {
      print("Nonlinearity")
    }
    nonlinearity <- lapply(data, function(x) {
      compare_nonlinear(x,
        verbose = TRUE,
        ...
      )
    })

    # Calculate aic ratios
    if (verbose == TRUE) {
      print("AIC ratios")
    }
    aic_ratios <- aic.ratios(nonlinearity)

    # Return everything
    l <- list(
      unlist(lapply(overall_entropy, function(x) x$entropy)),
      unlist(lapply(overall_entropy, function(x) x$normalized_entropy)),
      unlist(lapply(overall_entropy, function(x) x$determinant)),
      unlist(lapply(overall_entropy, function(x) x$n_phenotypes)),
      unlist(mutual_information),
      unlist(subsampled_entropy_slopes),
      unlist(aic_ratios),
      nonlinearity,
      subsampled_entropy
    )
    names(l) <- c(
      "entropy",
      "normalized_entropy",
      "determinant",
      "n_phenotypes",
      "mutual_information",
      "subsampled_entropy_slopes",
      "nonlinearity_aic_ratios",
      "nonlinearity_results",
      "subsampled_entropy_results"
    )
    return(l)
  } else {
    # Return everything
    l <- list(
      unlist(lapply(overall_entropy, function(x) x$entropy)),
      unlist(lapply(overall_entropy, function(x) x$normalized_entropy)),
      unlist(lapply(overall_entropy, function(x) x$determinant)),
      unlist(lapply(overall_entropy, function(x) x$n_phenotypes)),
      unlist(mutual_information),
      unlist(subsampled_entropy_slopes),
      subsampled_entropy
    )
    names(l) <- c(
      "entropy",
      "normalized_entropy",
      "determinant",
      "n_phenotypes",
      "mutual_information",
      "subsampled_entropy_slopes",
      "subsampled_entropy_results"
    )
    return(l)
  }
}
