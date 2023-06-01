#Set working directory
setwd('~/Documents/Research/github/accounting-for-nonlinear-phenotypes/')

#Source utility functions
source('01_code/R/nonlinear_phenotype_utils.R')

###################
#####Load data#####
###################
#Load phenotypes
phens = list()
files = list.files('02_output/synthetic-phenotypes/autoencoder_phenotypes/')
files = files[grep('train', files)]
for(i in 1:length(files)){
  tmp = pd$read_pickle(paste('02_output/synthetic-phenotypes/autoencoder_phenotypes/', files[i], sep = ''))
  phens[[files[i]]] = do.call(rbind, tmp)
}

#Load predictions
dat = read.delim('02_output/p_pleio_sweep_pred_out.txt',
                 header = FALSE)
dat5_20 = read.delim('02_output/p_pleio_sweep_pred_5_20.txt',
                     header = FALSE)
dat5_30 = read.delim('02_output/batch_out_5_30_2.txt',
                     header = FALSE)

#Replace 5->20 predictions in original file
#dat = rbind(dat[1:121,], dat[243:363,], dat5_20)
dat = rbind(dat[1:121,], dat[243:363,], dat5_20, dat5_30)

#Add colnames
colnames(dat) = c('p_pleio', 
                  'p_int',
                  'n_phens_predicted',
                  'n_phens_analyzed',
                  'average_mean_percent_error', 
                  'average_pearsonsr_real_pred',
                  'pearsonsr_values',
                  'mean_percent_error_values')

#Make matrix with all error values
all_error = split(dat, as.character(paste(dat$n_phens_predicted, 
                                          dat$n_phens_analyzed,
                                          dat$p_pleio,
                                          dat$p_int, 
                                          sep = '_')))
for(i in 1:length(all_error)){
  x = gsub('\\[', '', all_error[[i]]$mean_percent_error_values)
  x = gsub('\\]', '', x)
  x = unlist(strsplit(x, ', '))
  x = data.frame(p_pleio = rep(all_error[[i]]$p_pleio, length(x)),
                 p_int = rep(all_error[[i]]$p_int, length(x)),
                 n_phens_predicted = rep(all_error[[i]]$n_phens_predicted, length(x)),
                 n_phens_analyzed = rep(all_error[[i]]$n_phens_analyzed, length(x)),
                 mean_percent_error_values = as.numeric(x))
  all_error[[i]] = x
}
all_error = do.call(rbind, all_error)

#Split on n_phenos
dat_nphenos = split(dat, 
                    as.character(paste(dat$n_phens_predicted, 
                                       dat$n_phens_analyzed, sep = '_')))
all_error_nphenos = split(all_error, 
                          as.character(paste(all_error$n_phens_predicted, 
                                             all_error$n_phens_analyzed, sep = '_')))

#Reorder
dat_nphenos = c(dat_nphenos[4], dat_nphenos[1:3])
all_error_nphenos = c(all_error_nphenos[4], all_error_nphenos[1:3])

#Split on p_pleio and p_int
dat_pleio = split(all_error, 
                  as.character(all_error$p_pleio))
dat_int = split(all_error, 
                  as.character(all_error$p_int))

############################
#####Analyze phenotypes#####
############################
#Normalized entropy
toTest = c(5, 10, 20, 30)
ents = list()
for(i in 1:length(toTest)){
  tmp = c()
  for(j in 1:length(phens)){
    dat = phens[[j]][,1:toTest[i]]
    x = cov(dat)
    e = eigen(x)
    tmp = c(tmp, entropy::entropy.empirical(e$values, unit = 'log2')/log2(length(e$values)))
  }
  ents[[as.character(toTest[i])]] = tmp
}

#Split on pleiotropy
pleio = unlist(lapply(strsplit(names(phens), "_"), function(v){v[4]}))
e = lapply(ents, function(x){
  y = split(x, pleio)
  lapply(y, function(z) mean(z))
})

#Plot
plot(unlist(e[[1]]), 
     type = 'l',
     ylim = c(0, 1))
for(i in 2:4){
  lines(unlist(e[[i]]))
}

#Compute all stats
toTest = c(5, 10, 20, 30)
ents = list()
for(i in 1:length(toTest)){
  tmp = list()
  
  print(paste('n phenotypes =', toTest[i]))
  
  pb <- txtProgressBar(min = 1,
                       max = length(phens),
                       style = 3,
                       width = 50,
                       char = ".")
  for(j in 1:length(phens)){
    setTxtProgressBar(pb, j)
    dat = phens[[j]][sample(nrow(phens[[j]]), 600),1:toTest[i]]
    tmp[[j]] = compute.phenotype.stats(list(dat), run_nonlinear = TRUE)
  }
  ents[[as.character(toTest[i])]] = tmp
}

#Extract mean stats as a function of pleiotropy
stats = list()
for(i in 1:length(ents)){
  x = split(ents[[i]], pleio)
  tmp = list()
  for(j in 1:length(x)){
    subsampled_entropy_slopes = mean(unlist(lapply(x[[j]], function(x) unlist(x$subsampled_entropy_slopes))))
    nonlinearity_aic_ratios = mean(unlist(lapply(x[[j]], function(x) unlist(x$nonlinearity_aic_ratios))))
    normalized_entropy = mean(unlist(lapply(x[[j]], function(x) unlist(x$normalized_entropy))))
    mi = mean(unlist(lapply(x[[j]], function(x) unlist(x$mutual_information))))
    l = list(subsampled_entropy_slopes,
             nonlinearity_aic_ratios,
             normalized_entropy,
             mi)
    names(l) = c('subsampled_entropy_slopes',
                 'nonlinearity_aic_ratios',
                 'normalized_entropy',
                 'mutual_information')
    tmp[[names(x)[j]]] = l
  }
  stats[[names(ents)[i]]] = tmp
}

#Get colors
cols = RColorBrewer::brewer.pal(6, 'YlOrBr')[3:6]

#Plot entropy slopes
plot(seq(0, 1, 0.1),
     unlist(lapply(stats[[1]], function(x) x$subsampled_entropy_slopes)),
     type = 'l',
     col = cols[1],
     ylim = c(0, 0.6),
     cex.axis = 1.5,
     cex.lab = 1.5,
     lwd = 1.5,
     ylab = 'Entropy slope',
     xlab = 'P (pleiotropy)')
for(i in 2:4){
  lines(seq(0, 1, 0.1),
        unlist(lapply(stats[[i]], function(x) x$subsampled_entropy_slopes)),
        col = cols[i],
        lwd = 1.5)
}
text(1, 0.575, '5 -> 5', col = cols[1], adj = 1)
text(1, 0.55, '10 -> 5', col = cols[2], adj = 1)
text(1, 0.525, '20 -> 5', col = cols[3], adj = 1)
text(1, 0.5, '30 -> 5', col = cols[4], adj = 1)

#Plot nonlinearity
plot(seq(0, 1, 0.1),
     unlist(lapply(stats[[1]], function(x) x$nonlinearity_aic_ratios)),
     type = 'l',
     col = cols[1],
     ylim = c(0, 1),
     cex.axis = 1.5,
     cex.lab = 1.5,
     lwd = 1.5,
     ylab = 'Entropy slope',
     xlab = 'P (pleiotropy)')
for(i in 2:4){
  lines(seq(0, 1, 0.1),
        unlist(lapply(stats[[i]], function(x) x$nonlinearity_aic_ratios)),
        col = cols[i],
        lwd = 1.5)
}
text(1, 0.975, '5 -> 5', col = cols[1], adj = 1)
text(1, 0.95, '10 -> 5', col = cols[2], adj = 1)
text(1, 0.925, '20 -> 5', col = cols[3], adj = 1)
text(1, 0.9, '30 -> 5', col = cols[4], adj = 1)

#Plot entropy
plot(seq(0, 1, 0.1),
     unlist(lapply(stats[[1]], function(x) x$normalized_entropy)),
     type = 'l',
     col = cols[1],
     ylim = c(0, 1),
     cex.axis = 1.5,
     cex.lab = 1.5,
     lwd = 1.5,
     ylab = 'Normalized entropy',
     xlab = 'P (pleiotropy)')
for(i in 2:4){
  lines(seq(0, 1, 0.1),
        unlist(lapply(stats[[i]], function(x) x$normalized_entropy)),
        col = cols[i],
        lwd = 1.5)
}
text(1, 0.975, '5 -> 5', col = cols[1], adj = 1)
text(1, 0.95, '10 -> 5', col = cols[2], adj = 1)
text(1, 0.925, '20 -> 5', col = cols[3], adj = 1)
text(1, 0.9, '30 -> 5', col = cols[4], adj = 1)

###############################################
#####Compare error across phenotype number#####
###############################################
#Extract mimimum error
m = unlist(lapply(dat_nphenos, function(x) min(x$average_mean_percent_error)))

#Get colors
cols = RColorBrewer::brewer.pal(6, 'YlOrBr')[3:6]

#Plot
plot(m,
     pch = 21,
     bg = cols,
     col = darken_color(cols),
     cex = 2,
     ylim = c(0, 1),
     ylab = 'Minimum % error',
     cex.axis = 1.5,
     cex.lab = 1.5)

#Half violin plots
vioplot::vioplot(lapply(dat_nphenos, function(x) x$average_mean_percent_error),
                 col = cols,
                 side = "right",
                 border = darken_color(cols),
                 ylab = '% error',
                 xlab = '', 
                 font.main = 1,
                 cex.main = 1.5,
                 las = 2, 
                 ylim = c(0, 8),
                 names = names(dat_nphenos),
                 cex.axis = 1.5,
                 cex.lab = 1.5)
abline(h = median(dat_nphenos[[1]]$average_mean_percent_error))
abline(h = median(dat_nphenos[[2]]$average_mean_percent_error))
abline(h = median(dat_nphenos[[3]]$average_mean_percent_error))
abline(h = median(dat_nphenos[[4]]$average_mean_percent_error))

stripchart(lapply(dat_nphenos, function(x) x$average_mean_percent_error),
           col = cols,
           at = seq(0.8, (length(dat_nphenos)-1)+0.8, 1), 
           jitter = 0.1,
           method = "jitter", 
           vertical = TRUE, 
           cex = 1,
           pch = 20, 
           add = TRUE)

#Simple plot as percentage
vioplot::vioplot(lapply(dat_nphenos, function(x) 100-x$average_mean_percent_error),
                 col = cols,
                 side = "right",
                 border = darken_color(cols),
                 ylab = '% accuracy',
                 xlab = '', 
                 font.main = 1,
                 cex.main = 1.5,
                 las = 2, 
                 ylim = c(92, 100),
                 names = names(dat_nphenos),
                 cex.axis = 1.5,
                 cex.lab = 1.5)
stripchart(lapply(dat_nphenos, function(x) 100-x$average_mean_percent_error),
           col = cols,
           at = seq(0.8, (length(dat_nphenos)-1)+0.8, 1), 
           jitter = 0.1,
           method = "jitter", 
           vertical = TRUE, 
           cex = 1,
           pch = 20, 
           add = TRUE)


##########################
#####Plot error space#####
##########################
par(mfrow = c(2,2))
for(i in 1:length(all_error_nphenos)){
  x = all_error_nphenos[[i]]$p_pleio
  y = all_error_nphenos[[i]]$p_int
  h = as.numeric(all_error_nphenos[[i]]$mean_percent_error_values)
  p1 = all_error_nphenos[[i]]$n_phens_predicted[i]
  p2 = all_error_nphenos[[i]]$n_phens_analyzed[i]
  
  m = gam(h ~ s(x, y, k=10))
  
  fit <- predict(m)
  # vis.gam(m, 
  #         se=F,
  #         xlab = 'Prob pleiotropy',
  #         ylab = 'Prob gene interactions',
  #         plot.type = "contour", 
  #         cex.lab = 1.5,
  #         cex.axis = 1.5, 
  #         n.grid = 200,
  #         #color = 'terrain',
  #         zlim = c(0, max(as.numeric(all_error$mean_percent_error_values))),
  #         main = paste(p1, '->', p2),
  #         font.main = 1,
  #         cex.main = 1.5,
  #         type="response")
  
  vis.gam(m, 
          se=F,
          xlab = 'Prob pleiotropy',
          ylab = 'Prob gene interactions',
          plot.type = "persp", 
          zlim = c(0, max(as.numeric(all_error$mean_percent_error_values))),
          cex.lab = 1.5,
          cex.axis = 1.5,
          color = 'topo',
          main = paste(p1, '->', p2),
          zlab = '% non-linear',
          font.main = 1,
          cex.main = 1.5,
          theta = 160, 
          n.grid = 50,
          phi=40,
          axes = TRUE, 
          ticktype = "detailed")
}

###############################################
#####Plot error as function of pleiotropy######
###############################################
#Calculate error mean and variance as a function of pleiotropy
v = list()
means = list()
error = list()
for(i in 1:length(all_error_nphenos)){
  v[[i]] = unlist(lapply(split(all_error_nphenos[[i]], 
                               all_error_nphenos[[i]]$p_pleio), function(x) 
                                 var(x$mean_percent_error_values)/mean(x$mean_percent_error_values)))
  means[[i]] = unlist(lapply(split(all_error_nphenos[[i]], 
                                   all_error_nphenos[[i]]$p_pleio), function(x) 
                                     mean(x$mean_percent_error_values)))
  error[[i]] = unlist(lapply(split(all_error_nphenos[[i]], 
                                   all_error_nphenos[[i]]$p_pleio), function(x) 
                                     plotrix::std.error(x$mean_percent_error_values)))
}

#Get colors
cols = RColorBrewer::brewer.pal(6, 'YlOrBr')[3:6]

#Plot mean
plot(seq(0, 1, 0.1),
     means[[1]],
     type = 'l',
     col = cols[1],
     ylim = c(0, 7),
     cex.axis = 1.5,
     cex.lab = 1.5,
     lwd = 1.5,
     ylab = 'Mean % error',
     xlab = 'P (pleiotropy)')
for(i in 2:4){
  lines(seq(0, 1, 0.1),
        means[[i]],
        col = cols[i],
        lwd = 1.5)
}
text(1, 6.75, '5 -> 5', col = cols[1], adj = 1)
text(1, 6.5, '10 -> 5', col = cols[2], adj = 1)
text(1, 6.25, '20 -> 5', col = cols[3], adj = 1)
text(1, 6.0, '30 -> 5', col = cols[4], adj = 1)

#Plot as half violin plots
par(mfrow = c(2,2))
for(i in 1:4){
  
  tmp = split(all_error_nphenos[[i]]$mean_percent_error_values, 
              all_error_nphenos[[i]]$p_pleio)
  p1 = all_error_nphenos[[i]]$n_phens_predicted[1]
  p2 = all_error_nphenos[[i]]$n_phens_analyzed[1]
  
  vioplot::vioplot(tmp,
                   col = cols[i],
                   side = "right",
                   border = darken_color(cols[i]),
                   ylab = '% error',
                   xlab = '', 
                   main = paste(p2, '->', p1),
                   font.main = 1,
                   cex.main = 1.5,
                   las = 2, 
                   ylim = c(0, 10),
                   names = names(tmp),
                   cex.axis = 1.5,
                   cex.lab = 1.5)
  
  stripchart(tmp,
             col = cols[i],
             at = seq(0.8, (length(tmp)-1)+0.8, 1), 
             jitter = 0.1,
             method = "jitter", 
             vertical = TRUE, 
             cex = 1,
             pch = 20, 
             add = TRUE)
}

#Plot as mean with error
plot(seq(0, 1, 0.1),
     means[[1]],
     type = 'l',
     col = NULL,
     ylim = c(0, 7),
     cex.axis = 1.5,
     cex.lab = 1.5,
     lwd = 1.5,
     ylab = 'Mean % error',
     xlab = 'P (pleiotropy)')
for(i in 1:4){
  polygon(c(seq(0, 1, 0.1),rev(seq(0, 1, 0.1))),
          c(means[[i]]-error[[i]],rev(means[[i]]+error[[i]])),
          col = alpha(cols[i], 0.5), 
          border = FALSE)
  
  lines(seq(0, 1, 0.1),
        means[[i]],
        col = cols[i],
        lwd = 1.5)
}

text(1, 6.75, '5 -> 5', col = cols[1], adj = 1)
text(1, 6.5, '10 -> 5', col = cols[2], adj = 1)
text(1, 6.25, '20 -> 5', col = cols[3], adj = 1)
text(1, 6.0, '30 -> 5', col = cols[4], adj = 1)

#################################################
#####Plot error as a function of interaction#####
#################################################
#Calculate error mean and variance as a function of interaction
v = list()
means = list()
for(i in 1:length(all_error_nphenos)){
  v[[i]] = unlist(lapply(split(all_error_nphenos[[i]], 
                               all_error_nphenos[[i]]$p_int), function(x) 
                                 var(x$mean_percent_error_values)/mean(x$mean_percent_error_values)))
  means[[i]] = unlist(lapply(split(all_error_nphenos[[i]], 
                                   all_error_nphenos[[i]]$p_int), function(x) 
                                     mean(x$mean_percent_error_values)))
}

#Get colors
cols = RColorBrewer::brewer.pal(6, 'YlOrBr')[3:6]

#Plot mean
plot(seq(0, 1, 0.1),
     means[[1]],
     type = 'l',
     col = cols[1],
     ylim = c(0, 7),
     cex.axis = 1.5,
     cex.lab = 1.5,
     lwd = 1.5,
     ylab = 'Mean % error',
     xlab = 'P (interaction)')
for(i in 2:4){
  lines(seq(0, 1, 0.1),
        means[[i]],
        col = cols[i],
        lwd = 1.5)
}
text(1, 6.75, '5 -> 5', col = cols[1], adj = 1)
text(1, 6.5, '10 -> 5', col = cols[2], adj = 1)
text(1, 6.25, '20 -> 5', col = cols[3], adj = 1)
text(1, 6.0, '30 -> 5', col = cols[4], adj = 1)

############################################################
#####Compare n predictor phenotypes using linear models#####
############################################################
#Combine dat_pleio and dat_int
d = list(dat_pleio, dat_int)
names(d) = c('p_pleio', 'p_int')

#Plot
par(mfrow = c(1,2))
for(i in 1:2){
  #Calculate models
  mods = lapply(d[[i]], function(x) summary(lm(x$mean_percent_error_values~x$n_phens_predicted+x$n_phens_analyzed+x$p_int)))
  
  #Extract t-values
  t_values = lapply(mods, function(x) x$coefficients[2,3])
  
  #Plot
  cols = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(t_values)))
  
  plot(names(mods),
       t_values,
       pch = 21,
       bg = cols,
       cex = 1.5,
       ylim = c(-8, 5),
       ylab = 'Effect size (n phenos analyzed)',
       cex.axis = 1.5,
       cex.lab = 1.5,
       xlab = names(d)[i],
       col = darken_color(cols))
  for(j in 1:(length(t_values)-1)){
    segments(as.numeric(names(mods)[j]), 
             t_values[[j]], 
             as.numeric(names(mods)[j+1]), 
             t_values[[j+1]], col = cols[j])
  }
  abline(h = 0, lty = 'dashed', col = 'gray50', lwd = 1.5)
}

######################################################
#####Compare prediction error and phenotype stats#####
######################################################
#Calculate error mean and variance as a function of pleiotropy
v = list()
means = list()
for(i in 1:length(all_error_nphenos)){
  v[[i]] = unlist(lapply(split(all_error_nphenos[[i]], 
                               all_error_nphenos[[i]]$p_pleio), function(x) 
                                 var(x$mean_percent_error_values)/mean(x$mean_percent_error_values)))
  means[[i]] = unlist(lapply(split(all_error_nphenos[[i]], 
                                   all_error_nphenos[[i]]$p_pleio), function(x) 
                                     mean(x$mean_percent_error_values)))
}

#Compare error and entropy
plot(sort(unlist(lapply(stats, function(x) lapply(x, function(y) y$subsampled_entropy_slopes)))),
     sort(unlist(means)),
     xlab = 'Entropy (slope)',
     pch = 21,
     bg = 'gray70',
     col = 'gray50',
     cex = 2,
     ylim = c(0, 6.5),
     cex.axis = 1.5,
     cex.lab = 1.5,
     lwd = 1.5,
     ylab = 'Mean % error')
s = smooth.spline(sort(unlist(lapply(stats, function(x) lapply(x, function(y) y$subsampled_entropy_slopes)))),
                  sort(unlist(means)), spar = 1)
lines(s$x, s$y, lwd = 2, col = 'gray50', lty = 'dashed')

minerva::mine(unlist(lapply(stats, function(x) lapply(x, function(y) y$subsampled_entropy_slopes))),
              unlist(means))

#Zooming in on minimum error and values for each model
plot(unlist(lapply(stats, function(x) lapply(x, function(y) y$subsampled_entropy_slopes)[11])),
     unlist(lapply(means, function(x) x[11])),
     ylim = c(0, 2),
     pch = 21,
     bg = cols,
     col = darken_color(cols),
     cex = 2,
     xlab = 'Entropy (slope)',
     ylab = 'Mean % error',
     cex.lab = 1.5,
     cex.axis = 1.5)




#Model
mod = lm(dat$average_mean_percent_error~dat$p_pleio+dat$p_int+dat$n_phens_predicted:dat$n_phens_analyzed)

#Split on n_phenos
dat = split(dat, paste(dat$n_phens_predicted, dat$n_phens_analyzed, sep = '_'))

x = dat[[1]]$p_pleio
y = dat[[1]]$p_int
h = dat[[1]]$average_mean_percent_error

m = gam(h ~ s(x, y, k=30))
fit <- predict(m)
vis.gam(m, 
        se=F,
        xlab = 'Prob pleiotropy',
        ylab = 'Prob gene interactions',
        plot.type = "contour", 
        cex.lab = 1.5,
        cex.axis = 1.5, 
        #color = 'topo',
        main = 'Accuracy',
        font.main = 1,
        cex.main = 1.5,
        type="response")

#Split each on pleio_int combo
dat = lapply(dat, function(x) split(x, x$p_int))

#Plot
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(dat[[1]])))

plot(dat[[3]][[1]]$p_pleio[order(dat[[3]][[1]]$p_pleio)],
     dat[[3]][[1]]$average_mean_percent_error[order(dat[[3]][[1]]$p_pleio)],
     type = 'l',
     col = cols[1],
     ylab = '% error',
     xlab = 'Pleiotropy (probability)',
     cex.axis = 1.5,
     cex.lab = 1.5,
     ylim = c(0, 8))
for(i in 2:length(dat[[1]])){
  lines(dat[[3]][[i]]$p_pleio[order(dat[[3]][[i]]$p_pleio)],
        dat[[3]][[i]]$average_mean_percent_error[order(dat[[3]][[i]]$p_pleio)],
        col = cols[i])
}




