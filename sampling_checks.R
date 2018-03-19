your_directory<- "D:/code/proactive_ids"
setwd(your_directory)
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")
load("top_samples_LBA.RData")

#Checks convergence and mixing by comparing chains to one another
#SLOW for this many participants
# gelman.diag.dmc(top_samples)
#all <1.1


#We need to inspect the samples to make sure they look
#stationary and generally ok ("fuzzy caterpillars")
#h.RUN.dmc (and run.grid.dmc) do a bunch of automatic
#checks as well

#automaticlaly plot.dmc chooses participant 1
plot.dmc(top_samples)
#have a look at another random one
plot.dmc(top_samples[[50]])

#The trace plots look kind of ok

#Another thing to check is prior influence - 
#did the priors affect the values of the posteriors much
#will need to install the version of the code package that comes with
#dmc for these plots to work
plot.dmc(top_samples[[2]], p.prior=top_samples[[2]]$p.prior)
plot.dmc(top_samples[[4]], p.prior=top_samples[[4]]$p.prior)
plot.dmc(top_samples[[400]], p.prior=top_samples[[400]]$p.prior)

#plotting a few random participants seems to reveal some prior influence:
# the posterior densities (black lines) are not all that far from the 
# prior densities (red lines))
#pretty unsurprising given relatively low trial numbers/ error rates
#the hierarchical model hopefully will sort this out


#dmc sampler grows/shrinks sizes of samples
#in order to get good samples. Check 
# nmc across participants
unlist(lapply(top_samples, function(x) x$nmc))

fast.avgsamples <- function (hsamples) {
  nmcs<- sapply(hsamples, function(x) x$nmc)
  nmc <- min(nmcs)
#Different numbers of nmc for each participant... use the min number and then
  #for participants wtih more randomly sample out that many
  for (i in 1:length(hsamples)) if (nmcs[i] > nmc) hsamples[[i]]$theta <- 
    hsamples[[i]]$theta[,,sample(1:dim(hsamples[[i]]$theta)[3], nmc)]
  avg <- list()
  for (i in 1:length(hsamples)) avg [[i]] <- hsamples[[i]]$theta
  avg2 <- unlist(avg)
  dim3 <- c(dim(avg[[1]]), length(avg2)/prod(dim(avg[[1]])))
  dim(avg2) <- dim3
  out <- apply(avg2, c(1,2,3), mean)
  colnames(out) <- dimnames(avg[[1]])[[2]]
  out
}

avg_samples_top <- fast.avgsamples(h_top_samples2)
avg_samples_topt0 <- fast.avgsamples(top_samples_1t0)
apply(avg_samples_top, 2, mean)
apply(avg_samples_topt0 , 2, mean)
apply(avg_samples, 2, quantile, probs=c(0.1,0.5,0.9))

##expect slight differences due to randomly sampling the thetas down to min length
summary <- summary.dmc(top_samples)
save(summary, file="summary.RData")

# PP <- h.post.predict.dmc(top_samples, save.simulation=TRUE)
# save(PP, file= "top_lba_PP.RData")

load("top_lba_PP.RData")

GET.fitgglist.dmc <- function (
  PP, factors=NA, noR = FALSE,
  quantiles.to.get = c(0.1, 0.5, 0.9), CI= c(0.025, 0.975),
  acc.fun=function(x){as.numeric(x$S)==as.numeric(x$R)},
  correct.only=FALSE,error.only=FALSE
  
) {
  
  sim <- do.call(rbind, PP)
  # Do the same for the data
  data <- lapply(PP, function(x) attr(x, "data"))
  data <- do.call(rbind, data)
  get.fitgglist.dmc (sim,data, factors=factors, noR=noR, quantiles.to.get=quantiles.to.get,
                     CI=CI, acc.fun=acc.fun, correct.only=correct.only, error.only=
                       error.only)
  
}

fit <- GET.fitgglist.dmc(PP, noR=T)
ggplot.RT.dmc(fit[[3]])

fit <- GET.fitgglist.dmc(PP, correct.only=TRUE)
accs <- fit[[1]][fit[[1]]$R=="TRUE",]
accs <- accs[,-3]
ggplot.RP.dmc(accs, xaxis="cond")

#Fit not quite good. Fair enough as we did expect fixed-effects to go bad with these 
#trial numbers, high participant numbers may lead to good hierarchical fits.
