#Plan is for this file to have all the specification of the LBA models of the data
# 






rm(list=ls())
your_directory<- "~/proactive_ids"
setwd(your_directory)
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")
load("data/dats.RData")




h.RUN.dmc <- function(hsamples,
                      cores=1,
                      report=10,
                      p.migrate=.05,
                      h.p.migrate=.05,
                      max.try=100,
                      cut.unstuck=10,
                      cut.flat.location=.5,
                      cut.flat.scale=.5,
                      cut.converge=1.1,
                      split=TRUE,
                      minN=NA,
                      meanN=NA,
                      use.effectiveSize=TRUE,
                      n.add=NA,
                      force=FALSE,
                      thorough=TRUE,
                      verbose=FALSE,
                      gamma.mult=2.38,
                      h.gamma.mult=NA,
                      subjects.to.cores=TRUE,
                      slaveOutfile=NULL)

{

  if ( !is.null(hsamples$theta) ) {

    ### single participant ###

    stop("For a single subject use RUN.dmc")

  } else if (any(names(attributes(hsamples))=="hyper")) {

    ### multiple participants, truly hierarchical ###

    # functions for checking
    h.StuckTests <- function(samples,
                             thorough=TRUE,
                             verbose=FALSE,
                             cut) {

      if (verbose) cat("Stuck chains check\n")

      if (thorough) {
        # check both hyper & participant-level chains
        stucks <- lapply(samples,pick.stuck.dmc,cut=cut)
        ok <- c(hyper=length(pick.stuck.dmc(samples,hyper=TRUE,cut=cut)),
                unlist(lapply(stucks,length)))
      } else {
        # check only hyper chains
        ok <- c(hyper=length(pick.stuck.dmc(samples,hyper=TRUE,cut=cut)))
      }

      fail <- any( ok > 0 )
      if (verbose) {
        if (!fail) cat(": OK\n") else
          cat(paste(":",sum(ok),"\n"))
      }
      fail
    }

    h.FlatTests <- function(samples,
                            thorough=TRUE,
                            p1=1/3,
                            p2=1/3,
                            cut.location=0.25,
                            cut.scale=Inf,
                            verbose=FALSE) {

      gmcmc <- function(samples, thorough) {

        h <- attr(samples, "hyper")
        phi1 <- matrix(aperm(h$phi[[1]],c(1,3,2)),
                       ncol=dim(h$phi[[1]])[2],
                       dimnames=list(NULL, paste0(dimnames(h$phi[[1]])[[2]],
                                                  ".h1")))
        phi2 <- matrix(aperm(h$phi[[2]],c(1,3,2)),
                       ncol=dim(h$phi[[2]])[2],
                       dimnames=list(NULL, paste0(dimnames(h$phi[[2]])[[2]],
                                                  ".h2")))
        m <- cbind(phi1, phi2)

        if (thorough) {

          # add participant-level chains

          mindlist <- lapply(seq_along(samples), function(i) {
            m_i <- matrix(aperm(samples[[i]]$theta,c(1,3,2)),
                          ncol=dim(samples[[i]]$theta)[2],
                          dimnames=list(NULL,
                          paste0(i, ".", dimnames(samples[[i]]$theta)[[2]])))
          })

          mind <- do.call(cbind, mindlist)
          m <- cbind(m, mind)

        }

        mcmc(m)

      }

      mat <- gmcmc(samples, thorough = thorough)
      xlen <- round(dim(mat)[1] * p1)
      ylen <- round(dim(mat)[1] * p2)
      # change in mean relative to robst SD
      m.zs <- apply(mat,2,function(x){
        m1 <- median(x[1:xlen])
        m2 <- median(x[(length(x)-ylen):length(x)])
        abs(m1-m2)/IQR(x)
      })
      names(m.zs) <- paste("m",names(m.zs),sep="_")
      fail <- any(m.zs>cut.location)
      if (!fail) out <- "" else
        out <- paste(names(m.zs)[m.zs==max(m.zs)],"=",round(max(m.zs),2))
      if ( is.finite(cut.scale) ) {
        # Change in IQR relative to overall IQR
        s.zs <- apply(mat,2,function(x){
          m1 <- IQR(x[1:xlen])
          m2 <- IQR(x[(length(x)-ylen):length(x)])
          abs(m1-m2)/IQR(x)
        })
        names(s.zs) <- paste("s",names(s.zs),sep="_")
        if (out != "") out <- paste(out,", ",sep="")
        if (any(s.zs>cut.scale)) out <-
          paste(out,names(s.zs)[s.zs==max(s.zs)],"=",round(max(s.zs),2))
        fail <- fail | any(s.zs>cut.scale)
      }
      if (verbose) {
        cat("Flat check\n")
        print(round(m.zs,2))
        if ( is.finite(cut.scale) )
          print(round(s.zs,2)) else
            if (!fail) cat(": OK\n") else
              cat(paste(":",out,"\n"))
      }
      fail
    }

    h.MixTests <- function(samples,
                           thorough=TRUE,
                           verbose=FALSE,
                           cut,
                           split=TRUE) {

      tmp <- gelman.diag.dmc(samples, hyper=TRUE, split=split)
      names(tmp$mpsrf) <- "hyper_mpsrf"
      gds <- c(tmp$mpsrf,tmp$psrf[,1])

      if (thorough) {
        # add participant-level
        tmp2 <- gelman.diag.dmc(samples, split=split)
        tmp3 <- lapply(seq_along(tmp2), function(i) {
          psrf_i <- tmp2[[i]]$psrf[,1]
          names(psrf_i) <- paste0(i, ".", names(psrf_i))
          mpsrf_i <- tmp2[[i]]$mpsrf
          names(mpsrf_i) <- paste0(i, ".mpsrf")
          return(c(mpsrf_i, psrf_i))
        })
        gds <- c(gds, unlist(tmp3))
      }
      fail <- max(gds) > cut
      if (verbose) {
        cat("Mixing check\n")
        print(round(gds,2))
        if (!fail) cat(": OK\n") else {
          nam <- names(gds)[gds==max(gds)]
          cat(paste(":",nam,"=",round(max(gds),2),"\n"))
        }
      }
      fail
    }

    h.LengthTests <- function(samples,
                              minN,
                              nfun,
                              thorough=TRUE,
                              verbose=FALSE) {

      n <- do.call(nfun,list(h.get.size(samples, thorough = thorough)))

      fail <- n < minN
      if (verbose) {
        cat("Length check")
        if (!fail) cat(": OK\n") else cat(paste(":",n,"\n"))
      }
      fail
    }


    if (!verbose) report <- 1e8

    if (use.effectiveSize) {
      h.get.size <- function(hsamples, thorough) {
        neff <- effectiveSize.dmc(hsamples, hyper = TRUE)
        if (thorough) {
          neff <- c(neff, unlist(effectiveSize.dmc(hsamples)))
        }
        return(neff)
      }
    } else {
      h.get.size <- function(x, thorough){prod(dim(x[[1]]$theta)[-2])}
    }

    if ( !is.na(minN) & !is.na(meanN) ) {
      warning("Both minN and meanN specified, using minN")
      meanN <- NA
    }
    if ( !is.na(minN) ) nfun <- "min"
    if ( !is.na(meanN) ) {
      nfun <- "mean"
      minN <- meanN
    }

    if ( is.na(n.add) ) n <-  ceiling(hsamples[[1]]$nmc/3) else {
      if ( n.add >= floor(hsamples[[1]]$nmc/2) )
        stop(paste("n.flat.test to large, must be less than",
                   floor(hsamples[[1]]$nmc/2)))
      n <- n.add
    }

    if ( any(!is.finite(hsamples[[1]]$log_likelihoods[hsamples[[1]]$nmc,])) )
    { # New samples
      do.migrate=1
      if (verbose) {
        cat("\nGetting initial set of samples\n")
      }
      hsamples <- h.run.dmc(samples=hsamples,report=report,cores=cores,
                            gamma.mult=gamma.mult,h.gamma.mult=h.gamma.mult,
                            p.migrate=p.migrate,h.p.migrate=h.p.migrate)
    } else do.migrate=0 # Samples already good, just want more.

    if ( !force && !is.null(attr(hsamples,"auto")) &&
         !is.na(attr(hsamples,"auto")) && attr(hsamples,"auto")!="GRID FAIL" )
      return(hsamples) # If already sucessfully run then dont run again unless force=TRUE
    n.try <- 0
    repeat {

      # DO TESTS

      # Stuck check
      if ( h.StuckTests(hsamples,
                        thorough=thorough,
                        cut=cut.unstuck,
                        verbose=verbose) )
      { # Start again
        do.migrate <- 1
        get.new <- test.again <- TRUE
      } else get.new <- test.again <- FALSE

      # Stationarity check
      if ( !get.new && any(h.FlatTests(hsamples,
                                       verbose=verbose,
                                       thorough=thorough,
                                       cut.location=cut.flat.location,
                                       cut.scale=cut.flat.scale)) )
      { # remove first 1/3, add new 1/3
        if (verbose) cat("Removing initial 1/3 of samples\n")
        nshift <- ceiling(hsamples[[1]]$nmc/3)
        hsamples <- h.samples.dmc(samples=hsamples,remove=1:nshift,
                                  nmc=0,add=TRUE)
        hsamples <- h.samples.dmc(samples=hsamples,add=TRUE,nmc=n)
        test.again <- TRUE
        get.new <- FALSE
      }

      # Mixing check
      if ( (!get.new & !test.again) &&
           h.MixTests(hsamples,
                      thorough=thorough,
                      verbose=verbose,
                      cut=cut.converge,
                      split=split) ) {
        { # Make longer
          hsamples <- h.samples.dmc(samples=hsamples,nmc=n,add=TRUE)
          test.again <- TRUE
          get.new <- FALSE
        }
      }

      # Length check
      if ( (!get.new & !test.again) &&
           ( !is.na(minN) && h.LengthTests(hsamples,
                                           thorough=thorough,
                                           minN=minN,
                                           nfun=nfun,
                                           verbose=verbose)) ) { # Not long enough
        hsamples <- h.samples.dmc(samples=hsamples,nmc=n,add=TRUE)
        test.again <- TRUE
        get.new <- FALSE
      }

      # Sucess?
      if ( !get.new & !test.again ) {
        outcome <- "SUCESS"
        if (verbose) cat(paste("\nOUTCOME:",outcome,"AFTER TRY",n.try,"\n"))
        break
      }
      # GET MORE SAMPLES
      if ( get.new ) { # Stuck chains, start again
        if (verbose) cat("Getting new set of samples\n")
        hsamples <- h.run.dmc(h.samples.dmc(samples=hsamples,nmc=3*n,
                                            thin=hsamples[[1]]$thin),
                              report=report,cores=cores,
                              gamma.mult=gamma.mult,
                              h.gamma.mult=h.gamma.mult,
                              p.migrate=p.migrate,
                              h.p.migrate=h.p.migrate)
      } else {         # Test failed, update
        if (verbose) cat("Adding extra samples\n")
        hsamples <- h.run.dmc(hsamples,
                              report=report,cores=cores,
                              gamma.mult=gamma.mult,
                              h.gamma.mult=h.gamma.mult,
                              p.migrate=p.migrate*do.migrate)
      }

      # Give up?
      n.try <- n.try + 1
      if ( (n.try > max.try) ) {
        outcome <- "FAIL"
        if (verbose) cat(paste("\nOUTCOME:",outcome,"AFTER TRY",n.try,"\n"))
        break
      } else if (verbose) cat(paste("COMPLETED TRY",n.try,"\n\n"))
    }
    if (outcome=="FAIL") attr(hsamples,"auto") <- NA else
      attr(hsamples,"auto") <- n.try

  } else {

    ### multiple participants, not truly hierarchical ###

    if ( cores == 1 | !subjects.to.cores )
      hsamples <- lapply(hsamples,RUN.dmc,cores=cores,report=report,
                         p.migrate=p.migrate,max.try=max.try,cut.unstuck=cut.unstuck,
                         cut.flat.location=cut.flat.location,cut.flat.scale=cut.flat.scale,
                         cut.converge=cut.converge,split=split,
                         minN=minN,meanN=meanN,use.effectiveSize = use.effectiveSize,
                         verbose=verbose,gamma.mult=gamma.mult) else
    {
      os <- get.os()
      if ( os=="windows" ) {
        require(snowfall,quietly=TRUE)
        require(rlecuyer,quietly=TRUE)
        sfInit(parallel=TRUE, cpus=cores, type="SOCK",slaveOutfile=slaveOutfile)
        sfClusterSetupRNG()
        sfLibrary(msm)
        sfLibrary(rtdists)
        sfLibrary(statmod)
        sfLibrary(pracma)
        sfLibrary(coda)
        sfExportAll()
        if (subjects.to.cores==TRUE) hsamples <-
          sfLapply(hsamples,RUN.dmc,cores=1,report=report,
                   p.migrate=p.migrate,max.try=max.try,cut.unstuck=cut.unstuck,
                   cut.flat.location=cut.flat.location,cut.flat.scale=cut.flat.scale,
                   cut.converge=cut.converge,split=split,
                   minN=minN,meanN=meanN,use.effectiveSize = use.effectiveSize,
                   verbose=verbose,gamma.mult=gamma.mult)
        sfStop()
      } else {
        require(parallel, quietly=TRUE)
        hsamples <- mclapply(hsamples,RUN.dmc,cores=1,report=report,
                             p.migrate=p.migrate,max.try=max.try,cut.unstuck=cut.unstuck,
                             cut.flat.location=cut.flat.location,cut.flat.scale=cut.flat.scale,
                             cut.converge=cut.converge,split=split,
                             minN=minN,meanN=meanN,use.effectiveSize = use.effectiveSize,
                             verbose=verbose,gamma.mult=gamma.mult,mc.cores=cores)
      }
    }

    cat(paste("Number of trys (NA implies fails after max.try =",max.try,")\n"))
    print(unlist(lapply(hsamples,function(x){attr(x,"auto")})))

  }
  hsamples
}
# 
# 
model <- model.dmc(
    factors=list(S=c("nn", "ll"), cond= c("con", "pm")),
    responses=c("N", "L"),
    #Here is where we specify what can vary
    p.map=list(A="1",B=c("cond", "R"),t0=c("1"),mean_v=c("S", "cond", "R"),
               sd_v=c("M"),st0="1"),
    #match map scores true and false
    match.map=list(M=list(nn="N", ll="L")),
    #constants: fix an sd_v as a scaling parameter. fix variability in t0 at 0 (standard)
    constants=c(sd_v.false=1, st0=0))

p.vector <-   c(
  t0  = 0.3  ,
  A = 0.5  ,
  sd_v.true =1,
  B.con.N =1,         B.pm.N =1 ,
  B.con.L =1  ,      B.pm.L  =1,

  mean_v.nn.con.N =1, mean_v.ll.con.N =0,
  mean_v.nn.pm.N =1, mean_v.ll.pm.N =0, mean_v.nn.con.L=0,
  mean_v.ll.con.L =1, mean_v.nn.pm.L =0, mean_v.ll.pm.L =1)


check.p.vector(p.vector, model)


p.prior <-   prior.p.dmc(
  p1=p.vector,
  dists=  rep("tnorm", length(p.vector)),
  p2= c(rep(1,1),1,rep(1,1),rep(1,4),rep(2,8)),
  lower=c(rep(0.1,1),0,rep(0,1),rep(0,4), rep(NA, 8)),
  upper=c(rep(1,1),rep(NA, length(p.vector)-1))
  )


dm <- data.model.dmc(dats,model)
load("~/proactive_ids/fixedt0_samples_LBA.RData")
# # 
##Hierarchical priors
# use same priors for hierarchical means as the individual subject priors
# use gamma distributions for hierarchical sd priors
p1 <- get.p.vector(fixedt0_samples[[1]])[names(p.prior)]
#SDs of 2
p2 <- get.p.vector(fixedt0_samples[[1]])[names(p.prior)] + 1
s.prior <- prior.p.dmc(p1=p1,p2=p2,
dists=rep("tnorm",length(p1)),
lower= rep(0, length(p1)), upper= rep(Inf, length(p1))
)
pp.prior=list(p.prior,s.prior)
# hstart <- make.hstart(fixedt0_samples)
# theta1 <- make.theta1(fixedt0_samples)
# h_fixedt0_samples <- h.samples.dmc(nmc=180,p.prior,dm,pp.prior,
#    hstart.prior=hstart,theta1=theta1,thin=1)

##Request all the cores available on cl2
cores=80

# save(h_fixedt0_samples, file="startpoints_nosv.RData")
# load("startpoints_nosv.RData")

#Firstly start the sampling to look for 'stuck' chains- chains that run off
#and get stuck in unlikely regions
# #Note run.unstuck does not give you your valid final samples (p.migrate is on)
# #you have to run without it after
load("~/proactive_ids/h_nosv_fixedt0_LBA.RData")
h_fixedt0_samples  <- h.RUN.dmc(
  h.samples.dmc(samples=h_fixedt0_samples, nmc=60, thin=20, p.prior=p.prior, pp.prior=pp.prior)
  , cores = cores, max.try=15)
save(h_fixedt0_samples,file="h_nosv_fixedt0_LBA.RData")

# load("h_topunstuck_LBA.RData")
# # #Next run to convergence
# # h_fixedt0_samples1 <- h.run.converge.dmc(h.samples.dmc(nmc=60, samples=h_fixedt0_samples,thin=15), nmc=60,
# #   thorough=TRUE,cores=cores,finalrun=FALSE,max.try=3)
# # save(h_fixedt0_samples1,file="h_fixedt0_samples_LBA.RData")
# # 
# # # load("~/proactive_ids/h_fixedt0_samples_LBA.RData")
# # h_fixedt0_samples2 <- h.run.converge.dmc(h.samples.dmc(nmc=60, samples=h_fixedt0_samples1,thin=15), nmc=60,
# #   thorough=TRUE,cores=cores,max.try=10)
# # save(h_fixedt0_samples2,file="h_fixedt0_samples2_LBA.RData")
# 
# # load("h_fixedt0_samples2_LBA.RData")
# # h_fixedt0_samples3 <- h.run.converge.dmc(h.samples.dmc(nmc=60, samples=h_fixedt0_samples2,thin=1), nmc=60,
# #   thorough=TRUE,cores=cores,finalrun=TRUE, finalI=180, max.try=10)
# # save(h_fixedt0_samples3,file="h_fixedt0_samples3_LBA.RData")
# 
# # load("h_fixedt0_samples3_LBA.RData")
# # h_fixedt0_samples4 <- h.run.dmc(h.samples.dmc(nmc=90, samples=h_fixedt0_samples3,thin=20),cores=cores, report=5)
# # save(h_fixedt0_samples4,file="h_fixedt0_samples4_LBA.RData")
# 
# load("~/proactive_ids/h_fixedt0_samples4_LBA.RData")
# 
#  # cores=1,
#  #                      report=10,
#  #                      p.migrate=.05,
#  #                      h.p.migrate=.05,
#  #                      max.try=100,
#  #                      cut.unstuck=10,
#  #                      cut.flat.location=.5,
#  #                      cut.flat.scale=.5,
#  #                      cut.converge=1.1,
#  #                      split=TRUE,
#  #                      minN=NA,
#  #                      meanN=NA,
#  #                      use.effectiveSize=TRUE,
#  #                      n.add=NA,
#  #                      force=FALSE,
#  #                      thorough=TRUE,
#  #                      verbose=FALSE,
#  #                      gamma.mult=2.38,
#  #                      h.gamma.mult=NA,
#  #                      subjects.to.cores=TRUE,
#  #                      slaveOutfile=NUL
# 
# h_fixedt0_samples5 <- h.RUN.dmc(h.samples.dmc(nmc=60, samples=h_fixedt0_samples4,thin=20), 
#   thorough=TRUE,cores=cores, max.try=5)
# save(h_fixedt0_samples5,file="h_fixedt0_samples5_LBA.RData")
# 
# #t0 fixed over PM condition
# 
# #thresholds fixed over PM condition
# 
# #mean_v fixed over PM condition
# 
# #sd_v fixed over PM condition
# 
# 
# 
#   # 