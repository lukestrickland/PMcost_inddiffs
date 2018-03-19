#Plan is for this file to have all the specification of the LBA models of the data
# - i.e., the 'top' model and then the reduced models with parameters fixed

your_directory<- "~/proactive_ids"
setwd(your_directory)
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")
load("data/dats.RData")

###First I will specify a 'top' model that allows all 
#major parameters to vary by PM:
#B, v, sv, t0

#This model.dmc function creates a big array that will assign parameters
# to the appropriate cell of the design
#check the array (parameters x design cell x accumulator) to check it worked 
model <- model.dmc(
    factors=list(S=c("nn", "ll"), cond= c("con", "pm")),
    responses=c("N", "L"),
    #Here is where we specify what can vary
    p.map=list(A="1",B=c("cond", "R"),t0=c("cond"),mean_v=c("S", "cond", "R"),
               sd_v=c("cond", "M"),st0="1"),
    #match map scores true and false
    match.map=list(M=list(nn="N", ll="L")),
    #constants: fix an sd_v as a scaling parameter. fix variability in t0 at 0 (standard)
    constants=c(sd_v.con.false=1, st0=0))

#p.vector will be used to specify prior means
#pick some numbers that seem somewhat sensible
p.vector <-   c(
  t0.con  = 0.3  ,  t0.pm   = 0.3 ,    
  A = 0.5  ,
  sd_v.con.true =1,  sd_v.pm.true =1,   sd_v.pm.false  =1,
  B.con.N =1,         B.pm.N =1 ,       
  B.con.L =1  ,      B.pm.L  =1,            
  
  mean_v.nn.con.N =1, mean_v.ll.con.N =0,
  mean_v.nn.pm.N =1, mean_v.ll.pm.N =0, mean_v.nn.con.L=0,
  mean_v.ll.con.L =1, mean_v.nn.pm.L =0, mean_v.ll.pm.L =1)

#check.p.vector checks that the p.vector has all the parameters from the model
check.p.vector(p.vector, model)

#specify a prior with prior.p.dmc.
p.prior <-   prior.p.dmc(
#names and posterior means set by p1
  p1=p.vector,
#truncated normals  
  dists=  rep("tnorm", length(p.vector)),
#standard deviations of prior 
# param order :t0, A, sv, B, mv            
  p2= c(rep(1,2),1,rep(1,3),rep(1,4),rep(2,8)),
#lower bounds of prior. Same param order
#t0 lower bound 0.1 seconds
#starting point 0, sd_v 0, B (which is b-A) 0, rates no truncation
  lower=c(rep(0.1,2),0,rep(0,3),rep(0,4), rep(NA, 8)),
#t0 upper bound 1s, the rest inf
  upper=c(rep(1,2),rep(NA, length(p.vector)-2))
  )

#binds the data and the model together into one object
dm <- data.model.dmc(dats,model)
#generates start points for the sampling process out of the prior.
#We will save off samples here, then deploy them on the newcastle supercomputing
#grid (using grid_dispatch.R)
# top_samples <- h.samples.dmc(nmc = 180,p.prior,dm)

# #load fixed-effects fitting
# load("top_samples_LBA.RData")
# 
# ##Hierarchical priors
# # use same priors for hierarchical means as the individual subject priors
# # use gamma distributions for hierarchical sd priors
# p1 <- get.p.vector(top_samples[[1]])[names(p.prior)]
# s.prior <- prior.p.dmc(p1=p1,p2=p1,
# dists=rep("gamma",length(p1)))
# pp.prior=list(p.prior,s.prior)
# hstart <- make.hstart(top_samples)
# theta1 <- make.theta1(top_samples)
#h_top_samples <- h.samples.dmc(nmc=180,p.prior,dm,pp.prior,
 #    hstart.prior=hstart,theta1=theta1,thin=1)

##Request all the cores available on cl2
cores=64

#Firstly start the sampling to look for 'stuck' chains- chains that run off
#and get stuck in unlikely regions
#Note run.unstuck does not give you your valid final samples (p.migrate is on)
#you have to run without it after
# h_top_samples  <- h.run.unstuck.dmc(h_top_samples, p.migrate = .05, cores = cores)
# save(h_top_samples,file="h_topunstuck_LBA.RData")
# load("h_topunstuck_LBA.RData")
# #Next run to convergence
# h_top_samples1 <- h.run.converge.dmc(h.samples.dmc(nmc=60, samples=h_top_samples,thin=15), nmc=60,
#   thorough=TRUE,cores=cores,finalrun=FALSE,max.try=3)
# save(h_top_samples1,file="h_top_samples_LBA.RData")

# load("~/proactive_ids/h_top_samples_LBA.RData")
# h_top_samples2 <- h.run.converge.dmc(h.samples.dmc(nmc=60, samples=h_top_samples1,thin=15), nmc=60,
#   thorough=TRUE,cores=cores,max.try=10)
# save(h_top_samples2,file="h_top_samples2_LBA.RData")
# 
# h_top_samples3 <- h.run.converge.dmc(h.samples.dmc(nmc=60, samples=h_top_samples2,thin=15), nmc=60,
#   thorough=TRUE,cores=cores,finalrun=TRUE, finalI=180, max.try=10)
# save(h_top_samples3,file="h_top_samples3_LBA.RData")



p1 <- get.p.vector(top_samples[[1]])[names(p.prior)]
s.prior <- prior.p.dmc(p1=p1,p2=p1,
dists=rep("gamma",length(p1)))
pp.prior=list(p.prior,s.prior)
hstart <- make.hstart(top_samples)
theta1 <- make.theta1(top_samples)
h_top_samples <- h.samples.dmc(nmc=180,p.prior,dm,pp.prior,
   hstart.prior=hstart,theta1=theta1,thin=1)

#t0 fixed over PM condition

#thresholds fixed over PM condition

#mean_v fixed over PM condition

#sd_v fixed over PM condition



  # 