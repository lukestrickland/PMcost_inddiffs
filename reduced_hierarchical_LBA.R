#Plan is for this file to have all the specification of the LBA models of the data
# 

your_directory<- "~/proactive_ids"
setwd(your_directory)
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")
load("data/dats.RData")


model <- model.dmc(
    factors=list(S=c("nn", "ll"), cond= c("con", "pm")),
    responses=c("N", "L"),
    #Here is where we specify what can vary
    p.map=list(A="1",B=c("cond", "R"),t0=c("1"),mean_v=c("S", "R"),
               sd_v=c("M"),st0="1"),
    #match map scores true and false
    match.map=list(M=list(nn="N", ll="L")),
    #constants: fix an sd_v as a scaling parameter. fix variability in t0 at 0 (standard)
    constants=c(sd_v.false=1, st0=0))

#p.vector will be used to specify prior means
#pick some numbers that seem somewhat sensible
p.vector <-   c(
  t0  = 0.3  ,   
  A = 0.5  ,
  sd_v.true =1,  
  B.con.N =1,         B.pm.N =1 ,       
  B.con.L =1  ,      B.pm.L  =1,            
  
  mean_v.nn.N =1, mean_v.ll.N =0,
  mean_v.nn.L=0,mean_v.ll.L =1)

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
  p2= c(rep(1,1),1,rep(1,1),rep(1,4),rep(1,4)),
#lower bounds of prior. Same param order
#t0 lower bound 0.1 seconds
#starting point 0, sd_v 0, B (which is b-A) 0, rates no truncation
  lower=c(rep(0.1,1),0,rep(0,1),rep(0,4), rep(NA, 4)),
#t0 upper bound 1s, the rest inf
  upper=c(rep(1,1),rep(NA, length(p.vector)-1))
  )

#binds the data and the model together into one object
dm <- data.model.dmc(dats,model)


load("~/proactive_ids/samples/B_samples_LBA.RData")
# # 
##Hierarchical priors
# use same priors for hierarchical means as the individual subject priors
# use gamma distributions for hierarchical sd priors
p1 <- get.p.vector(B_samples[[1]])[names(p.prior)]
#SDs of 2
p2 <- get.p.vector(B_samples[[1]])[names(p.prior)] + 1
s.prior <- prior.p.dmc(p1=p1,p2=p2,
dists=rep("tnorm",length(p1)),
lower= rep(0, length(p1)), upper= rep(Inf, length(p1))
)
pp.prior=list(p.prior,s.prior)
hstart <- make.hstart(B_samples)
theta1 <- make.theta1(B_samples)
h_B_samples <- h.samples.dmc(nmc=60,p.prior,dm,pp.prior,
   hstart.prior=hstart,theta1=theta1,thin=20)

##Request all the cores available on cl2
cores=80

save(h_B_samples, file="samples/startpoints_h_B_LBA.RData")
load("samples/startpoints_h_B_LBA.RData")

#Firstly start the sampling to look for 'stuck' chains- chains that run off
#and get stuck in unlikely regions
# #Note run.unstuck does not give you your valid final samples (p.migrate is on)
# #you have to run without it after
h_B_samples  <- h.RUN.dmc(h_B_samples, cores = cores, max.try=15)
save(h_B_samples,file="samples/h_B_LBA.RData")
