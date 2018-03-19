setwd("D:/Dropbox/Dropbox/shared/proactive_ids")
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")
# load("D:/Dropbox/Dropbox/shared/proactive_ids/h_top_samples3_LBA.RData")
load("h_top_samples3_LBA.RData")
subject.cv <- read.csv("data/covariates.csv")
names(subject.cv)[1] <- "s"
#check that participant ids are ordered the same
any(!(subject.cv$s == as.numeric(names(h_top_samples3))))
colnames(h_top_samples3[[1]]$theta)

# [1] "A"               "B.con.N"         "B.pm.N"         
# [4] "B.con.L"         "B.pm.L"          "t0.con"         
# [7] "t0.pm"           "mean_v.nn.con.N" "mean_v.ll.con.N"
# [10] "mean_v.nn.pm.N"  "mean_v.ll.pm.N"  "mean_v.nn.con.L"
# [13] "mean_v.ll.con.L" "mean_v.nn.pm.L"  "mean_v.ll.pm.L" 
# [16] "sd_v.con.true"   "sd_v.pm.true"    "sd_v.pm.false" 

names(subject.cv)[3] <- "proactive"
names(subject.cv)[2] <- "PM"

thres_control = function(x) (x["B.pm.N"] - x["B.con.N"] + x["B.pm.L"] - x["B.con.L"]) /2
cor_rates_cost =  function(x) (x["mean_v.nn.con.N"] - x["mean_v.nn.pm.N"] + 
                               x["mean_v.ll.con.L"] - x["mean_v.ll.pm.L"]    )/2

##Error rate capacity would presumably go the other way - higher error drift = less capacity
err_rates_cost =  function(x) (x["mean_v.nn.pm.L"] - x["mean_v.nn.con.L"] + 
                                 x["mean_v.ll.pm.N"] - x["mean_v.ll.con.N"]    )/2

t0_cost = function(x) (x["t0.pm"] - x["t0.con"])

sdv_cost <- function(x)   x["sd_v.pm.true"] - x["sd_v.con.true"] 

sdvF_cost <- function(x)    x["sd_v.pm.false"]


cor.test <- cor.plausible(h_top_samples3,p.name="proactive",cv=subject.cv,plot=TRUE,
                          fun = thres_control)
mean(cor.test)

cor.test <- cor.plausible(h_top_samples3,p.name="proactive",cv=subject.cv,plot=TRUE,
                          fun = cor_rates_cost)
mean(cor.test)

cor.test <- cor.plausible(h_top_samples3,p.name="proactive",cv=subject.cv,plot=TRUE,
                          fun = err_rates_cost)
mean(cor.test)

cor.test <- cor.plausible(h_top_samples3,p.name="proactive",cv=subject.cv,plot=TRUE,
                          fun = t0_cost)
mean(cor.test)

cor.test <- cor.plausible(h_top_samples3,p.name="proactive",cv=subject.cv,plot=TRUE,
                          fun = sdv_cost)
mean(cor.test)

cor.test <- cor.plausible(h_top_samples3,p.name="proactive",cv=subject.cv,plot=TRUE,
                          fun = sdvF_cost)
mean(cor.test)


#Threshold cost is the only thing supported at the sample level. Does this remain 
#when we generalise to the population?
cor.r = cor.plausible(h_top_samples3,p.name="proactive",cv=subject.cv,plot=TRUE,
                          fun = thres_control)

kappa <- 1
dens.r <- postRav(n=100, r=cor.r,spacing=.01,kappa=kappa) 

postRav.mean(dens.r)
postRav.p(dens.r,upper=0)


cor.test <- cor.plausible(h_top_samples3,p.name="PM",cv=subject.cv,plot=TRUE,
                          fun = thres_control)
mean(cor.test)

cor.test <- cor.plausible(h_top_samples3,p.name="PM",cv=subject.cv,plot=TRUE,
                          fun = cor_rates_cost)
mean(cor.test)

cor.test <- cor.plausible(h_top_samples3,p.name="PM",cv=subject.cv,plot=TRUE,
                          fun = err_rates_cost)
mean(cor.test)

cor.test <- cor.plausible(h_top_samples3,p.name="PM",cv=subject.cv,plot=TRUE,
                          fun = t0_cost)
mean(cor.test)

cor.test <- cor.plausible(h_top_samples3,p.name="PM",cv=subject.cv,plot=TRUE,
                          fun = sdvF_cost)
mean(cor.test)

#Threshold cost is the only thing supported at the sample level. Does this remain 
#when we generalise to the population?
cor.r = cor.plausible(h_top_samples3,p.name="PM",cv=subject.cv,plot=TRUE,
                      fun = thres_control)

kappa <- 1
dens.r <- postRav(n=100, r=cor.r,spacing=.01,kappa=kappa) 

postRav.mean(dens.r)
postRav.p(dens.r,upper=0)
