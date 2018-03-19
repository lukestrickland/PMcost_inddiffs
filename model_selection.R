your_directory<- "D:/code/proactive_ids"
setwd(your_directory)
setwd("~/proactive_ids")
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")

load("top_samples_LBA.RData")
load("top_samplest_LBA.RData")
load("fixedthreshold_samples_LBA.RData")
load("sdvfixed_samples_LBA.RData")
load("meanvfixed_samples_LBA.RData")

h.IC.dmc(top_samples, DIC=TRUE)
h.IC.dmc(top_samples_1t0, DIC=TRUE)
h.IC.dmc(sdvfixed_samples, DIC=TRUE)
h.IC.dmc(meanvfixed_samples, DIC=TRUE)
h.IC.dmc(fixedthreshold_samples, DIC=TRUE)