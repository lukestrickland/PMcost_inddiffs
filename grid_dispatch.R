rm(list=ls())
setwd("~/proactive_ids")
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")

run.grid.dmc("B_samples_LBA",model.dir ="LBA", 
             model.file="lba_B.R",user="ljs392",n.add=60, wall.hours = 300,
             GB = 1, max.try=1000)

rm(list=ls())
setwd("~/proactive_ids")
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")

run.grid.dmc("V_samples_LBA",model.dir ="LBA", 
             model.file="lba_B.R",user="ljs392",n.add=60, wall.hours = 300,
             GB = 1, max.try=1000)

rm(list=ls())
setwd("~/proactive_ids")
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")

run.grid.dmc("t0_samples_LBA",model.dir ="LBA", 
             model.file="lba_B.R",user="ljs392",n.add=60, wall.hours = 300,
             GB = 1, max.try=1000)

rm(list=ls())
setwd("~/proactive_ids")
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")

run.grid.dmc("sv_samples_LBA",model.dir ="LBA", 
             model.file="lba_B.R",user="ljs392",n.add=60, wall.hours = 300,
             GB = 1, max.try=1000)

