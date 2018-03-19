# Initial coding:
# accuracy: 0 = incorrect, 1 = correct (which you can delete)
# RT (in seconds)
# block: 0 = control block, 1 = PM block
# stim type (correct response): 0 = nonliving, 1 = living
# response (actual response): 0 = nonliving, 1 = living

your_directory<- "~/proactive_ids"
setwd(your_directory)

dats <- read.csv("data/data.csv")

#Easier levels
dats$stimtype <- factor(dats$stimtype, labels= c("nn", "ll"),
                        levels= c(0, 1))
dats$RESPONSE <- factor(dats$RESPONSE, labels= c("N", "L"),
                        levels= c(0, 1))
dats$Block <- factor(dats$Block, labels= c("con", "pm"),
                        levels= c(0, 1))
#drop acc column
dats <- dats[-2]

#column order expected by dmc (R/RT last two columns, s first)
dats <- dats[,c(1,3,4,5,2)]

#easier colnames/colnames expected by dmc (s,..other factors..,S,R,RT)
names(dats) <- c("s", "cond", "S", "R", "RT")
dats$s <- factor(dats$s)
save(dats, file="data/dats.RData")
