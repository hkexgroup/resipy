# %  Implementation of linear mixed effect error model in R for pyR2
# %  cite: Tso C.-H.M. et al. (2017) Journal of Applied Geophysics (doi:10.1016/j.jappgeo.2017.09.009)
# %  Michael Tso, Lancaster University, March 2019
# %
# %  This script fits to LME model to data and write to protocal.dat for R2
# %  
# %  Note 1: needs statstical and machine learning toolbox
# %  Note 2: the "fitlme" functions requires data to be put in the "table"
# %  format (type on command window: help table)
# %  Note 3: we have assumed negligible forward modeling errors. They can be
# %  added to prescribed error levels (column 7 of protocol.dat)

library(lme4) # install in terminal: e.g. conda install r-lme4

#script.dir <- dirname(sys.frame(0)$ofile)
#setwd(script.dir)
setwd("resipy")
print(getwd())


## fitting
fname = file.path("invdir","protocol-lmeInRecip.dat")
fid <- file(fname,"r")
nmeas <- strtoi(readLines(fid,n=1))
R2data <- scan(fname, skip=1, list(i=0, p1=0, p2=0, c1=0, c2=0, avgR=0, Err=0))
close(fid)
R2data$avgR = abs(R2data$avgR)
R2data$Err = abs(R2data$Err)
R2data <- as.data.frame(R2data)


lme4 <- lmer(Err ~ avgR + (1| p1)  + (1 | p2)  + (1 | c1)  + (1 | c2), data=R2data, REML=F) #fit lme model
R2data$pred <- predict(lme4)
summary(lme4)
plot(R2data$Err,R2data$pred,pch=16)

#print(R2data)

fname=file.path("invdir","protocol-lmeOutRecip.dat")
write.table(nmeas,file=fname,col.names=FALSE,row.names=FALSE,append=FALSE)
write.table( cbind(R2data$i,R2data$p1,R2data$p2,R2data$c1,R2data$c2,R2data$avgR,R2data$pred) ,file=fname, 
             sep = "\t", col.names=FALSE,row.names=FALSE, append = TRUE)


## prediction
fname = file.path("invdir","protocol-lmeIn.dat")
fid <- file(fname,"r")
nmeas <- strtoi(readLines(fid,n=1))
R2data <- scan(fname, skip=1, list(i=0, p1=0, p2=0, c1=0, c2=0, avgR=0))
close(fid)
R2data$avgR = abs(R2data$avgR)
R2data <- as.data.frame(R2data)

R2data$pred <- abs(predict(lme4,newdata = R2data, allow.new.levels=T))
fname=file.path("invdir","protocol-lmeOut.dat")
write.table(nmeas,file=fname,col.names=FALSE,row.names=FALSE,append=FALSE)
write.table( cbind(R2data$i,R2data$p1,R2data$p2,R2data$c1,R2data$c2,R2data$avgR,abs(R2data$pred)) ,file=fname, 
             sep = "\t", col.names=FALSE,row.names=FALSE, append = TRUE)
print('end R')



