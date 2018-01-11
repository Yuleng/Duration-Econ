## This is the Rcode for testing the duration paper
## Note 2018 Deadline March

setwd("E:/OneDrive2nd/OneDrive - 广厚设计学校/GIT/Duration-Econ")
library(survival)
library(survminer)
## load time dependent covariates data
load("DurationTimeCov_v2.RData")
tdData <- durB # change to durG for robustness check

tdData$joint_demo <- ifelse(tdData$demo1>5 & tdData$demo2>5, 1, 0)
## run preliminary test
fit1 <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+demo1+demo2+mindist+alliance+powerratio, data=tdData);summary(fit1)
## mindist can be rerun by conttype, contbinary
## alliance can be rerun by igo, affinity
## powerratio can be rerun by major1+major2

## check time dependent coefficient
cox.zph(fit1) ## indicate demo1 demo2 mindist should have time varying impacts

## plot estimate survival time
ggsurvplot(survfit(fit1), fun = "pct", ggtheme = theme_minimal())


###################
## Plot the impact of tradeshare1
temp <- as.data.frame(model.matrix(fit1))
newdt1 <- cbind(rep(quantile(temp[,1],na.rm=TRUE,probs=c(.001,.99)), each=dim(temp)[1]),
                rbind(temp[,-1],temp[,-1]))
names(newdt1) <- names(temp)
ggcoxadjustedcurves(fit1, data = newdt1, fun = "cumhaz", xlim=c(0,quantile(tdData$dur, na.rm=TRUE, probs=.75)),
                    legend.title="trade", variable = newdt1[,1], ylab="Prob of Quiting")

## a different version
temp <- as.data.frame(model.matrix(fit1))
newdt1 <- cbind(rep(c(0.001,0.1,0.3), each=dim(temp)[1]),
                rbind(temp[,-1],temp[,-1],temp[,-1]))
names(newdt1) <- names(temp)
ggcoxadjustedcurves(fit1, data = newdt1, fun = "pct", xlim=c(0,quantile(tdData$dur, na.rm=TRUE, probs=.75)),
                    legend.title="trade", legend.labs=c("L","M","H"), variable = newdt1[,1])

## add time dep effect by demos
fit1t <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+demo1+demo2+tt(demo1)+tt(demo2)+mindist+alliance+tt(alliance)+powerratio, data=tdData, tt=function(x, t, ...) {x*t});summary(fit1t)
## check again
cox.zph(fit1t)




## use stratification instead of interacting with time given I have no interest in demo
fit1t <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+strata(demo1)+strata(demo2)+mindist+alliance+powerratio, data=tdData);summary(fit1t)
## check again
cox.zph(fit1t)

## using the suggestion by
## https://www.r-bloggers.com/dealing-with-non-proportional-hazards-in-r/
tdData$demo1_time = tdData$demo1*tdData$tstart; tdData$demo2_time = tdData$demo2*tdData$tstart
tdData$alliance_time = tdData$alliance*tdData$tstart;
fit2t <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+demo1+demo1_time+demo2+demo2_time+mindist+alliance+alliance_time+powerratio, data=tdData);summary(fit2t)
## check again
cox.zph(fit2t);dim(model.matrix(fit2t))


temp <- as.data.frame(model.matrix(fit2t))
newdt1 <- cbind(rep(c(0.001,0.1,0.3), each=dim(temp)[1]),
                rbind(temp[,-1],temp[,-1],temp[,-1]))
names(newdt1) <- names(temp)
ggcoxadjustedcurves(fit2t, data = newdt1, fun = "pct", xlim=c(1,quantile(tdData$dur, na.rm=TRUE, probs=.75)),
                    legend.title="trade", variable = newdt1[,1])


ggforest(fit2t)
ggcoxzph(cox.zph(fit2t))
ggcoxdiagnostics(fit2t)


## testing influential outliers
outliertest <- ggcoxdiagnostics(fit1, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())

## reestimate ofter deleting outliers deviance greater than 1
tdDataNoOutlier <- na.omit(tdData[,c(names(tdData)[1:6],colnames(model.matrix(fit1)))])[outliertest$data$res<1,]

fit2 <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+demo1+demo2+mindist+alliance+powerratio, data=tdDataNoOutlier);summary(fit2)
## check again
cox.zph(fit2)

## testing nonlinearty
ggcoxfunctional(Surv(tstart, tstop, quit) ~ tradeshare1 + log(tradeshare1) + sqrt(tradeshare1), data = tdData)

###################
## Plot the impact of tradeshare1
newdt1 <- as.data.frame(model.matrix(fit1t))
tradeshare1 <- rep(c(0.00001,0.5), each=dim(newdt1)[1])
newdt1 <- cbind(tradeshare1, rbind(subset(newdt1,select=-tradeshare1),subset(newdt1,select=-tradeshare1)))
ggcoxadjustedcurves(fit1t, data = newdt1, fun = "pct", variable = newdt1[,"tradeshare1"])




## Add splines for tradeshare 1 because the model indicates its impacts may vary
fit1s <- coxph(Surv(tstart, tstop, quit) ~ pspline(tradeshare1)+pspline(tradeshare2)+demo1+demo2+contbinary+alliance+cinc1+cinc2,data=tdData)
summary(fit1s)
termplot(fit1s,term=1,se=TRUE, ylim=c(-2.5,0.5)) # se=TRUE
termplot(fit1s,term=1,se=TRUE, ylim=c(-0.5,1.5))

fit2s <- coxph(Surv(tstart, tstop, quit) ~ pspline(tradeshare1)+tradeshare2+demo1+demo2+contbinary+alliance+cinc1+cinc2,data=na.omit(tdData))
summary(fit2s)
termplot(fit2s,term=1,se=TRUE, ylim=c(-10,10))


ptemp <- termplot(fit1s, term=1, se=TRUE, plot=F)
ptemp1 <- ptemp$tradeshare1
ptemp1 <- ptemp1[ptemp1$x<0.2,]
h <- ggplot(ptemp1, aes(x=x, y=y))
h+geom_ribbon(aes(ymin=y-1.96*se,ymax=y+1.96*se),fill="grey70")+geom_line(aes(y=y))
