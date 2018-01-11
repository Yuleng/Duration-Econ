## This is an Rcode for analyzing the duration econ data

load("C:/Users/YULENG/OneDrive/Documents/R Programming/DurationPaper/Duration.RData")

# incorporate hhi 
# load("C:/Users/YULENG/OneDrive/Documents/RESEARCH/HHI Paper/HHI_Data.RData")
#hhi <- dat[,c("ccode1", "ccode2", "year", "hhi_ex1", "hhi_ex2")] # get only the export hhi
#durB <- merge(durB,hhi,by=c("ccode1","ccode2", "year"),all.x = TRUE)
#durB$tradeshare1hhi <- durB$tradeshare1*durB$hhi_ex1
#durB$tradeshare2hhi <- durB$tradeshare2*durB$hhi_ex2
#durG <- merge(durG,hhi,by=c("ccode1","ccode2", "year"),all.x = TRUE)
#durG$tradeshare1hhi <- durG$tradeshare1*durG$hhi_ex1
#durG$tradeshare2hhi <- durG$tradeshare2*durG$hhi_ex2

library(survival)
## cox proportional hazard
mod1 <- coxph(Surv(maxdur,cens)~demo1+demo2+alliance+IGOTally+contiguity+mindist+powerratio+tradeshare1+tradeshare2,durB)
summary(mod1)

mod2 <- coxph(Surv(maxdur)~demo1+alliance+contiguity+powerratio+tradedepend1+tradedepend2,durB)

mod3 <- coxph(Surv(maxdur)~demo1+alliance+contiguity+powerratio+tradeshare1+tradeshare2,durG)

mod4 <- coxph(Surv(maxdur)~demo1+alliance+contiguity+powerratio+tradedepend1+tradedepend2,durG)

# Validating Cox PH Assumptions
validate_coxph = cox.zph(mod1, transform = "km")
validate_coxph

## Accelerated failure-time models
## note that while PH model models hazard
## AFT models survival time
## check http://myweb.uiowa.edu/pbreheny/7210/f15/notes/10-15.pdf
mod1 <- survreg(Surv(maxdur,cens)~demo1+alliance+IGOTally+contiguity+mindist+powerratio+tradeshare1+tradeshare2,durB, dist="weibull")
summary(mod1)


## Kaplan-Meier estimator looks like not doable
kmfit <- survfit(Surv(maxdur, cens) ~ tradeshare1+tradeshare2, 
                     data=durB, conf.type="log-log")
summary(kmfit)

## Nelson-Altschuler estimate looks like not doable
nafit <- survfit(Surv(maxdur, cens) ~ demo1+alliance+IGOTally+contiguity+mindist+powerratio+tradeshare1+tradeshare2, 
                 data=durB, conf.type="log-log", type="fh")
summary(nafit)


## load time dependent covariates
load("C:/Users/YULENG/OneDrive/Documents/R Programming/DurationPaper/DurationTimeCov_v1.RData")
tdData <- durB
fit1 <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+demo1+demo2+mindist+alliance+powerratio, data=tdData)
fit2 <- coxph(Surv(tstart, tstop, quit) ~ tradedepend1+tradedepend2+demo1+demo2+mindist+alliance+powerratio, data=tdData)
fit3 <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+demo1+demo2+contiguity+alliance+powerratio, data=tdData)
fit4 <- coxph(Surv(tstart, tstop, quit) ~ tradedepend1+tradedepend2+demo1+demo2+contiguity+alliance+powerratio, data=tdData)
fit5 <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+demo1+demo2+mindist+igo+powerratio, data=tdData)
fit6 <- coxph(Surv(tstart, tstop, quit) ~ tradedepend1+tradedepend2+demo1+demo2+mindist+igo+powerratio, data=tdData)

#####
## consider time dependent coefficients
fit1 <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+demo1+demo2+mindist+alliance+igo+powerratio, data=tdData);summary(fit1)
## check time dependent coefficient
cox.zph(fit1) ## indicate demo1 demo2 mindist should have time varying impacts
fit1t <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+demo1+demo2+tt(demo1)+tt(demo2)+mindist+tt(mindist)+alliance+igo+powerratio, data=tdData, tt=function(x, t, ...) {x*t});summary(fit1t)
## check again
cox.zph(fit1t) ## still some demo1 demo2 but gengerally no time varying impacts for all other variables

fit1 <- coxph(Surv(tstart, tstop, quit) ~ tradedepend1+tradedepend2+demo1+demo2+mindist+alliance+igo+powerratio, data=tdData);summary(fit1)
## check time dependent coefficient
cox.zph(fit1) ## indicate demo1 demo2 mindist should have time varying impacts
fit1t <- coxph(Surv(tstart, tstop, quit) ~ tradedepend1+tradedepend2+demo1+demo2+tt(demo1)+tt(demo2)+mindist+tt(mindist)+alliance+igo+powerratio, data=tdData, tt=function(x, t, ...) {x*t});summary(fit1t)
## check again
cox.zph(fit1t) ## still some demo1 demo2 but gengerally no time varying impacts for all other variables

## althought no indication for trade, give it a try
fit1t <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+tt(tradeshare1)+tt(tradeshare2)+demo1+demo2+tt(demo1)+tt(demo2)+mindist+tt(mindist)+alliance+igo+powerratio, data=tdData, tt=function(x, t, ...) {x*t});summary(fit1t)
## check again
cox.zph(fit1t) ## still some demo1 demo2 but gengerally no time varying impacts for all other variables


fit1t <- coxph(Surv(tstart, tstop, quit) ~ tradedepend1+tradedepend2+tt(tradedepend1)+tt(tradedepend2)+demo1+demo2+tt(demo1)+tt(demo2)+mindist+tt(mindist)+alliance+igo+powerratio, data=tdData, tt=function(x, t, ...) {x*t});summary(fit1t)
## check again
cox.zph(fit1t) ## still some demo1 demo2 but gengerally no time varying impacts for all other variables
## seems to turn trade into nonsignificant


####
## Now COMPETING RISKS
####################