## This is the Rcode for testing the duration paper
## Note 2018 Deadline March
## Rethink the data first

setwd("E:/OneDrive2nd/OneDrive - 广厚设计学校/GIT/Duration-Econ")
library(survival)
library(survminer)

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
## load time dependent covariates data
load("DurationTimeCov_v2.RData")
tdData <- durB # change to durG for robustness check

tdData$joint_demo <- ifelse(tdData$demo1>5 & tdData$demo2>5, 1, 0)
tdData$lossratio <- tdData$trooploss1/(tdData$trooploss2+1)
#tdData$defense <- tdData$defense1/(tdData$defense2+2) ##defense 2 has 23 -1 values, figure it out later
tdData$traderatio <- tdData$tradeshare1/(tdData$tradeshare2+1)
tdData$resolve <- ifelse(tdData$revtype11==1,1,0)

## km plot
km_fit <- survfit(Surv(tstart,tstop, quit) ~ 1, data=tdData)
ggsurvplot(km_fit)

## run preliminary model
# fit1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+resolve+joint_demo+lossratio+contbinary+powerratio, data=tdData);summary(fit1)
# fit1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+resolve+traderatio:resolve+joint_demo+lossratio+contbinary+powerratio, data=tdData);summary(fit1)
# fit1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+resolve+traderatio:resolve+joint_demo+lossratio+contbinary+powerratio+cluster(dispnum3), data=tdData);summary(fit1)

fit1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+resolve+traderatio:resolve+joint_demo+lossratio+contbinary+powerratio, data=tdData);summary(fit1)
###################
## Plot the impact of traderatio
temp <- as.data.frame(model.matrix(fit1))
## check the distribution of traderatio quantile(temp[,1])
newdt1 <- data.frame(traderatio=rep(c(0.1,1,2,5,10), each=dim(temp)[1]),
                rbind(temp[,-1],temp[,-1],temp[,-1],temp[,-1],temp[,-1]))
newdt1$`traderatio:resolve`<- newdt1$traderatio * newdt1$resolve
ggadjustedcurves(fit1, data = newdt1, fun = "pct", method="average",
                 legend.title="Impact of Trade Ratio", 
                 variable ="traderatio", ylab="Prob of Surv") # add xlim if needed xlim=c(0,quantile(tdData$dur, na.rm=TRUE, probs=.75)),

## Plot the impact of traderatio for resolved types
newdt2 <- newdt1; newdt2$resolve=1
newdt2$`traderatio:resolve`<- newdt2$traderatio * newdt2$resolve
ggadjustedcurves(fit1, data = newdt2, fun = "pct", method="average",
                 legend.title="Impact of Trade Ratio", 
                 variable ="traderatio", ylab="Prob of Surv") # add xlim if needed xlim=c(0,quantile(tdData$dur, na.rm=TRUE, probs=.75)),

## check the impact of resolve: conditions the impact of traderatio
newdt3 <- rbind(temp,temp);newdt3$resolve <- rep(c(0,1),each=dim(temp)[1])
newdt3$`traderatio:resolve`<- newdt3$traderatio * newdt3$resolve
ggadjustedcurves(fit1, data = newdt3, fun = "pct", method="average",
                 legend.title="Impact of Resolve", 
                 variable ="resolve", ylab="Prob of Surv") # add xlim if needed xlim=c(0,quantile(tdData$dur, na.rm=TRUE, probs=.75)),

## This part follows STHDA
## http://www.sthda.com/english/wiki/cox-model-assumptions
##
## testing proportional hazard
## using the default instead of identity or log
## see Terry T.'s response here
## http://r.789695.n4.nabble.com/Why-is-transform-quot-km-quot-the-default-for-cox-zph-td797535.html
test.ph <- cox.zph(fit1)
test.ph
ggcoxzph(test.ph)[1]
plot(test.ph[2])
abline(h=0,col=2)
abline(h=fit1$coef[2], col=3, lwd=2, lty=2)
## zoom in
plot(test.ph[2], ylim=c(-2,2))
abline(h=0,col=2)
abline(h=fit1$coef[2], col=3, lwd=2, lty=2)
## zoom in
plot(test.ph[1], ylim=c(-1,1))
abline(h=0,col=2)
abline(h=fit1$coef[1], col=3, lwd=2, lty=2)


## deal with ph violation
#####################################3
## follow Patrick Breheny
## http://myweb.uiowa.edu/pbreheny/7210/f17/notes.html
## to deal with ph violation, search for best fit, and plot
## use his function for plotting effect
effectPlot <- function(fit, t, fun, term, ...) {
  b <- coef(fit)
  y <- lwr <- upr <- numeric(length(t))
  for (i in 1:length(t)) {
    delta <- rep(0, length(b))
    names(delta) <- names(b)
    delta[term] <- 1
    delta[paste0('tt(', term, ')')] <- fun(t[i])
    y[i] <- delta %*% coef(fit)
    v <- delta %*% vcov(fit) %*% delta
    lwr[i] <- y[i] + qnorm(.025)*sqrt(v)
    upr[i] <- y[i] + qnorm(.975)*sqrt(v)
  }
  plot(t, y, ylim=range(c(lwr, upr)), type='l', las=1, bty='n', ...)
  polygon(c(t, rev(t)), c(lwr, rev(upr)), col='gray85', border=NA)
  lines(t, y, lwd=3)
  abline(h=0, col='gray', lwd=2, lty=2)
}

effectPlot(mod1, 0:1000, as.numeric, 'traderatio', xlab='Time (Days)', ylab='Treatment effect')
## still need to tweak it to plot the survival rate
## consider park and hendry 2015, licht2011

## select the best fit
## using log likelihood
logLik(mod1)
lam <- seq(0.1, 10, len=99)
l <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(Surv(tstart, tstop, quit) ~ traderatio + resolve + 
                 traderatio:resolve + joint_demo + lossratio + contbinary + 
                 powerratio + tt(resolve) + tt(traderesolve),
               data=tdData1, tt=function(x, t, ...) {x*exp(-t/lam[i])})
  l[i] <- logLik(fit)
}
par(mar=c(5,6,1,1))
plot(lam, l, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
mtext('Log likelihood', 2, line=4)
## now use the best lam to refit the model
lamhat <- lam[which.max(l)]
fit <- coxph(Surv(tstart, tstop, quit) ~ traderatio + resolve + 
               traderatio:resolve + joint_demo + lossratio + contbinary + 
               powerratio + tt(resolve) + tt(traderesolve),
             data=tdData1,tt=function(x, t, ...) {x*exp(-t/lamhat)})
effectPlot(fit, seq(0, 1.5, len=99), function(t) exp(-t/lamhat), 'resolve', xlab='Time (Days)', ylab='Treatment effect')




## People seem to disagree on time interaction
## here it is interacted with start time
## https://www.r-bloggers.com/dealing-with-non-proportional-hazards-in-r/
## fox interacted with stop time
## https://socserv.socsci.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Cox-Regression.pdf
## Therneau urged the use of tt function
## because interacting with time breaks the rule of
## not looking into the future
## in the document for survival package p20
tdData1 <- tdData; tdData1$traderatio <- tdData1$traderatio+0.0001 #add this to exclude 0 so that my plots to get around tt function can work
tdData1$traderesolve <- tdData1$traderatio*tdData1$resolve
mod1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+resolve+traderatio:resolve+joint_demo+lossratio+contbinary+powerratio+
              tt(resolve)+tt(traderesolve),
              data=tdData1,tt = function(x, t, ...) x * log(t+20));summary(mod1)

## cannot be plotted because survfit cannot process tt terms
## possible solutions
## https://stackoverflow.com/questions/31105216/plotting-estimated-hr-from-coxph-object-with-time-dependent-coefficient-and-spli/31316057#31316057
## https://cran.r-project.org/web/packages/Greg/vignettes/timeSplitter.html
## my thought: extract the model.matrx then refit the model

## Plot the impact of traderatio
temp <- as.data.frame(model.matrix(mod1))
weight <- c(0.1/temp$traderatio, 10/temp$traderatio)
newdt2 <- rbind(temp,temp)
newdt2$traderatio <- newdt2$traderatio*weight
newdt2$`traderatio:resolve` <- newdt2$`traderatio:resolve`*weight
newdt2$`tt(traderesolve)` <- newdt2$`tt(traderesolve)`*weight
names(newdt2)[6:8] <- c("tresolve","ttraderesolve","traderesolve")
newdt2 <- cbind(tdData1[ceiling(as.numeric(rownames(temp))),c("tstart","tstop","quit")],newdt2)
rmod1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+resolve+traderesolve+joint_demo+lossratio+contbinary+
                 tresolve+ttraderesolve,
               data=newdt2);summary(rmod1)
ggadjustedcurves(mod1, data = newdt2, fun = "pct", method="average",
                 legend.title="Impact of Trade Ratio", 
                 variable ="traderatio", ylab="Prob of Surv")

## try these suggested methods
## https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
cox_fit1 <- survfit(fit1)
autoplot(cox_fit1)

aa_fit <-aareg(Surv(tstart, tstop, quit) ~ traderatio+resolve+joint_demo+lossratio+contbinary+powerratio, data=tdData);summary(fit1)
autoplot(aa_fit)


## run preliminary test
fit1 <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+demo1+demo2+trooploss1+trooploss2+mindist+defense1+defense2+powerratio, data=tdData);summary(fit1)
fit1 <- coxph(Surv(tstart, tstop, quit) ~ tradeshare1+tradeshare2+joint_demo+lossratio+mindist+defense+powerratio, data=tdData);summary(fit1)
## mindist can be rerun by conttype, contbinary
## alliance can be rerun by igo, affinity
## powerratio can be rerun by major1+major2

## check time dependent coefficient
cox.zph(fit1) ## indicate demo1 demo2 mindist should have time varying impacts

## plot estimate survival time
ggsurvplot(survfit(fit1), fun = "pct", ggtheme = theme_minimal())


###################
## Plot the impact of traderatio
temp <- as.data.frame(model.matrix(fit1))
## check the distribution of traderatio quantile(temp[,1])
newdt1 <- cbind(rep(c(0.1,1,2,5,10), each=dim(temp)[1]),
                rbind(temp[,-1],temp[,-1],temp[,-1],temp[,-1],temp[,-1]))
names(newdt1) <- names(temp)
ggadjustedcurves(fit1, data = newdt1, fun = "pct", method="average",
                 legend.title="Impact of Trade Ratio", 
                variable ="traderatio", ylab="Prob of Surv") # add xlim if needed xlim=c(0,quantile(tdData$dur, na.rm=TRUE, probs=.75)),

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
