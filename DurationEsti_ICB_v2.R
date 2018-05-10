setwd("E:/OneDrive2nd/OneDrive - 广厚设计学校/GIT/Duration-Econ")
library(survival)
library(survminer)

library(survival)
library(ranger)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggfortify)
library(grid)
library(xtable)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## load time dependent covariates data
load("DurationTimeCov_ICB.RData")
tdData <- durB # change to durG for robustness check

tdData$joint_demo <- ifelse(tdData$demo1>5 & tdData$demo2>5, 1, 0)
tdData$traderatio <- tdData$tradeshare1/(tdData$tradeshare2+tdData$tradeshare1+1) ## not like power ratio, because the latter produce too many NA values
## need to think more about viol and salience
## because viol may affect termination
## Note 4/23/2018
## the choice of gravity matters a lot; think more about justification
tdData$issuesalience <- factor(ifelse(tdData$gravty %in% c(5,6),2,ifelse(tdData$gravty %in% c(2,3,4),1,0)), levels=c(0,1,2))
temp <- sapply(levels(tdData$issuesalience), function(x) as.integer(x == tdData$issuesalience))
colnames(temp) <- paste0("issuesalience",colnames(temp))
tdData <- cbind(tdData,temp)
## instead of factor; maybe transform into binary with three variables can be better
tdData$traderatio.net <- tdData$tradeshare1.net/(tdData$tradeshare2.net+tdData$tradeshare1.net)
tdData$defenseratio <- tdData$defense1/(1+tdData$defense2+tdData$defense1)
tdData$tradeissuesalience1 <- tdData$traderatio*tdData$issuesalience1
tdData$tradeissuesalience2 <- tdData$traderatio*tdData$issuesalience2

fit1b <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+tradeissuesalience1+tradeissuesalience2+joint_demo+contbinary+powerratio+defenseratio, data=tdData);summary(fit1b)
fit1bc <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+tradeissuesalience1+tradeissuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(fit1bc)
fit1bf<- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+tradeissuesalience1+tradeissuesalience2+joint_demo+contbinary+powerratio+defenseratio+frailty(crisno), data=tdData);summary(fit1bf)

test.ph <- cox.zph(fit1b,transform='rank')#given the proportion of censoring
## I need to use time transform rank or km suggested by Park2015
test.ph 


######################
## Main Model Clustered Model
#####################

####################
## tt without traderatio
####################
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(update.formula(fit1bc, ~. + tt(issuesalience1) +tt(tradeissuesalience1)
                              + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
               data=tdData1, tt=function(x, t, ...) {x*log(t+lam[i])})
  logL[i] <- logLik(fit)
  aic[i] <- AIC(fit)
  bic[i] <- BIC(fit)
}
#par(mar=c(5,6,1,1))
plot(lam, logL, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
plot(lam, aic, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
plot(lam, bic, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
mtext('Log likelihood', 2, line=4)
## now use the best lam to refit the model
lamhat <- lam[which.min(bic)]
best_fitc_main <- coxph(update.formula(fit1bc, ~. + tt(issuesalience1) +tt(tradeissuesalience1)
                                  + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
                   data=tdData1,tt=function(x, t, ...) {x*log(t+lamhat)})
summary(best_fitc_main)

##########
## First Difference of Trade Ratio
## tt without traderatio
###################
par_traderatio <- summary(best_fitc_main)$coefficients["traderatio",c(1,4)]
par_tradeissuesalience1 <- summary(best_fitc_main)$coefficients["tradeissuesalience1",c(1,4)]
par_ttradeissuesalience1 <- summary(best_fitc_main)$coefficients["tt(tradeissuesalience1)",c(1,4)]
par_tradeissuesalience2 <- summary(best_fitc_main)$coefficients["tradeissuesalience2",c(1,4)]
par_ttradeissuesalience2 <- summary(best_fitc_main)$coefficients["tt(tradeissuesalience2)",c(1,4)]
#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par01, par2, par3) (exp(0.01*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par01, par2, par3) (exp(0.01*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par01) (exp(0.01*(par01))-1)*100
t <- seq(0,999, length=1000) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = log(t[i]+lamhat)
  impact0=impact1=impact2=numeric()
  for (j in 1:sim) {
    par01 = rnorm(1, mean=par_traderatio[1], sd=par_traderatio[2])
    par12 = rnorm(1, mean=par_tradeissuesalience1[1], sd=par_tradeissuesalience1[2])
    par13 = rnorm(1, mean=par_ttradeissuesalience1[1], sd=par_ttradeissuesalience1[2])
    par22 = rnorm(1, mean=par_tradeissuesalience2[1], sd=par_tradeissuesalience2[2])
    par23 = rnorm(1, mean=par_ttradeissuesalience2[1], sd=par_ttradeissuesalience2[2])
    impact2[j] = hrtrade_issuesalience2(x, par01,par22, par23)
    impact1[j] = hrtrade_issuesalience1(x, par01, par12, par13)
    impact0[j] = hrtrade_issuesalience0(x, par01)
  }
  temp = data.frame (time=i, firstdiff1=impact1, firstdiff2=impact2, firstdiff0=impact0)
  plot_dat = rbind(plot_dat, temp)
}
plot_dat = reshape(plot_dat, varying=c("firstdiff1","firstdiff2", "firstdiff0"), v.names="impact", timevar="type",times=c("firstdiff1","firstdiff2","firstdiff0"), direction="long", new.row.names = NULL)
## Plot the impact of trade ratio
tradeimpactplot_c_main <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.001)+
  geom_smooth(colour="black",aes(linetype=type))+
  xlab("Days")+ylab(expression("First Difference"))+
  theme_bw() + theme(legend.position="bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_linetype_manual(values=c("twodash", "dotted","solid"),
                        name="Salience", 
                        labels = c("Low","Median","High"))


####################
## tt for traderatio
####################
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(update.formula(fit1bc, ~. +tt(traderatio)+ tt(issuesalience1) +tt(tradeissuesalience1)
                              + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
               data=tdData1, tt=function(x, t, ...) {x*log(t+lam[i])})
  logL[i] <- logLik(fit)
  aic[i] <- AIC(fit)
  bic[i] <- BIC(fit)
}
#par(mar=c(5,6,1,1))
plot(lam, logL, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
plot(lam, aic, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
plot(lam, bic, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
mtext('Log likelihood', 2, line=4)
## now use the best lam to refit the model
lamhat <- lam[which.min(bic)]
best_fitc <- coxph(update.formula(fit1bc, ~. +tt(traderatio)+ tt(issuesalience1) +tt(tradeissuesalience1)
                                  + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
                   data=tdData1,tt=function(x, t, ...) {x*log(t+lamhat)})
summary(best_fitc)

##########
## First Difference of Trade Ratio
## tt for traderatio
###################
par_traderatio <- summary(best_fitc)$coefficients["traderatio",c(1,4)]
par_ttraderatio <- summary(best_fitc)$coefficients["tt(traderatio)",c(1,4)]
par_tradeissuesalience1 <- summary(best_fitc)$coefficients["tradeissuesalience1",c(1,4)]
par_ttradeissuesalience1 <- summary(best_fitc)$coefficients["tt(tradeissuesalience1)",c(1,4)]
par_tradeissuesalience2 <- summary(best_fitc)$coefficients["tradeissuesalience2",c(1,4)]
par_ttradeissuesalience2 <- summary(best_fitc)$coefficients["tt(tradeissuesalience2)",c(1,4)]
#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par01, par2, par3) (exp(0.01*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par01, par2, par3) (exp(0.01*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par01,par02) (exp(0.01*(par01+par02*x))-1)*100
t <- seq(0,999, length=1000) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = log(t[i]+lamhat)
  impact0=impact1=impact2=numeric()
  for (j in 1:sim) {
    par01 = rnorm(1, mean=par_traderatio[1], sd=par_traderatio[2])
    par02 = rnorm(1, mean=par_ttraderatio[1], sd=par_ttraderatio[2])
    par12 = rnorm(1, mean=par_tradeissuesalience1[1], sd=par_tradeissuesalience1[2])
    par13 = rnorm(1, mean=par_ttradeissuesalience1[1], sd=par_ttradeissuesalience1[2])
    par22 = rnorm(1, mean=par_tradeissuesalience2[1], sd=par_tradeissuesalience2[2])
    par23 = rnorm(1, mean=par_ttradeissuesalience2[1], sd=par_ttradeissuesalience2[2])
    impact2[j] = hrtrade_issuesalience2(x, par01,par22, par23)
    impact1[j] = hrtrade_issuesalience1(x, par01, par12, par13)
    impact0[j] = hrtrade_issuesalience0(x, par01,par02)
  }
  temp = data.frame (time=i, firstdiff1=impact1, firstdiff2=impact2, firstdiff0=impact0)
  plot_dat = rbind(plot_dat, temp)
}
plot_dat = reshape(plot_dat, varying=c("firstdiff1","firstdiff2", "firstdiff0"), v.names="impact", timevar="type",times=c("firstdiff1","firstdiff2","firstdiff0"), direction="long", new.row.names = NULL)
## Plot the impact of trade ratio
tradeimpactplot_c <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.001)+
  geom_smooth(colour="black",aes(linetype=type))+
  xlab("Days")+ylab(expression("First Difference"))+
  theme_bw() + theme(legend.position="bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_linetype_manual(values=c("twodash", "dotted","solid"),
                        name="Salience", 
                        labels = c("Low","Median","High"))
#########################
## The main graph to show
#########################
tradeimpactplot_c+scale_y_continuous(breaks =c(-20,0,20),limits = c(-40,40))







#########################
## Frailty Model
#########################

## select the best fit
## using log likelihood
tdData1 <- tdData
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(update.formula(fit1bf, ~. +tt(traderatio)+ tt(issuesalience1) +tt(tradeissuesalience1)
                              + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
               data=tdData1, tt=function(x, t, ...) {x*log(t+lam[i])})
  logL[i] <- logLik(fit)
  aic[i] <- AIC(fit)
  bic[i] <- BIC(fit)
}
#par(mar=c(5,6,1,1))
plot(lam, logL, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
plot(lam, aic, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
plot(lam, bic, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
mtext('Log likelihood', 2, line=4)
## now use the best lam to refit the model
lamhat <- lam[which.min(bic)]
best_fitf <- coxph(update.formula(fit1bf, ~.+tt(traderatio) + tt(issuesalience1) +tt(tradeissuesalience1)
                                 + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
                  data=tdData1,tt=function(x, t, ...) {x*log(t+lamhat)})
summary(best_fitf)
cox.zph(best_fitf,transform='rank')

##########
## First Difference of Trade Ratio
## tt for traderatio
###################
par_traderatio <- summary(best_fitf)$coefficients["traderatio",c(1,3)]
par_ttraderatio <- summary(best_fitf)$coefficients["tt(traderatio)",c(1,3)]
par_tradeissuesalience1 <- summary(best_fitf)$coefficients["tradeissuesalience1",c(1,3)]
par_ttradeissuesalience1 <- summary(best_fitf)$coefficients["tt(tradeissuesalience1)",c(1,3)]
par_tradeissuesalience2 <- summary(best_fitf)$coefficients["tradeissuesalience2",c(1,3)]
par_ttradeissuesalience2 <- summary(best_fitf)$coefficients["tt(tradeissuesalience2)",c(1,3)]
#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par01, par2, par3) (exp(0.01*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par01, par2, par3) (exp(0.01*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par01,par02) (exp(0.01*(par01+par02*x))-1)*100
t <- seq(0,999, length=1000) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = log(t[i]+lamhat)
  impact0=impact1=impact2=numeric()
  for (j in 1:sim) {
    par01 = rnorm(1, mean=par_traderatio[1], sd=par_traderatio[2])
    par02 = rnorm(1, mean=par_ttraderatio[1], sd=par_ttraderatio[2])
    par12 = rnorm(1, mean=par_tradeissuesalience1[1], sd=par_tradeissuesalience1[2])
    par13 = rnorm(1, mean=par_ttradeissuesalience1[1], sd=par_ttradeissuesalience1[2])
    par22 = rnorm(1, mean=par_tradeissuesalience2[1], sd=par_tradeissuesalience2[2])
    par23 = rnorm(1, mean=par_ttradeissuesalience2[1], sd=par_ttradeissuesalience2[2])
    impact2[j] = hrtrade_issuesalience2(x, par01,par22, par23)
    impact1[j] = hrtrade_issuesalience1(x, par01, par12, par13)
    impact0[j] = hrtrade_issuesalience0(x, par01,par02)
  }
  temp = data.frame (time=i, firstdiff1=impact1, firstdiff2=impact2, firstdiff0=impact0)
  plot_dat = rbind(plot_dat, temp)
}
plot_dat = reshape(plot_dat, varying=c("firstdiff1","firstdiff2", "firstdiff0"), v.names="impact", timevar="type",times=c("firstdiff1","firstdiff2","firstdiff0"), direction="long", new.row.names = NULL)
## Plot the impact of trade ratio
tradeimpactplot_f <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.001)+
  geom_smooth(colour="black",aes(linetype=type))+
  xlab("Days")+ylab(expression("First Difference"))+
  theme_bw() + theme(legend.position="bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_linetype_manual(values=c("twodash", "dotted","solid"),
                        name="Salience", 
                        labels = c("Low","Median","High"))
tradeimpactplot_f+geom_hline(yintercept = 0)
tradeimpactplot_f+ylim(-15,15)+geom_hline(yintercept = 0)

## Get around way to deal with time interaction for plotting survival curve
## needs more work





################################
## Model without cluster
###############################
## select the best fit
## using log likelihood
tdData1 <- tdData
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(update.formula(fit1b, ~. +tt(traderatio)+ tt(issuesalience1) +tt(tradeissuesalience1)
                              + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
               data=tdData1, tt=function(x, t, ...) {x*log(t+lam[i])})
  logL[i] <- logLik(fit)
  aic[i] <- AIC(fit)
  bic[i] <- BIC(fit)
}
#par(mar=c(5,6,1,1))
plot(lam, logL, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
plot(lam, aic, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
plot(lam, bic, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
mtext('Log likelihood', 2, line=4)
## now use the best lam to refit the model
lamhat <- lam[which.min(bic)]
best_fit <- coxph(update.formula(fit1b, ~.+tt(traderatio) + tt(issuesalience1) +tt(tradeissuesalience1)
                                 + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
                  data=tdData1,tt=function(x, t, ...) {x*log(t+lamhat)})
summary(best_fit)
cox.zph(best_fit,transform='rank')

##################
## First Difference of Salience
###################
## Check 

## Median Salience
par_salience1 <- summary(best_fit)$coefficients["issuesalience1",c(1,3)]
par_tradeissuesalience1 <- summary(best_fit)$coefficients["tradeissuesalience1",c(1,3)]
par_ttradeissuesalience1 <- summary(best_fit)$coefficients["tt(tradeissuesalience1)",c(1,3)]
salienceimpact <- function(x, trade,par1,par2,par3) ((exp(par1+par2*trade+par3*trade*x)))

t <- seq(0,999, length=1000) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = log(t[i]+lamhat)
  impactl=impacth=numeric()
  for (j in 1:sim) {
    par1 = rnorm(1, mean=par_salience1[1], sd=par_salience1[2])
    par2 = rnorm(1, mean=par_tradeissuesalience1[1], sd=par_tradeissuesalience1[2])
    par3 = rnorm(1, mean=par_ttradeissuesalience1[1], sd=par_ttradeissuesalience1[2])
    impactl[j]=salienceimpact(x,quantile(tdData$traderatio,.25,na.rm=TRUE), par1, par2, par3)
    impacth[j]=salienceimpact(x,quantile(tdData$traderatio,.75,na.rm=TRUE), par1, par2, par3)
  }
  temp = data.frame (time=i, firstdiffl=impactl, firstdiffh=impacth)
  plot_dat = rbind(plot_dat, temp)
}
plot_dat = reshape(plot_dat, varying=c("firstdiffl","firstdiffh"), v.names="impact", timevar="type",times=c("firstdiffl","firstdiffh"), direction="long", new.row.names = NULL)
## Plot the impact of trade ratio
salience1impactplot <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.01)+
  geom_smooth(colour="black",aes(linetype=type))+ylim(0,1)+ xlim(0,50)+
  xlab("Days")+ylab(expression("Relative Hazard"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_linetype_discrete(name="Trade Ratio", 
                          labels = c("High", "Low"))

## High Salience
par_salience2 <- summary(best_fit)$coefficients["issuesalience2",c(1,3)]
par_tradeissuesalience2 <- summary(best_fit)$coefficients["tradeissuesalience2",c(1,3)]
par_ttradeissuesalience2 <- summary(best_fit)$coefficients["tt(tradeissuesalience2)",c(1,3)]
salienceimpact <- function(x, trade,par1,par2,par3) ((exp(par1+par2*trade+par3*trade*x)))

t <- seq(0,999, length=1000) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = log(t[i]+lamhat)
  impactl=impacth=numeric()
  for (j in 1:sim) {
    par1 = rnorm(1, mean=par_salience2[1], sd=par_salience2[2])
    par2 = rnorm(1, mean=par_tradeissuesalience2[1], sd=par_tradeissuesalience2[2])
    par3 = rnorm(1, mean=par_ttradeissuesalience2[1], sd=par_ttradeissuesalience2[2])
    impactl[j]=salienceimpact(x,quantile(tdData$traderatio,.25,na.rm=TRUE), par1, par2, par3)
    impacth[j]=salienceimpact(x,quantile(tdData$traderatio,.75,na.rm=TRUE), par1, par2, par3)
  }
  temp = data.frame (time=i, firstdiffl=impactl, firstdiffh=impacth)
  plot_dat = rbind(plot_dat, temp)
}
plot_dat = reshape(plot_dat, varying=c("firstdiffl","firstdiffh"), v.names="impact", timevar="type",times=c("firstdiffl","firstdiffh"), direction="long", new.row.names = NULL)
## Plot the impact of trade ratio
salience2impactplot <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.01)+
  geom_smooth(colour="black",aes(linetype=type))+ylim(0,1)+ xlim(0,50)+
  xlab("Days")+ylab(expression("Relative Hazard"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_linetype_discrete(name="Trade Ratio", 
                          labels = c("High", "Low"))


##########
## First Difference of Trade Ratio
###################
par_traderatio <- summary(best_fit)$coefficients["traderatio",c(1,3)]
par_tradeissuesalience1 <- summary(best_fit)$coefficients["tradeissuesalience1",c(1,3)]
par_ttradeissuesalience1 <- summary(best_fit)$coefficients["tt(tradeissuesalience1)",c(1,3)]
par_tradeissuesalience2 <- summary(best_fit)$coefficients["tradeissuesalience2",c(1,3)]
par_ttradeissuesalience2 <- summary(best_fit)$coefficients["tt(tradeissuesalience2)",c(1,3)]
#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par1, par2, par3) (exp(0.01*(par1+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par1, par2, par3) (exp(0.01*(par1+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par1) (exp(0.01*(par1))-1)*100
t <- seq(0,999, length=1000) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = log(t[i]+lamhat)
  impact0=impact1=impact2=numeric()
  for (j in 1:sim) {
    par1 = rnorm(1, mean=par_traderatio[1], sd=par_traderatio[2])
    par12 = rnorm(1, mean=par_tradeissuesalience1[1], sd=par_tradeissuesalience1[2])
    par13 = rnorm(1, mean=par_ttradeissuesalience1[1], sd=par_ttradeissuesalience1[2])
    par22 = rnorm(1, mean=par_tradeissuesalience2[1], sd=par_tradeissuesalience2[2])
    par23 = rnorm(1, mean=par_ttradeissuesalience2[1], sd=par_ttradeissuesalience2[2])
    impact2[j] = hrtrade_issuesalience2(x, par1, par22, par23)
    impact1[j] = hrtrade_issuesalience1(x, par1, par12, par13)
    impact0[j] = hrtrade_issuesalience0(x, par1)
  }
  temp = data.frame (time=i, firstdiff1=impact1, firstdiff2=impact2, firstdiff0=impact0)
  plot_dat = rbind(plot_dat, temp)
}
plot_dat = reshape(plot_dat, varying=c("firstdiff1","firstdiff2", "firstdiff0"), v.names="impact", timevar="type",times=c("firstdiff1","firstdiff2","firstdiff0"), direction="long", new.row.names = NULL)
## Plot the impact of trade ratio
tradeimpactplot <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.001)+
  geom_smooth(colour="black",aes(linetype=type))+ylim(-5,5) +
  xlab("Days")+ylab(expression("First Difference"))+
  theme_bw() + theme(legend.position="bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_linetype_manual(values=c("twodash", "dotted","solid"),
                        name="Salience", 
                        labels = c("Low","Median","High"))
  

##########
## First Difference of Trade Ratio
## tt for traderatio
## 5/8 sth still not right,need more work
###################
par_traderatio <- summary(best_fit)$coefficients["traderatio",c(1,3)]
par_ttraderatio <- summary(best_fit)$coefficients["tt(traderatio)",c(1,3)]
par_tradeissuesalience1 <- summary(best_fit)$coefficients["tradeissuesalience1",c(1,3)]
par_ttradeissuesalience1 <- summary(best_fit)$coefficients["tt(tradeissuesalience1)",c(1,3)]
par_tradeissuesalience2 <- summary(best_fit)$coefficients["tradeissuesalience2",c(1,3)]
par_ttradeissuesalience2 <- summary(best_fit)$coefficients["tt(tradeissuesalience2)",c(1,3)]
#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par01, par2, par3) (exp(0.01*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par01, par2, par3) (exp(0.01*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par01,par02) (exp(0.01*(par01+par02*x))-1)*100
t <- seq(0,999, length=1000) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = log(t[i]+lamhat)
  impact0=impact1=impact2=numeric()
  for (j in 1:sim) {
    par01 = rnorm(1, mean=par_traderatio[1], sd=par_traderatio[2])
    par02 = rnorm(1, mean=par_ttraderatio[1], sd=par_ttraderatio[2])
    par12 = rnorm(1, mean=par_tradeissuesalience1[1], sd=par_tradeissuesalience1[2])
    par13 = rnorm(1, mean=par_ttradeissuesalience1[1], sd=par_ttradeissuesalience1[2])
    par22 = rnorm(1, mean=par_tradeissuesalience2[1], sd=par_tradeissuesalience2[2])
    par23 = rnorm(1, mean=par_ttradeissuesalience2[1], sd=par_ttradeissuesalience2[2])
    impact2[j] = hrtrade_issuesalience2(x, par01,par22, par23)
    impact1[j] = hrtrade_issuesalience1(x, par01, par12, par13)
    impact0[j] = hrtrade_issuesalience0(x, par01,par02)
  }
  temp = data.frame (time=i, firstdiff1=impact1, firstdiff2=impact2, firstdiff0=impact0)
  plot_dat = rbind(plot_dat, temp)
}
plot_dat = reshape(plot_dat, varying=c("firstdiff1","firstdiff2", "firstdiff0"), v.names="impact", timevar="type",times=c("firstdiff1","firstdiff2","firstdiff0"), direction="long", new.row.names = NULL)
## Plot the impact of trade ratio
tradeimpactplot <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.001)+
  geom_smooth(colour="black",aes(linetype=type))+ylim(-5,5) +
  xlab("Days")+ylab(expression("First Difference"))+
  theme_bw() + theme(legend.position=c(.8,.8), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_linetype_manual(values=c("twodash", "dotted","solid"),
                        name="Salience", 
                        labels = c("Low","Median","High"))

######################
## Clustered Model
#####################

####################
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(update.formula(fit1bc, ~. +tt(traderatio)+ tt(issuesalience1) +tt(tradeissuesalience1)
                              + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
               data=tdData1, tt=function(x, t, ...) {x*log(t+lam[i])})
  logL[i] <- logLik(fit)
  aic[i] <- AIC(fit)
  bic[i] <- BIC(fit)
}
#par(mar=c(5,6,1,1))
plot(lam, logL, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
plot(lam, aic, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
plot(lam, bic, type='l', bty='n', las=1, lwd=2,
     xlab=expression(lambda), ylab='')
mtext('Log likelihood', 2, line=4)
## now use the best lam to refit the model
lamhat <- lam[which.min(bic)]
best_fitc <- coxph(update.formula(fit1bc, ~. +tt(traderatio)+ tt(issuesalience1) +tt(tradeissuesalience1)
                                 + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
                  data=tdData1,tt=function(x, t, ...) {x*log(t+lamhat)})
summary(best_fitc)

## robust se use column 4
par_traderatio <- summary(best_fitc)$coefficients["tt(traderatio)",c(1,4)]
par_tradeissuesalience1 <- summary(best_fitc)$coefficients["tradeissuesalience1",c(1,4)]
par_ttradeissuesalience1 <- summary(best_fitc)$coefficients["tt(tradeissuesalience1)",c(1,4)]
par_tradeissuesalience2 <- summary(best_fitc)$coefficients["tradeissuesalience2",c(1,4)]
par_ttradeissuesalience2 <- summary(best_fitc)$coefficients["tt(tradeissuesalience2)",c(1,4)]
#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par1, par2, par3) (exp(0.01*(par1*x+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par1, par2, par3) (exp(0.01*(par1*x+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par1) (exp(0.01*(par1*x))-1)*100
t <- seq(0,999, length=1000) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = t[i]
  impact0=impact1=impact2=numeric()
  for (j in 1:sim) {
    par1 = rnorm(1, mean=par_traderatio[1], sd=par_traderatio[2])
    par12 = rnorm(1, mean=par_tradeissuesalience1[1], sd=par_tradeissuesalience1[2])
    par13 = rnorm(1, mean=par_ttradeissuesalience1[1], sd=par_ttradeissuesalience1[2])
    par22 = rnorm(1, mean=par_tradeissuesalience2[1], sd=par_tradeissuesalience2[2])
    par23 = rnorm(1, mean=par_ttradeissuesalience2[1], sd=par_ttradeissuesalience2[2])
    impact2[j] = hrtrade_issuesalience2(x, par1, par22, par23)
    impact1[j] = hrtrade_issuesalience1(x, par1, par12, par13)
    impact0[j] = hrtrade_issuesalience0(x, par1)
  }
  temp = data.frame (time=x, firstdiff1=impact1, firstdiff2=impact2, firstdiff0=impact0)
  plot_dat = rbind(plot_dat, temp)
}
plot_dat = reshape(plot_dat, varying=c("firstdiff1","firstdiff2", "firstdiff0"), v.names="impact", timevar="type",times=c("firstdiff1","firstdiff2","firstdiff0"), direction="long", new.row.names = NULL)
## Plot the impact of trade ratio
tradeimpactplotc <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0)+
  geom_smooth(colour="black")+ylim(-100,250) +
  xlab("Days")+ylab(expression("First Difference"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", 750, 120, label="Low Salience")+ 
  annotate("text",500,0,label="Median Salience")+
  annotate("text",100,-75,label="High Salience")

plot_dat2 <- ddply(plot_dat, .(type,time), summarize, 
                   mean=mean(impact), sd=sd(impact),
                   lower=mean(impact)-sd(impact)*1.96,
                   upper=mean(impact)+sd(impact)*1.96)

ggplot(plot_dat2,aes(x=time, y=mean, group=type)) + geom_line()+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=time, fill = "band"), alpha = 0.3)



best_fitc_se <- summary(best_fitc)$coefficients[,4]
fit1bc_se <- summary(fit1bc)$coefficients[,4]

stargazer(best_fitc, best_fit, fit1bc, fit1b, ci=TRUE,digits=2,
          se=list(best_fitc_se,NULL,fit1bc_se,NULL),
          title="Cox Regression Results", dep.var.labels.include = FALSE,
          dep.var.caption  = "Days before Quitting",
          column.labels   = c("Adjusted Models", "Original Models"),
          column.separate = c(2, 2),
          covariate.labels=c("Cost Ratio","Median Salience", "High Salience",
                             "Cost*MSalience", "Cost*HSalience",
                             "Joint Democracy","Contiguity","Power Ratio",
                             "Defense Ratio", "tt(Cost Ratio)","tt(Median Salience)",
                             "tt(Cost*MSalience)", "tt(High Salience)",
                             "tt(Cost*HSalience)", "tt(Contiguity)"
                             ),
          add.lines = list(c("Clustered?","Yes","No","Yes","No"),
                           c("BIC", round(BIC(best_fitc),2),round(BIC(best_fit),2), round(BIC(fit1bc),2),round(BIC(fit1b),2))),
          keep.stat=c("wald","ll", "n"),
          no.space=TRUE)