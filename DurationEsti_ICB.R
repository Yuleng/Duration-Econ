## This is the Rcode for testing the duration paper
## Using ICB data
## Note that the data on Github has not been updated automatically
## since I only save data on my local files

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
tdData$traderatio <- tdData$tradeshare1/(tdData$tradeshare2+tdData$tradeshare1) ## not like power ratio, because the latter produce too many NA values
## need to think more about viol and salience
## because viol may affect termination

## Note 4/23/2018
## the choice of gravity matters a lot; think more about justification
tdData$issuesalience <- factor(ifelse(tdData$gravty %in% c(5,6),2,ifelse(tdData$gravty %in% c(2,3,4),1,0)), levels=c(0,1,2))
temp <- sapply(levels(tdData$issuesalience), function(x) as.integer(x == tdData$issuesalience))
colnames(temp) <- paste0("issuesalience",colnames(temp))
tdData <- cbind(tdData,temp)
## instead of factor; maybe transform into binary with three variables can be better
tdData$traderatio.net <- tdData$tradeshare1.net/(tdData$tradeshare2.net+1)
tdData$defenseratio <- tdData$defense1/(tdData$defense1+tdData$defense2)

## km plot
km_fit <- survfit(Surv(tstart,tstop, quit) ~ 1, data=tdData)
ggsurvplot(km_fit, palette="grey",legend="none")

## run preliminary model
# fit1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience+joint_demo+lossratio+contbinary+powerratio, data=tdData);summary(fit1)
# fit1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience+traderatio:issuesalience+joint_demo+lossratio+contbinary+powerratio, data=tdData);summary(fit1)
# fit1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience+traderatio:issuesalience+joint_demo+lossratio+contbinary+powerratio+cluster(dispnum3), data=tdData);summary(fit1)

#########################################
## First, fit the model with coxph, jump to the second step as the plots are misguiding and meaningless
## if ph is violated
#########################################
fit1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience+traderatio:issuesalience+joint_demo+contbinary+powerratio+defense1+defense2, data=tdData);summary(fit1)
fit1b <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio, data=tdData);summary(fit1b)
###################
## Plot the impact of traderatio
temp <- as.data.frame(model.matrix(fit1))
## check the distribution of traderatio quantile(temp[,1])
newdt1 <- data.frame(traderatio=rep(c(0.1,1,2,5,10), each=dim(temp)[1]),
                     rbind(temp[,-1],temp[,-1],temp[,-1],temp[,-1],temp[,-1]))
newdt1$`traderatio:issuesalience`<- newdt1$traderatio * newdt1$issuesalience
ggadjustedcurves(fit1, data = newdt1, fun = "pct", method="average",
                 legend.title="Impact of Trade Ratio", 
                 variable ="traderatio", ylab="Prob of Surv") # add xlim if needed xlim=c(0,quantile(tdData$dur, na.rm=TRUE, probs=.75)),

## Plot the impact of traderatio for issuesalienced types
newdt2 <- newdt1; newdt2$issuesalience=1
newdt2$`traderatio:issuesalience`<- newdt2$traderatio * newdt2$issuesalience
ggadjustedcurves(fit1, data = newdt2, fun = "pct", method="average",
                 legend.title="Impact of Trade Ratio", 
                 variable ="traderatio", ylab="Prob of Surv") # add xlim if needed xlim=c(0,quantile(tdData$dur, na.rm=TRUE, probs=.75)),

## check the impact of issuesalience: conditions the impact of traderatio
newdt3 <- rbind(temp,temp);newdt3$issuesalience <- rep(c(0,1),each=dim(temp)[1])
newdt3$`traderatio:issuesalience`<- newdt3$traderatio * newdt3$issuesalience
ggadjustedcurves(fit1, data = newdt3, fun = "pct", method="average",
                 legend.title="Impact of issuesalience", 
                 variable ="issuesalience", ylab="Prob of Surv") # add xlim if needed xlim=c(0,quantile(tdData$dur, na.rm=TRUE, probs=.75)),

## This part follows STHDA
## http://www.sthda.com/english/wiki/cox-model-assumptions

#################################
## Second, testing proportional hazard
#################################
## using the default instead of identity or log
## see Terry T.'s response here
## http://r.789695.n4.nabble.com/Why-is-transform-quot-km-quot-the-default-for-cox-zph-td797535.html

test.ph <- cox.zph(fit1,transform='rank')#given the proportion of censoring
## I need to use time transform rank or km suggested by Park2015
test.ph ## note, also try km to see if results change much
## graph the level of censoring
# naive way because missing data are ignored
cen_plot1 <- ggplot(as.data.frame(table(tdData$cens)),aes(x=Var1,y=Freq))+
  geom_bar(stat="identity")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("(a) Original Data") 
# proper way
cen_plot2 <- ggplot(as.data.frame(table(tdData[rownames(as.data.frame(model.matrix(fit1))), "cens"])),aes(x=Var1,y=Freq))+
  geom_bar(stat="identity")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("(b) Data in the Model") 
# plot the two together
multiplot(cen_plot1,cen_plot2, cols=2)

## indentity graph shows the impact of outliers
ggcoxzph(cox.zph(fit1,transform='identity'), point.col="grey")[2] #Showcase the need of time transformation Park2015
## graph using rank transformation
ggcoxzph(test.ph, point.col="grey")[2]
## zoom in shows the relation clearer
p1 <- ggcoxzph(test.ph, ylim=c(-4,4), xlab=FALSE, ylab=FALSE, point.col="grey")[2]
p2 <- ggcoxzph(test.ph, ylim=c(-4,4), xlab=FALSE, ylab=FALSE, point.col="grey")[3]
multiplot(p1, p2, cols=2)
## an alternative way to plot
## plot(test.ph[2], ylim=c(-2,2))
## abline(h=0,col=2)
## abline(h=fit1$coef[2], col=3, lwd=2, lty=2)

################################
## Now select the best fit
###############################
## select the best fit
## using log likelihood
tdData1 <- tdData
tdData1$tradeissuesalience1 <- tdData1$traderatio*tdData1$issuesalience1 #create the interaction term for tt function
tdData1$tradeissuesalience2 <- tdData1$traderatio*tdData1$issuesalience2 #create the interaction term for tt function
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(update.formula(fit1b, ~. + tt(issuesalience1) +tt(tradeissuesalience1)
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
best_fit <- coxph(update.formula(fit1b, ~. + tt(issuesalience1) +tt(tradeissuesalience1)
                 + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
                  data=tdData1,tt=function(x, t, ...) {x*log(t+lamhat)})
cox.zph(best_fit,transform='rank')
stargazer(best_fit, fit1b, 
          title="Regression Results", align=TRUE, dep.var.labels=c("Corrected Cox Model","Original Cox Model"), 
          covariate.labels=c("Cost Ratio","Median Salience", "High Salience",
                             "Joint Democracy","Contiguity","Power Ratio",
                             "Defense Ratio", "tt(Median Salience)",
                             "tt(Cost Ratio*Median Salience)", "tt(High Salience)",
                             "tt(Cost Ratio*High Salience)", "tt(Contiguity)",
                             "Cost Ratio*Median Salience", "Cost Ratio*High Salience"), 
          keep.stat=c("ll", "rsq", "n"), no.space=TRUE,
          add.lines=list(c("BIC", round(BIC(best_fit),1), round(BIC(fit1b),1))))


## Now Plot the first difference (percentage change in hazard rate) as
## suggested by Licht_2011 and by Box-Steffensmeier and Jones (2004, 60)
## A 1OO percent increase of trade ratio can be unrealistic
## hence, use 1 percent
par_traderatio <- summary(best_fit)$coefficients["traderatio",c(1,3)]
par_tradeissuesalience1 <- summary(best_fit)$coefficients["traderatio:issuesalience1",c(1,3)]
par_ttradeissuesalience1 <- summary(best_fit)$coefficients["tt(tradeissuesalience1)",c(1,3)]
par_tradeissuesalience2 <- summary(best_fit)$coefficients["traderatio:issuesalience2",c(1,3)]
par_ttradeissuesalience2 <- summary(best_fit)$coefficients["tt(tradeissuesalience2)",c(1,3)]
#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par1, par2, par3) (exp(0.01*(par1+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par1, par2, par3) (exp(0.01*(par1+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par1, par2, par3) (exp(0.01*(par1+par2*0+par3*x*0))-1)*100
t <- seq(0,999, length=1000) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = t[i]
  impact1=impac2=numeric()
  for (j in 1:sim) {
    par1 = rnorm(1, mean=par_traderatio[1], sd=par_traderatio[2])
    par12 = rnorm(1, mean=par_tradeissuesalience1[1], sd=par_tradeissuesalience1[2])
    par13 = rnorm(1, mean=par_ttradeissuesalience1[1], sd=par_ttradeissuesalience1[2])
    par22 = rnorm(1, mean=par_tradeissuesalience2[1], sd=par_tradeissuesalience2[2])
    par23 = rnorm(1, mean=par_ttradeissuesalience2[1], sd=par_ttradeissuesalience2[2])
    impact2= hrtrade_issuesalience2(x, par1, par22, par23)
    impact1 = hrtrade_issuesalience1(x, par1, par12, par13)
  }
  temp = data.frame (time=x, firstdiff1=impact1, firstdiff2=impact2)
  plot_dat = rbind(plot_dat, temp)
}
plot_dat = reshape(plot_dat, varying=c("firstdiff1","firstdiff2"), v.names="impact", timevar="type",times=c("firstdiff1","firstdiff2"), direction="long", new.row.names = NULL)
## Plot the impact of trade ratio
tradeimpactplot <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.05)+
  geom_smooth(colour="black")+ylim(-100,10) +
  xlab("Days")+ylab(expression("Hazard Rate"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

## Zoom in plot and incorporate the impact of issuesalience0 type 
zoominplot <- tradeimpactplot+
  geom_hline(linetype=2,yintercept=hrtrade_issuesalience0(0,par_traderatio[1],par_tradeissuesalience2[1],par_ttradeissuesalience2[1]))+
  annotate("text", 500, hrtrade_issuesalience0(0,par_traderatio[1],par_tradeissuesalience2[1],par_ttradeissuesalience2[1]), vjust = -1, label = paste("Low Salience:",round(hrtrade_issuesalience0(0,par_traderatio[1],par_tradeissuesalience2[1],par_ttradeissuesalience2[1]),digits=2)))+
  annotate("text",400,-40,label="Median Salience")+
  annotate("text",100,-75,label="High Salience")

## come back to Licht2011 when writing about interpretation


#########################
## Finally, if Hazard Rate change sign
## According to Ruhe2018 should plot survival rate
#########################
## tt function cannot handle prediction now
## second best choice interaction with time
## now there seems to be different takes on interacting with start or stop time
## John Fox interact it with stop time
## which appears to be close to the best fit model
tdData2 <- tdData1
#tdData2$traderatio <- tdData2$traderatio+0.0001 #add this to exclude 0 so that my plots to get around tt function can work
tdData2$time <- tdData2$tstop
tdData2$tissuesalience <- tdData2$issuesalience*log(tdData1$tstop)
tdData2$ttradeissuesalience <- tdData2$tradeissuesalience*log(tdData1$tstop)
secondbest_fit <- coxph(Surv(tstart, tstop, quit) ~ traderatio + issuesalience + 
                          tradeissuesalience + joint_demo + lossratio + contbinary + 
                          powerratio + tissuesalience+ttradeissuesalience,
                        data=tdData2); summary(secondbest_fit)
## check the difference
## summary(best_fit)$coefficients[,1];summary(secondbest_fit)$coefficients[,1]
## now used survminer's plot
## Plot the impact of traderatio

dat_fun <- function (tstart, tstop, traderatio, issuesalience=1, joint_demo=0, lossratio=1, contbinary=1, powerratio=1){
  num = length(traderatio)
  len = length(tstart)
  temp = data.frame(tstart=rep(tstart,num), tstop=rep(tstop,num), quit=rep(0,len*num),
                    curve=rep(1:num,each=len),
                    traderatio=rep(traderatio, each=len),
                    issuesalience=issuesalience, tradeissuesalience=rep(traderatio,each=len)*issuesalience,
                    joint_demo=joint_demo, lossratio=lossratio, contbinary=contbinary, powerratio=powerratio,
                    tissuesalience=log(rep(tstop,2))*issuesalience,
                    ttradeissuesalience=as.vector(outer(log(tstop),traderatio))*issuesalience)
  return(temp)
}
end <- 200; step <- 1
t1 <- seq(from=0,to=end-step,by=step)
t2 <- seq(from=step, to=end, by=step)
traderatio <- c(0.5,1)
temp <- dat_fun(t1,t2,traderatio,issuesalience=1,lossratio=0.1)


ggadjustedcurves(secondbest_fit,data=temp,variable="traderatio",method="average")

###############################################################
## The parts below can be deleted
#########################################################
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
tdData1$tradeissuesalience <- tdData1$traderatio*tdData1$issuesalience
mod1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience+traderatio:issuesalience+joint_demo+lossratio+contbinary+powerratio+
                tt(issuesalience)+tt(tradeissuesalience),
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
newdt2$`traderatio:issuesalience` <- newdt2$`traderatio:issuesalience`*weight
newdt2$`tt(tradeissuesalience)` <- newdt2$`tt(tradeissuesalience)`*weight
names(newdt2)[6:8] <- c("tissuesalience","ttradeissuesalience","tradeissuesalience")
newdt2 <- cbind(tdData1[ceiling(as.numeric(rownames(temp))),c("tstart","tstop","quit")],newdt2)
rmod1 <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience+tradeissuesalience+joint_demo+lossratio+contbinary+
                 tissuesalience+ttradeissuesalience,
               data=newdt2);summary(rmod1)
ggadjustedcurves(mod1, data = newdt2, fun = "pct", method="average",
                 legend.title="Impact of Trade Ratio", 
                 variable ="traderatio", ylab="Prob of Surv")

## try these suggested methods
## https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
cox_fit1 <- survfit(fit1)
autoplot(cox_fit1)

aa_fit <-aareg(Surv(tstart, tstop, quit) ~ traderatio+issuesalience+joint_demo+lossratio+contbinary+powerratio, data=tdData);summary(fit1)
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
