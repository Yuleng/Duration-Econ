##############
## This is R code for my attition paper: Using MID
#############
## 0. Preliminary stuff (data preparation, tests)
## 1. Main model (clustered and timeinteraction): First Difference Plot; Proxy survival Plot
## 2. tt function Model
## 4. Latex Table

##############
## libraries and functions
#############

setwd("E:/OneDrive2nd/OneDrive - 广厚设计学校/GIT/Duration-Econ")
library(MASS)
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
library(gridExtra)
library(coxme)
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
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}


############
## 0. Data
############
## load time dependent covariates data
load("DurationTimeCov_v2.RData")
tdData <- durB # change to durG for robustness check

#use network trade
tdData$joint_demo <- ifelse(tdData$demo1>5 & tdData$demo2>5, 1, 0)
## contest success function provides probability of winning
## cost ratio is closer to what I seek to capture
tdData$traderatio <- log(tdData$tradeshare1/(tdData$tradeshare2+.Machine$double.eps) +1)
tdData$issuesalience <- as.integer(ifelse(tdData$revtype11==1,1,0))
tdData$traderatio.net <- log(tdData$tradeshare1.net/(tdData$tradeshare2.net+.Machine$double.eps)+1)
tdData$defenseratio <- log(tdData$defense1/(tdData$defense2+.Machine$double.eps)+1)
tdData$lossratio <- log(tdData$trooploss1/(tdData$trooploss2+.Machine$double.eps)+1)
## alternative measurement
## tdData$traderatio <- tdData$tradedepend1/(tdData$tradedepend2+.Machine$double.eps)
## tdData$traderatio.net <- tdData$tradedepend1.net/(tdData$tradedepend2.net+.Machine$double.eps)

## change to using net measurement
tdData$traderatio <- tdData$traderatio.net ## comment out if not use

## construct the interaction term after deciding on trade net or not
tdData$tradeissuesalience <- tdData$traderatio*tdData$issuesalience

## file to save outputs; depend on data choice
file <- paste0("C:/Users/YULENG/OneDrive/Documents/RESEARCH/Manuscript/Duration/",
               "MID_Tradeshare_Net/")#tweak to different files

##############
## Basic Models and Tests

fit1b <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience+traderatio:issuesalience+joint_demo+contbinary+powerratio+defenseratio+lossratio, data=tdData);summary(fit1b)
fit1bc <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience+traderatio:issuesalience+joint_demo+contbinary+powerratio+defenseratio+lossratio+cluster(dispnum3), data=tdData);summary(fit1bc)
cox.zph(fit1b,transform='rank');cox.zph(fit1bc,transform='rank')
cox.zph(fit1b,transform='km');cox.zph(fit1bc,transform='km')

## graph the level of censoring
# naive way because missing data are ignored
cen_plot1 <- ggplot(as.data.frame(table(tdData$cens)),aes(x=Var1,y=Freq))+
  geom_bar(stat="identity")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("(a) Original Data") 
# proper way
cen_plot2 <- ggplot(as.data.frame(table(tdData[rownames(as.data.frame(model.matrix(fit1b))), "cens"])),aes(x=Var1,y=Freq))+
  geom_bar(stat="identity")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("(b) Data in the Model") 
# plot the two together
CensorPlot <- grid.arrange(cen_plot1,cen_plot2, nrow = 1)
ggsave(filename=paste0(file,"Censor.pdf"),
       plot=CensorPlot)

## outliers
## indentity graph shows the impact of outliers
OutlierPlot <- ggcoxzph(cox.zph(fit1bc,transform='identity'), point.col="grey")[2] 
#ggcoxzph(cox.zph(fit1bc,transform='identity'), point.col="grey")[3] 
ggsave(filename=paste0(file,"Outlier.pdf"))

###########
## COX.zph test indicate trade salience condiguity
phtest_rank <- cox.zph(fit1bc,transform='rank')$table
phtest_km <- cox.zph(fit1bc,transform='km')$table
print(xtable(phtest_rank, type = "latex",digits=4), file = paste0(file,"PHTestRank.tex"))
print(xtable(phtest_km, type = "latex",digits=4), file = paste0(file,"PHTestKM.tex"))

######################
## 1. Main Model Clustered Model
#####################

## conventional way using time interaction
## prep
tdData$time <- tdData$tstop ## using fox
cmod <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience+traderatio:issuesalience+joint_demo+contbinary+powerratio+defenseratio+lossratio+cluster(dispnum3), data=tdData);summary(cmod)
cox.zph(cmod,transform='rank')

## new interaciton formula
## Using only variables that violate the assumption most
cmodform <- update.formula(cmod, ~.+tcontbinary+tissuesalience+traderatio:tissuesalience)
## select the best model
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  #tdData$ttraderatio <- tdData$traderatio*log(tdData$time+lam)
  lam0 <- lam[i]
  tdData$tcontbinary <- tdData$contbinary*log(tdData$time+lam0)
  tdData$tpowerratio <- tdData$powerratio*log(tdData$time+lam0)
  tdData$tissuesalience <- tdData$issuesalience*log(tdData$time+lam0)
  fit <- coxph(cmodform, data=tdData)
  logL[i] <- logLik(fit)
  aic[i] <- AIC(fit)
  bic[i] <- BIC(fit)
}
## By BIC; others give similar results
lamhat <- lam[which.min(bic)]
#tdData$ttraderatio <- tdData$traderatio*log(tdData$time+lamhat)
tdData$tcontbinary <- tdData$contbinary*log(tdData$time+lamhat)
tdData$tpowerratio <- tdData$powerratio*log(tdData$time+lamhat)
tdData$tissuesalience <- tdData$issuesalience*log(tdData$time+lamhat)
cmod2 <- coxph(cmodform, data=tdData);summary(cmod2)
cox.zph(cmod2,transform='rank')

##########
## First Difference of Trade Ratio

##########
## check the distribution of traderatio
library(car)
densityPlot(tdData$traderatio,na.rm=T);
densityPlot(tdData$traderatio[tdData$traderatio<quantile(tdData$traderatio,.75,na.rm=T)],na.rm=T)

## proportional of effect
## prop <- quantile(tdData$traderatio,.75,na.rm=T)-quantile(tdData$traderatio,.25,na.rm=T)
prop = 4.6-0.01
## get the coef and vcov
coef_name <- c("traderatio","traderatio:issuesalience","traderatio:tissuesalience")
coef_mean <- coef(cmod2)[coef_name]
coef_vcov <- vcov(cmod2)[coef_name,coef_name]

#first difference function from Licht2011
hrtrade_issuesalience1 <- function(x, par01, par2, par3) (exp(prop*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par01) (exp(prop*(par01))-1)*100
t <- seq(1,5000, length=5000) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = log(t[i]+lamhat)
  impact0=impact1=numeric()
  for (j in 1:sim) {
    par = mvrnorm(1, mu=coef_mean, Sigma=coef_vcov)
    par01 = par[1]
    par12 = par[2]
    par13 = par[3]
    impact1[j] = hrtrade_issuesalience1(x, par01,par12, par13)
    impact0[j] = hrtrade_issuesalience0(x, par01)
  }
  temp = data.frame (time=i, firstdiff1=impact1, firstdiff0=impact0)
  plot_dat = rbind(plot_dat, temp)
}
plot_dat = reshape(plot_dat, varying=c("firstdiff1", "firstdiff0"), v.names="impact", timevar="type",times=c("firstdiff1","firstdiff0"), direction="long", new.row.names = NULL)
## Plot the impact of trade ratio
tradeimpactplot_c_main <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.001)+
  geom_smooth(colour="black",aes(linetype=type))+
  xlab("Days")+ylab(expression("First Difference"))+
  theme_bw() + theme(legend.position="bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_linetype_manual(values=c("dotted","solid"),
                        name="Salience", 
                        labels = c("Low","High"))
main_fd_c <-tradeimpactplot_c_main+coord_cartesian(ylim=c(-50,50))
ggsave(filename=paste0(file,"FDPlot_Main.pdf"),
       plot=main_fd_c)

## Now get the confidence interval
salienceplot_dat <- ddply(plot_dat[plot_dat$type=="firstdiff1",], .(time), summarise, mean=mean(impact,na.rm = TRUE), sd=sd(impact,na.rm = TRUE))
salienceplot <- ggplot(salienceplot_dat, aes(time,mean))+geom_line()+
  geom_ribbon(aes(ymin=mean-1.96*sd,ymax=mean+1.96*sd),, fill = "grey70", alpha=0.3)+
  theme_bw() + theme(legend.position="bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Days")+ylab(expression("First Difference"))+geom_hline(yintercept=0)
ggsave(filename=paste0(file,"FDPlot_sal2.pdf"),
       plot=salienceplot)

allin1FD <- grid.arrange(tradeimpactplot_c_main+coord_cartesian(ylim=c(-50,50))+theme(legend.position = c(.1,.1))
                         ,salienceplot+coord_cartesian(ylim=c(-50,50)),nrow=1)
ggsave(filename=paste0(file,"FDPlot_allin1.pdf"),
       plot=allin1FD, width=12, height=8) 


## Plot the survival rate
## Using the marginal approach
###############
## Impact of issuesalience and traderatio
summary(tdData$traderatio) #summary(tdData$powerratio); summary(tdData$defenseratio)
ran <- c(.01,4.6) #range to plot
## need to demonstrate what it means when traderatio is 1000
pdata <- data.frame(traderatio=rep(ran,times=length(1:5000)), #justify on contrasting that of the opponent
                    powerratio=median(tdData$powerratio,na.rm=TRUE),
                    defenseratio=median(tdData$defenseratio,na.rm=TRUE),
                    lossratio=median(tdData$lossratio,na.rm=TRUE),
                    time=rep(1:5000,times=2),
                    dispnum3=median(tdData$dispnum3,na.rm=TRUE),
                    joint_demo=1,contbinary=0)
## low salience
pdata2 <- pdata
pdata2$issuesalience=0
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time+lamhat)
pdata2$tpowerratio <- pdata2$powerratio*log(pdata2$time+lamhat)
pdata2$tissuesalience <- pdata2$issuesalience*log(pdata2$time+lamhat)
p1 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",palette="lancet",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(a) Low Salience", method="marginal",legend.title="Cost Ratio")

## High salience
pdata2 <- pdata
pdata2$issuesalience=1
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time+lamhat)
pdata2$tpowerratio <- pdata2$powerratio*log(pdata2$time+lamhat)
pdata2$tissuesalience <- pdata2$issuesalience*log(pdata2$time+lamhat)
p2 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",palette="lancet",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(b) High Salience", method="marginal",legend.title="Cost Ratio")


main_surv <- grid_arrange_shared_legend(p1,p2,nrow=1)


ggsave(filename=paste0(file,"SurvivalPlot.pdf"),
       plot=main_surv, width=10, height=6.6)  

####################
## 2. tt function model
####################
clusterform <- update.formula(fit1bc, ~. +tt(issuesalience) +tt(tradeissuesalience)+ 
                                tt(contbinary))

lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(clusterform, data=tdData, tt=function(x, t, ...) {x*log(t+lam[i])})
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
lamhat <- lamhat_maincluster <- lam[which.min(bic)] 
fit2c_tt <- coxph(clusterform, data=tdData,tt=function(x, t, ...) {x*log(t+lamhat)})
summary(fit2c_tt)$coefficients[,c(1,6)]

####################
## Latex table for clustered models
###################
cmod2_se <- summary(cmod2)$coefficients[,4]
fit2c_tt_se <- summary(fit2c_tt)$coefficients[,4]
fit1bc_se <- summary(fit1bc)$coefficients[,4]
library(stargazer)
cluster_tex <- stargazer(cmod2, fit2c_tt, fit1bc, single.row = FALSE,
                         se=list(cmod2_se,fit2c_tt_se,fit1bc_se), ci=TRUE,
                         title="Cox Regression Results", 
                         label="CoxTable",dep.var.labels.include = FALSE,
                         dep.var.caption  = "Days before Quitting",
                         column.labels   = c("Adjusted Models", "Original Models"),
                         column.separate = c(2, 1),
                         covariate.labels=c("Cost Ratio","High Salience",
                                            "Cost*HSalience",
                                            "Joint Democracy","Contiguity","Power Ratio",
                                            "Defense Ratio", "Troop Loss Ratio",
                                            "Time*Contiguity", "Time*HSalience","Time*Cost*HSalience",
                                            "tt(Contiguity)", "tt(High Salience)",
                                            "tt(Cost*HSalience)"
                         ),
                         order=c(1,2,13,3:9,14,12,10,11),
                         add.lines = list(c("Agjusted Method","Time Interaction","tt","No"),
                                          c("BIC", round(BIC(cmod2),2), round(BIC(fit2c_tt),2),round(BIC(fit1bc),2))),
                         keep.stat=c("wald","ll", "n"),
                         no.space=TRUE,
                         out=paste0(file,"ClusterTable.tex"), out.header = FALSE)

