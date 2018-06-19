##############
## This is R code for my attition paper
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
load("DurationTimeCov_ICB.RData")
tdData <- durB # change to durG for robustness check

#use network trade
tdData$joint_demo <- ifelse(tdData$demo1>5 & tdData$demo2>5, 1, 0)
## contest success function provides probability of winning
## cost ratio is closer to what I seek to capture
tdData$traderatio <- tdData$tradeshare1/(tdData$tradeshare2+.Machine$double.eps) 
tdData$issuesalience <- factor(ifelse(tdData$gravty %in% c(5,6),2,ifelse(tdData$gravty %in% c(2,3,4),1,0)), levels=c(0,1,2))
temp <- sapply(levels(tdData$issuesalience), function(x) as.integer(x == tdData$issuesalience))
colnames(temp) <- paste0("issuesalience",colnames(temp))
tdData <- cbind(tdData,temp)
tdData$traderatio.net <- tdData$tradeshare1.net/(tdData$tradeshare2.net+.Machine$double.eps)
tdData$defenseratio <- log(tdData$defense1/(tdData$defense2+.Machine$double.eps)+1)

## alternative measurement
## tdData$traderatio <- tdData$tradedepend1/(tdData$tradedepend2+.Machine$double.eps)
## tdData$traderatio.net <- tdData$tradedepend1.net/(tdData$tradedepend2.net+.Machine$double.eps)

## change to using net measurement
tdData$traderatio <- tdData$traderatio.net ## comment out if not use

## construct the interaction term after deciding on trade net or not
tdData$tradeissuesalience1 <- tdData$traderatio*tdData$issuesalience1
tdData$tradeissuesalience2 <- tdData$traderatio*tdData$issuesalience2

## file to save outputs; depend on data choice
file <- paste0("C:/Users/YULENG/OneDrive/Documents/RESEARCH/Manuscript/Duration/",
               "TradeshareNet/")#tweak to different files

## codes for getting the table of trade ratio >1000
## table1<- (tdData[which(tdData$traderatio>1000),c("crisno","tradeshare1","tradeshare2","tradeshare1.net","tradeshare2.net","ccode1","ccode2")])
## table1$ccode1 <- countrycode(table1$ccode1,'cown','country.name')
## table1$ccode2 <- countrycode(table1$ccode2,'cown','country.name')
## also for long table
## print(xtable(table1[,c("crisno","ccode1","ccode2","tradeshare1.net","tradeshare2.net","tradeshare1","tradeshare2")], type = "latex",digits=8, caption = "Trade Ratio Larger Than 1000."), file = paste0(file,"TradeRatio_hundred.tex"), include.rownames = F,floating=FALSE,tabular.environment='longtable',add.to.row = list(pos = list(0), command = "\\hline \\endhead "))
##############
## Basic Models and Tests

fit1b <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio, data=tdData);summary(fit1b)
fit1bc <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(fit1bc)
cox.zph(fit1b,transform='rank');cox.zph(fit1bc,transform='rank')
cox.zph(fit1b,transform='km');cox.zph(fit1bc,transform='km')

#############
## Graph the level of censoring and outliers
## which justify the use of time transformation is cox.zph

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
cmod <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(cmod)
cox.zph(cmod,transform='rank')

## new interaciton formula
## Using only variables that violate the assumption most
cmodform <- update.formula(cmod, ~.+tcontbinary+tissuesalience2+traderatio:tissuesalience2)
## select the best model
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  #tdData$ttraderatio <- tdData$traderatio*log(tdData$time+lam)
  lam0 <- lam[i]
  tdData$tcontbinary <- tdData$contbinary*log(tdData$time+lam0)
  tdData$tpowerratio <- tdData$powerratio*log(tdData$time+lam0)
  tdData$tissuesalience1 <- tdData$issuesalience1*log(tdData$time+lam0)
  tdData$tissuesalience2 <- tdData$issuesalience2*log(tdData$time+lam0)
  
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
tdData$tissuesalience1 <- tdData$issuesalience1*log(tdData$time+lamhat)
tdData$tissuesalience2 <- tdData$issuesalience2*log(tdData$time+lamhat)
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
##prop <- quantile(tdData$traderatio,.75,na.rm=T)-quantile(tdData$traderatio,.25,na.rm=T)
prop <- 100-0.01

## get the coef and vcov
coef_name <- c("traderatio","traderatio:issuesalience1","traderatio:issuesalience2","traderatio:tissuesalience2")
coef_mean <- coef(cmod2)[coef_name]
coef_vcov <- vcov(cmod2)[coef_name,coef_name]

#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par01, par2, par3) (exp(prop*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par01, par2) (exp(prop*(par01+par2))-1)*100
hrtrade_issuesalience0 <- function(x, par01) (exp(prop*(par01))-1)*100
t <- seq(1,1500, length=1500) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = log(t[i]+lamhat)
  impact0=impact1=impact2=numeric()
  for (j in 1:sim) {
    par = mvrnorm(1, mu=coef_mean, Sigma=coef_vcov)
    par01 = par[1]
    par12 = par[2]
    par22 = par[3]
    par23 = par[4]
    impact2[j] = hrtrade_issuesalience2(x, par01,par22, par23)
    impact1[j] = hrtrade_issuesalience1(x, par01, par12)
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
##main_fd_c <- tradeimpactplot_c_main+scale_y_continuous(breaks =c(-1,0,1),limits = c(-2,1))
main_fd_c <-tradeimpactplot_c_main+coord_cartesian(ylim=c(-20,20))
ggsave(filename=paste0(file,"FDPlot_Main.pdf"),
       plot=tradeimpactplot_c_main)

## Now get the confidence interval
salience2plot_dat <- ddply(plot_dat[plot_dat$type=="firstdiff2",], .(time), summarise, mean=mean(impact,na.rm = TRUE), sd=sd(impact,na.rm = TRUE))
salience2plot <- ggplot(salience2plot_dat, aes(time,mean))+geom_line()+
              geom_ribbon(aes(ymin=mean-1.96*sd,ymax=mean+1.96*sd), fill = "grey70", alpha=0.3)+
              theme_bw() + theme(legend.position="bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
              xlab("Days")+ylab(expression("First Difference"))+geom_hline(yintercept=0)
ggsave(filename=paste0(file,"FDPlot_sal2.pdf"),
       plot=salience2plot)

allin1FD <- grid.arrange(tradeimpactplot_c_main+coord_cartesian(ylim=c(-20,20))+theme(legend.position = c(.1,.1))
  ,salience2plot+coord_cartesian(ylim=c(-20,20)),nrow=1)
ggsave(filename=paste0(file,"FDPlot_allin1.pdf"),
       plot=allin1FD, width=12, height=8)    
  
###############
## Plot the survival rate
## Using the marginal approach
###############
## Impact of issuesalience and traderatio
summary(tdData$traderatio) #summary(tdData$powerratio); summary(tdData$defenseratio)
ran <- c(.01,100) #range to plot
## need to demonstrate what it means when traderatio is 1000
pdata <- data.frame(traderatio=rep(ran,times=length(1:1500)), #justify on contrasting that of the opponent
                    powerratio=median(tdData$powerratio,na.rm=TRUE),
                    defenseratio=median(tdData$defenseratio,na.rm=TRUE),
                    time=rep(1:1500,times=2),
                    crisno=median(tdData$crisno,na.rm=TRUE),
                    joint_demo=1,contbinary=0)
## low salience
pdata2 <- pdata
pdata2$issuesalience1=0
pdata2$issuesalience2=0
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time+lamhat)
pdata2$tpowerratio <- pdata2$powerratio*log(pdata2$time+lamhat)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time+lamhat)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time+lamhat)
p1 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",palette="lancet",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(a1) Low Salience", method="marginal",legend.title="Cost Ratio")

## Median salience
pdata2 <- pdata
pdata2$issuesalience1=1
pdata2$issuesalience2=0
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time+lamhat)
pdata2$tpowerratio <- pdata2$powerratio*log(pdata2$time+lamhat)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time+lamhat)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time+lamhat)
p2 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",palette="lancet",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(b1) Median Salience", method="marginal",legend.title="Cost Ratio")

## High salience
pdata2 <- pdata
pdata2$issuesalience1=0
pdata2$issuesalience2=1
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time+lamhat)
pdata2$tpowerratio <- pdata2$powerratio*log(pdata2$time+lamhat)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time+lamhat)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time+lamhat)
p3 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",palette="lancet",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(c1) High Salience", method="marginal",legend.title="Cost Ratio")


main_surv <- grid_arrange_shared_legend(p1,p2,p3,nrow=1)


p1z <- p1+coord_cartesian(xlim = c(500, 1500),ylim=c(0,.25))+labs(title="(a2) Low Salience")
p2z <- p2+coord_cartesian(xlim = c(500, 1500),ylim=c(0,.25))+labs(title="(b2) Median Salience")
p3z <- p3+coord_cartesian(xlim = c(500, 1500),ylim=c(0,.25))+labs(title="(c2) High Salience")
zoomin_sur <- grid_arrange_shared_legend(p1,p2,p3,
                                         p1z,p2z,p3z,nrow=2,ncol=3)

ggsave(filename=paste0(file,"SurvivalPlot.pdf"),
       plot=main_surv, width=12, height=5)  
ggsave(filename=paste0(file,"SurvivalPlotZoomin.pdf"),
       plot=zoomin_sur, width=12, height=10)  

####################
## 2. tt function model
####################
clusterform <- update.formula(fit1bc, ~. +tt(issuesalience2) +tt(tradeissuesalience2)+ 
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
                         covariate.labels=c("Cost Ratio","Median Salience", "High Salience",
                                            "Cost*MSalience","Cost*HSalience",
                                            "Joint Democracy","Contiguity","Power Ratio",
                                            "Defense Ratio", 
                                            "Time*Contiguity", "Time*HSalience","Time*Cost*HSalience",
                                            "tt(Contiguity)", "tt(High Salience)",
                                            "tt(Cost*HSalience)"
                         ),
                         order=c(1,2,3,13,14,4,5,6,7,8,9,15,13,11,12), 
                         add.lines = list(c("Agjusted Method","Time Interaction","tt","No"),
                                          c("BIC", round(BIC(cmod2),2), round(BIC(fit2c_tt),2),round(BIC(fit1bc),2))),
                         keep.stat=c("wald","ll", "n"),
                         no.space=TRUE,
                         out=paste0(file,"ClusterTable.tex"), out.header = FALSE)


#############################################
#############################################
## Additional code for Robustness check
## 0. using only the interaction suggested by the nph test
## 1. using more interaction
## 2. using the frailty model 
## 3. using the trade dependence measurement
## 4. using log for tradeshare measurement
#############################################
#############################################

#####################
## 0. using only the interaction suggested by the nph test
file = "C:/Users/YULENG/OneDrive/Documents/RESEARCH/Manuscript/Duration/InteractAllNPHVariables/"

fit1b <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio, data=tdData);summary(fit1b)
fit1bc <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(fit1bc)
cox.zph(fit1b,transform='rank');cox.zph(fit1bc,transform='rank')
cox.zph(fit1b,transform='km');cox.zph(fit1bc,transform='km')

#############
## Graph the level of censoring and outliers
## which justify the use of time transformation is cox.zph

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
## 0.1 Main Model Clustered Model
#####################

## conventional way using time interaction
## prep
tdData$time <- tdData$tstop ## using fox
cmod <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(cmod)
cox.zph(cmod,transform='rank')

## new interaciton formula
## Using only variables that violate the assumption most
cmodform <- update.formula(cmod, ~.+tcontbinary+tissuesalience1+traderatio:tissuesalience1+tissuesalience2+traderatio:tissuesalience2)
## select the best model
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  lam0 <- lam[i]
  tdData$tcontbinary <- tdData$contbinary*log(tdData$time+lam0)
  tdData$tpowerratio <- tdData$powerratio*log(tdData$time+lam0)
  tdData$tissuesalience1 <- tdData$issuesalience1*log(tdData$time+lam0)
  tdData$tissuesalience2 <- tdData$issuesalience2*log(tdData$time+lam0)
  
  fit <- coxph(cmodform, data=tdData)
  logL[i] <- logLik(fit)
  aic[i] <- AIC(fit)
  bic[i] <- BIC(fit)
}
## By BIC; others give similar results
lamhat <- lam[which.min(bic)]
tdData$ttraderatio <- tdData$traderatio*log(tdData$time+lamhat)
tdData$tcontbinary <- tdData$contbinary*log(tdData$time+lamhat)
tdData$tpowerratio <- tdData$powerratio*log(tdData$time+lamhat)
tdData$tissuesalience1 <- tdData$issuesalience1*log(tdData$time+lamhat)
tdData$tissuesalience2 <- tdData$issuesalience2*log(tdData$time+lamhat)
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
##prop <- quantile(tdData$traderatio,.75,na.rm=T)-quantile(tdData$traderatio,.25,na.rm=T)
prop <- 100-0.01

## get the coef and vcov
coef_name <- c("traderatio","traderatio:issuesalience1","traderatio:tissuesalience1","traderatio:issuesalience2","traderatio:tissuesalience2")
coef_mean <- coef(cmod2)[coef_name]
coef_vcov <- vcov(cmod2)[coef_name,coef_name]

#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par01, par2, par3) (exp(prop*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par01, par2, par3) (exp(prop*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par01) (exp(prop*(par01))-1)*100
t <- seq(1,1500, length=1500) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = log(t[i]+lamhat)
  impact0=impact1=impact2=numeric()
  for (j in 1:sim) {
    par = mvrnorm(1, mu=coef_mean, Sigma=coef_vcov)
    par01 = par[1]
    par12 = par[2]
    par13 = par[3]
    par22 = par[4]
    par23 = par[5]
    impact2[j] = hrtrade_issuesalience2(x, par01, par22, par23)
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
##main_fd_c <- tradeimpactplot_c_main+scale_y_continuous(breaks =c(-1,0,1),limits = c(-2,1))
main_fd_c <-tradeimpactplot_c_main+coord_cartesian(ylim=c(-50,50))
ggsave(filename=paste0(file,"FDPlot_Main.pdf"),
       plot=tradeimpactplot_c_main)

## Now get the confidence interval
salience2plot_dat <- ddply(plot_dat[plot_dat$type=="firstdiff2",], .(time), summarise, mean=mean(impact,na.rm = TRUE), sd=sd(impact,na.rm = TRUE))
salience2plot <- ggplot(salience2plot_dat, aes(time,mean))+geom_line()+
  geom_ribbon(aes(ymin=mean-1.96*sd,ymax=mean+1.96*sd), fill = "grey70", alpha=0.3)+
  theme_bw() + theme(legend.position="bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Days")+ylab(expression("First Difference"))+geom_hline(yintercept=0)
ggsave(filename=paste0(file,"FDPlot_sal2.pdf"),
       plot=salience2plot)

allin1FD <- grid.arrange(tradeimpactplot_c_main+coord_cartesian(ylim=c(-20,20))+theme(legend.position = c(.1,.1))
                         ,salience2plot+coord_cartesian(ylim=c(-20,20)),nrow=1)
ggsave(filename=paste0(file,"FDPlot_allin1.pdf"),
       plot=allin1FD, width=12, height=8)    

###############
## Plot the survival rate
## Using the marginal approach
###############
## Impact of issuesalience and traderatio
summary(tdData$traderatio) #summary(tdData$powerratio); summary(tdData$defenseratio)
ran <- c(.01,100) #range to plot
## need to demonstrate what it means when traderatio is 1000
pdata <- data.frame(traderatio=rep(ran,times=length(1:1500)), #justify on contrasting that of the opponent
                    powerratio=median(tdData$powerratio,na.rm=TRUE),
                    defenseratio=median(tdData$defenseratio,na.rm=TRUE),
                    time=rep(1:1500,times=2),
                    crisno=median(tdData$crisno,na.rm=TRUE),
                    joint_demo=1,contbinary=0)
## low salience
pdata2 <- pdata
pdata2$issuesalience1=0
pdata2$issuesalience2=0
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tpowerratio <- pdata2$powerratio*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p1 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",palette="lancet",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(a1) Low Salience", method="marginal",legend.title="Cost Ratio")

## Median salience
pdata2 <- pdata
pdata2$issuesalience1=1
pdata2$issuesalience2=0
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tpowerratio <- pdata2$powerratio*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p2 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",palette="lancet",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(b1) Median Salience", method="marginal",legend.title="Cost Ratio")

## High salience
pdata2 <- pdata
pdata2$issuesalience1=0
pdata2$issuesalience2=1
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tpowerratio <- pdata2$powerratio*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p3 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",palette="lancet",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(c1) High Salience", method="marginal",legend.title="Cost Ratio")


main_surv <- grid_arrange_shared_legend(p1,p2,p3,nrow=1)


p1z <- p1+coord_cartesian(xlim = c(500, 1500),ylim=c(0,.1))+labs(title="(a2) Low Salience")
p2z <- p2+coord_cartesian(xlim = c(500, 1500),ylim=c(0.5,.7))+labs(title="(b2) Median Salience")
p3z <- p3+coord_cartesian(xlim = c(500, 1500),ylim=c(.4,.65))+labs(title="(c2) High Salience")
zoomin_sur <- grid_arrange_shared_legend(p1,p2,p3,
                                         p1z,p2z,p3z,nrow=2,ncol=3)

ggsave(filename=paste0(file,"SurvivalPlot.pdf"),
       plot=main_surv, width=12, height=5)  
ggsave(filename=paste0(file,"SurvivalPlotZoomin.pdf"),
       plot=zoomin_sur, width=12, height=10)  

####################
## 0.2 tt function model
####################
clusterform <- update.formula(fit1bc, ~. +tt(issuesalience1) +tt(tradeissuesalience1)+tt(issuesalience2) +tt(tradeissuesalience2)+ 
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
cluster_tex <- stargazer(cmod2, fit2c_tt, fit1bc, single.row = TRUE,
                         se=list(cmod2_se,fit2c_tt_se,fit1bc_se), ci=TRUE,
                         title="Cox Regression Results", 
                         label="CoxTable",dep.var.labels.include = FALSE,
                         dep.var.caption  = "Days before Quitting",
                         column.labels   = c("Adjusted Models", "Original Models"),
                         column.separate = c(2, 1),
                         covariate.labels=c("Cost Ratio","Median Salience", "High Salience",
                                            "Cost*MSalience","Cost*HSalience",
                                            "Joint Democracy","Contiguity","Power Ratio",
                                            "Defense Ratio", 
                                            "Time*Contiguity", "Time*MSalience","Time*HSalience",
                                            "Time*Cost*MSalience","Time*Cost*HSalience",
                                            "tt(Median Salience)","tt(Cost*MSalience)",
                                            "tt(High Salience)",
                                            "tt(Cost*HSalience)","tt(Contiguity)"
                         ),
                         order=c(1:3,16,17,4:10,18,19,11:15), 
                         add.lines = list(c("Agjusted Method","Time Interaction","tt","No"),
                                          c("BIC", round(BIC(cmod2),2), round(BIC(fit2c_tt),2),round(BIC(fit1bc),2))),
                         keep.stat=c("wald","ll", "n"),
                         no.space=TRUE,
                         out=paste0(file,"ClusterTable.tex"), out.header = FALSE)




###################
## 1. using more interaction
##############
## Basic Models and Tests
file = "C:/Users/YULENG/OneDrive/Documents/RESEARCH/Manuscript/Duration/MoreInteraction/"

fit1b <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio, data=tdData);summary(fit1b)
fit1bc <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(fit1bc)
cox.zph(fit1b,transform='rank');cox.zph(fit1bc,transform='rank')
cox.zph(fit1b,transform='km');cox.zph(fit1bc,transform='km')

#############
## Graph the level of censoring and outliers
## which justify the use of time transformation is cox.zph

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
## 1.1 Main Model Clustered Model
#####################

## conventional way using time interaction
## prep
tdData$time <- tdData$tstop ## using fox
cmod <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(cmod)
cox.zph(cmod,transform='rank')

## new interaciton formula
## Using only variables that violate the assumption most
cmodform <- update.formula(cmod, ~.+ttraderatio+tcontbinary+tissuesalience1+traderatio:tissuesalience1+tissuesalience2+traderatio:tissuesalience2)
## select the best model
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  lam0 <- lam[i]
  tdData$ttraderatio <- tdData$traderatio*log(tdData$time+lam0)
  tdData$tcontbinary <- tdData$contbinary*log(tdData$time+lam0)
  tdData$tpowerratio <- tdData$powerratio*log(tdData$time+lam0)
  tdData$tissuesalience1 <- tdData$issuesalience1*log(tdData$time+lam0)
  tdData$tissuesalience2 <- tdData$issuesalience2*log(tdData$time+lam0)
  
  fit <- coxph(cmodform, data=tdData)
  logL[i] <- logLik(fit)
  aic[i] <- AIC(fit)
  bic[i] <- BIC(fit)
}
## By BIC; others give similar results
lamhat <- lam[which.min(bic)]
tdData$ttraderatio <- tdData$traderatio*log(tdData$time+lamhat)
tdData$tcontbinary <- tdData$contbinary*log(tdData$time+lamhat)
tdData$tpowerratio <- tdData$powerratio*log(tdData$time+lamhat)
tdData$tissuesalience1 <- tdData$issuesalience1*log(tdData$time+lamhat)
tdData$tissuesalience2 <- tdData$issuesalience2*log(tdData$time+lamhat)
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
##prop <- quantile(tdData$traderatio,.75,na.rm=T)-quantile(tdData$traderatio,.25,na.rm=T)
prop <- 100-0.01

## get the coef and vcov
coef_name <- c("traderatio","ttraderatio","traderatio:issuesalience1","traderatio:tissuesalience1","traderatio:issuesalience2","traderatio:tissuesalience2")
coef_mean <- coef(cmod2)[coef_name]
coef_vcov <- vcov(cmod2)[coef_name,coef_name]

#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par01, par02, par2, par3) (exp(prop*(par01+par02*x+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par01, par02, par2, par3) (exp(prop*(par01+par02*x+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par01, par02) (exp(prop*(par01+par02*x))-1)*100
t <- seq(1,1500, length=1500) # time length by day
sim <- 100 # number of simulation
set.seed(11) 
plot_dat <- data.frame()
for (i in 1:length(t)){
  x = log(t[i]+lamhat)
  impact0=impact1=impact2=numeric()
  for (j in 1:sim) {
    par = mvrnorm(1, mu=coef_mean, Sigma=coef_vcov)
    par01 = par[1]
    par02 = par[2]
    par12 = par[3]
    par13 = par[4]
    par22 = par[5]
    par23 = par[6]
    impact2[j] = hrtrade_issuesalience2(x, par01, par02, par22, par23)
    impact1[j] = hrtrade_issuesalience1(x, par01, par02, par12, par13)
    impact0[j] = hrtrade_issuesalience0(x, par01, par02)
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
##main_fd_c <- tradeimpactplot_c_main+scale_y_continuous(breaks =c(-1,0,1),limits = c(-2,1))
main_fd_c <-tradeimpactplot_c_main+coord_cartesian(ylim=c(-50,50))
ggsave(filename=paste0(file,"FDPlot_Main.pdf"),
       plot=main_fd_c)

## Now get the confidence interval
salience2plot_dat <- ddply(plot_dat[plot_dat$type=="firstdiff2",], .(time), summarise, mean=mean(impact,na.rm = TRUE), sd=sd(impact,na.rm = TRUE))
salience2plot <- ggplot(salience2plot_dat, aes(time,mean))+geom_line()+
  geom_ribbon(aes(ymin=mean-1.96*sd,ymax=mean+1.96*sd),, fill = "grey70", alpha=0.3)+
  theme_bw() + theme(legend.position="bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Days")+ylab(expression("First Difference"))+geom_hline(yintercept=0)
ggsave(filename=paste0(file,"FDPlot_sal2.pdf"),
       plot=salience2plot)

allin1FD <- grid.arrange(tradeimpactplot_c_main+coord_cartesian(ylim=c(-20,20))+theme(legend.position = c(.1,.1))
                         ,salience2plot+coord_cartesian(ylim=c(-20,20)),nrow=1)
ggsave(filename=paste0(file,"FDPlot_allin1.pdf"),
       plot=allin1FD, width=12, height=8)    

###############
## Plot the survival rate
## Using the marginal approach
###############
## Impact of issuesalience and traderatio
summary(tdData$traderatio) #summary(tdData$powerratio); summary(tdData$defenseratio)
ran <- c(.01,100) #range to plot
## need to demonstrate what it means when traderatio is 1000
pdata <- data.frame(traderatio=rep(ran,times=length(1:1500)), #justify on contrasting that of the opponent
                    powerratio=median(tdData$powerratio,na.rm=TRUE),
                    defenseratio=median(tdData$defenseratio,na.rm=TRUE),
                    time=rep(1:1500,times=2),
                    crisno=median(tdData$crisno,na.rm=TRUE),
                    joint_demo=1,contbinary=0)
## low salience
pdata2 <- pdata
pdata2$issuesalience1=0
pdata2$issuesalience2=0
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tpowerratio <- pdata2$powerratio*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p1 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",palette="lancet",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(a1) Low Salience", method="marginal",legend.title="Cost Ratio")

## Median salience
pdata2 <- pdata
pdata2$issuesalience1=1
pdata2$issuesalience2=0
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tpowerratio <- pdata2$powerratio*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p2 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",palette="lancet",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(b1) Median Salience", method="marginal",legend.title="Cost Ratio")

## High salience
pdata2 <- pdata
pdata2$issuesalience1=0
pdata2$issuesalience2=1
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tpowerratio <- pdata2$powerratio*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p3 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",palette="lancet",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(c1) High Salience", method="marginal",legend.title="Cost Ratio")


main_surv <- grid_arrange_shared_legend(p1,p2,p3,nrow=1)


p1z <- p1+coord_cartesian(xlim = c(500, 1500),ylim=c(0,.25))+labs(title="(a2) Low Salience")
p2z <- p2+coord_cartesian(xlim = c(500, 1500),ylim=c(0,.25))+labs(title="(b2) Median Salience")
p3z <- p3+coord_cartesian(xlim = c(500, 1500),ylim=c(0,.25))+labs(title="(c2) High Salience")
zoomin_sur <- grid_arrange_shared_legend(p1,p2,p3,
                                         p1z,p2z,p3z,nrow=2,ncol=3)

ggsave(filename=paste0(file,"SurvivalPlot.pdf"),
       plot=main_surv, width=12, height=5)  
ggsave(filename=paste0(file,"SurvivalPlotZoomin.pdf"),
       plot=zoomin_sur, width=12, height=10)  

####################
## 1.2 tt function model
####################
clusterform <- update.formula(fit1bc, ~. +tt(traderatio)+tt(issuesalience1) +tt(tradeissuesalience1)+tt(issuesalience2) +tt(tradeissuesalience2)+ 
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
cluster_tex <- stargazer(cmod2, fit2c_tt, fit1bc, single.row = TRUE,
                         se=list(cmod2_se,fit2c_tt_se,fit1bc_se), ci=TRUE,
                         title="Cox Regression Results", 
                         label="CoxTable",dep.var.labels.include = FALSE,
                         dep.var.caption  = "Days before Quitting",
                         column.labels   = c("Adjusted Models", "Original Models"),
                         column.separate = c(2, 1),
                         covariate.labels=c("Cost Ratio","Median Salience", "High Salience",
                                            "Cost*MSalience","Cost*HSalience",
                                            "Joint Democracy","Contiguity","Power Ratio",
                                            "Defense Ratio", "Time*Traderatio",
                                            "Time*Contiguity", "Time*MSalience","Time*HSalience",
                                            "Time*Cost*MSalience","Time*Cost*HSalience",
                                            "tt(Cost Ratio)","tt(Median Salience)","tt(Cost*MSalience)",
                                             "tt(High Salience)",
                                            "tt(Cost*HSalience)","tt(Contiguity)"
                         ),
                         order=c(1:3,18,19,4:11,20,21,12:17), 
                         add.lines = list(c("Agjusted Method","Time Interaction","tt","No"),
                                          c("BIC", round(BIC(cmod2),2), round(BIC(fit2c_tt),2),round(BIC(fit1bc),2))),
                         keep.stat=c("wald","ll", "n"),
                         no.space=TRUE,
                         out=paste0(file,"ClusterTable.tex"), out.header = FALSE)




######
## 2. using the frailty model 
fit1bf <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio+frailty(crisno), data=tdData);summary(fit1bf)
fit1bf_me <- coxme(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio+(1|crisno), data=tdData);summary(fit1bf_me)

## Using traditional interaction results in error
## Error in if (y[iter] == max(y) && x[iter] == max(x)) newtheta <- 2 * max(x) else newtheta <- frailty.brent(sqrt(x),  : 
## missing value where TRUE/FALSE needed
frailtyform <- update.formula(fit1bf, ~. +tt(issuesalience2) +tt(tradeissuesalience2)+ 
                                tt(contbinary))

lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(frailtyform, data=tdData, tt=function(x, t, ...) {x*log(t+lam[i])})
  logL[i] <- logLik(fit)
  aic[i] <- AIC(fit)
  bic[i] <- BIC(fit)
}

## now use the best lam to refit the model
lamhat <- lam[which.min(bic)] 
fit2f_tt <- coxph(frailtyform, data=tdData,tt=function(x, t, ...) {x*log(t+lamhat)})
summary(fit2f_tt)$coefficients[,c(1,6)]

####################
## Latex table for frailty models
###################

## Extract statistics
fit2f_tt_coef <- coef(fit2f_tt)
fit1bf_coef <- coef(fit1bf)
fit1bf_me_coef <- coef(fit1bf_me)

fit2f_tt_se <- summary(fit2f_tt)$coefficients[-8,4]
fit1bf_se <- summary(fit1bf)$coefficients[-8,4]
fit1bf_me_se <- sqrt(diag(vcov(fit1bf_me)))

## Now fake models without the frailty term for tables
fit2f_tt_fake <- coxph(update.formula(fit2f_tt, ~. -frailty(crisno)), data=tdData)
fit1bf_fake <- coxph(update.formula(fit1bf, ~. -frailty(crisno)), data=tdData)
fit1bf_me_fake <- coxph(update.formula(fit1bf, ~. -frailty(crisno)), data=tdData)


library(stargazer)
frailty_tex <- stargazer(fit2f_tt_fake,fit1bf_fake,fit1bf_me_fake, single.row = FALSE,
                         coef=list(fit2f_tt_coef,fit1bf_coef,fit1bf_me_coef),
                         se=list(fit2f_tt_se,fit1bf_se,fit1bf_me_se), ci=TRUE,
                         title="Cox Regression Results for Frailty Models", 
                         label="CoxTable",dep.var.labels.include = FALSE,
                         dep.var.caption  = "Days before Quitting",
                         column.labels   = c("Adjusted Models", "Frailty Models"),
                         column.separate = c(1, 2),
                         covariate.labels=c("Cost Ratio","Median Salience", "High Salience",
                                            "Joint Democracy","Contiguity","Power Ratio",
                                            "Defense Ratio",  "tt(High Salience)",
                                            "tt(Cost*HSalience)", 
                                            "tt(Contiguity)",
                                            "Cost*MSalience","Cost*HSalience"
                         ),
                         add.lines = list(c("Agjusted Method","tt","No","No"),
                                          c("BIC", round(BIC(fit2f_tt),2), round(BIC(fit1bf),2),round(BIC(fit1bf_me),2))
                                          ),
                         keep.stat=c("n"),
                         no.space=TRUE,
                         out=paste0(file,"FrailtyTable.tex"), out.header = FALSE)

###########
## 3. using the trade dependence measurement
tdData$traderatio.net <- tdData$tradedepend1.net/(tdData$tradedepend2.net+.Machine$double.eps)

## change to using net measurement
tdData$traderatio <- tdData$traderatio.net ## comment out if not use

## construct the interaction term after deciding on trade net or not
tdData$tradeissuesalience1 <- tdData$traderatio*tdData$issuesalience1
tdData$tradeissuesalience2 <- tdData$traderatio*tdData$issuesalience2
#######
## Now use similar codes as 
## .......
#######

## a number of things to take note
## results too large, need to rescale to observe
tradeimpactplot_c_main <- tradeimpactplot_c_main+scale_y_continuous(limits = c(-50,50))

allin1FD <- grid.arrange(tradeimpactplot_c_main+scale_y_continuous(limits = c(-50,50))+theme(legend.position = c(.1,.1))
                         ,salience2plot+scale_y_continuous(limits = c(-50,50)),nrow=1)


###########
## 4. using the trade dependence measurement
tdData$traderatio <- tdData$tradeshare1/(tdData$tradeshare2+.Machine$double.eps) 
tdData$traderatio.net <- tdData$tradeshare1.net/(tdData$tradeshare2.net+.Machine$double.eps)
## change to using net measurement
tdData$traderatio <- tdData$traderatio.net ## comment out if not use

## construct the interaction term after deciding on trade net or not
tdData$tradeissuesalience1 <- tdData$traderatio*tdData$issuesalience1
tdData$tradeissuesalience2 <- tdData$traderatio*tdData$issuesalience2

## file to save outputs; depend on data choice
file <- paste0("C:/Users/YULENG/OneDrive/Documents/RESEARCH/Manuscript/Duration/",
               "Tradeshare_Net_log/")#tweak to different files

prop <- 4.6-0.01
ran <- c(0.01,4.6)

allin1FD <- grid.arrange(tradeimpactplot_c_main+coord_cartesian(ylim=c(-100,100))+theme(legend.position = c(.1,.1))
                                                   ,salience2plot+coord_cartesian(ylim=c(-100,100)),nrow=1)
