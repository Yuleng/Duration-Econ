##############
## This is R code for my attition paper
#############
## 0. Preliminary stuff (data preparation, tests)
## 1. Main model (clustered and tt): First Difference Plot; Proxy survival Plot
## 2. Proxy clustered without tt
## 3. Frailty Model with tt
## 4. Original Model without time interaction

##############
## libraries and functions
#############

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
## Data
############
## load time dependent covariates data
load("DurationTimeCov_ICB.RData")
tdData <- durB # change to durG for robustness check

#use network trade
tdData$joint_demo <- ifelse(tdData$demo1>5 & tdData$demo2>5, 1, 0)
## contest success function provides probability of winning
## cost ratio is closer to what I seek to capture
tdData$traderatio <- tdData$tradeshare1/(tdData$tradeshare2+.Machine$double.eps) ## not like power ratio, because the latter produce too many NA values
tdData$issuesalience <- factor(ifelse(tdData$gravty %in% c(5,6),2,ifelse(tdData$gravty %in% c(2,3,4),1,0)), levels=c(0,1,2))
temp <- sapply(levels(tdData$issuesalience), function(x) as.integer(x == tdData$issuesalience))
colnames(temp) <- paste0("issuesalience",colnames(temp))
tdData <- cbind(tdData,temp)
tdData$traderatio.net <- tdData$tradeshare1.net/(tdData$tradeshare2.net+.Machine$double.eps)
tdData$defenseratio <- log(tdData$defense1/(tdData$defense2+.Machine$double.eps)+1)
tdData$tradeissuesalience1 <- tdData$traderatio*tdData$issuesalience1
tdData$tradeissuesalience2 <- tdData$traderatio*tdData$issuesalience2
file <- paste0("C:/Users/YULENG/OneDrive/Documents/RESEARCH/Manuscript/Duration/ISQ/",
          "TradeshareNet/")#tweak to different files
## tdData$traderatio <- tdData$tradedepend1/(tdData$tradedepend2+.Machine$double.eps) ## not like power ratio, because the latter produce too many NA values
## tdData$traderatio.net <- tdData$tradedepend1.net/(tdData$tradedepend2.net+.Machine$double.eps)

## tdData$traderatio <- tdData$traderatio.net

##############
## Basic Models and Tests
##############
fit1b <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+tradeissuesalience1+tradeissuesalience2+joint_demo+contbinary+powerratio+defenseratio, data=tdData);summary(fit1b)
fit1bc <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+tradeissuesalience1+tradeissuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(fit1bc)
fit1bf<- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+tradeissuesalience1+tradeissuesalience2+joint_demo+contbinary+powerratio+defenseratio+frailty(crisno), data=tdData);summary(fit1bf)
fit1bf_me<- coxme(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+tradeissuesalience1+tradeissuesalience2+joint_demo+contbinary+powerratio+defenseratio+(1|crisno), data=tdData);summary(fit1bf_me)
## coxme is prefered in frailty model

cox.zph(fit1b,transform='rank');cox.zph(fit1bc,transform='rank');cox.zph(fit1bf,transform='rank')
cox.zph(fit1b,transform='km');cox.zph(fit1bc,transform='km');cox.zph(fit1bf,transform='km')


#############
## Graph the level of censoring and outliers
## which justify the use of time transformation is cox.zph
#############
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

####################
## tt without traderatio
####################
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(update.formula(fit1bc, ~. + tt(issuesalience1) +tt(tradeissuesalience1)
                              + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)+tt(powerratio)),
               data=tdData, tt=function(x, t, ...) {x*log(t+lam[i])})
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
best_fitc_main <- coxph(update.formula(fit1bc, ~. + tt(issuesalience1) +tt(tradeissuesalience1)
                                       + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)+tt(powerratio)),
                        data=tdData,tt=function(x, t, ...) {x*log(t+lamhat)})
summary(best_fitc_main)$coefficients[,c(1,6)]

##########
## First Difference of Trade Ratio
## tt without traderatio
##########
## check the distribution of traderatio
library(car)
densityPlot(tdData$traderatio,na.rm=T);
densityPlot(tdData$traderatio[tdData$traderatio<quantile(tdData$traderatio,.75,na.rm=T)],na.rm=T)

## proportional of effect
prop <- quantile(tdData$traderatio,.75,na.rm=T)-quantile(tdData$traderatio,.25,na.rm=T)
par_traderatio <- summary(best_fitc_main)$coefficients["traderatio",c(1,4)]
par_tradeissuesalience1 <- summary(best_fitc_main)$coefficients["tradeissuesalience1",c(1,4)]
par_ttradeissuesalience1 <- summary(best_fitc_main)$coefficients["tt(tradeissuesalience1)",c(1,4)]
par_tradeissuesalience2 <- summary(best_fitc_main)$coefficients["tradeissuesalience2",c(1,4)]
par_ttradeissuesalience2 <- summary(best_fitc_main)$coefficients["tt(tradeissuesalience2)",c(1,4)]
#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par01, par2, par3) (exp(prop*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par01, par2, par3) (exp(prop*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par01) (exp(prop*(par01))-1)*100
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
main_fd_c <- tradeimpactplot_c_main+scale_y_continuous(breaks =c(-1,0,1),limits = c(-2,1))
ggsave(filename=paste0(file,"FDPlot_Main.pdf"),
       plot=main_fd_c)

#######################
## Include tt for traderatio
###########
## COX.zph test indicate trade salience condiguity
cox.zph(fit1bc,transform='rank');cox.zph(fit1bc,transform='km')

## select the best model
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(update.formula(fit1bc, ~. +tt(traderatio)+ tt(issuesalience1) +tt(tradeissuesalience1)
                              + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)+tt(powerratio)),
               data=tdData, tt=function(x, t, ...) {x*log(t+lam[i])})
  logL[i] <- logLik(fit)
  aic[i] <- AIC(fit)
  bic[i] <- BIC(fit)
}
## By BIC; others give similar results
lamhat <- lam[which.min(bic)]
best_fitc <- coxph(update.formula(fit1bc, ~. +tt(traderatio)+ tt(issuesalience1) +tt(tradeissuesalience1)
                                  + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)+tt(powerratio)),
                   data=tdData,tt=function(x, t, ...) {x*log(t+lamhat)})
summary(best_fitc)
summary(best_fitc)$coefficients[,c(1,6)]
cox.zph(best_fitc,transform='rank');cox.zph(best_fitc,transform='km')

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
hrtrade_issuesalience2 <- function(x, par01, par02, par2, par3) (exp(prop*(par01+par02*x+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par01, par02, par2, par3) (exp(prop*(par01+par02*x+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par01, par02) (exp(prop*(par01+par02*x))-1)*100
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
    impact2[j] = hrtrade_issuesalience2(x, par01,par02,par22, par23)
    impact1[j] = hrtrade_issuesalience1(x, par01, par02,par12, par13)
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

main_fd_tttrade <- tradeimpactplot_c+scale_y_continuous(breaks =c(-5,0,5),limits = c(-10,5))
## tradeshare withou net
## use this
## main_fd_tttrade <- tradeimpactplot_c+scale_y_continuous(breaks =c(-5,0,5),limits = c(-20,20))
ggsave(filename=paste0(file,"FDPlotTTtrade.pdf"),
       plot=main_fd_tttrade)

#grid_arrange_shared_legend(main_fd_c,main_fd_tttrade,nrow=1)


########################
## 2. Proxy clustered model and survival plot
########################
tdData$time <- tdData$tstop ## using fox
cmod <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(cmod)
cox.zph(cmod,transform='rank')

## select the best model
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  #tdData$ttraderatio <- tdData$traderatio*log(tdData$time+lam)
  tdData$tcontbinary <- tdData$contbinary*log(tdData$time+lam)
  tdData$tpowerratio <- tdData$powerratio*log(tdData$time+lam)
  tdData$tissuesalience1 <- tdData$issuesalience1*log(tdData$time+lam)
  tdData$tissuesalience2 <- tdData$issuesalience2*log(tdData$time+lam)
  
  fit <- coxph(update.formula(cmod, ~.+tcontbinary+tpowerratio+tissuesalience1+tissuesalience2+
                                traderatio:tissuesalience1+traderatio:tissuesalience2), data=tdData)
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
cmod2 <- coxph(update.formula(cmod, ~.+tcontbinary+tpowerratio+tissuesalience1+tissuesalience2+
                                traderatio:tissuesalience1+traderatio:tissuesalience2), data=tdData);summary(cmod2)



## using the lamhat from the main model




tdData$tcontbinary <- tdData$contbinary*log(tdData$time+lamhat_maincluster)
tdData$tpowerratio <- tdData$powerratio*log(tdData$time+lamhat_maincluster)
tdData$tissuesalience1 <- tdData$issuesalience1*log(tdData$time+lamhat_maincluster)
tdData$tissuesalience2 <- tdData$issuesalience2*log(tdData$time+lamhat_maincluster)
cmod2 <- coxph(update.formula(cmod, ~.+tcontbinary+tpowerratio+tissuesalience1+tissuesalience2+
                                traderatio:tissuesalience1+traderatio:tissuesalience2), data=tdData);summary(cmod2)
summary(cmod2)$coefficients[,c(1,6)]
summary(best_fitc_main)$coefficients[,c(1,6)]



###############
## Plot the proxy survival rate
## Using the marginal approach
###############
## Impact of issuesalience and traderatio
summary(tdData$traderatio)
ran <- c(.01,100) #range to plot
pdata <- data.frame(traderatio=rep(ran,times=length(1:1500)), #justify on contrasting that of the opponent
                    powerratio=median(tdData$powerratio,na.rm=TRUE),
                    defenseratio=median(tdData$defenseratio,na.rm=TRUE),
                    time=rep(1:1500,times=2),
                    crisno=median(tdData$crisno,na.rm=TRUE),
                    joint_demo=0,contbinary=1)
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
                       main="(a) Low Salience", method="marginal",legend.title="Cost Ratio")+
  annotate("text", x = 700, y = .7, label = "Low Cost")

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
                       main="(b) Median Salience", method="marginal",legend.title="Cost Ratio")+
  annotate("text", x = 700, y = .7, label = " ")

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
                       main="(c) High Salience", method="marginal",legend.title="Cost Ratio")+
  annotate("text", x = 700, y = .7, label = "High Cost")


main_surv <- grid_arrange_shared_legend(p1,p2,p3,nrow=1)
ggsave(filename=paste0(file,"SurvivalPlot.pdf"),
       plot=main_surv, width=12, height=5)

####################
## ## 3. Frailty Model with tt
####################
## select the best fit
## tt(powerratio) cannot be included causing frailty model to collapse
lam <- seq(0, 99, len=100)
logL <- aic <- bic <- numeric(length(lam))
for (i in 1:length(lam)) {
  fit <- coxph(update.formula(fit1bf, ~. + tt(issuesalience1) +tt(tradeissuesalience1)
                              + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
               data=tdData, tt=function(x, t, ...) {x*log(t+lam[i])})
  logL[i] <- logLik(fit)
  aic[i] <- AIC(fit)
  bic[i] <- BIC(fit)
}
## now use the best lam to refit the model
lamhat <- lam[which.min(bic)]
best_fitf <- coxph(update.formula(fit1bf, ~. + tt(issuesalience1) +tt(tradeissuesalience1)
                                  + tt(issuesalience2) +tt(tradeissuesalience2)+ tt(contbinary)),
                   data=tdData,tt=function(x, t, ...) {x*log(t+lamhat)})
summary(best_fitf)
cox.zph(best_fitf,transform='rank')

##########
## First Difference of Trade Ratio
## tt without traderatio
###################
par_traderatio <- summary(best_fitf)$coefficients["traderatio",c(1,3)]
par_tradeissuesalience1 <- summary(best_fitf)$coefficients["tradeissuesalience1",c(1,3)]
par_ttradeissuesalience1 <- summary(best_fitf)$coefficients["tt(tradeissuesalience1)",c(1,3)]
par_tradeissuesalience2 <- summary(best_fitf)$coefficients["tradeissuesalience2",c(1,3)]
par_ttradeissuesalience2 <- summary(best_fitf)$coefficients["tt(tradeissuesalience2)",c(1,3)]
#first difference function from Licht2011
hrtrade_issuesalience2 <- function(x, par01, par2, par3) (exp(prop*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience1 <- function(x, par01, par2, par3) (exp(prop*(par01+par2+par3*x))-1)*100
hrtrade_issuesalience0 <- function(x, par01) (exp(prop*(par01))-1)*100
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
tradeimpactplot_f <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.001)+
  geom_smooth(colour="black",aes(linetype=type))+
  xlab("Days")+ylab(expression("First Difference"))+
  theme_bw() + theme(legend.position="bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_linetype_manual(values=c("twodash", "dotted","solid"),
                        name="Salience", 
                        labels = c("Low","Median","High"))

main_fd_f <- tradeimpactplot_f+scale_y_continuous(breaks =c(-5,0,5),limits = c(-10,5))
ggsave(filename=paste0(file,"FDPlot_Frailty.pdf"),
       plot=main_fd_f)

####################
## Latex table
#####################
## clustered models
best_fitc_main_se <- summary(best_fitc_main)$coefficients[,4]
best_fitc_se <- summary(best_fitc)$coefficients[,4]
cmod2_se <- summary(cmod2)$coefficients[,4]
fit1bc_se <- summary(fit1bc)$coefficients[,4]
library(stargazer)
cluster_tex <- stargazer(best_fitc_main, best_fitc, cmod2, fit1bc, single.row = TRUE,
          se=list(best_fitc_main_se,best_fitc_se,cmod2_se,fit1bc_se),
          title="Cox Regression Results", dep.var.labels.include = FALSE,
          dep.var.caption  = "Days before Quitting",
          column.labels   = c("Adjusted Models", "Original Models"),
          column.separate = c(3, 1),
          covariate.labels=c("Cost Ratio","Median Salience", "High Salience",
                             "Cost*MSalience", "Cost*HSalience",
                             "Joint Democracy","Contiguity","Power Ratio",
                             "Defense Ratio", "tt(Cost Ratio)","tt(Median Salience)",
                             "tt(Cost*MSalience)", "tt(High Salience)",
                             "tt(Cost*HSalience)", "tt(Contiguity)","tt(Power Ratio)",
                             "T(Contiguity)","T(Power Ratio)",
                             "T(Median Salience)","T(High Salience)",
                             "Cost*MSalience", "Cost*HSalience",
                             "T(Cost*MSalience)", "T(Cost*HSalience)"
          ),
          add.lines = list(c("Agjusted Method","tt","tt","Time Interaction","No"),
                           c("BIC", round(BIC(best_fitc_main),2),round(BIC(best_fitc),2), round(BIC(cmod2),2),round(BIC(fit1bc),2))),
          keep.stat=c("wald","ll", "n"),
          no.space=TRUE,
          out=paste0(file,"ClusterTable.tex"), out.header = FALSE)



## frailty models
## get around way to trick stargazer
vec.coef <- list()
vec.se <- list()
vec.t <- list()
vec.ci <- list()
vec.coef[[1]] <- best_fitf$coefficients 
vec.se[[1]] <- diag(sqrt(vcov(best_fitf)))
vec.ci[[1]] <- as.matrix(cbind(vec.coef[[1]]-1.96*vec.se[[1]],vec.coef[[1]]+1.96*vec.se[[1]])) 
vec.t[[1]] <- vec.coef[[1]]/vec.se[[1]]

vec.coef[[2]] <- fit1bf$coefficients 
vec.se[[2]] <- diag(sqrt(vcov(fit1bf)))
vec.ci[[2]] <- as.matrix(cbind(vec.coef[[2]]-1.96*vec.se[[2]],vec.coef[[2]]+1.96*vec.se[[2]])) 
vec.t[[2]] <- vec.coef[[2]]/vec.se[[2]]


vec.coef[[3]] <- fit1bf_me$coefficients 
vec.se[[3]] <- diag(sqrt(vcov(fit1bf_me)))
vec.ci[[3]] <- as.matrix(cbind(vec.coef[[3]]-1.96*vec.se[[3]],vec.coef[[3]]+1.96*vec.se[[3]])) 
vec.t[[3]] <- vec.coef[[3]]/vec.se[[3]]

## report yearly fixed effect
vec.names <- c("Cost Ratio","Median Salience", "High Salience",
               "Cost*MSalience", "Cost*HSalience",
               "Joint Democracy","Contiguity","Power Ratio",
               "Defense Ratio", "tt(Median Salience)",
               "tt(Cost*MSalience)", "tt(High Salience)",
               "tt(Cost*HSalience)", "tt(Contiguity)","tt(Power Ratio)"
)
stargazer(best_fitc_main,fit1bc,fit1bc, single.row = FALSE, omit.stat="all", dep.var.caption  = "Dependent Variables",
          dep.var.labels= c("Adjusted","Original","Original (coxme)"), 
          covariate.labels=vec.names,
          coef = vec.coef,
          ci.custom = vec.ci,
          t = vec.t, 
          add.lines = list(c("Observations", summary(best_fitf)$n,summary(fit1bf)$n, fit1bf_me$n[2]),
                           c("AIC",AIC(best_fitf),AIC(fit1bf),AIC(fit1bf_me)),
                           c("BIC",BIC(best_fitf),BIC(fit1bf),BIC(fit1bf_me))
                           ),
          out=paste0(file,"FrailtyTable.tex"), out.header = FALSE)

