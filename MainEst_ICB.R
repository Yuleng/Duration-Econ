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

tdData$joint_demo <- ifelse(tdData$demo1>5 & tdData$demo2>5, 1, 0)
tdData$traderatio <- tdData$tradeshare1/(tdData$tradeshare2+tdData$tradeshare1+1) ## not like power ratio, because the latter produce too many NA values
tdData$issuesalience <- factor(ifelse(tdData$gravty %in% c(5,6),2,ifelse(tdData$gravty %in% c(2,3,4),1,0)), levels=c(0,1,2))
temp <- sapply(levels(tdData$issuesalience), function(x) as.integer(x == tdData$issuesalience))
colnames(temp) <- paste0("issuesalience",colnames(temp))
tdData <- cbind(tdData,temp)
tdData$traderatio.net <- tdData$tradeshare1.net/(tdData$tradeshare2.net+tdData$tradeshare1.net+1)
tdData$defenseratio <- tdData$defense1/(1+tdData$defense2+tdData$defense1)
tdData$tradeissuesalience1 <- tdData$traderatio*tdData$issuesalience1
tdData$tradeissuesalience2 <- tdData$traderatio*tdData$issuesalience2

##############
## Basic Models and Tests
##############
fit1b <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+tradeissuesalience1+tradeissuesalience2+joint_demo+contbinary+powerratio+defenseratio, data=tdData);summary(fit1b)
fit1bc <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+tradeissuesalience1+tradeissuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(fit1bc)
fit1bf<- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+tradeissuesalience1+tradeissuesalience2+joint_demo+contbinary+powerratio+defenseratio+frailty(crisno), data=tdData);summary(fit1bf)

cox.zph(fit1b,transform='rank');cox.zph(fit1bc,transform='rank');cox.zph(fit1bf,transform='rank')

#############
## Graph the level of censoring and outliers
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
multiplot(cen_plot1,cen_plot2, cols=2)

## outliers
## indentity graph shows the impact of outliers
ggcoxzph(cox.zph(fit1bc,transform='identity'), point.col="grey")[2] 
ggcoxzph(cox.zph(fit1bc,transform='identity'), point.col="grey")[3] 

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
main_fd_c <- tradeimpactplot_c_main+scale_y_continuous(breaks =c(-10,0,10),limits = c(-20,10))
ggsave(filename="C:/Users/YULENG/OneDrive/Documents/RESEARCH/Manuscript/Duration/ISQ/FDPlot_Main.pdf",
       plot=main_fd_c)


##############
## with tt for traderatio
##############
## select the best model
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
## By BIC; others give similar results
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
## The other main graph to show
#########################
main_fd <- tradeimpactplot_c+scale_y_continuous(breaks =c(-20,0,20),limits = c(-40,40))
ggsave(filename="C:/Users/YULENG/OneDrive/Documents/RESEARCH/Manuscript/Duration/ISQ/FDPlot.pdf",
       plot=main_fd)

########################
## 2. Proxy clustered model and survival plot
########################
tdData$time <- tdData$tstop ## using fox
cmod <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(cmod)
cox.zph(cmod,transform='rank')
## Get around way to deal with time interaction for plotting survival curve
tdData$ttraderatio <- tdData$traderatio*log(tdData$time)
tdData$tcontbinary <- tdData$contbinary*log(tdData$time)
tdData$tissuesalience1 <- tdData$issuesalience1*log(tdData$time)
tdData$tissuesalience2 <- tdData$issuesalience2*log(tdData$time)
cmod2 <- coxph(update.formula(cmod, ~.+tcontbinary+tissuesalience1+tissuesalience2+
                                traderatio:tissuesalience1+traderatio:tissuesalience2), data=tdData);summary(cmod2)
summary(cmod2)$coefficients[,c(1,6)]
summary(best_fitc_main)$coefficients[,c(1,6)]

###############
## Plot the proxy survival rate
###############
## Impact of issuesalience and traderatio
summary(tdData$traderatio)
pdata <- data.frame(traderatio=c(.01,.1),
                    powerratio=median(tdData$powerratio,na.rm=TRUE),
                    defenseratio=median(tdData$defenseratio,na.rm=TRUE),
                    time=365,#one year
                    crisno=median(tdData$crisno,na.rm=TRUE),
                    joint_demo=1,contbinary=1)
## low salience
pdata2 <- pdata
pdata2$issuesalience1=0
pdata2$issuesalience2=0
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p1 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(a) Low Salience", method="marginal",legend.title="Cost Ratio")+
  annotate("text", x = 600, y = .7, label = "Low Cost")

## Median salience
pdata2 <- pdata
pdata2$issuesalience1=1
pdata2$issuesalience2=0
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p2 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(b) Median Salience", method="marginal",legend.title="Cost Ratio")+
  annotate("text", x = 600, y = .7, label = " ")

## High salience
pdata2 <- pdata
pdata2$issuesalience1=0
pdata2$issuesalience2=1
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p3 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(c) High Salience", method="marginal",legend.title="Cost Ratio")+
  annotate("text", x = 600, y = .7, label = "High Cost")


main_surv <- grid_arrange_shared_legend(p1,p2,p3,nrow=1)
ggsave(filename="C:/Users/YULENG/OneDrive/Documents/RESEARCH/Manuscript/Duration/ISQ/SurvivalPlot.pdf",
       plot=main_surv, width=12, height=5)

####################
## ## 3. Frailty Model with tt
####################
## select the best fit
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
tradeimpactplot_f <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.001)+
  geom_smooth(colour="black",aes(linetype=type))+
  xlab("Days")+ylab(expression("First Difference"))+
  theme_bw() + theme(legend.position="bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_linetype_manual(values=c("twodash", "dotted","solid"),
                        name="Salience", 
                        labels = c("Low","Median","High"))

main_fd_f <- tradeimpactplot_f+scale_y_continuous(breaks =c(-10,0,10),limits = c(-20,10))
ggsave(filename="C:/Users/YULENG/OneDrive/Documents/RESEARCH/Manuscript/Duration/ISQ/FDPlot_Frailty.pdf",
       plot=main_fd_f)