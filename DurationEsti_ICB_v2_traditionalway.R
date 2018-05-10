## This is an R code using traditional interaction with log time solution

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

library(gridExtra)
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
tdData$traderatio.net <- tdData$tradeshare1.net/(tdData$tradeshare2.net+tdData$tradeshare1.net)
tdData$defenseratio <- tdData$defense1/(1+tdData$defense2) #because we have a large number of 0 defense pacts, n
## tdData$defenseratio <- tdData$defense1/(tdData$defense1+tdData$defense2);tdData$traderatio <- tdData$tradeshare1/(tdData$tradeshare2+tdData$tradeshare1) 
tdData$tradeissuesalience1 <- tdData$traderatio*tdData$issuesalience1
tdData$tradeissuesalience2 <- tdData$traderatio*tdData$issuesalience2
tdData$type <- as.factor(tdData$major1*10+tdData$major2)

###############################
## Models
################################
## justify not using frailty because it is not about different groups having 
## different intercepts
## But cluster, because it is about the correlation within a dyad
## In medical studies, like twin studies or blindness on both eyes
## Though there is no clear guidance on which model to use
## I may check the model fit though

tdData$time <- tdData$tstop ## using fox

fmod <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+tradeissuesalience1+tradeissuesalience2+joint_demo+contbinary+powerratio+defenseratio+frailty(crisno), data=tdData);summary(fmod)
cox.zph(fmod,transform='rank')

fmod1 <- coxph(update.formula(fmod, ~. +contbinary:log(time)+tradeissuesalience1:log(time)+tradeissuesalience2:log(time)), data=tdData);summary(fmod1)
cox.zph(fmod1,transform='rank')


cmod <- coxph(Surv(tstart, tstop, quit) ~ traderatio+issuesalience1+issuesalience2+traderatio:issuesalience1+traderatio:issuesalience2+joint_demo+contbinary+powerratio+defenseratio+cluster(crisno), data=tdData);summary(cmod)
cox.zph(cmod,transform='rank')

cmod1 <- coxph(update.formula(cmod, ~. +contbinary:log(time)+issuesalience1:log(time)+issuesalience2:log(time)+
                 traderatio:issuesalience1:log(time)+traderatio:issuesalience2:log(time)), data=tdData);summary(cmod1)
cox.zph(cmod1,transform='rank')

## Get around way to deal with time interaction for plotting survival curve
tdData$ttraderatio <- tdData$traderatio*log(tdData$time)
tdData$tcontbinary <- tdData$contbinary*log(tdData$time)
tdData$tissuesalience1 <- tdData$issuesalience1*log(tdData$time)
tdData$tissuesalience2 <- tdData$issuesalience2*log(tdData$time)
cmod2 <- coxph(update.formula(cmod, ~.+tcontbinary+tissuesalience1+tissuesalience2+
                                traderatio:tissuesalience1+traderatio:tissuesalience2), data=tdData);summary(cmod2)
cox.zph(cmod2,transform='rank')

fmod2 <- coxph(update.formula(fmod, ~.+tcontbinary+tissuesalience1+tissuesalience2+
                                traderatio:tissuesalience1+traderatio:tissuesalience2), data=tdData);summary(fmod2)


########################
## Hazard Ratio and Survival curve
#############################

## Impact of issuesalience1
## Plot the impact of traderatio
pdata <- data.frame(traderatio=c(0.2,0.9),
                   powerratio=median(tdData$powerratio,na.rm=TRUE),
                   defenseratio=median(tdData$defenseratio,na.rm=TRUE),
                   time=365,#one year
                   crisno=median(tdData$crisno,na.rm=TRUE),
                   joint_demo=1,contbinary=1,
                   issuesalience1=0,
                   issuesalience2=0
                   )
#ggadjustedcurves(cmod,data=pdata,variable="traderatio",method="marginal",xlim=c(0,300))
pdata2 <- pdata
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p1 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(a) Low Salience", method="marginal",legend.title="Trade Ratio")+
      annotate("text", x = 600, y = .7, label = "Low Cost")

pdata <- data.frame(traderatio=c(0.2,0.9),
                    powerratio=median(tdData$powerratio,na.rm=TRUE),
                    defenseratio=median(tdData$defenseratio,na.rm=TRUE),
                    time=365,#one year
                    crisno=median(tdData$crisno,na.rm=TRUE),
                    joint_demo=1,contbinary=1,
                    issuesalience1=1,
                    issuesalience2=0
)
#ggadjustedcurves(cmod,data=pdata,variable="traderatio",method="marginal",xlim=c(0,300))
pdata2 <- pdata
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p2 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(b) Median Salience", method="marginal",legend.title="Trade Ratio")+
      annotate("text", x = 600, y = .7, label = "High Cost")

pdata <- data.frame(traderatio=c(0.2,0.9),
                    powerratio=median(tdData$powerratio,na.rm=TRUE),
                    defenseratio=median(tdData$defenseratio,na.rm=TRUE),
                    time=365,#one year
                    crisno=median(tdData$crisno,na.rm=TRUE),
                    joint_demo=1,contbinary=1,
                    issuesalience1=0,
                    issuesalience2=1
)
#ggadjustedcurves(cmod,data=pdata,variable="traderatio",method="marginal",xlim=c(0,300))
pdata2 <- pdata
pdata2$tcontbinary <- pdata2$contbinary*log(pdata2$time)
pdata2$tissuesalience1 <- pdata2$issuesalience1*log(pdata2$time)
pdata2$tissuesalience2 <- pdata2$issuesalience2*log(pdata2$time)
p3 <- ggadjustedcurves(cmod2,data=pdata2,variable="traderatio",
                       xlab=FALSE, ylab="", font.tickslab=10,
                       main="(c) High Salience", method="marginal",legend.title="Trade Ratio")+
      annotate("text", x = 600, y = .7, label = "High Cost")


grid_arrange_shared_legend(p1,p2,p3,nrow=1)

############################
## Plot the first difference
###########################
par_traderatio <- summary(cmod2)$coefficients["ttraderatio",c(1,4)]
par_tradeissuesalience1 <- summary(cmod2)$coefficients["issuesalience1:traderatio",c(1,4)]
par_ttradeissuesalience1 <- summary(cmod2)$coefficients["tissuesalience1:traderatio",c(1,4)]
par_tradeissuesalience2 <- summary(cmod2)$coefficients["issuesalience2:traderatio",c(1,4)]
par_ttradeissuesalience2 <- summary(cmod2)$coefficients["tissuesalience2:traderatio",c(1,4)]
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
tradeimpactplot <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.01)+
  geom_smooth(colour="black")+ylim(-100,250) +
  xlab("Days")+ylab(expression("First Difference"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", 750, 120, label="Low Salience")+ 
  annotate("text",500,0,label="Median Salience")+
  annotate("text",100,-75,label="High Salience")

## Zoom in plot and incorporate the impact of issuesalience0 type 
zoominplot <- ggplot(plot_dat, aes(time, impact, group=type))+
  geom_point(colour="grey",alpha=0.01)+
  geom_smooth(colour="black")+ylim(-100,250) +xlim(0,730)+
  xlab("Days")+ylab(expression("First Difference"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", 600, 120, label="Low Salience")+ 
  annotate("text",400,0,label="Median Salience")+
  annotate("text",100,-75,label="High Salience")
