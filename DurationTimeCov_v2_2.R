# This is an R code to transform the duration data into having time dependent covariates
# This is a revised version on v_1

setwd("E:/OneDrive2nd/OneDrive - 广厚设计学校")

library(survival)
## read mid data from gml 1816-2010
mid <- read.csv("./DATA/COW/MID/GML-MID/gml-ndy.csv")[,c("dispnum3","ccode1","ccode2","year","mindur","maxdur","outcome","orig1","orig2","fatalpre1","fatalpre2")]
mid <- na.omit(mid)## deleting na values
## focus on counting quitters
## In MID outcome 3 yield by A, 4 yield by B, 6 both side yield
## Released 7 can be counted as yield
## check Jones1996 "whenever the seizure of material or personnel culminates with 
## their release from captivity"
## As such, I exclude 8 Unclear and -9 Missing
mid <- mid[mid$outcome %in% c(1:7,9),] 
mid$outcome[mid$outcome==7&mid$orig1==1&mid$orig2==0] <- 7.1 # side 1 as the quitter
mid$outcome[mid$outcome==7&mid$orig1==0&mid$orig2==1] <- 7.2 # side 2 as the quitter
mid$outcome[mid$outcome==7&mid$orig1==1&mid$orig2==1] <- 7.3 # side 2 as the quitter
## to avoid double counting mids lasting longer than 1 year
mid <- mid[!duplicated(mid[,1:3]),]; mid <- subset(mid,select=-dispnum3)
## Now add causualty NA values
mid$fatalpre1[mid$fatalpre1==-9] <- NA; mid$fatalpre2[mid$fatalpre2==-9] <- NA
## Cannot find annual battle death data; resort to using a time-invariant depend variable
## Fatality should be weighted by army size
## USE the cinc data
## National Power v5.0 1816-2012
milipersonnel <- read.csv("./DATA/COW/NMC_v4_0/NMC_5_0.csv",na.strings=-9)[,c("ccode","year","milper")]
milipersonnel$ccode1 <- milipersonnel$ccode2 <- milipersonnel$ccode;
milipersonnel$milper1 <- milipersonnel$milper2 <- milipersonnel$milper
mid <- merge(mid, milipersonnel[,c("ccode1", "year", "milper1")], by=c("ccode1","year"), all.x=TRUE)
mid <- merge(mid, milipersonnel[,c("ccode2", "year", "milper2")], by=c("ccode2","year"), all.x=TRUE)
mid$trooploss1 <- mid$fatalpre1/(mid$milper1+1)
mid$trooploss2 <- mid$fatalpre1/(mid$milper2+1)
library(dplyr)
mid <- mid %>% select(ccode1, everything())

## Krustev only counts mid duration; however, the theory predicts on states rationally choosing a time to quit
temp1 <- mid
temp1$cens <- ifelse(temp1$outcome %in% c(3,6,7.1,7.3), 1, 0)# censor as 0
temp2 <- mid; names(temp2)[1:2] <- c("ccode2","ccode1"); temp2 <- temp2[,names(mid)]
temp2$cens <- ifelse(temp1$outcome %in% c(4,6,7.2,7.3), 1, 0)
temp <- rbind(temp1,temp2)

## create the basic frame for duration data
dur <- subset(temp, select=-c(outcome,orig1,orig2));rm(temp1,temp2,temp)
dur$id <- paste(dur$ccode1,dur$ccode2,dur$year)

## creating a dataset to hold covariates from the duration dataset
## use "maxdur"; I don't expect the results to differ
pool <- dur[,c("id","ccode1","ccode2","year","maxdur")]
hold <- data.frame()
for (i in 1:nrow(pool)){
  temp = pool[i,]
  temp = cbind(temp[,c("id","ccode1","ccode2")], 
               dur=temp$maxdur,
               year=temp$year+0:floor(temp$maxdur/365.25),
               day=365.25*(0:floor(temp$maxdur/365.25)))
  hold = rbind(hold,temp); rm(temp)
}
## NOW merge all ind variables
## Read the controls
load("./DATA/COW/COWPolityControls2017.RData") # need to specify the working directory first, chinese character doesn't work for reading RData

##########################
## democracy
##########################
## polity 1800-2015
polity$ccode2 <- polity$ccode1 <- polity$ccode; polity$demo1 <- polity$demo2 <- polity$polity2
## MERGE with hold
hold <- merge(hold, polity[,c("ccode1","year","demo1")], by=c("ccode1","year"), all.x=TRUE)
hold <- merge(hold, polity[,c("ccode2","year","demo2")], by=c("ccode2","year"), all.x=TRUE)

##########################
## inter-state relation
##########################
## alliance 1816-2012
## outside defense packs
hold <- merge(hold, allydyad, by=c("ccode1","ccode2","year"), all.x=TRUE)
allystate$ccode2 <- allystate$ccode1; allystate$defense1 <- allystate$defense2 <- allystate$defense
hold <- merge(hold, allystate[,c("ccode1","year","defense1")], by=c("ccode1","year"), all.x=TRUE)
hold <- merge(hold, allystate[,c("ccode2","year","defense2")], by=c("ccode2","year"), all.x=TRUE)
hold$defense[is.na(hold$defense)] <- 0 
hold$defense1[is.na(hold$defense1)] <- 0 
hold$defense2[is.na(hold$defense2)] <- 0 
hold$defense1 <- hold$defense1-hold$defense
hold$defense2 <- hold$defense2-hold$defense


#######################
## distance
#######################
## contiguity 1816-2016
hold <- merge(hold, contiguity[,c("ccode1","ccode2","year","conttype")], by=c("ccode1","ccode2","year"), all.x=TRUE)
hold[hold$year>1815 & hold$year<2017 & is.na(hold$conttype),]$conttype <- 6 # following Slantchev code noncontiguous as 6
hold$contbinary <- ifelse(hold$conttype==6, 0, 1)

## distance 1946-2015
hold <- merge(hold, capdist, by=c("ccode1","ccode2","year"), all.x=TRUE)
hold <- merge(hold, centdist, by=c("ccode1","ccode2","year"), all.x=TRUE)
hold <- merge(hold, mindist, by=c("ccode1","ccode2","year"), all.x=TRUE)

#########################
## power
########################
## national capability 1816-2012
hold <- merge(hold, cinc[,c("ccode1","year","cinc1")], by=c("ccode1","year"), all.x=TRUE)
hold <- merge(hold, cinc[,c("ccode2","year","cinc2")], by=c("ccode2","year"), all.x=TRUE)

library(plyr)
## curate the data a bit: calculate powerratio
hold$powerratio <- log(hold$cinc1/hold$cinc2+1)

## major power status
library(countrycode)
hold$major1 <- ifelse((hold$ccode1==2 & hold$year >1897 & hold$year < 2017) | 
                        (hold$ccode1==200 & hold$year >1815 & hold$year < 2017) |  
                        (hold$ccode1==220 & hold$year >1815 & hold$year < 1941) |
                        (hold$ccode1==220 & hold$year >1944 & hold$year < 2017) |
                        (hold$ccode1==255 & hold$year >1815 & hold$year < 1919) |
                        (hold$ccode1==255 & hold$year >1924 & hold$year < 1946) |
                        (hold$ccode1==255 & hold$year >1990 & hold$year < 2017) |
                        (hold$ccode1==300 & hold$year >1815 & hold$year < 1919) |
                        (hold$ccode1==325 & hold$year >1859 & hold$year < 1944) |
                        (hold$ccode1==365 & hold$year >1815 & hold$year < 1918) |
                        (hold$ccode1==365 & hold$year >1921 & hold$year < 2017) |
                        (hold$ccode1==710 & hold$year >1949 & hold$year < 2017) |
                        (hold$ccode1==740 & hold$year >1894 & hold$year < 1946) |
                        (hold$ccode1==740 & hold$year >1990 & hold$year < 2017),
                      1, 0)
hold$major2 <- ifelse((hold$ccode2==2 & hold$year >1897 & hold$year < 2017) | 
                        (hold$ccode2==200 & hold$year >1815 & hold$year < 2017) |  
                        (hold$ccode2==220 & hold$year >1815 & hold$year < 1941) |
                        (hold$ccode2==220 & hold$year >1944 & hold$year < 2017) |
                        (hold$ccode2==255 & hold$year >1815 & hold$year < 1919) |
                        (hold$ccode2==255 & hold$year >1924 & hold$year < 1946) |
                        (hold$ccode2==255 & hold$year >1990 & hold$year < 2017) |
                        (hold$ccode2==300 & hold$year >1815 & hold$year < 1919) |
                        (hold$ccode2==325 & hold$year >1859 & hold$year < 1944) |
                        (hold$ccode2==365 & hold$year >1815 & hold$year < 1918) |
                        (hold$ccode2==365 & hold$year >1921 & hold$year < 2017) |
                        (hold$ccode2==710 & hold$year >1949 & hold$year < 2017) |
                        (hold$ccode2==740 & hold$year >1894 & hold$year < 1946) |
                        (hold$ccode2==740 & hold$year >1990 & hold$year < 2017),
                      1, 0)


## Read trade data
load("./DATA/TRADE/COW/trade2017.RData")
load("./DATA/TRADE/tradeNet_constUSd.RData")


holdB <- merge(hold, tradeB[,c("ccode1","ccode2","year","tradeshare1","tradeshare2","tradedepend1","tradedepend2")],by=c("ccode1","ccode2","year"), all.x=TRUE)
holdG <- merge(hold, tradeG[,c("ccode1","ccode2","year","tradeshare1","tradeshare2","tradedepend1","tradedepend2")],by=c("ccode1","ccode2","year"), all.x=TRUE)
rm(tradeB, tradeG) # remove the raw data

trade.closeness.constusd$ccode1 <- trade.closeness.constusd$ccode2 <- trade.closeness.constusd$node
trade.closeness.constusd$integration1 <- trade.closeness.constusd$integration2 <- trade.closeness.constusd$total
holdB <- merge(holdB, trade.closeness.constusd[,c("ccode1","year","integration1")],by=c("ccode1","year"), all.x=TRUE)
holdB <- merge(holdB, trade.closeness.constusd[,c("ccode2","year","integration2")],by=c("ccode2","year"), all.x=TRUE)
holdG <- merge(holdG, trade.closeness.constusd[,c("ccode1","year","integration1")],by=c("ccode1","year"), all.x=TRUE)
holdG <- merge(holdG, trade.closeness.constusd[,c("ccode2","year","integration2")],by=c("ccode2","year"), all.x=TRUE)

holdB$tradedepend1.net <- holdB$tradedepend1*exp(-holdB$integration1)
holdB$tradedepend2.net <- holdB$tradedepend2*exp(-holdB$integration2)
holdB$tradeshare1.net <- holdB$tradeshare1*exp(-holdB$integration1)
holdB$tradeshare2.net <- holdB$tradeshare2*exp(-holdB$integration2)

holdG$tradedepend1.net <- holdG$tradedepend1*exp(-holdG$integration1)
holdG$tradedepend2.net <- holdG$tradedepend2*exp(-holdG$integration2)
holdG$tradeshare1.net <- holdG$tradeshare1*exp(-holdG$integration1)
holdG$tradeshare2.net <- holdG$tradeshare2*exp(-holdG$integration2)


## now assemble the data into having time varying covariates
dur$dur <- dur$maxdur
dur <- dur[,c("id","dur","cens","trooploss1","trooploss2")] 
dur <- tmerge(dur, dur, id=id, quit=event(dur,cens)) # set the range
durB <- tmerge(dur, holdB, id=id, 
               tradedepend1=tdc(day,tradedepend1), tradedepend2=tdc(day,tradedepend2),
               tradeshare1=tdc(day,tradeshare1), tradeshare2=tdc(day,tradeshare2),
               tradedepend1.net=tdc(day,tradedepend1.net), tradedepend2=tdc(day,tradedepend2),
               tradeshare1.net=tdc(day,tradeshare1.net), tradeshare2=tdc(day,tradeshare2),
               demo1=tdc(day,demo1), demo2=tdc(day,demo2),
               major1=tdc(day,major1), major2=tdc(day,major2),
               cinc1=tdc(day,cinc1), cinc2=tdc(day,cinc2),
               mindist=tdc(day,mindist), conttype=tdc(day,conttype), contbinary=tdc(day,contbinary),
               defense1=tdc(day,defense1), defense2=tdc(day, defense2),
               powerratio=tdc(day,powerratio))
durG <- tmerge(dur, holdG, id=id, 
               tradedepend1=tdc(day,tradedepend1), tradedepend2=tdc(day,tradedepend2),
               tradeshare1=tdc(day,tradeshare1), tradeshare2=tdc(day,tradeshare2),
               tradedepend1.net=tdc(day,tradedepend1.net), tradedepend2=tdc(day,tradedepend2),
               tradeshare1.net=tdc(day,tradeshare1.net), tradeshare2=tdc(day,tradeshare2),
               demo1=tdc(day,demo1), demo2=tdc(day,demo2),
               major1=tdc(day,major1), major2=tdc(day,major2),
               cinc1=tdc(day,cinc1), cinc2=tdc(day,cinc2),
               mindist=tdc(day,mindist), conttype=tdc(day,conttype), contbinary=tdc(day,contbinary),
               defense1=tdc(day,defense1), defense2=tdc(day, defense2),
               powerratio=tdc(day,powerratio))

## save the data
save(durB,durG, file="./GIT/Duration-Econ/DurationTimeCov_v2.RData")

