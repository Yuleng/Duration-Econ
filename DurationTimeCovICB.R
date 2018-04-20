## This is an Rcode using ICB dataset

setwd("E:/OneDrive2nd/OneDrive - 广厚设计学校")

library(survival)
library(countrycode)
## read icb data from icb 1918-2015
icb <- read.csv("./DATA/icb/icbdy_v12.csv")[,c("crisno","year","statea","stateb")]
icb2 <- read.csv("./DATA/icb/icb2v12.csv")[,c("crisno","cracid","trgterra","resterra","outfor","viol","issue","gravty")]
icb2$statea <- icb2$cracid
## I choose gravity over issue because the former varies within a dyad, the latter does no
## use crisno 445 Russian-Georgia to illustrate
## Russia saw it as 1 limited military threat; Georgia saw it as 3 territorial threat
## the issue variable cannot capture this given it assign both sides as dealing with the same issue
## in this case 1 military-security issue
## first, rearrange icb into directed dyad
## also, eliminate double counting by year
icb <- icb[!duplicated(subset(icb,select=-year)),]
temp <- icb; names(temp)[3:4] <- c("stateb","statea"); temp <- temp[,names(icb)]
icb <- rbind(icb,temp); rm(temp)
## now merge the two datasets
icbdy <- merge(icb, subset(icb2,select=-cracid), by=c("crisno","statea"), all.x=TRUE)
icbdy <- na.omit(icbdy)
## consult the icb dataviewer
## http://www.icb.umd.edu/dataviewer/?crisno=1
## decide if it is voluntary agreement or tacit understanding (1,2,3) and compliance (7), it counts as quitting
## rest as censored
icbdy$cens <- ifelse(icbdy$outfor %in% c(1,2,3,7), 1, 0)# censor as 0
icbdy$ccode1 <- icbdy$statea; icbdy$ccode2 <- icbdy$stateb

## create the basic frame for duration data
dur <- icbdy[,c("ccode1","ccode2","year","trgterra","resterra","cens","viol","gravty","crisno")]
dur$id <- 1:length(dur$crisno)

## creating a dataset to hold covariates from the duration dataset
## use "maxdur"; I don't expect the results to differ
pool <- dur[,c("id","ccode1","ccode2","year","trgterra")]
hold <- data.frame()
for (i in 1:nrow(pool)){
  temp = pool[i,]
  temp = cbind(temp[,c("id","ccode1","ccode2")], 
               dur=temp$trgterra,
               year=temp$year+0:floor(temp$trgterra/365.25),
               day=365.25*(0:floor(temp$trgterra/365.25)))
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
## four -1 values for defense2
## three of them concerns Turkey and Cyprus (Treaty of Alliance (1960),concerning the independence of Cyprus from UK )
## the correlates of war dataset code the dyad inconsistently
## Turkey-Cyprus coded as 1 defense pack; yet Cyprus-Turkey as 0 defense pack
## For my purpose, this outside defense pack can be coded as 0 (for Cyprus)
## the last case concerns Russia and Lithuania 1939
## again, due to inconsistent coding
## Russian-L has 1 defense pack; yet L-R has 0 defense pack
## again, code the latter as 0 
hold$defense2[hold$defense2==-1] <- 0


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
hold$powerratio <- log(hold$cinc1/hold$cinc2+hold$ccode1)

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
dur$dur <- dur$trgterra; dur0<-dur
dur <- dur[,c("id","dur","cens","crisno","gravty","viol")] 
dur <- tmerge(dur, dur, id=id, quit=event(dur,cens)) # set the range
durB <- tmerge(dur, holdB, id=id, 
               tradedepend1=tdc(day,tradedepend1), tradedepend2=tdc(day,tradedepend2),
               tradeshare1=tdc(day,tradeshare1), tradeshare2=tdc(day,tradeshare2),
               tradedepend1.net=tdc(day,tradedepend1.net), tradedepend2.net=tdc(day,tradedepend2.net),
               tradeshare1.net=tdc(day,tradeshare1.net), tradeshare2.net=tdc(day,tradeshare2.net),
               demo1=tdc(day,demo1), demo2=tdc(day,demo2),
               major1=tdc(day,major1), major2=tdc(day,major2),
               cinc1=tdc(day,cinc1), cinc2=tdc(day,cinc2),
               mindist=tdc(day,mindist), conttype=tdc(day,conttype), contbinary=tdc(day,contbinary),
               defense1=tdc(day,defense1), defense2=tdc(day, defense2),
               powerratio=tdc(day,powerratio))
durG <- tmerge(dur, holdG, id=id, 
               tradedepend1=tdc(day,tradedepend1), tradedepend2=tdc(day,tradedepend2),
               tradeshare1=tdc(day,tradeshare1), tradeshare2=tdc(day,tradeshare2),
               tradedepend1.net=tdc(day,tradedepend1.net), tradedepend2.net=tdc(day,tradedepend2.net),
               tradeshare1.net=tdc(day,tradeshare1.net), tradeshare2.net=tdc(day,tradeshare2.net),
               demo1=tdc(day,demo1), demo2=tdc(day,demo2),
               major1=tdc(day,major1), major2=tdc(day,major2),
               cinc1=tdc(day,cinc1), cinc2=tdc(day,cinc2),
               mindist=tdc(day,mindist), conttype=tdc(day,conttype), contbinary=tdc(day,contbinary),
               defense1=tdc(day,defense1), defense2=tdc(day, defense2),
               powerratio=tdc(day,powerratio))

## save the data
save(durB,durG, file="./GIT/Duration-Econ/DurationTimeCov_ICB.RData")

