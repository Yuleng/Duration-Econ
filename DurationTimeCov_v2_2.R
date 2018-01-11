# This is an R code to transform the duration data into having time dependent covariates
# This is a revised version on v_1

setwd("E:/OneDrive2nd/OneDrive - 广厚设计学校/GIT/Duration-Econ")

library(survival)
## read mid data from gml 1816-2010
mid <- read.csv("E:/OneDrive2nd/OneDrive - 广厚设计学校/DATA/COW/MID/GML-MID/gml-ndy.csv")[,c("dispnum3","ccode1","ccode2","year","mindur","maxdur","outcome")]
mid <- na.omit(mid)## deleting na values
## focus on counting quitters
mid <- mid[mid$outcome %in% c(3,4,6),] 
## to avoid double counting mids lasting longer than 1 year
mid <- mid[!duplicated(mid[,1:3]),]; mid <- subset(mid,select=-dispnum3)


## Krustev only counts mid duration; however, the theory predicts on states rationally choosing a time to quit
temp1 <- mid
temp1$cens <- ifelse(temp1$outcome %in% c(3,6), 1, 0)
temp2 <- mid; names(temp2)[1:2] <- c("ccode2","ccode1"); temp2 <- temp2[,names(mid)]
temp2$cens <- ifelse(temp1$outcome %in% c(4,6), 1, 0)
temp <- rbind(temp1,temp2)

## create the basic frame for duration data
dur <- subset(temp, select=-outcome);rm(temp1,temp2,temp)
library(countrycode)
dur$ccode1 <- countrycode(dur$ccode1, "cown", "country.name")
dur$ccode2 <- countrycode(dur$ccode2, "cown", "country.name")
dur$id <- paste(dur$ccode1,dur$ccode2,dur$year)

## creating a dataset to hold covariates from the duration dataset
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
load("D:/DATA/COW/COWPolityControls2017.RData")

##########################
## democracy
##########################
## polity 1800-2015
polity$ccode2 <- polity$ccode1; polity$demo1 <- polity$polity2; polity$demo2 <- polity$polity2
## MERGE with hold
hold <- merge(hold, polity[,c("ccode1","year","demo1")], by=c("ccode1","year"), all.x=TRUE)
hold <- merge(hold, polity[,c("ccode2","year","demo2")], by=c("ccode2","year"), all.x=TRUE)

##########################
## inter-state relation
##########################
## alliance 1816-2012
hold <- merge(hold, allydyad, by=c("ccode1","ccode2","year"), all.x=TRUE)

## igo 1820-2005
hold <- merge(hold, igodyad, by=c("ccode1","ccode2","year"), all.x=TRUE)

## affinity 1946-2008
library(foreign)
affin <- read.dta("D:/DATA/POLICY PREFERENCE/affinity_01242010.dta")[,c("s3un4608","year","ccodea","ccodeb")]
affin$ccode1 <- countrycode(affin$ccodea,"cown","country.name")
affin$ccode2 <- countrycode(affin$ccodeb,"cown","country.name")
affin$affinity <-affin$s3un4608
affin <- subset(affin, select=c(ccode1, ccode2, year, affinity))
hold <- merge(hold, affin, by=c("ccode1","ccode2","year"), all.x=TRUE)

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
cinc$ccode2 <- cinc$ccode1; cinc$cinc1 <- cinc$cinc; cinc$cinc2 <- cinc$cinc
hold <- merge(hold, cinc[,c("ccode1","year","cinc1")], by=c("ccode1","year"), all.x=TRUE)
hold <- merge(hold, cinc[,c("ccode2","year","cinc2")], by=c("ccode2","year"), all.x=TRUE)

library(plyr)
## curate the data a bit: adding alliance total; calculate powerratio
hold$alliance <- rowSums(hold[,c("defense","nonaggression","neutrality")],na.rm=TRUE) # excluding entente because it is of much less significance
hold$powerratio <- log(hold$cinc1/hold$cinc2+1)

## major power status
library(countrycode)
hold$major1 <- ifelse((hold$ccode1==countrycode(2,"cown","country.name") & hold$year >1897 & hold$year < 2017) | 
                        (hold$ccode1==countrycode(200,"cown","country.name") & hold$year >1815 & hold$year < 2017) |  
                        (hold$ccode1==countrycode(220,"cown","country.name") & hold$year >1815 & hold$year < 1941) |
                        (hold$ccode1==countrycode(220,"cown","country.name") & hold$year >1944 & hold$year < 2017) |
                        (hold$ccode1==countrycode(255,"cown","country.name") & hold$year >1815 & hold$year < 1919) |
                        (hold$ccode1==countrycode(255,"cown","country.name") & hold$year >1924 & hold$year < 1946) |
                        (hold$ccode1==countrycode(255,"cown","country.name") & hold$year >1990 & hold$year < 2017) |
                        (hold$ccode1==countrycode(300,"cown","country.name") & hold$year >1815 & hold$year < 1919) |
                        (hold$ccode1==countrycode(325,"cown","country.name") & hold$year >1859 & hold$year < 1944) |
                        (hold$ccode1==countrycode(365,"cown","country.name") & hold$year >1815 & hold$year < 1918) |
                        (hold$ccode1==countrycode(365,"cown","country.name") & hold$year >1921 & hold$year < 2017) |
                        (hold$ccode1==countrycode(710,"cown","country.name") & hold$year >1949 & hold$year < 2017) |
                        (hold$ccode1==countrycode(740,"cown","country.name") & hold$year >1894 & hold$year < 1946) |
                        (hold$ccode1==countrycode(740,"cown","country.name") & hold$year >1990 & hold$year < 2017),
                      1, 0)
hold$major2 <- ifelse((hold$ccode2==countrycode(2,"cown","country.name") & hold$year >1897 & hold$year < 2017) | 
                        (hold$ccode2==countrycode(200,"cown","country.name") & hold$year >1815 & hold$year < 2017) |  
                        (hold$ccode2==countrycode(220,"cown","country.name") & hold$year >1815 & hold$year < 1941) |
                        (hold$ccode2==countrycode(220,"cown","country.name") & hold$year >1944 & hold$year < 2017) |
                        (hold$ccode2==countrycode(255,"cown","country.name") & hold$year >1815 & hold$year < 1919) |
                        (hold$ccode2==countrycode(255,"cown","country.name") & hold$year >1924 & hold$year < 1946) |
                        (hold$ccode2==countrycode(255,"cown","country.name") & hold$year >1990 & hold$year < 2017) |
                        (hold$ccode2==countrycode(300,"cown","country.name") & hold$year >1815 & hold$year < 1919) |
                        (hold$ccode2==countrycode(325,"cown","country.name") & hold$year >1859 & hold$year < 1944) |
                        (hold$ccode2==countrycode(365,"cown","country.name") & hold$year >1815 & hold$year < 1918) |
                        (hold$ccode2==countrycode(365,"cown","country.name") & hold$year >1921 & hold$year < 2017) |
                        (hold$ccode2==countrycode(710,"cown","country.name") & hold$year >1949 & hold$year < 2017) |
                        (hold$ccode2==countrycode(740,"cown","country.name") & hold$year >1894 & hold$year < 1946) |
                        (hold$ccode2==countrycode(740,"cown","country.name") & hold$year >1990 & hold$year < 2017),
                      1, 0)


## Read trade data
load("D:/DATA/TRADE/COW/trade.RData")
holdB <- merge(hold, tradeB2014,by=c("ccode1","ccode2","year"), all.x=TRUE)
holdG <- merge(hold, tradeG,by=c("ccode1","ccode2","year"), all.x=TRUE)
rm(tradeB2014, tradeG) # remove the raw data

## now assemble the data into having time varying covariates
dur$dur <- dur$maxdur
dur <- dur[,c("id","dur","cens")] ## if doing competing risks may change this to status?
dur <- tmerge(dur, dur, id=id, quit=event(dur,cens)) # set the range
durB <- tmerge(dur, holdB, id=id, 
               tradedepend1=tdc(day,tradedepend1), tradedepend2=tdc(day,tradedepend2),
               tradeshare1=tdc(day,tradeshare1), tradeshare2=tdc(day,tradeshare2),
               tradeopen1=tdc(day,tradeopen1), tradeopen2=tdc(day,tradeopen2),
               demo1=tdc(day,demo1), demo2=tdc(day,demo2),
               major1=tdc(day,major1), major2=tdc(day,major2),
               cinc1=tdc(day,cinc1), cinc2=tdc(day,cinc2),
               mindist=tdc(day,mindist), conttype=tdc(day,conttype), contbinary=tdc(day,contbinary),
               alliance=tdc(day,alliance),igo=tdc(day, jointmemtotal), affinity=tdc(day,affinity),
               powerratio=tdc(day,powerratio))
durG <- tmerge(dur, holdG, id=id, 
               tradedepend1=tdc(day,tradedepend1), tradedepend2=tdc(day,tradedepend2),
               tradeshare1=tdc(day,tradeshare1), tradeshare2=tdc(day,tradeshare2),
               tradeopen1=tdc(day,tradeopen1), tradeopen2=tdc(day,tradeopen2),
               demo1=tdc(day,demo1), demo2=tdc(day,demo2),
               major1=tdc(day,major1), major2=tdc(day,major2),
               cinc1=tdc(day,cinc1), cinc2=tdc(day,cinc2),
               mindist=tdc(day,mindist), conttype=tdc(day,conttype), contbinary=tdc(day,contbinary),
               alliance=tdc(day,alliance),igo=tdc(day,jointmemtotal), affinity=tdc(day,affinity),
               powerratio=tdc(day,powerratio))

## save the data
save(durB,durG, file="C:/Users/YULENG/OneDrive/Documents/R Programming/DurationPaper/DurationTimeCov_v2.2.RData")

