# This is an R code to transform the duration data into having time dependent covariates

library(survival)
## read mid data from gml 1816-2010
mid <- read.csv("D:/DATA/COW/MID/GML-MID/gml-ddy.csv")[,c("ccode1","ccode2","year","mindur","maxdur","outcome")]
mid <- na.omit(mid)
## simply count side A quitting outcome 3 as not sensoring the rest as censoring
mid$cens <- ifelse(mid$outcome==3, 1, 0) # censor as 0
mid$status <- 0
for (i in 1:length(mid$outcome)) {
    if (mid$outcome[i] == 3){ # side A quits
      mid$status[i] == 1
  } else if (mid$outcome[i] ==  2){  # side B wins and A quits as a result
    mid$status[i] = 2
  } else if (mid$outcome[i] == 1){ # side A wins
    mid$status[i] = 3 
  } else if (mid$outcome[i] == 4){ # side B quits
    mid$status[i] = 4
  } else if (mid$outcome[i] == 6){ # compromise, both quits
    mid$status[i] = 5
  } else {mid$status[i] == 0}
}
dur <- mid[,-6]
library(countrycode)
dur$ccode1 <- countrycode(dur$ccode1, "cown", "country.name")
dur$ccode2 <- countrycode(dur$ccode2, "cown", "country.name")
dur$id <- paste(dur$ccode1,dur$ccode2,dur$year)

## Read trade data
load("D:/DATA/TRADE/COW/trade.RData")
## creating a dataset to hold covariates from the duration dataset
pool <- dur[,c("id","ccode1","ccode2","year","maxdur")]
hold <- data.frame()
for (i in 1:nrow(pool)){
  temp = pool[i,]
  temp = cbind(temp[,c("id","ccode1","ccode2")], 
               dur=temp$maxdur,
               year=temp$year+0:floor(temp$maxdur/365.25),
               day=365.25*(0:floor(temp$maxdur/365.25))+1)
  hold = rbind(hold,temp); rm(temp)
}
## NOW merge all ind variables
## Read the controls
## polity
library(readxl)
library(countrycode)
polity <- read_excel("D:/DATA/POLITY/p4v2015.xls")[,c("ccode","year","polity2")]
polity$ccode1 <- countrycode(polity$ccode,"p4_ccode","country.name")
## These codes are not matched unambiguously
## 89, 99, 260, 265, 324, 345, 364, 530, 564, 625, 730, 769, 816, 817
polity$ccode1[polity$ccode==89] <- "United Province CA"
polity$ccode1[polity$ccode==99] <- "Colombia" ## Gran Colombia
polity$ccode1[polity$ccode==260] <- "Federal Republic of Germany"
polity$ccode1[polity$ccode==265] <- "German Democratic Republic"
polity$ccode1[polity$ccode==324] <- "Sardinia"
polity$ccode1[polity$ccode==345] <- "Yugoslavia"
polity$ccode1[polity$ccode==364] <- "Russian Federation"
polity$ccode1[polity$ccode==530] <- "Ethiopia"
polity$ccode1[polity$ccode==564] <- "Orange Free State"
polity$ccode1[polity$ccode==625] <- "Sudan"
polity$ccode1[polity$ccode==730] <- "Korea"
polity$ccode1[polity$ccode==769] <- "Pakistan"
polity$ccode1[polity$ccode==816] <- "Viet Nam"
polity$ccode1[polity$ccode==817] <- "Republic of Vietnam"


polity$ccode2 <- polity$ccode1
polity$demo1 <- polity$polity2
polity$demo2 <- polity$polity2
polity$t1 <- paste(polity$ccode1,polity$year)
## note that there are still 6 duplications
## year polity2             ccode1
## 1993       1           Ethiopia
## 1832      -5           Colombia
## 2011      -2              Sudan
## 1922      -7 Russian Federation
## 1976      -7           Viet Nam
## 1991      -5         Yugoslavia
## except for the last two, all have differnt polity scores
## loos like it happen holding change of regimes
## in that case, I will simply delete the duplicated ones
polity <- polity[!duplicated(polity$t1),]; polity$t1 <- NULL
## MERGE with hold
hold <- merge(hold, polity[,c("ccode1","year","demo1")], by=c("ccode1","year"), all.x=TRUE)
hold <- merge(hold, polity[,c("ccode2","year","demo2")], by=c("ccode2","year"), all.x=TRUE)

######################################
## now start merging COWControls
load("D:/DATA/COW/COWControls.RData")
## alliance up to 2012
alliance <- Alliance[,c(1,3,5:9)]
alliance$ccode1 <- countrycode(alliance$ccode1,"cown","country.name")
alliance$ccode2 <- countrycode(alliance$ccode2,"cown","country.name")

library(plyr)
alliance <- ddply(alliance, .(ccode1, ccode2, year), summarize, 
                  defense=sum(defense,na.rm=TRUE), 
                  neutrality=sum(neutrality,na.rm=TRUE),
                  nonaggression=sum(nonaggression,na.rm=TRUE),
                  entente=sum(entente,na.rm=TRUE))
## check if duplicated results
## alliance$t1 <- paste(alliance$ccode1,alliance$ccode2,alliance$year)
## alliance <- alliance[!duplicated(alliance$t1),]; alliance$t1 <- NULL

hold <- merge(hold, alliance, by=c("ccode1","ccode2","year"), all.x=TRUE)
hold[hold$year>min(alliance$year)-1 & hold$year<max(alliance$year)+1 & is.na(hold$defense),"defense"] <- 0
hold[hold$year>min(alliance$year)-1 & hold$year<max(alliance$year)+1 & is.na(hold$neutrality),"neutrality"] <- 0
hold[hold$year>min(alliance$year)-1 & hold$year<max(alliance$year)+1 & is.na(hold$nonaggression),"nonaggression"] <- 0
hold[hold$year>min(alliance$year)-1 & hold$year<max(alliance$year)+1 & is.na(hold$entente),"entente"] <- 0

## contiguity NOTE that it is up too 2006, can expand, maybe by not including "year" in the merge by
## looks like countrycode 348 has no match
contiguity <- Continuity[,c(2,4,6,7)]
contiguity$state1no[contiguity$state1no==348] <- 341
contiguity$state2no[contiguity$state2no==348] <- 341
contiguity$ccode1 <- countrycode(contiguity$state1no,"cown","country.name")
contiguity$ccode2 <- countrycode(contiguity$state2no,"cown","country.name")

contiguity$t1 <- paste(contiguity$ccode1,contiguity$ccode2,contiguity$year)
## 58 duplicated cases; for now, jsut delete the duplicated cases
contiguity <- contiguity[!duplicated(contiguity$t1),]; contiguity$t1 <- NULL

hold <- merge(hold, contiguity[,c("ccode1","ccode2","year","conttype")], by=c("ccode1","ccode2","year"), all.x=TRUE)
hold[hold$year>min(contiguity$year)-1 & hold$year<max(contiguity$year)+1 & is.na(hold$conttype),"conttype"] <- 0
## IGO Note that it is up to 2005
igo <- IGO[,c(1,3,5,6)]
igo$ccode1 <- countrycode(igo$ccode1,"cown","country.name")
igo$ccode2 <- countrycode(igo$ccode2,"cown","country.name")
## check for duplication
## igo$t1 <- paste(igo$ccode1,igo$ccode2,igo$year)
##  dim(igo[!duplicated(igo$t1),])
hold <- merge(hold, igo, by=c("ccode1","ccode2","year"), all.x=TRUE)
hold[hold$year>min(igo$year)-1 & hold$year<max(igo$year)+1 & is.na(hold$IGOTally),"IGOTally"] <- 0
## national capability up to 2007
cinc <- NationalCapability[,-1]
cinc[cinc$cinc==-9,"cinc"] <- NA 
cinc$cinc1 <- cinc$cinc; cinc$cinc2 <- cinc$cinc
cinc$ccode1 <- countrycode(cinc$ccode,"cown","country.name")
cinc$ccode2 <- countrycode(cinc$ccode,"cown","country.name")
## check for duplication
## cinc$t1 <- paste(cinc$ccode1,cinc$ccode2,cinc$year)
## dim(cinc[!duplicated(cinc$t1),])
hold <- merge(hold, cinc[,c("ccode1","year","cinc1")], by=c("ccode1","year"), all.x=TRUE)
hold <- merge(hold, cinc[,c("ccode2","year","cinc2")], by=c("ccode2","year"), all.x=TRUE)
library(plyr)
## curate the data a bit
hold$alliance <- rowSums(hold[,c("defense","nonaggression","neutrality","entente")],na.rm=TRUE) 
hold <- ddply(hold, .(ccode1,year), mutate, defense1=sum(defense,na.rm=TRUE)-defense)
hold <- ddply(hold, .(ccode2,year), mutate, defense2=sum(defense,na.rm=TRUE)-defense)
hold$powerratio <- log(hold$cinc1/hold$cinc2+1)
hold$contiguity <- ifelse(hold$conttype>0, 1, 0)

## add major power
#############################
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
## Read Distance Data
load("D:/DATA/DISTANCE/dis.RData")
capdist$ccode1 <- countrycode(capdist$ccode1, "cown", "country.name")
capdist$ccode2 <- countrycode(capdist$ccode2, "cown", "country.name")
centdist$ccode1 <- countrycode(centdist$ccode1, "cown", "country.name")
centdist$ccode2 <- countrycode(centdist$ccode2, "cown", "country.name")
mindist$ccode1 <- countrycode(mindist$ccode1, "cown", "country.name")
mindist$ccode2 <- countrycode(mindist$ccode2, "cown", "country.name")

hold <- merge(hold, capdist, by=c("ccode1","ccode2","year"), all.x=TRUE)
hold <- merge(hold, centdist, by=c("ccode1","ccode2","year"), all.x=TRUE)
hold <- merge(hold, mindist, by=c("ccode1","ccode2","year"), all.x=TRUE)

holdB <- merge(hold, tradeB,by=c("ccode1","ccode2","year"), all.x=TRUE)
holdG <- merge(hold, tradeG,by=c("ccode1","ccode2","year"), all.x=TRUE)
rm(tradeB, tradeG) # remove the raw data

## now assemble the data
dur$dur <- dur$maxdur
dur <- dur[,c("id","dur","cens")]
dur <- tmerge(dur, dur, id=id, quit=event(dur,cens)) # set the range
durB <- tmerge(dur, holdB, id=id, 
               tradedepend1=tdc(day,tradedepend1), tradedepend2=tdc(day,tradedepend2),
               tradeshare1=tdc(day,tradeshare1), tradeshare2=tdc(day,tradeshare2),
               tradeopen1=tdc(day,tradeopen1), tradeopen2=tdc(day,tradeopen2),
               demo1=tdc(day,demo1), demo2=tdc(day,demo2),
               major1=tdc(day,major1), major2=tdc(day,major2),
               mindist=tdc(day,mindist), contiguity=tdc(day,contiguity),
               alliance=tdc(day,alliance),igo=tdc(IGOTally),
               powerratio=tdc(day,powerratio))
durG <- tmerge(dur, holdG, id=id, 
               tradedepend1=tdc(day,tradedepend1), tradedepend2=tdc(day,tradedepend2),
               tradeshare1=tdc(day,tradeshare1), tradeshare2=tdc(day,tradeshare2),
               tradeopen1=tdc(day,tradeopen1), tradeopen2=tdc(day,tradeopen2),
               demo1=tdc(day,demo1), demo2=tdc(day,demo2),
               major1=tdc(day,major1), major2=tdc(day,major2),
               mindist=tdc(day,mindist), contiguity=tdc(day,contiguity),
               alliance=tdc(day,alliance),igo=tdc(IGOTally),
               powerratio=tdc(day,powerratio))

## save the data
save(durB,durG, file="C:/Users/YULENG/OneDrive/Documents/R Programming/DurationPaper/DurationTimeCov.RData")

