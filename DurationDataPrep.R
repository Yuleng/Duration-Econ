## This is an R code for preparing the dataset for duration and econ

## Read in the mids data and take out duration (in days)
## mindur and maxdur caused by missing day; estimate from month
#####################
## old codes where I only use ndy and data where one side quits
####################
## mid <- read.csv("D:/DATA/COW/MID/GML-MID/gml-ndy.csv")[,c("ccode1","ccode2","year","mindur","maxdur","outcome")]
## mid <- na.omit(mid)
## temp1 <- mid[mid$outcome==3,] #yield by side a
## temp2 <- mid[mid$outcome==4,] #yield by side b 
## names(temp2)[1:2] <- c("ccode2","ccode1")
## have a dataset that has ccode 1 being the side that yields
## dur <- merge(temp1, temp2, by = c("ccode1","ccode2", "year","mindur","maxdur"), all=TRUE)[,-c(6:7)]
## rm(mid, temp1, temp2) # remove the raw data
## note that the majority of mindur and maxdur are the same
## hist(dur$maxdur-dur$mindur); sum((dur$maxdur-dur$mindur)>0)
#####################################
## Think more about how to count censoring ASAUR p6 informative vs. noninfo censoring
#####################################
mid <- read.csv("D:/DATA/COW/MID/GML-MID/gml-ddy.csv")[,c("ccode1","ccode2","year","mindur","maxdur","outcome")]
mid <- na.omit(mid)
## simply count side A quitting outcome 3 as not sensoring the rest as censoring
mid$cens <- ifelse(mid$outcome==3, 1, 0) # censor as 0
dur <- mid[,-6]
library(countrycode)
dur$ccode1 <- countrycode(dur$ccode1, "cown", "country.name")
dur$ccode2 <- countrycode(dur$ccode2, "cown", "country.name")

## Read the controls
## now get all the controls
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
## loos like it happen during change of regimes
## in that case, I will simply delete the duplicated ones
polity <- polity[!duplicated(polity$t1),]; polity$t1 <- NULL
## MERGE with dur
dur <- merge(dur, polity[,c("ccode1","year","demo1")], by=c("ccode1","year"), all.x=TRUE)
dur <- merge(dur, polity[,c("ccode2","year","demo2")], by=c("ccode2","year"), all.x=TRUE)

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

dur <- merge(dur, alliance, by=c("ccode1","ccode2","year"), all.x=TRUE)
dur[dur$year>min(alliance$year)-1 & dur$year<max(alliance$year)+1 & is.na(dur$defense),"defense"] <- 0
dur[dur$year>min(alliance$year)-1 & dur$year<max(alliance$year)+1 & is.na(dur$neutrality),"neutrality"] <- 0
dur[dur$year>min(alliance$year)-1 & dur$year<max(alliance$year)+1 & is.na(dur$nonaggression),"nonaggression"] <- 0
dur[dur$year>min(alliance$year)-1 & dur$year<max(alliance$year)+1 & is.na(dur$entente),"entente"] <- 0

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

dur <- merge(dur, contiguity[,c("ccode1","ccode2","year","conttype")], by=c("ccode1","ccode2","year"), all.x=TRUE)
dur[dur$year>min(contiguity$year)-1 & dur$year<max(contiguity$year)+1 & is.na(dur$conttype),"conttype"] <- 0
## IGO Note that it is up to 2005
igo <- IGO[,c(1,3,5,6)]
igo$ccode1 <- countrycode(igo$ccode1,"cown","country.name")
igo$ccode2 <- countrycode(igo$ccode2,"cown","country.name")
## check for duplication
## igo$t1 <- paste(igo$ccode1,igo$ccode2,igo$year)
##  dim(igo[!duplicated(igo$t1),])
dur <- merge(dur, igo, by=c("ccode1","ccode2","year"), all.x=TRUE)
dur[dur$year>min(igo$year)-1 & dur$year<max(igo$year)+1 & is.na(dur$IGOTally),"IGOTally"] <- 0
## national capability up to 2007
cinc <- NationalCapability[,-1]
cinc[cinc$cinc==-9,"cinc"] <- NA 
cinc$cinc1 <- cinc$cinc; cinc$cinc2 <- cinc$cinc
cinc$ccode1 <- countrycode(cinc$ccode,"cown","country.name")
cinc$ccode2 <- countrycode(cinc$ccode,"cown","country.name")
## check for duplication
## cinc$t1 <- paste(cinc$ccode1,cinc$ccode2,cinc$year)
## dim(cinc[!duplicated(cinc$t1),])
dur <- merge(dur, cinc[,c("ccode1","year","cinc1")], by=c("ccode1","year"), all.x=TRUE)
dur <- merge(dur, cinc[,c("ccode2","year","cinc2")], by=c("ccode2","year"), all.x=TRUE)
library(plyr)
## curate the data a bit
dur$alliance <- rowSums(dur[,c("defense","nonaggression","neutrality","entente")],na.rm=TRUE) 
dur <- ddply(dur, .(ccode1,year), mutate, defense1=sum(defense,na.rm=TRUE)-defense)
dur <- ddply(dur, .(ccode2,year), mutate, defense2=sum(defense,na.rm=TRUE)-defense)
dur$powerratio <- log(dur$cinc1/dur$cinc2+1)
dur$contiguity <- ifelse(dur$conttype>0, 1, 0)

## add major power
#############################
library(countrycode)
dur$major1 <- ifelse((dur$ccode1==countrycode(2,"cown","country.name") & dur$year >1897 & dur$year < 2017) | 
                           (dur$ccode1==countrycode(200,"cown","country.name") & dur$year >1815 & dur$year < 2017) |  
                           (dur$ccode1==countrycode(220,"cown","country.name") & dur$year >1815 & dur$year < 1941) |
                           (dur$ccode1==countrycode(220,"cown","country.name") & dur$year >1944 & dur$year < 2017) |
                           (dur$ccode1==countrycode(255,"cown","country.name") & dur$year >1815 & dur$year < 1919) |
                           (dur$ccode1==countrycode(255,"cown","country.name") & dur$year >1924 & dur$year < 1946) |
                           (dur$ccode1==countrycode(255,"cown","country.name") & dur$year >1990 & dur$year < 2017) |
                           (dur$ccode1==countrycode(300,"cown","country.name") & dur$year >1815 & dur$year < 1919) |
                           (dur$ccode1==countrycode(325,"cown","country.name") & dur$year >1859 & dur$year < 1944) |
                           (dur$ccode1==countrycode(365,"cown","country.name") & dur$year >1815 & dur$year < 1918) |
                           (dur$ccode1==countrycode(365,"cown","country.name") & dur$year >1921 & dur$year < 2017) |
                           (dur$ccode1==countrycode(710,"cown","country.name") & dur$year >1949 & dur$year < 2017) |
                           (dur$ccode1==countrycode(740,"cown","country.name") & dur$year >1894 & dur$year < 1946) |
                           (dur$ccode1==countrycode(740,"cown","country.name") & dur$year >1990 & dur$year < 2017),
                         1, 0)
dur$major2 <- ifelse((dur$ccode2==countrycode(2,"cown","country.name") & dur$year >1897 & dur$year < 2017) | 
                           (dur$ccode2==countrycode(200,"cown","country.name") & dur$year >1815 & dur$year < 2017) |  
                           (dur$ccode2==countrycode(220,"cown","country.name") & dur$year >1815 & dur$year < 1941) |
                           (dur$ccode2==countrycode(220,"cown","country.name") & dur$year >1944 & dur$year < 2017) |
                           (dur$ccode2==countrycode(255,"cown","country.name") & dur$year >1815 & dur$year < 1919) |
                           (dur$ccode2==countrycode(255,"cown","country.name") & dur$year >1924 & dur$year < 1946) |
                           (dur$ccode2==countrycode(255,"cown","country.name") & dur$year >1990 & dur$year < 2017) |
                           (dur$ccode2==countrycode(300,"cown","country.name") & dur$year >1815 & dur$year < 1919) |
                           (dur$ccode2==countrycode(325,"cown","country.name") & dur$year >1859 & dur$year < 1944) |
                           (dur$ccode2==countrycode(365,"cown","country.name") & dur$year >1815 & dur$year < 1918) |
                           (dur$ccode2==countrycode(365,"cown","country.name") & dur$year >1921 & dur$year < 2017) |
                           (dur$ccode2==countrycode(710,"cown","country.name") & dur$year >1949 & dur$year < 2017) |
                           (dur$ccode2==countrycode(740,"cown","country.name") & dur$year >1894 & dur$year < 1946) |
                           (dur$ccode2==countrycode(740,"cown","country.name") & dur$year >1990 & dur$year < 2017),
                         1, 0)
## Read Distance Data
load("D:/DATA/DISTANCE/dis.RData")
capdist$ccode1 <- countrycode(capdist$ccode1, "cown", "country.name")
capdist$ccode2 <- countrycode(capdist$ccode2, "cown", "country.name")
centdist$ccode1 <- countrycode(centdist$ccode1, "cown", "country.name")
centdist$ccode2 <- countrycode(centdist$ccode2, "cown", "country.name")
mindist$ccode1 <- countrycode(mindist$ccode1, "cown", "country.name")
mindist$ccode2 <- countrycode(mindist$ccode2, "cown", "country.name")

dur <- merge(dur, capdist, by=c("ccode1","ccode2","year"), all.x=TRUE)
dur <- merge(dur, centdist, by=c("ccode1","ccode2","year"), all.x=TRUE)
dur <- merge(dur, mindist, by=c("ccode1","ccode2","year"), all.x=TRUE)

## Read trade data
load("D:/DATA/TRADE/COW/trade.RData")
durB <- merge(dur, tradeB,by=c("ccode1","ccode2","year"), all.x=TRUE)
durG <- merge(dur, tradeG,by=c("ccode1","ccode2","year"), all.x=TRUE)
rm(tradeB, tradeG) # remove the raw data

## save the data
save(durB,durG, file="C:/Users/YULENG/OneDrive/Documents/R Programming/DurationPaper/Duration.RData")

