#2015 variance partitioning

######get SPRING data and libraries#####
library(vegan)
library(PerformanceAnalytics)

setwd("~/Documents/Research/PCAP/PCAPdata/Dataframes")
#raw community pattern data:
spplot.patterns.a<-read.csv("spplot.patterns.allforms.csv",header=TRUE,sep=',')
spplot.patterns.a<-spplot.patterns.a[,-1]
spplot.patterns.f<-read.csv("spplot.patterns.forbs.csv",header=TRUE,sep=',')
spplot.patterns.f<-spplot.patterns.f[,-1]
spplot.patterns.o<-read.csv("spplot.patterns.otherforms.csv",header=TRUE,sep=',')
spplot.patterns.o<-spplot.patterns.o[,-1]

#transformed/normalized abiotic data:
spplot.ab.t.z<-read.csv("spplot.ab.t.z.csv",header=TRUE,sep=',')
spplot.ab.t.z<-spplot.ab.t.z[,-1]

#categorical data (reservation and community type):
setwd("~/Documents/Research/PCAP/PCAPdata")
cats<-read.csv("2015 categorical data.csv",header=TRUE,sep=',')
categories<-matrix(as.numeric(as.factor(as.matrix(cats[,3:4]))),ncol=2)
colnames(categories)<-c("comm","res")

########################################################################
#SPRING GAMMA RICHNESS
########################################################################
#make the data frame
concatenate<-data.frame(spplot.ab.t.z,categories,spplot.patterns.f)
#remove plot 10 because it's an outlier in relationship (2), and I want to be 
#able to compare models via AICc.
concatenate<-concatenate[c(1:9,11:29),]

#####partition spring variables EXCLUDING trends along abiotic factors that change with forest####
#mean Ca, cv P, and cv K
chart.Correlation(cbind(concatenate[,c(10)],concatenate[,c(19,25)]))
plot(varpart(concatenate$spgamma.rich,concatenate[,c(10)],concatenate[,c(19,25)],data=concatenate))
#output values: 28, 25, 8, 39

#####partition spring variables INCLUDING trends along abiotic factors that change with forest####

#within Beech-Maple forests:
concatenateBM<-concatenate[concatenate$comm==1,]
chart.Correlation(cbind(concatenateBM[,c(10,12,17)],concatenateBM[,c(19,25)]))

#Choose pair that minimizes unexplained variance:

#plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(10)],concatenateBM[,c(19)],data=concatenateBM))
#res=1.43
#plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(10)],concatenateBM[,c(25)],data=concatenateBM))
#res=0.63
#plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(12)],concatenateBM[,c(19)],data=concatenateBM))
#res=1.25
plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(12)],concatenateBM[,c(25)],data=concatenateBM))
#res=0.18 -> mean calcium and cv potassium
#plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(17)],concatenateBM[,c(19)],data=concatenateBM))
#res=1.63
#plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(17)],concatenateBM[,c(25)],data=concatenateBM))
#res=0.62

#within Floodplain forests:
concatenateF<-concatenate[concatenate$comm==3,]
chart.Correlation(cbind(concatenateF[,c(10,12,17)],concatenateF[,c(19,25)]))

#Choose pair that minimizes unexplained variance:

#plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(10)],concatenateF[,c(19)],data=concatenateF))
#res=0.53
#plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(10)],concatenateF[,c(25)],data=concatenateF))
#res=0.38
#plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(12)],concatenateF[,c(19)],data=concatenateF))
#res=0.61
#plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(12)],concatenateF[,c(25)],data=concatenateF))
#res=0.46
#plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(17)],concatenateF[,c(19)],data=concatenateF))
#res=0.52
plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(17)],concatenateF[,c(25)],data=concatenateF))
#res=0.35-> mean restrictive layer depth and cv potassium

#within Mixed forests
concatenateM<-concatenate[concatenate$comm==5,]
chart.Correlation(cbind(concatenateM[,c(10,12,17)],concatenateM[,c(19,25)]))

#Choose pair that minimizes unexplained variance:

#plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(10)],concatenateM[,c(19)],data=concatenateM))
#res=1.67
#plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(10)],concatenateM[,c(25)],data=concatenateM))
#res=1.44
plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(12)],concatenateM[,c(19)],data=concatenateM))
#res=0.15 -> mean nitrogen and cv phosphorus
#plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(12)],concatenateM[,c(25)],data=concatenateM))
#res=0.01, mean>1
#plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(17)],concatenateM[,c(19)],data=concatenateM))
#res=0.21
#plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(17)],concatenateM[,c(25)],data=concatenateM))
#res=0.0, but mean >1

#within Oak forests:
concatenateO<-concatenate[concatenate$comm==8,]
chart.Correlation(cbind(concatenateO[,c(10,12,17)],concatenateO[,c(19,25)]))

#Choose pair that minimizes unexplained variance:

plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(10)],concatenateO[,c(19)],data=concatenateO))
#res=0.4 -> mean calcium and cv phosphorus
#plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(10)],concatenateO[,c(25)],data=concatenateO))
#res=0.63
#plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(12)],concatenateO[,c(19)],data=concatenateO))
#res=0.50
#plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(12)],concatenateO[,c(25)],data=concatenateO))
#res=0.54
#plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(17)],concatenateO[,c(19)],data=concatenateO))
#res=0.53
#plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(17)],concatenateO[,c(25)],data=concatenateO))
#res=0.58

#####plot within-forest SPRING results######
sp.bm<-c(18,80,29,0)
sp.fp<-c(35,29,5,32)
sp.m<-c(15,0,83,10)
sp.oak<-c(40,15,27,18)
sp.varpartition<-as.matrix(cbind(sp.bm,sp.fp,sp.m,sp.oak))
barplot(sp.varpartition,col=c('black','slategrey','lightgrey','white'))#,legend=c('unexplained','heterogeneity','availability','both'))
 

######get SUMMER data#####
setwd("~/Documents/Research/PCAP/PCAPdata/Dataframes")

#call raw community pattern data
summplot.patterns.a<-read.csv("summplot.patterns.allforms.csv",header=TRUE,sep=',')
summplot.patterns.a<-summplot.patterns.a[,-1]
summplot.patterns.f<-read.csv("summplot.patterns.forbs.csv",header=TRUE,sep=',')
summplot.patterns.f<-summplot.patterns.f[,-1]
summplot.patterns.o<-read.csv("summplot.patterns.otherforms.csv",header=TRUE,sep=',')
summplot.patterns.o<-summplot.patterns.o[,-1]

#call the transformed/normalized abiotic data
summplot.ab.t.z<-read.csv("summplot.ab.t.z.csv",header=TRUE,sep=',')
summplot.ab.t.z<-summplot.ab.t.z[,-1]

#call categorical data reservation and community
setwd("~/Documents/Research/PCAP/PCAPdata")
cats<-read.csv("2015 categorical data.csv",header=TRUE,sep=',')
categories<-matrix(as.numeric(as.factor(as.matrix(cats[,3:4]))),ncol=2)
colnames(categories)<-c("comm","res")
########################################################################
#SUMMER GAMMA RICHNESS
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
#remove 6 so we can compare to (2) and combination models
concatenate<-concatenate[c(1:5,7:29),] 


#####partition summer variables EXCLUDING trends along abiotic factors that change with forest#####
chart.Correlation(cbind(concatenate[,c(2,8,10)],concatenate[,c(23,29)]))
plot(varpart(concatenate$summgamma.rich,concatenate[,c(2,8,10)],concatenate[,c(23,29)],data=concatenate))

#####partition summer variables INCLUDING trends along abiotic factors that change with forest#####
concatenateBM<-concatenate[concatenate$comm==1,]
chart.Correlation(cbind(concatenateBM[,c(2,8,10)],concatenateBM[,c(23,24,29)]))


plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(2)],concatenateBM[,c(23)]))
#.59
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(2)],concatenateBM[,c(24)],data=concatenateBM))
#1.44
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(2)],concatenateBM[,c(29)],data=concatenateBM))
#.96
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(8)],concatenateBM[,c(23)],data=concatenateBM))
#.65
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(8)],concatenateBM[,c(24)],data=concatenateBM))
#1.22
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(8)],concatenateBM[,c(29)],data=concatenateBM))
#.93
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(10)],concatenateBM[,c(23)],data=concatenateBM))
#.84
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(10)],concatenateBM[,c(24)],data=concatenateBM))
#1.28
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(10)],concatenateBM[,c(29)],data=concatenateBM))
#1.05

concatenateF<-concatenate[concatenate$comm==3,]
chart.Correlation(cbind(concatenateF[,c(2,8,10)],concatenateF[,c(23,24,29)]))

plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(2)],concatenateF[,c(23)],data=concatenateF))
#.81
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(2)],concatenateF[,c(24)],data=concatenateF))
#.69
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(2)],concatenateF[,c(29)],data=concatenateF))
#.65
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(8)],concatenateF[,c(23)],data=concatenateF))
#.40
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(8)],concatenateF[,c(24)],data=concatenateF))
#.67
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(8)],concatenateF[,c(29)],data=concatenateF))
#.70
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(10)],concatenateF[,c(23)],data=concatenateF))
#.72
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(10)],concatenateF[,c(24)],data=concatenateF))
#1.02
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(10)],concatenateF[,c(29)],data=concatenateF))
#.86


concatenateM<-concatenate[concatenate$comm==5,]
chart.Correlation(cbind(concatenateM[,c(2,8,10)],concatenateM[,c(23,24,29)]))

plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(2)],concatenateM[,c(23)],data=concatenateM))
#1.51
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(2)],concatenateM[,c(24)],data=concatenateM))
#.86
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(2)],concatenateM[,c(29)],data=concatenateM))
#.56
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(8)],concatenateM[,c(23)],data=concatenateM))
#1.34
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(8)],concatenateM[,c(24)],data=concatenateM))
#.87
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(8)],concatenateM[,c(29)],data=concatenateM))
#.59
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(10)],concatenateM[,c(23)],data=concatenateM))
#1.43
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(10)],concatenateM[,c(24)],data=concatenateM))
#.36
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(10)],concatenateM[,c(29)],data=concatenateM))
#.57

concatenateO<-concatenate[concatenate$comm==8,]
chart.Correlation(cbind(concatenateO[,c(2,8,10)],concatenateO[,c(23,24,29)]))

plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(2)],concatenateO[,c(23)],data=concatenateO))
#.32
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(2)],concatenateO[,c(24)],data=concatenateO))
#1.29
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(2)],concatenateO[,c(29)],data=concatenateO))
#1.08
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(8)],concatenateO[,c(23)],data=concatenateO))
#.30
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(8)],concatenateO[,c(24)],data=concatenateO))
#1.30
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(8)],concatenateO[,c(29)],data=concatenateO))
#1.08
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(10)],concatenateO[,c(23)],data=concatenateO))
#.32
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(10)],concatenateO[,c(24)],data=concatenateO))
#.69
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(10)],concatenateO[,c(29)],data=concatenateO))
#.63

######plot within-forest SUMMER results#####
summ.bm<-c(59,49,4,0)
summ.fp<-c(40,22,31,7)
summ.m<-c(36,96,26,0)
summ.oak<-c(30,82,0,0)
summ.varpartition<-as.matrix(cbind(summ.bm,summ.fp,summ.m,summ.oak))
barplot(summ.varpartition,col=c('black','slategrey','lightgrey','white'))#,legend=c('unexplained','heterogeneity','availability','both'))


#####plot SPRING AND SUMMER results by season (EXCLUDING trends along abiotic factors that change across forests)#####
spring<-c(39,8,28,25)
summer<-c(39,0,51,11)
varpartition<-as.matrix(cbind(spring,summer))
barplot(varpartition,col=c('black','slategrey','lightgrey','white'),legend=c('unexplained','heterogeneity','availability','both'))
